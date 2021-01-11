#include <ostream>
#include "Grid.h"
#include "LocalMatrix.h"
#include <iostream>
#include <iomanip>

// Симметричная матрица в профильном формате
template<typename T>
class Matrix
{
private:
    int dim = 0;
    std::vector<T> al, di;
    std::vector<int> ia;
    // В симметричной матрице они одинаковы.
    std::vector<T>& au = al;


    // Инициализировать матрицу для решения МКЭ.
    // Размерность матрицы зависит от количества конечных элементов и базиса.
    void init(int count_elems, int basis, std::vector<T>& x)
    {
        int size = 0;
        size = basis * count_elems + 1;
        x.resize(basis + 1);

        ia.resize(size + 1);
        di.resize(size);
        dim = size;

        // Сколько будет элементов вне диагонали.
        int count = 0;
        for (int i = 1; i < basis + 1; i++)
            count += i;

        al.resize(count * count_elems);
    }

    // Получить указатель на элемент по индексации плотной матрицы
    T* getElemPtr(int i, int j)
    {
        if (i == j)
            return &di[i];
        if (i > j)
            if (j < i - ia[i + 1] + ia[i])
                return nullptr;
            else
                return &al[j - i + ia[i + 1]];
        else
            if (i < j - ia[j + 1] + ia[j])
                return nullptr;
            else
                return &au[i - j + ia[j + 1]];
    }

    // Прямой ход
    // Для решения СЛАУ использовать метод solve_matrix(vector<T>, vector<T>)
    void forward(std::vector<T>& y, const std::vector<T>& b)
    {
        for (int i = 0; i < dim; i++)
        {
            T elem = b[i];
            int i0 = ia[i];
            int i1 = ia[i + 1];

            for (int j = ia[i], k = i1 - i0 < i ? i - (i1 - i0) : 0; j < i1; j++, k++)
                elem -= al[j] * b[k];

            elem /= di[i];
            y[i] = elem;
        }
    }

    // Обратный ход
    // Для решения СЛАУ использовать метод solve_matrix(vector<T>, vector<T>)
    void backward(std::vector<T>& x, const std::vector<T>& y)
    {
        for (int i = dim - 1; i >= 0; i--)
        {
            T sum = 0;
            T xi = x[i] = y[i] / di[i];

            for (int j = dim - 1; j > i; j--)
            {
                int bias = ia[j + 1] - ia[j] - j + i;
                if (bias >= 0)
                    sum += al[ia[j] + bias] * y[j];
            }
            sum /= di[i];
            x[i] -= sum;
        }
    }

    // Собрать глобальную матрицу.
    void global_matrix(grid_in& in, ILocalMatrix<T>& localMatrix, ILocalVector<T>& localVector, std::vector<T>& b)
    {
        // Размер глобальной матрицы определяется базисными функциями
        // И количеством конечных элементов.
        // Нужно для каждого элемента построить локальную матрицу
        // И внести её в глобальную, учесть краевые и, 
        // если матрица получилась несимметричной, то сделать её симметричной.
        std::vector<T> x;

        init(in.count_elements, in.basis, x);
        b.resize(this->dim);

        // Сборка глобальной матрицы
        for (int k = 0; k < in.count_elements; k++)
        {
            int num_material = in.elems[k];
            T x0 = in.nodes[k], h = (in.nodes[k + 1] - in.nodes[k]) / in.basis;

            // Вычисление иксов 
            for (int i = 0; i <= in.basis; i++)
                x[i] = x0 + i * h;

            if (x[0] == 0) x[0] += 1e-14;
            else x[0] += pow(10, int(log10(x[0])) - 14);

            if (x[in.basis] == 0) x[in.basis] -= 1e-14;
            else x[in.basis] -= pow(10, int(log10(x[0])) - 14);

            // Внесение локальной матрицы в глобальную
            std::vector<std::vector<T>>* l_m = localMatrix.get_matrix(x, in.materials[num_material]);
            insert_local(*l_m, k);

            // Внесение локального вектора в глобальный
            std::vector<T>* l_v = localVector.get_vector(x);
            int size = in.basis + 1;
            for (int i = 0; i < size; i++)
                b[k * in.basis + i] += l_v->at(i);
        }
    }

    // Учесть краевые условия всех родов.
    void conditions(grid_in& in, std::vector<double>& b)
    {
        int first_border = std::get<0>(in.r_cond);
        int second_border = std::get<1>(in.r_cond);
        // Индекс в массиве с коэффициентами для краевых условий.
        int index = 0;

        // Правильно учитывать сначала вторые и третьи(естественные), только потом главные(первые).
        if (first_border == 2)
            b[0] += in.conditions[index++];
        else if (first_border == 3)
        {
            this->di[0] += in.conditions[index++];  // beta
            b[0] += in.conditions[index++];         // ubeta
        }

        // Учитываем вторые или третьи на правой границе
        if (second_border == 2)
            b[dim - 1] += in.conditions[index++];
        else if (second_border == 3)
        {
            double beta = in.conditions[index++];
            this->di[dim - 1] += beta;
            b[dim - 1] += beta * in.conditions[index++];
        }

        if (first_border == 1)
        {
            // Зануляем первую строку, приводим к симметричному виду
            this->di[0] = 1;
            b[0] = in.conditions[index++];

            for (int i = 0; i < in.basis; i++)
            {
                b[i + 1] += -this->al[this->ia[i + 1]] * b[0];
                // Не будем менять профиль, чтобы не делать n операций
                this->al[this->ia[i + 1]] = 0;
            }
        }

        if (second_border == 1)
        {
            this->di[dim - 1] = 1;
            b[dim - 1] = in.conditions[index];

            // Смотрим элементы на последней строке
            for (int i = 0; i < in.basis; i++)
                b[dim - i - 2] += -this->al[this->ia[dim - 1] + in.basis - 1 - i] * b[dim - 1];

            this->ia[dim] -= in.basis;
        }
    }

    // Вставить локальную матрицу в глобальную на позицию k-го конечного элемента.
    void insert_local(std::vector<std::vector<T>>& l_m, int k)
    {
        // Индекс массива al для заполнения глобальной матрицы.
        static int stop = 0;

        // Размерность локальной матрицы.
        int size = l_m.size();

        for (int i = 0; i < size; i++)
        {
            if (k == 0)
            {
                ia[i + 1] = ia[i] + i;
                di[i] += l_m[i][i];
            }
            else
            {
                // Если никаких внедиагональных элементов не вносится, то профиль не меняется.
                if (i != 0)
                    ia[k * size + i - k + 1] = ia[k * size + i - k] + i;

                di[k * size + i - k] += l_m[i][i];
            }

            for (int j = 0; j < i; j++)
                al[stop++] = l_m[i][j];
        }
    }


public:
    ~Matrix()
    {
        al.clear();
        di.clear();
        ia.clear();
    }

    // LLT - измененная матрица.
    // Метод факторизует матрицу LLT по разложению Холецкого.
    // Не меняет исходную матрицу.
    // Для решения СЛАУ использовать метод solve_matrix(vector<T>, vector<T>)
    void factorization(Matrix<T>& LLT)
    {
        for (int i = 0; i < dim; i++)
        {
            // i - index of row
            // j - index of column
            int i0 = ia[i];
            int i1 = ia[i + 1];
            int j = i - (i1 - i0);
            T sum_di = 0;

            for (int k = i0; k < i1; k++, j++)
            {
                int j0 = ia[j];
                int j1 = ia[j + 1];

                int ki = i0;
                int kj = j0;

                int kui = k - i0;
                int kuj = j1 - j0;
                int kur = kui - kuj;

                if (kur > 0)
                    ki += kur;
                else
                    kj -= kur;

                T sum = 0;
                for (; ki < k; ki++, kj++)
                    sum += al[ki] * al[kj];

                LLT.al[k] = (al[k] - sum) / LLT.di[j];
                sum_di += LLT.al[k] * LLT.al[k];
            }
            LLT.di[i] = sqrt(di[i] - sum_di);
        }
    }

    // Вывести матрицу в плотном формате
    void display(std::ostream& out)
    {
        for (int i = 0; i < dim; i++)
        {
            for (int j = 0; j < dim; j++)
                out << getElem(i, j) << " ";
            out << std::endl;
        }
    }

    // Получить элемент по индексации плотной матрицы
    T getElem(int i, int j)
    {
        if (i == j)
            return di[i];
        if (i > j)
            if (j < i - ia[i + 1] + ia[i])
                return 0;
            else
                return al[j - i + ia[i + 1]];
        else
            if (i < j - ia[j + 1] + ia[j])
                return 0;
            else
                return au[i - j + ia[j + 1]];
    }

    // Решить уравнение вида Ax = b
    // х - полученное решение той же размерности, что и b.
    void solve_matrix(Matrix<T>& A, std::vector<T>& b, std::vector<T>& x)
    {
        factorization(A);
        forward(x, b);
        backward(x, x);
    }

    // Решение задачи МКЭ для уравнения вида 
    // -div(lambda(x)*grad(u(x))) + gamma * u(x) = f(x)
    // q - полученное решение
    // FEM - Finite Element Method
    void solve_FEM(grid_in& in, IInputFunctions<T>& Functions, std::vector<T>& q)
    {
        // Собираем глобальную матрицу
        if (in.basis == 2)
            global_matrix(in, *(new LocalMatrix2_lambda<T>(Functions)), *(new LocalVector2<T>(Functions)), q);
        else if (in.basis == 3)
            global_matrix(in, *(new LocalMatrix3_lambda<T>(Functions)), *(new LocalVector3<T>(Functions)), q);
        else
            throw new std::invalid_argument("Invalid basis in input");

        // Учитываем краевые условия
        conditions(in, q);

        // Решаем СЛАУ
        solve_matrix(*this, q, q);
    }
};