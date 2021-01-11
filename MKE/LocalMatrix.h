#pragma once
#include <vector>
#include <tuple>
#include "Functions.h"

/*
    ILocalMatrix и ILocalVector - интерфейсы, необходимые для решения МКЭ.
    В этом файле реализованы основные классы для квадратичных и кубических базисов 
    с разложением коэффициента лямбда и без него.
    Чтобы использовать другую локальную матрицу - необходимо унаследовать новый класс от ILocaxMatrix и 
    реализовать конструктор с заполнением коэффициентов.
    Аналогично с локальным вектором.
    
    При реализации локального вектора необходимо в конструкторе задать объект класса, 
    унаследованного от IInputFunctions для задания правой части.
*/

template<typename T>
class ILocalMatrix
{
public:
    ILocalMatrix<T>() {}
    // Локальная матрица
    virtual std::vector<std::vector<T>>* get_matrix(std::vector<T>& x, std::tuple<double, double>& mat) = 0;
};

template<typename T>
class ILocalVector
{
public:
    // Локальный вектор
    virtual std::vector<T>* get_vector(std::vector<T>& x) = 0;
};

#pragma region Локальные матрицы и вектора без разложения базисных функций

template<typename T>
class LocalMatrix2 : public ILocalMatrix<T>
{
public:
    std::vector<std::vector<T>> M;
    std::vector<std::vector<T>> G;
    std::vector<std::vector<T>> Res;

    LocalMatrix2<T>()
    {
        M.resize(3);
        G.resize(3);
        Res.resize(3);
        Res[0].resize(3);
        Res[1].resize(3);
        Res[2].resize(3);

        M[0] = { 4, 2, -1 };
        M[1] = { 2, 16, 2 };
        M[2] = { -1, 2, 4 };

        G[0] = { 7, -8, 1 };
        G[1] = { -8, 16, -8 };
        G[2] = { 1, -8, 7 };
    }

    // Локальная матрица. mat - лямбда и гамма
    virtual std::vector<std::vector<T>>* get_matrix(std::vector<T>& x, std::tuple<double, double>& mat)
    {
        T h = x[x.size() - 1] - x[0];
        T coefG = std::get<0>(mat) / (3 * h);
        T coefM = std::get<1>(mat) * h / 30;

        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                Res[i][j] = G[i][j] * coefG + M[i][j] * coefM;

        return &Res;
    }
};

template<typename T>
class LocalMatrix3 : public ILocalMatrix<T>
{
public:
    std::vector<std::vector<T>> M;
    std::vector<std::vector<T>> G;
    std::vector<std::vector<T>> Res;

    LocalMatrix3<T>()
    {
        M.resize(4);
        G.resize(4);
        Res.resize(4);
        Res[0].resize(4);
        Res[1].resize(4);
        Res[2].resize(4);
        Res[3].resize(4);

        M[0] = { 128, 99, -36, 19 };
        M[1] = { 99, 648, -81, -36 };
        M[2] = { -36, -81, 648, 99 };
        M[3] = { 19, -36, 99, 128 };

        G[0] = { 148, -189, 54, -13 };
        G[1] = { -189, 432, -297, 54 };
        G[2] = { 54, -297, 432, -189 };
        G[3] = { -13, 54, -189, 148 };
    }

    // Локальная матрица. mat - лямбда и гамма
    virtual std::vector<std::vector<T>>* get_matrix(std::vector<T>& x, std::tuple<double, double>& mat)
    {
        T h = x[x.size() - 1] - x[0];
        T coefG = std::get<0>(mat) / (40 * h);
        T coefM = std::get<1>(mat) * h / 1680;

        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                Res[i][j] = G[i][j] * coefG + M[i][j] * coefM;

        return &Res;
    }
};

template<typename T>
class LocalVector2 : public ILocalVector<T>
{
public:
    std::vector<T> v;
    IInputFunctions<T>* Functions;

    LocalVector2<T>(IInputFunctions<T>& Functions)
    {
        v.resize(3);
        this->Functions = &Functions;
    }

    virtual std::vector<T>* get_vector(std::vector<T>& x)
    {
        T h = x[x.size() - 1] - x[0];
        T f1 = Functions->f(x[0]), f2 = Functions->f(x[1]), f3 = Functions->f(x[2]);

        v[0] = h / 30 * (4 * f1 + 2 * f2 - f3);
        v[1] = h / 30 * (2 * f1 + 16 * f2 + 2 * f3);
        v[2] = h / 30 * (-f1 + 2 * f2 + 4 * f3);

        return &v;
    }
};

template<typename T>
class LocalVector3 : public ILocalVector<T>
{
public:
    std::vector<T> v;
    IInputFunctions<T>* Functions;
    LocalVector3<T>(IInputFunctions<T>& Functions)
    {
        v.resize(4);
        this->Functions = &Functions;
    }

    virtual std::vector<T>* get_vector(std::vector<T>& x)
    {
        T h = x[x.size() - 1] - x[0];
        T f1 = Functions->f(x[0]), f2 = Functions->f(x[1]), f3 = Functions->f(x[2]), f4 = Functions->f(x[3]);

        v[0] = h / 1680 * (128 * f1 + 99 * f2 - 36 * f3 + 19 * f4);
        v[1] = h / 1680 * (99 * f1 + 648 * f2 - 81 * f3 - 36 * f4);
        v[2] = h / 1680 * (-36 * f1 - 81 * f2 + 648 * f3 + 99 * f4);
        v[3] = h / 1680 * (19 * f1 - 36 * f2 + 99 * f3 + 128 * f4);

        return &v;
    }
};

#pragma endregion

#pragma region Локальные матрицы с разложением коэффициента лямбда по базису
// Для них вектора правой части не меняются

// По квадратичному базису
template<typename T>
class LocalMatrix2_lambda : public ILocalMatrix<T>
{
public:
    std::vector<std::vector<T>> M;
    std::vector<std::vector<T>> G;
    IInputFunctions<T>* Functions;

    LocalMatrix2_lambda<T>(IInputFunctions<T>& Functions)
    {
        M.resize(3);
        G.resize(3);

        M[0] = { 4, 2, -1 };
        M[1] = { 2, 16, 2 };
        M[2] = { -1, 2, 4 };

        G[0].resize(3);
        G[1].resize(3);
        G[2].resize(3);

        this->Functions = &Functions;
    }

    // Локальная матрица. mat - лямбда и гамма
    virtual std::vector<std::vector<T>>* get_matrix(std::vector<T>& x, std::tuple<double, double>& mat)
    {
        T h = x[x.size() - 1] - x[0];
        T coefG = 1 / h;
        double L1 = Functions->lambda(x[0]), L2 = Functions->lambda(x[1]), L3 = Functions->lambda(x[2]);
        G[0][0] = coefG * (-0.1 * L3 + 1.2 * L2 + 1.2333333333333333333 * L1);
        G[0][1] = coefG * (-1.4666666666666666667 * L1 - 1.0666666666666666667 * L2 - 0.13333333333333333333 * L3);
        G[0][2] = coefG * (0.23333333333333333333 * L1 + 0.23333333333333333333 * L3 - 0.13333333333333333333 * L2);
        G[1][0] = G[0][1];
        G[1][1] = coefG * (1.6 * L1 + 1.6 * L3 + 2.1333333333333333333 * L2);
        G[1][2] = coefG * (-1.4666666666666666667 * L3 + -0.13333333333333333333 * L1 + -1.0666666666666666667 * L2);
        G[2][0] = G[0][2];
        G[2][1] = G[1][2];
        G[2][2] = coefG * (-0.1 * L1 + 1.2 * L2 + 1.2333333333333333333 * L3);

        T coefM = std::get<1>(mat) * h / 30;

        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                G[i][j] += M[i][j] * coefM;

        return &G;
    }
};

// По кубическому
template<typename T>
class LocalMatrix3_lambda : public ILocalMatrix<T>
{
public:
    std::vector<std::vector<T>> M;
    std::vector<std::vector<T>> G;
    IInputFunctions<T>* Functions;

    LocalMatrix3_lambda<T>(IInputFunctions<T>& Functions)
    {
        M.resize(4);
        G.resize(4);

        M[0] = { 128, 99, -36, 19 };
        M[1] = { 99, 648, -81, -36 };
        M[2] = { -36, -81, 648, 99 };
        M[3] = { 19, -36, 99, 128 };

        G[0].resize(4);
        G[1].resize(4);
        G[2].resize(4);
        G[3].resize(4);

        this->Functions = &Functions;
    }

    // Локальная матрица. 
    // mat - лямбда и гамма соответственно.
    virtual std::vector<std::vector<T>>* get_matrix(std::vector<T>& x, std::tuple<double, double>& mat)
    {
        T h = x[x.size() - 1] - x[0];
        T coefG = 1 / h; 
        double L1 = Functions->lambda(x[0]), L2 = Functions->lambda(x[1]), L3 = Functions->lambda(x[2]), L4 = Functions->lambda(x[3]);

        G[0][0] = coefG * (0.18258928571428571429 * L4 + 2.0263392857142857143 * L2 + 2.140625 * L1 + -0.64955357142857142857 * L3);
        G[0][1] = coefG * (-3.0147321428571428571 * L1 - 1.8441964285714285714 * L2 + 0.47008928571428571429 * L3 - 0.33616071428571428571 * L4);
        G[0][2] = coefG * (0.38705357142857142857 * L4 + 0.10848214285714285714 * L3 - 0.253125 * L2 + 1.1075892857142857143 * L1);
        G[0][3] = coefG * (-0.23348214285714285714 * L1 - 0.23348214285714285714 * L4 + 0.070982142857142857143 * L2 + 0.070982142857142857143 * L3);
        G[1][0] = G[0][1];
        G[1][1] = coefG * (0.87991071428571428571 * L4 + 1.8441964285714285714 * L3 + 3.796875 * L2 + 4.2790178571428571429 * L1);
        G[1][2] = coefG * (-1.6513392857142857143 * L1 - 1.6513392857142857143 * L4 - 2.0611607142857142857 * L2 + -2.0611607142857142857 * L3);
        G[1][3] = coefG * (0.38705357142857142857 * L1 + 0.10848214285714285714 * L2 - 0.253125 * L3 + 1.1075892857142857143 * L4);
        G[2][0] = G[0][2];
        G[2][1] = G[1][2];
        G[2][2] = coefG * (0.87991071428571428571 * L1 + 1.8441964285714285714 * L2 + 3.796875 * L3 + 4.2790178571428571429 * L4);
        G[2][3] = coefG * (-3.0147321428571428571 * L4 - 1.8441964285714285714 * L3 + 0.47008928571428571429 * L2 - 0.33616071428571428571 * L1);
        G[3][0] = G[0][3];
        G[3][1] = G[1][3];
        G[3][2] = G[2][3];
        G[3][3] = coefG * (0.18258928571428571429 * L1 + 2.0263392857142857143 * L3 + 2.140625 * L4 - 0.64955357142857142857 * L2);

        T coefM = std::get<1>(mat) * h / 1680;

        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                G[i][j] += M[i][j] * coefM;

        return &G;
    }
};

#pragma endregion