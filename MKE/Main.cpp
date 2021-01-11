#include "Matrix.cpp"
#include <iostream>
#include <fstream>
#include <iomanip>

// Задается тест, в том же классе можно задать аналитическую функцию для исследований
typedef test2<double> caseName;

#pragma region Вспомогательные функции

// Получить решение одномерной задачи
void run(std::vector<double>& q, grid_in& in)
{
	caseName Functions;
	input(Functions.ToString(), in);

	Matrix<double> m;
	m.solve_FEM(in, Functions, q);
	m.~Matrix();
}

// Получить решение в произвольной точке (в интервале, заданном в краевых условиях)
double get_solve(double x, std::vector<double>& q, grid_in& in)
{
	int k = 0;
	for (int i = 0; i < in.count_nodes; i++)
		if (x >= in.nodes[i])
			k = i;

	// В k элементе находится узел
	// Необходимо найти линейную комбинацию базисных функций
	if (in.basis == 2)
	{
		double ksi = (x - in.nodes[k]) / (in.nodes[k + 1] - in.nodes[k]);
		double fi1 = 2 * (ksi - 0.5)*(ksi - 1),
			fi2 = -4 * ksi*(ksi - 1),
			fi3 = 2*ksi*(ksi - 0.5);
		return fi1 * q[2*k] + fi2 * q[2*k + 1] + fi3 * q[2*k + 2];
	}
	else if (in.basis == 3)
	{
		double ksi = (x - in.nodes[k]) / (in.nodes[k + 1] - in.nodes[k]);
		double fi1 = -4.5 * (ksi - 0.3333333333333333333) * (ksi - 1.) * (ksi - 0.66666666666666666666666),
			fi2 = 13.5 * ksi * (ksi - 0.6666666666666666666)*(ksi - 1.),
			fi3 = -13.5 * ksi * (ksi - 0.333333333333333333333) * (ksi - 1.),
			fi4 = 4.5 * ksi * (ksi - 0.3333333333333333333333) * (ksi - 0.666666666666666666666);
		return fi1 * q[3*k] + fi2 * q[3 * k + 1] + fi3 * q[3 * k + 2] + fi4 * q[3 * k + 3];
	}
	else
		throw new std::invalid_argument("Invalid basis in input");
}

// Вывести погрешность полученного решения
void print_solve_accuracy(std::vector<double>& w, std::vector<double>& q, grid_in& in, std::ostream& out)
{
	caseName Functions;
	for (int i = 0; i < w.size(); i++)
		out << get_solve(w[i], q, in) - Functions.u(w[i]) << " ";
	out << std::endl;
}

// Получить вектор погрешности полученного решения
void get_solve_accuracy(std::vector<double>& w, std::vector<double>& q, grid_in& in, std::string path, std::ostream& out, std::vector<double>& result)
{
	caseName Functions;
	for (int i = 0; i < w.size(); i++)
		result[i] = get_solve(w[i], q, in) - Functions.u(w[i]);
}

// Вывести погрешность полученного решения
void print_solve_nodes(std::vector<double>& w, std::vector<double>& q, grid_in& in, std::ostream& out)
{
	for (int i = 0; i < w.size(); i++)
		out << get_solve(w[i], q, in) << " ";
	out << std::endl;
}

#pragma endregion

int main()
{
	grid_in in;
	std::vector<double> q;

	run(q, in);
	
	setlocale(LC_ALL, "Russian");
	std::cout << "1) Получить решение" << std::endl
			  << "2) Получить решение и погрешность решения в точках" << std::endl;

	int command = 0;
	std::cin >> command;

	switch (command)
	{
	case 1:
	{
		for (int i = 0; i < q.size(); i++)
			std::cout << q[i] << " ";
		std::cout << std::endl;
		q.clear();
		break;
	}
	case 2:
	{
		int count = 0;
		std::cout << "Задайте количество узлов, в которых хотите посчитать погрешность: ";
		std::cin >> count;

		std::vector<double> w;
		for (int i = 0; i < count; i++)
		{
			double node = 0;
			std::cin >> node;
			w.push_back(node);
		}

		std::cout << "Погрешность решения в узлах:" << std::endl;
		print_solve_accuracy(w, q, in, std::cout);

		std::cout << "Полученное решение в узлах:" << std::endl;
		print_solve_nodes(w, q, in, std::cout);

		q.clear();
		w.clear();
		break;
	}
	default:
		std::cout << "Введена неверная команда" << std::endl;
	}

	return 0;
}
