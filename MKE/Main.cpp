#include "Matrix.cpp"
#include <iostream>
#include <fstream>
#include <iomanip>

// �������� ����, � ��� �� ������ ����� ������ ������������� ������� ��� ������������
typedef test2<double> caseName;

#pragma region ��������������� �������

// �������� ������� ���������� ������
void run(std::vector<double>& q, grid_in& in)
{
	caseName Functions;
	input(Functions.ToString(), in);

	Matrix<double> m;
	m.solve_FEM(in, Functions, q);
	m.~Matrix();
}

// �������� ������� � ������������ ����� (� ���������, �������� � ������� ��������)
double get_solve(double x, std::vector<double>& q, grid_in& in)
{
	int k = 0;
	for (int i = 0; i < in.count_nodes; i++)
		if (x >= in.nodes[i])
			k = i;

	// � k �������� ��������� ����
	// ���������� ����� �������� ���������� �������� �������
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

// ������� ����������� ����������� �������
void print_solve_accuracy(std::vector<double>& w, std::vector<double>& q, grid_in& in, std::ostream& out)
{
	caseName Functions;
	for (int i = 0; i < w.size(); i++)
		out << get_solve(w[i], q, in) - Functions.u(w[i]) << " ";
	out << std::endl;
}

// �������� ������ ����������� ����������� �������
void get_solve_accuracy(std::vector<double>& w, std::vector<double>& q, grid_in& in, std::string path, std::ostream& out, std::vector<double>& result)
{
	caseName Functions;
	for (int i = 0; i < w.size(); i++)
		result[i] = get_solve(w[i], q, in) - Functions.u(w[i]);
}

// ������� ����������� ����������� �������
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
	std::cout << "1) �������� �������" << std::endl
			  << "2) �������� ������� � ����������� ������� � ������" << std::endl;

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
		std::cout << "������� ���������� �����, � ������� ������ ��������� �����������: ";
		std::cin >> count;

		std::vector<double> w;
		for (int i = 0; i < count; i++)
		{
			double node = 0;
			std::cin >> node;
			w.push_back(node);
		}

		std::cout << "����������� ������� � �����:" << std::endl;
		print_solve_accuracy(w, q, in, std::cout);

		std::cout << "���������� ������� � �����:" << std::endl;
		print_solve_nodes(w, q, in, std::cout);

		q.clear();
		w.clear();
		break;
	}
	default:
		std::cout << "������� �������� �������" << std::endl;
	}

	return 0;
}
