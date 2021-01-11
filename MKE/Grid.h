#pragma once
#include <vector>
#include <string>
#include <tuple>

struct grid_in
{
	int count_nodes, count_elements, count_materials, basis;
	std::vector<double> nodes;

	// ����� ��������� ��� ���������������� ��������.
	std::vector<int> elems;

	// ������ ����� - ����� �������. ������ - ������.
	// �������� ������ ������� �������.
	std::tuple<int, int> r_cond;

	// �������� �����, ����������� ��� ����� ������� �������
	std::vector<double> conditions;

	// ��������� - get<����� �������>(����������).
	// ������, �����.
	std::vector<std::tuple<double, double>> materials;
};

void input(std::string path, grid_in& out);