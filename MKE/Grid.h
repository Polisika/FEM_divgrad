#pragma once
#include <vector>
#include <string>
#include <tuple>

struct grid_in
{
	int count_nodes, count_elements, count_materials, basis;
	std::vector<double> nodes;

	// Номер материала для соответствующего элемента.
	std::vector<int> elems;

	// Первое число - левая граница. Второе - правая.
	// Написаны номера краевых условий.
	std::tuple<int, int> r_cond;

	// Записаны числа, необходимые для учёта краевых условий
	std::vector<double> conditions;

	// Обращение - get<номер позиции>(переменная).
	// Лямбда, гамма.
	std::vector<std::tuple<double, double>> materials;
};

void input(std::string path, grid_in& out);