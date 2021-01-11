#include "Grid.h"
#include <fstream>

// ѕуть указываетс€ до директории.  оса€ черта в конце не нужна.
// info.txt - информаци€ с конечных элементов , количеством узлов, материалов,
// информаци€ по базисам (использовать квадратичный (2) или кубический (3))
// номерами краевых условий на левой и правой границах.
// ћатериалы нумеруютс€ с нул€.
// ћассив с узлами должен быть монотонно возрастаюсщим.
void input(std::string path, grid_in& out)
{
	std::ifstream info(path + "/info.txt");

	info >> out.count_elements >> out.count_nodes >> out.count_materials;
	info >> out.basis;

	if (out.count_elements + 1 != out.count_nodes)
		throw new std::invalid_argument("Count nodes have to be equals count elements plus one. For example: 4 5 ...");

	out.elems.resize(out.count_elements);
	out.nodes.resize(out.count_nodes);
	out.materials.resize(out.count_materials);

	int condleft, condright;
	info >> condleft >> condright;
	if (condleft < 1 || condleft > 3)
		throw new std::invalid_argument("Left condition have to be in range [1, 3]");
	if (condright < 1 || condright > 3)
		throw new std::invalid_argument("RIght condition have to be in range [1, 3]");

	out.r_cond = std::make_tuple(condleft, condright);
	info.close();

	// ¬вод значений по краевым условий.
	std::ifstream conditions(path + "/conditions.txt");
	out.conditions.resize(4);
	int index = 0;

	// —начала числа дл€ вторых и третьих, потом дл€ первых краевых условий.
	if (condleft == 2)
		conditions >> out.conditions[index++];
	else if (condleft == 3)
	{
		conditions >> out.conditions[index++];
		conditions >> out.conditions[index++];
	}

	if (condright == 2)
		conditions >> out.conditions[index++];
	else if (condright == 3)
	{
		conditions >> out.conditions[index++];
		conditions >> out.conditions[index++];
	}

	if (condleft == 1)
		conditions >> out.conditions[index++];

	if (condright == 1)
		conditions >> out.conditions[index++];

	conditions.close();

	std::ifstream nodes(path + "/nodes.txt");
	for (int i = 0; i < out.count_nodes; i++)
		nodes >> out.nodes[i];
	nodes.close();

	std::ifstream elements(path + "/elements.txt");
	for (int i = 0; i < out.count_elements; i++)
		elements >> out.elems[i];

	std::ifstream mat(path + "/materials.txt");
	for (int i = 0; i < out.count_materials; i++)
	{
		double lambda, gamma;
		mat >> lambda >> gamma;
		out.materials[i] = std::make_tuple(lambda, gamma);
	}
}