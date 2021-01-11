#pragma once
//#include <ostream>
//#include "Grid.h"
//#include "LocalMatrices.h"
//
//template<typename T>
//class Matrix
//{
//public:
//	int dim = 0;
//	std::vector<T> al, di;
//	std::vector<T>& au = al;
//	std::vector<int> ia;
//	Matrix();
//	~Matrix();
//
//	// Размер в зависимости от количества конечных элементов и базиса
//	void init(int count_elems, int basis, std::vector<double>& x);
//	void display(std::ostream& out);
//	T getElem(int i, int j);
//	T* getElemPtr(int i, int j);
//	void factorization(Matrix<T>& LU);
//	void forward(std::vector<T>& y, std::vector<T>& b);
//	void backward(std::vector<T>& x, const std::vector<T>& y);
//	void global_matrix(grid_in& in, ILocalMatrix& localMatrix);
//	void insert_local(std::vector<std::vector<double>>& l_m, int k);
//};