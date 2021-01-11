#pragma once
#define M_PI 3.14159265358979323846

/*
	В данном файле реализован интерфейс для создания класса с входными данными для задачи.
	Желательно назвать каждый класс так же, как и папку, в которой лежат входные файлы.
	Классы наследуем от интерфейса, подаем на вход функции solve_FEM.
*/

template<typename T>
class IInputFunctions
{
public:
	virtual T f(T& x) = 0;
	virtual T lambda(T& x) = 0;
	virtual std::string ToString() = 0;
};

#pragma region Реализованные классы для тестов

// Соловейчик Рояк Персова, страница 140
template<typename T>
class test1 : public IInputFunctions<T>
{
public:
	virtual T f(T& x)
	{
		if (x <= 2)
			return 0;
		else if (x <= 3)
			return 1;
		else return 0.25;
	}

	virtual T lambda(T& x)
	{
		return 1;
	}

	virtual std::string ToString()
	{
		return "test1";
	}

	T u(T& x)
	{
		if (x <= 2)
			return x + 1;
		else if (x <= 3)
			return -(x - 3) * (x - 3) / 2 + 3.5;
		else
			return -(x - 3) * (x - 3) / 8 + 3.5;
	}
};

// Соловейчик Рояк Персова, страница 188
template<typename T>
class test2 : public IInputFunctions<T>
{
public:
	virtual T f(T& x)
	{
		if (x <= 1)
			return -6 * x + 2 * (x * x * x + 7 * x);
		else
			return -20 + (x * x - x + 8);
	}

	virtual T lambda(T& x)
	{
		if (x <= 1)
			return 1;
		else 
			return 10;
	}

	virtual std::string ToString()
	{
		return "test2";
	}

	T u(T& x)
	{
		if (x <= 1)
			return x * x * x + 7 * x;
		else
			return x * x - x + 8;
	}
};

#pragma endregion