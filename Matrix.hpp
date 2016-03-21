#ifndef Matrix_hpp_
#define Matrix_hpp_

#include <vector>
#include <algorithm>
#include <initializer_list>
#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <random>

template<typename T>
class Matrix
{
public:
	// constructs an empty 0x0 matrix 
	Matrix() : data_(nullptr), row_(0), col_(0) {}
	// constructs a matrix of given size with default value in each element
	Matrix(const int r, const int c) : data_((T*)calloc(r * c, sizeof(T))), row_(r), col_(c) {}
	// constructs a matrix of given size with given value in each element
	Matrix(const int r, const int c, const T& def) : data_((T*)malloc(r*c*sizeof(T))), row_(r), col_(c) { for (int i = 0; i != row_ * col_; ++i) data_[i] = def; }
	// constructs a matrix from initializer list of initializer lists
	Matrix(const std::initializer_list<std::initializer_list<T> >& il) : row_(il.size()), col_(il.begin()->size())
	{
		data_ = (T*)malloc(row_ * col_ * sizeof(T));
		int i = -1;

		for (auto it = il.begin(); it != il.end(); it++)
			for (auto it2 = it->begin(); it2 != it->end(); it2++)
				data_[++i] = *it2;
	}
	// constructs a matrix from input file, expected format is: number of rows, number of columns, rows*cols elements of type T
	explicit Matrix(const std::string& input)
	{
		std::ifstream iFile(input);
		iFile >> row_ >> col_;
		data_ = (T*)malloc(row_ * col_ * sizeof(T));

		for (int index = 0; index < row_ * col_; ++index)
			iFile >> data_[index];

		iFile.close();
	}
	// constructs a matrix of size given in r and c from input file, expected format is: rows*cols elements of type T
	Matrix(const int r, const int c, const std::string& input) : data_((T*)malloc(r * c * sizeof(T))), row_(r), col_(c)
	{
		std::ifstream iFile(input);

		for (int index = 0; index < row_ * col_; ++index)
			iFile >> data_[index];

		iFile.close();
	}

	Matrix(const Matrix<T>& m) : data_((T*)malloc(m.row_ * m.col_ * sizeof(T))), row_(m.row_), col_(m.col_) { for (int i = 0; i != row_ * col_; ++i) data_[i] = m.data_[i]; }
	Matrix(Matrix<T>&& m) : data_(m.data_), row_(m.row_), col_(m.col_) { m.data_ = nullptr; /* I've stolen the data */ }
	Matrix& operator=(const Matrix<T>& m)
	{
		row_ = m.row_;
		col_ = m.col_;
		data_ = (T*)malloc(row_ * col_ * sizeof(T));

		for (int i = 0; i != row_ * col_; ++i)
			data_[i] = m.data_[i];

		return *this;
	}
	Matrix& operator=(Matrix<T>&& m)
	{
		row_ = m.row_;
		col_ = m.col_;
		data_ = m.data_;
		m.data_ = nullptr;
		return *this;
	}
	~Matrix() { free(data_); data_ = nullptr; /* to be sure */ }

	// element access, non-const and const versions
	T& at(const int i, const int j)
	{
		if (i >= row_ || j >= col_)
			throw std::out_of_range("Index out of bounds.");

		return data_[i * col_ + j];
	}
	const T& at(const int i, const int j) const
	{
		if (i >= row_ || j >= col_)
			throw std::out_of_range("Index out of bounds.");

		return data_[i * col_ + j];
	}
	T& at(const std::pair<int, int>& index)
	{
		if (index.first >= row_ || index.second >= col_)
			throw std::out_of_range("Index out of bounds.");

		return data_[index.first * col_ + index.second];
	}
	const T& at(const std::pair<int, int>& index) const
	{
		if (index.first >= row_ || index.second >= col_)
			throw std::out_of_range("Index out of bounds.");

		return data_[index.first * col_ + index.second];
	}
	T& flip(const int i, const int j)
	{
		if (i >= row_ || j >= col_)
			throw std::out_of_range("Index out of bounds.");

		data_[i * col_ + j] = data_[i * col_ + j] ? false : true;
		return data_[i * col_ + j];
	}

	// returns the number of rows
	int getRow() const	{ return row_; }
	// returns the number of columns
	int getCol() const	{ return col_; }
	// returns the number of one-entrie
	int getOnes() const
	{
		int ones = 0;

		for (int i = 0; i != row_ * col_; ++i)
			if (data_[i])
				++ones;

		return ones;
	}
	// returns a formatted string containing the matrix
	std::string Print() const
	{
		std::stringstream out;
		out.precision(3);

		for (int i = 0; i < row_; i++)
		{
			for (int j = 0; j < col_; j++)
			{
				if (j != 0)
					out << " ";

				out << data_[i * col_ + j];
			}

			out << "\n";
		}

		return out.str();
	}
private:
	T* data_;
	int row_, col_;
};

#endif