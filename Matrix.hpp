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

template< typename T>
class Matrix
{
public:
	// constructs an empty 0x0 matrix 
	Matrix() : row_(0), col_(0) {}
	// constructs a matrix of given size with default value in each element
	Matrix(size_t r, size_t c) : data_(r*c), row_(r), col_(c) {}
	// constructs a matrix of given size with given value in each element
	Matrix(size_t r, size_t c, const T& def) : data_(r*c, def), row_(r), col_(c) {}
	// constructs a matrix from initializer list of initializer lists
	Matrix(const std::initializer_list<std::initializer_list<T> >& il) : row_(il.size()), col_(il.begin()->size())
	{
		data_.reserve(row_ * col_);

		for (auto it = il.begin(); it != il.end(); it++)
			for (auto it2 = it->begin(); it2 != it->end(); it2++)
				data_.push_back(*it2);
	}
	// constructs a matrix from input file, expected format is: number of rows, number of columns, rows*cols elements of type T
	explicit Matrix(const std::string& input)
	{
		std::ifstream iFile(input);
		iFile >> row_ >> col_;
		data_.resize(row_ * col_);

		for (auto &it : data_)
			iFile >> it;

		iFile.close();
	}
	// constructs a matrix of size given in r and c from input file, expected format is: rows*cols elements of type T
	Matrix(size_t r, size_t c, const std::string& input) : data_(r*c), row_(r), col_(c)
	{
		std::ifstream iFile(input);

		for (auto &it : data_)
			iFile >> it;

		iFile.close();
	}

	Matrix(const Matrix<T>& m) : data_(m.data_.begin(), m.data_.end()), row_(m.row_), col_(m.col_) {}
	Matrix(Matrix<T>&& m) : row_(m.row_), col_(m.col_) { data_ = std::move(m.data_); }
	Matrix& operator=(const Matrix<T>& m)
	{
		row_ = m.row_;
		col_ = m.col_;
		data_.resize(row_ * col_);
		std::copy(m.data_.begin(), m.data_.end(), data_.begin());
		return *this;
	}
	Matrix& operator=(Matrix<T>&& m)
	{
		row_ = m.row_;
		col_ = m.col_;
		data_ = std::move(m.data_);
		return *this;
	}

	// element access, non-const and const versions
	T& at(const size_t i, const size_t j)
	{
		if (i >= row_ || j >= col_)
			throw std::out_of_range("Index out of bounds.");

		return data_[i * col_ + j];
	}
	const T& at(const size_t i, const size_t j) const
	{
		if (i >= row_ || j >= col_)
			throw std::out_of_range("Index out of bounds.");

		return data_[i * col_ + j];
	}
	T& at(const std::pair<size_t, size_t>& index)
	{
		if (index.first >= row_ || index.second >= col_)
			throw std::out_of_range("Index out of bounds.");

		return data_[index.first * col_ + index.second];
	}
	const T& at(const std::pair<size_t, size_t>& index) const
	{
		if (index.first >= row_ || index.second >= col_)
			throw std::out_of_range("Index out of bounds.");

		return data_[index.first * col_ + index.second];
	}
	bool flip(const size_t i, const size_t j)
	{
		if (i >= row_ || j >= col_)
			throw std::out_of_range("Index out of bounds.");

		data_[i * col_ + j] = data_[i * col_ + j] ? 0 : 1;
		return data_[i * col_ + j];
	}

	// returns the number of rows
	size_t getRow() const	{ return row_; }
	// returns the number of columns
	size_t getCol()	const	{ return col_; }
	// returns the number of one-entrie
	size_t getOnes() const
	{
		size_t ones = 0;

		for (const auto& entry : data_)
			if (entry != 0)
				++ones;

		return ones;
	}
	// returns a formatted string containing the matrix
	std::string Print() const
	{
		std::stringstream out;
		out.precision(3);

		for (size_t i = 0; i < row_; i++)
		{
			for (size_t j = 0; j < col_; j++)
			{
				if (j != 0)
					out << " ";

				out << data_[i * col_ + j];
			}

			out << "\n";
		}

		return out.str();
	}
	// returns a random binary matrix
	static Matrix<T> random_bin_matrix(size_t r, size_t c)
	{
		// random generator from uniform distribution [0, n-1]
		std::random_device rd;
		std::mt19937 rng(rd());
		std::uniform_int_distribution<size_t> uni(0, 1);

		Matrix<T> ret(r, c);

		for (auto& it : ret.data_)
			it = uni(rng);

		return ret;
	}
private:
	std::vector<T> data_;
	size_t row_, col_;
};

/*
template<>
bool Matrix<bool>::flip(const size_t i, const size_t j)
{
	if (i >= row_ || j >= col_)
		throw std::out_of_range("Index out of bounds.");

	data_[i * col_ + j].flip();
	return data_[i * col_ + j];
}*/

#endif