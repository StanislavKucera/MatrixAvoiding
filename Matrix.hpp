#ifndef Matrix_hpp_
#define Matrix_hpp_

#include <vector>
#include <algorithm>
#include <initializer_list>
#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>

template< typename T>
class matrix
{
public:
	// constructs a matrix of given size with default value in each element
	matrix(size_t r, size_t c) : row_(r), col_(c), data_(r*c) {}
	// constructs a matrix of given size with given value in each element
	matrix(size_t r, size_t c, const T& def) : row_(r), col_(c), data_(r*c, def) {}
	// constructs a matrix from initializer list of initializer lists
	matrix(const std::initializer_list<std::initializer_list<T> >& il)
		: row_(il.size()), col_(il.begin()->size())
	{
		data_.reserve(row_ * col_);
		for (auto it = il.begin(); it != il.end(); it++)
		{
			for (auto it2 = it->begin(); it2 != it->end(); it2++)
			{
				data_.push_back(*it2);
			}
		}
	}
	// constructs a matrix from input file, expected format is: number of rows, number of columns, rows*cols elements of type T
	explicit matrix(const std::string& input)
	{
		std::ifstream iFile(input);
		iFile >> row_ >> col_;
		data_.resize(row_ * col_);
		for (auto &it : data_)
			iFile >> it;
		iFile.close();
	}

	matrix(const matrix<T>& m) : row_(m.row_), col_(m.col_), data_(m.data_.begin(), m.data_.end()) {}
	matrix(matrix<T>&& m)
		: row_(m.row_), col_(m.col_)
	{
		data_ = std::move(m.data_);
	}
	matrix& operator=(const matrix<T>& m)
	{
		row_ = m.row_;
		col_ = m.col_;
		data_.resize(row_ * col_);
		std::copy(m.data_.begin(), m.data_.end(), data_.begin());
		return *this;
	}
	matrix& operator=(matrix&& m)
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
	// returns number of rows
	size_t getRow() const	{ return row_; }
	// returns number of columns
	size_t getCol()	const	{ return col_; }
	// returns string containing formatted matrix
	std::string Print()
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
private:
	size_t row_, col_;
	std::vector<T> data_;
};

#endif