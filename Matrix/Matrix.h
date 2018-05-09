#pragma once
#include "stdafx.h"

namespace mango
{
	// the axis or direction of the matrix
	enum class Axis
	{
		Row,
		Col,
		Page
	};

	class Matrix
	{
	private:
		unsigned	rows_ = 0;
		unsigned	cols_ = 0;
		unsigned	pages_ = 0;

		float*		data_ = nullptr;


	public:

		Matrix(unsigned rows, unsigned cols, unsigned pages = 1);
		Matrix(const Matrix& m);
		Matrix(Matrix&& m);
		~Matrix();

		inline unsigned Rows()  const { return rows_; };
		inline unsigned Cols()  const { return cols_; };
		inline unsigned Pages() const { return pages_; };

		inline float& operator() (unsigned row, unsigned col, unsigned page = 0)
		{
			if (row >= rows_ || col >= cols_ || page >= pages_)
				throw std::out_of_range("Index of range");
			return data_[page*rows_*cols_ + row * cols_ + col];
		}

		inline float operator() (unsigned row, unsigned col, unsigned page = 0) const
		{
			if (row >= rows_ || col >= cols_ || page >= pages_)
				throw std::out_of_range("Index of range");
			return data_[page*rows_*cols_ + row * cols_ + col];
		}

		Matrix& operator= (const Matrix& m);
		Matrix& operator= (Matrix&& m);
		Matrix& operator+= (const float val);
		Matrix& operator-= (const float val);
		Matrix& operator*= (const float val);
		Matrix& operator/= (const float val);

		//----------------------------------------------------------------------------

		// Reshape the matrix, does not delete or add elements.
		Matrix& Reshape(unsigned rows, unsigned cols, unsigned pages = 1);

		// Sum the matrix along the specified axis, the length along this axis will become one.
		Matrix& Sum(Axis axis);

		// Sum the matrix along the specified axis (range from start to end), the length along this axis will become one.
		Matrix& Sum(Axis axis, unsigned start, unsigned end);

		// Take the average along the specified axis, the length along this axis will become one.
		Matrix& Average(Axis axis);

		// Take the average along the specified axis (range from start to end), the length along this axis will become one.
		Matrix& Average(Axis axis, unsigned start, unsigned end);

		// Rebin the matrix laong the specified axis, the length along this axis will reduce by a factor of rebin size.
		Matrix& Rebin(unsigned rebinSize, Axis axis, bool ignoreTail=false);

		// For any element less than threshold in the matrix, make the value equal to threshold.
		Matrix& AllNoLessThan(float threshold);

		// Take the logarithm of each element.
		Matrix& Log();

		// Save to raw file.
		bool SaveRawFile(const char* filename);

		// Append to raw file.
		bool AppendRawFile(const char* filename);


		//////////////////////////////////////////////////////////////////////////////////////
		// static functions
		//////////////////////////////////////////////////////////////////////////////////////

		// Take the logarithm of each element.
		static Matrix Log(const Matrix& m);

		// Sum the matrix along the specified axis, and return the result new matrix.
		static Matrix Sum(const Matrix& m, Axis axis);

		// Sum the matrix along the specified axis (range from start to end), and return the result new matrix.
		static Matrix Sum(const Matrix& m, Axis axis, unsigned start, unsigned end);

		// Take the average along the specified axis, and return the result new matrix.
		static Matrix Average(const Matrix& m, Axis axis);

		// Take the average along the specified axis (range from start to end), and return the result new matrix.
		static Matrix Average(const Matrix& m, Axis axis, unsigned start, unsigned end);

		// Read data from raw file.
		static Matrix ReadRawFile(const char* filename, const unsigned& rows, const unsigned& cols, const unsigned& pages = 1);

	};


	Matrix operator+(const Matrix& m1, const Matrix& m2);
	Matrix operator-(const Matrix& m1, const Matrix& m2);
	Matrix operator*(const Matrix& m1, const Matrix& m2);
	Matrix operator/(const Matrix& m1, const Matrix& m2);

	Matrix operator+(const Matrix& m, const float& val);
	Matrix operator-(const Matrix& m, const float& val);
	Matrix operator*(const Matrix& m, const float& val);
	Matrix operator/(const Matrix& m, const float& val);

	Matrix operator+(const float& val, const Matrix& m);
	Matrix operator-(const float& val, const Matrix& m);
	Matrix operator*(const float& val, const Matrix& m);
	Matrix operator/(const float& val, const Matrix& m);

}


