// Authors: Joe Pitt-Francis, Nicholas Barton

#include <iostream>
#include "Vector.hpp"

using namespace std;


// constructor that creates vector of given size with
// double precision entries all initially set to zero
Vector::Vector(int sizeVal)
{
  mData=new double[sizeVal];
  mSize = sizeVal;
  for (int i = 0; i < mSize; i++)
  {
     mData[i] = 0.0;
  }
}


// copy constructor - creates vector with the same entries as v1
Vector::Vector(const Vector& v1)
{
  mSize = v1.mSize;
  mData = new double[mSize];
  for (int i = 0; i < v1.mSize; i++)
  {
    mData[i] = v1.mData[i];
  }
}


// destructor - deletes pointer
Vector::~Vector()
{
  delete[] mData;
}


// definition of + between two vectors
Vector operator+(const Vector& v1, 
                        const Vector& v2)
{
  // dimension checking
  if (v1.mSize != v2.mSize)
  {
    throw Exception("dimension mismatch",
      "Attempt to add vectors of different length.");
  }

  Vector w(v1.mSize);

  for (int i = 0; i < v1.mSize; i++)
  {
    w.mData[i] = v1.mData[i] + v2.mData[i];
  }

  return w;
}


// definition of - between two vectors
// take unary minus of second argument and add
Vector operator-(const Vector& v1, const Vector& v2)
{
  return v1 + (-v2);
}


// definition of scalar product between two vectors
double operator*(const Vector& v1, const Vector& v2)
{
  // dimension checking
  if (v1.mSize != v2.mSize)
  {
    throw Exception("dimension mismatch",
      "Attempt to add vectors of different length.");
  }

  //  compute scalar product
  double dp = 0.0;

  for (int i = 0; i < v1.mSize; i++)
    {
      dp += v1.mData[i] * v2.mData[i];
    }

  return dp;

}


// definition of multiplication between a vector and a scalar
Vector operator*(const Vector& v, const double& a)
{
  //  create a vector of the same length as v with entries equal to a*v
  Vector w(v.mSize);

  for (int i = 0; i < v.mSize; i++)
    {
      w.mData[i] = a * v.mData[i];
    }

  return w;
}


// definition of multiplication between a scalar and a vector
Vector operator*(const double& a, const Vector& v)
{
  // don't duplicate code
  return v * a;
}


// definition of division of a vector by a scalar
Vector operator/(const Vector& v, const double& a)
{
  // parameter checks
  if (a == 0.0)
  {
     throw Exception("div 0", "Attempt to divide by zero.");
  }

  double arecip = 1.0 / a;

  // don't duplicate code
  return v * arecip;
}


// definition of the unary operator -
Vector operator-(const Vector& v)
{
  Vector w(v.mSize);

  for (int i = 0; i < v.mSize; i++)
    {
      w.mData[i] = -v.mData[i];
    }

  return w;
}


// definition of vector operator =
Vector& Vector::operator=(const Vector& v)
{
  // dimension checking
  if (v.mSize != mSize)
    {
      throw Exception("dimension mismatch",
		  "Attempt to assign to vector of incompatible length.");
    }
  for (int i = 0; i < mSize; i++)
	{
	  mData[i] = v.mData[i];
  }

  return *this;
}


// definition of vector operator ()
// allows v.mData[i] to be written as v(i+1), as in Matlab and FORTRAN
double& Vector::operator()(int i)
{
  if (i < 1)
    {
      throw Exception("out of range",
		  "accessing vector through () - index too small");
    }
  else if (i > mSize)
    {
      throw Exception("length mismatch",
		  "accessing vector through () - index too high");
    }

  return mData[i - 1];
}


// to print the vector in a consistent format
ostream& operator<<(ostream& output, const Vector& v)
{
  output << "(";
  for (int i = 0; i < v.mSize; i++)
    { 
      output << v.mData[i];
      if (i != v.mSize - 1)
      {
	       output << ", ";
      }
      else
      {
	       output << ")";
      }
    }

  return output;  // for multiple << operators.
}


// to write a vector as a string -- good for error checking!
string to_string(const Vector& v)
{
    string output = "(";

    for (int i = 0; i < v.mSize; i++)
     { 
        output += to_string(v.mData[i]);
        if (i != v.mSize - 1)
        {
          output += ", ";
        }
        else
        {
          output += ")";
        }
      }

  return output;
}


// slice a vector like in MATLAB
Vector slice(const Vector& v, int start, int end)
{
    // return a Vector of the implicitly specified length
    Vector output(end - start + 1);

    for (int i = 0; i < end - start + 1; i++)
    {
      output.mData[i] = v.mData[start + i - 1];
    }

    return output;
}

// Friend function
double norm(Vector& v, int p)
{
  // don't duplicate code
  return v.norm(p);
}


// Member method
// calculate p-norm of a vector v
// default value for p is 2
double Vector::norm(int p) const
{
  double temp, norm_val;

  norm_val = 0.0;
  for (int i = 0; i < mSize; i++)
    {
      temp = fabs(mData[i]);
      norm_val += pow(temp, p);
    }

  return pow(norm_val, 1.0 / ((double) (p)));
}


// return length of a vector
int length(const Vector& v)
{
  return v.mSize;
}
