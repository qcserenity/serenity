/**
 * @file Eigen3HDF5_test.cpp
 *
 * @date 27 Jul 2019
 * @author Moritz Bensberg
 *
 * Taken and modified from:
 * https://github.com/garrison/eigen3-hdf5/commits/master
 * Tag: 2c782414251e75a2de9b0441c349f5f18fe929a2
 *
 * @copyright
 *  SPECIAL LICENSE FOR THIS FILE ONLY
 *
 *  The MIT License (MIT)
 *
 *  Copyright (c) 2013 James R. Garrison
 *
 *  Permission is hereby granted, free of charge, to any person obtaining a copy
 *  of this software and associated documentation files (the "Software"), to deal
 *  in the Software without restriction, including without limitation the rights
 *  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 *  copies of the Software, and to permit persons to whom the Software is
 *  furnished to do so, subject to the following conditions:
 *
 *  The above copyright notice and this permission notice shall be included in
 *  all copies or substantial portions of the Software.
 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 *  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 *  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 *  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 *  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 *  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 *  THE SOFTWARE.
 */

/* Include Serenity Internal Headers */
#include "io/Eigen3HDF5.h"
/* Include Std and External Headers */
#include <H5Cpp.h>
#include <gtest/gtest.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <complex>
#include <iostream>

namespace Eigen {
// C++ and/or gtest require that these two methods, which are used by calling
// ASSERT_PRED_FORMAT2, be in the namespace of its argument.

// utility function to print an eigen object to an ostream; gtest will use this when
// it outputs a matrix used in a failed assertion. Without this function, gtest seems
// to dump some kind of byte representation of an eigen matrix, which is not very
// helpful.
template<class Derived>
void PrintTo(const Eigen::EigenBase<Derived>& mat, ::std::ostream* os) {
  (*os) << mat.derived() << "\n";
}

// utility function for gtest to use to check if two eigen objects are identical.
// returns assertion success when they are identical; returns assertion failure along
// with a nicely formatted message with the matrix contents when they are not
// identical.
//
// I put this function in this matrix test cpp file for a few reasons: 1) there is
// not already a header file to put common test code for eigen3-hdf5, and 2) because
// I needed it to help me debug test failures as I implemented the no copy read and
// write functions. I really think that providing a header for common test code
// should be addressed at some point, and then this function (and its companion
// PrintTo) should be moved there.
//
// Usage:
//
// ASSERT_PRED_FORMAT2(assert_same, mat, mat2);
template<class DerivedExp, class DerivedAct>
::testing::AssertionResult assert_same(const char* exp_expr, const char* act_expr,
                                       const Eigen::EigenBase<DerivedExp>& exp, const Eigen::EigenBase<DerivedAct>& act) {
  if (exp.rows() == act.rows() && exp.cols() == act.cols() && exp.derived() == act.derived()) {
    return ::testing::AssertionSuccess();
  }

  // if eigen did not define the == operator, you could use
  // exp.derived().cwiseEqual(act.derived()).all();

  ::testing::AssertionResult result = ::testing::AssertionFailure()
                                      << "Eigen objects are not the same: (" << exp_expr << ", " << act_expr << ")\n"
                                      << exp_expr << ":\n"
                                      << ::testing::PrintToString(exp) << "\n---and\n"
                                      << act_expr << ":\n"
                                      << ::testing::PrintToString(act) << "\n---are not equal!\n";

  return result;
}
} // namespace Eigen

namespace Serenity {
TEST(Eigen3HDF5_SparseMatrix, Double) {
  Eigen::SparseMatrix<double> mat(3, 3), mat2;
  mat.insert(0, 1) = 2.7;
  mat.insert(2, 0) = 82;
  {
    H5::H5File file("test_SparseMatrix_Double.h5", H5F_ACC_TRUNC);
    EigenHDF5::save_sparse(file, "mat", mat);
  }
  {
    H5::H5File file("test_SparseMatrix_Double.h5", H5F_ACC_RDONLY);
    EigenHDF5::load_sparse(file, "mat", mat2);
  }
#ifdef LOGGING
  std::cout << mat2 << std::endl;
#endif
  ASSERT_EQ(Eigen::MatrixXd(mat), Eigen::MatrixXd(mat2));
  std::remove("test_SparseMatrix_Double.h5");
}

TEST(Eigen3HDF5_SparseMatrix, Complex) {
  Eigen::SparseMatrix<std::complex<double>> mat(4, 4), mat2;
  mat.insert(0, 1) = std::complex<double>(2, 4.5);
  mat.insert(1, 2) = std::complex<double>(82, 1);
  {
    H5::H5File file("test_SparseMatrix_Complex.h5", H5F_ACC_TRUNC);
    EigenHDF5::save_sparse(file, "mat", mat);
  }
  {
    H5::H5File file("test_SparseMatrix_Complex.h5", H5F_ACC_RDONLY);
    EigenHDF5::load_sparse(file, "mat", mat2);
  }
#ifdef LOGGING
  std::cout << mat2 << std::endl;
#endif
  ASSERT_EQ(Eigen::MatrixXcd(mat), Eigen::MatrixXcd(mat2));
  std::remove("test_SparseMatrix_Complex.h5");
}

TEST(Eigen3HDF5_MatrixRoundTrip, Double) {
  Eigen::MatrixXd mat(3, 4), mat2;
  mat << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12;
#ifdef LOGGING
  std::cout << mat << std::endl;
#endif
  {
    H5::H5File file("test_MatrixRoundTrip_Double.h5", H5F_ACC_TRUNC);
    EigenHDF5::save(file, "double_matrix", mat);
  }
  {
    H5::H5File file("test_MatrixRoundTrip_Double.h5", H5F_ACC_RDONLY);
    EigenHDF5::load(file, "double_matrix", mat2);
  }
  ASSERT_PRED_FORMAT2(assert_same, mat, mat2);
  std::remove("test_MatrixRoundTrip_Double.h5");
}

TEST(Eigen3HDF5_MatrixRoundTrip, LongDouble) {
  Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic> mat(3, 4), mat2;
  mat << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12;
#ifdef LOGGING
  std::cout << mat << std::endl;
#endif
  {
    H5::H5File file("test_MatrixRoundTrip_LongDouble.h5", H5F_ACC_TRUNC);
    EigenHDF5::save(file, "longdouble_matrix", mat);
  }
  {
    H5::H5File file("test_MatrixRoundTrip_LongDouble.h5", H5F_ACC_RDONLY);
    EigenHDF5::load(file, "longdouble_matrix", mat2);
  }
  ASSERT_EQ(mat, mat2);
  std::remove("test_MatrixRoundTrip_LongDouble.h5");
}

TEST(Eigen3HDF5_MatrixRoundTrip, Int) {
  Eigen::MatrixXi mat(3, 4), mat2;
  mat << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12;
#ifdef LOGGING
  std::cout << mat << std::endl;
#endif
  {
    H5::H5File file("test_MatrixRoundTrip_Int.h5", H5F_ACC_TRUNC);
    EigenHDF5::save(file, "int_matrix", mat);
  }
  {
    H5::H5File file("test_MatrixRoundTrip_Int.h5", H5F_ACC_RDONLY);
    EigenHDF5::load(file, "int_matrix", mat2);
  }
  ASSERT_EQ(mat, mat2);
  std::remove("test_MatrixRoundTrip_Int.h5");
}

TEST(Eigen3HDF5_MatrixRoundTrip, ULongLong) {
  Eigen::Matrix<unsigned long long, Eigen::Dynamic, Eigen::Dynamic> mat(3, 4), mat2;
  mat << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12;
#ifdef LOGGING
  std::cout << mat << std::endl;
#endif
  {
    H5::H5File file("test_MatrixRoundTrip_ULongLong.h5", H5F_ACC_TRUNC);
    EigenHDF5::save(file, "ull_matrix", mat);
  }
  {
    H5::H5File file("test_MatrixRoundTrip_ULongLong.h5", H5F_ACC_RDONLY);
    EigenHDF5::load(file, "ull_matrix", mat2);
  }
  ASSERT_EQ(mat, mat2);
  std::remove("test_MatrixRoundTrip_ULongLong.h5");
}

TEST(Eigen3HDF5_MatrixRoundTrip, ComplexDouble) {
  Eigen::MatrixXcd mat(3, 4), mat2;
  mat << 1, std::complex<double>(0, 2), 3, 4, 5, 6, 7, 8, 9, 10, 11, 12;
#ifdef LOGGING
  std::cout << mat << std::endl;
#endif
  {
    H5::H5File file("test_MatrixRoundTrip_ComplexDouble.h5", H5F_ACC_TRUNC);
    EigenHDF5::save(file, "complex_matrix", mat);
  }
  {
    H5::H5File file("test_MatrixRoundTrip_ComplexDouble.h5", H5F_ACC_RDONLY);
    EigenHDF5::load(file, "complex_matrix", mat2);
  }
  ASSERT_EQ(mat, mat2);
  std::remove("test_MatrixRoundTrip_ComplexDouble.h5");
}

TEST(Eigen3HDF5_MatrixRoundTrip, IntBlock) {
  Eigen::Matrix4i mat(Eigen::Matrix4i::Zero());
  Eigen::Matrix4i mat2(Eigen::Matrix4i::Zero());
  mat(0, 0) = 1;
  mat(0, 1) = 2;
  mat(1, 0) = 3;
  mat(1, 1) = 4;
  mat(2, 2) = 5;
  mat2(2, 2) = 5;
#ifdef LOGGING
  std::cout << mat << std::endl;
#endif
  {
    H5::H5File file("test_MatrixRoundTrip_IntBlock.h5", H5F_ACC_TRUNC);
    EigenHDF5::save(file, "int_block", mat.block(0, 0, 2, 2));
  }
  {
    H5::H5File file("test_MatrixRoundTrip_IntBlock.h5", H5F_ACC_RDONLY);
    EigenHDF5::load(file, "int_block", mat2.block(0, 0, 2, 2));
  }
  ASSERT_PRED_FORMAT2(assert_same, mat, mat2);
  std::remove("test_MatrixRoundTrip_IntBlock.h5");
}

TEST(Eigen3HDF5_MatrixRoundTrip, IntBlockRowMajor) {
  typedef Eigen::Matrix<int, 4, 4, Eigen::RowMajor> Matrix4RowMajor;
  Matrix4RowMajor mat(Eigen::Matrix4i::Zero());
  Matrix4RowMajor mat2(Eigen::Matrix4i::Zero());
  mat << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16;

  mat2(0, 2) = 3;
  mat2(0, 3) = 4;
  mat2(1, 2) = 7;
  mat2(1, 3) = 8;
  mat2(2, 0) = 9;
  mat2(2, 1) = 10;
  mat2(2, 2) = 11;
  mat2(2, 3) = 12;
  mat2(3, 0) = 13;
  mat2(3, 1) = 14;
  mat2(3, 2) = 15;
  mat2(3, 3) = 16;
#ifdef LOGGING
  std::cout << mat << std::endl;
#endif
  {
    H5::H5File file("test_MatrixRoundTrip_IntBlockRowMajor.h5", H5F_ACC_TRUNC);
    EigenHDF5::save(file, "int_block", mat.block(0, 0, 2, 2));
  }
  {
    H5::H5File file("test_MatrixRoundTrip_IntBlockRowMajor.h5", H5F_ACC_RDONLY);
    EigenHDF5::load(file, "int_block", mat2.block(0, 0, 2, 2));
  }
  ASSERT_PRED_FORMAT2(assert_same, mat, mat2);
  std::remove("test_MatrixRoundTrip_IntBlockRowMajor.h5");
}

TEST(Eigen3HDF5_MatrixRoundTrip, DoubleSkipInternalCopy) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> mat(3, 4); // , mat2;
  Eigen::MatrixXd mat2;
  mat << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12;
#ifdef LOGGING
  std::cout << mat << std::endl;
#endif
  {
    H5::H5File file("test_MatrixRoundTrip_DoubleSkipInternalCopy.h5", H5F_ACC_TRUNC);
    EigenHDF5::save(file, "double_matrix", mat);
  }
  {
    H5::H5File file("test_MatrixRoundTrip_DoubleSkipInternalCopy.h5", H5F_ACC_RDONLY);
    EigenHDF5::load(file, "double_matrix", mat2);
  }
  ASSERT_PRED_FORMAT2(assert_same, mat, mat2);
  std::remove("test_MatrixRoundTrip_DoubleSkipInternalCopy.h5");
}

TEST(Eigen3HDF5_MatrixRoundTrip, DoubleSkipInternalCopyBlock) {
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MyMatrixXdRowMajor;
  MyMatrixXdRowMajor mat(3, 4);
  Eigen::Block<MyMatrixXdRowMajor> matblock = mat.block(1, 1, 2, 3);
  Eigen::MatrixXd mat2;
  mat << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12;
#ifdef LOGGING
  std::cout << mat << std::endl;
#endif
  {
    H5::H5File file("test_MatrixRoundTrip_DoubleSkipInternalCopy.h5", H5F_ACC_TRUNC);
    EigenHDF5::save(file, "double_matrix", matblock);
  }
  {
    H5::H5File file("test_MatrixRoundTrip_DoubleSkipInternalCopy.h5", H5F_ACC_RDONLY);
    EigenHDF5::load(file, "double_matrix", mat2);
  }
  ASSERT_PRED_FORMAT2(assert_same, matblock, mat2);
  std::remove("test_MatrixRoundTrip_DoubleSkipInternalCopy.h5");
}

TEST(Eigen3HDF5_MatrixRoundTrip, DoubleFixedRow) {
  typedef Eigen::Matrix<double, 4, 6, Eigen::RowMajor> MyMatrix46RowMajor;
  MyMatrix46RowMajor mat;
  Eigen::Block<MyMatrix46RowMajor> matblock = mat.block(1, 2, 2, 3);
  MyMatrix46RowMajor fmat2;
  Eigen::Matrix<double, 2, 3, Eigen::RowMajor> fmatblock2;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> dmat2, dmatblock2;

  mat << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24;
#ifdef LOGGING
  std::cout << mat << std::endl;
#endif
  {
    H5::H5File file("test_MatrixRoundTrip_DoubleFixedRow.h5", H5F_ACC_TRUNC);
    EigenHDF5::save(file, "double_matrix", mat);
    EigenHDF5::save(file, "matrix_block", matblock);
  }
  {
    H5::H5File file("test_MatrixRoundTrip_DoubleFixedRow.h5", H5F_ACC_RDONLY);
    // read into a dynamic sized matrix and then copy into fixed size
    EigenHDF5::load(file, "double_matrix", dmat2);
    EigenHDF5::load(file, "matrix_block", dmatblock2);
    fmat2 = dmat2;
    fmatblock2 = dmatblock2;
  }
  ASSERT_PRED_FORMAT2(assert_same, mat, fmat2);
  ASSERT_PRED_FORMAT2(assert_same, matblock, fmatblock2);
  std::remove("test_MatrixRoundTrip_DoubleFixedRow.h5");
}

TEST(Eigen3HDF5_MatrixRoundTrip, DoubleFixedCol) {
  typedef Eigen::Matrix<double, 4, 6, Eigen::ColMajor> MyMatrix46ColMajor;
  MyMatrix46ColMajor mat;
  Eigen::Block<MyMatrix46ColMajor> matblock = mat.block(1, 2, 2, 3);
  MyMatrix46ColMajor fmat2;
  Eigen::Matrix<double, 2, 3, Eigen::RowMajor> fmatblock2;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> dmat2, dmatblock2;

  mat << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24;
#ifdef LOGGING
  std::cout << mat << std::endl;
#endif
  {
    H5::H5File file("test_MatrixRoundTrip_DoubleFixedRow.h5", H5F_ACC_TRUNC);
    EigenHDF5::save(file, "double_matrix", mat);
    EigenHDF5::save(file, "matrix_block", matblock);
  }
  {
    H5::H5File file("test_MatrixRoundTrip_DoubleFixedRow.h5", H5F_ACC_RDONLY);
    EigenHDF5::load(file, "double_matrix", fmat2);
    EigenHDF5::load(file, "matrix_block", fmatblock2);
  }
  ASSERT_PRED_FORMAT2(assert_same, mat, fmat2);
  ASSERT_PRED_FORMAT2(assert_same, matblock, fmatblock2);
  std::remove("test_MatrixRoundTrip_DoubleFixedRow.h5");
}

TEST(Eigen3HDF5_Attribute, Matrix) {
  Eigen::Matrix<double, 2, 3, Eigen::RowMajor> rmat1, rmat2;
  Eigen::Matrix<double, 2, 3, Eigen::ColMajor> cmat1, cmat2;
  rmat1 << 1, 2, 3, 4, 5, 6;
  cmat1 << 1, 2, 3, 4, 5, 6;
  {
    H5::H5File file("test_Attribute_Matrix.h5", H5F_ACC_TRUNC);
    EigenHDF5::save_attribute(file, "rowmat", rmat1);
    EigenHDF5::save_attribute(file, "colmat", cmat1);
  }
  {
    H5::H5File file("test_Attribute_Matrix.h5", H5F_ACC_RDONLY);
    EigenHDF5::load_attribute(file, "rowmat", rmat2);
    EigenHDF5::load_attribute(file, "colmat", cmat2);
  }
  ASSERT_PRED_FORMAT2(assert_same, rmat1, rmat2);
  ASSERT_PRED_FORMAT2(assert_same, cmat1, cmat2);
  ASSERT_PRED_FORMAT2(assert_same, rmat2, cmat2);
  std::remove("test_Attribute_Matrix.h5");
}

TEST(Eigen3HDF5_Attribute, Integer) {
  H5::H5File file("test_Attribute_Integer.h5", H5F_ACC_TRUNC);
  EigenHDF5::save_scalar_attribute(file, "integer", 23);
  std::remove("test_Attribute_Integer.h5");
}

TEST(Eigen3HDF5_Attribute, Double) {
  H5::H5File file("test_Attribute_Double.h5", H5F_ACC_TRUNC);
  EigenHDF5::save_scalar_attribute(file, "double", 23.7);
  std::remove("test_Attribute_Double.h5");
}

TEST(Eigen3HDF5_Attribute, String) {
  H5::H5File file("test_Attribute_String.h5", H5F_ACC_TRUNC);
  EigenHDF5::save_scalar_attribute(file, "str1", std::string("hello"));
  EigenHDF5::save_scalar_attribute(file, "str2", "goodbye");
  const char* s = "again";
  EigenHDF5::save_scalar_attribute(file, "str3", s);
  std::remove("test_Attribute_String.h5");
}

TEST(Eigen3HDF5_VectorRoundTrip, Double) {
  Eigen::Vector4d mat, mat2;
  mat << 1, 2, 3, 4;
#ifdef LOGGING
  std::cout << mat << std::endl;
#endif
  {
    H5::H5File file("test_VectorRoundTrip_Double.h5", H5F_ACC_TRUNC);
    EigenHDF5::save(file, "double_vector", mat);
  }
  {
    H5::H5File file("test_VectorRoundTrip_Double.h5", H5F_ACC_RDONLY);
    EigenHDF5::load(file, "double_vector", mat2);
  }
  ASSERT_EQ(mat, mat2);
  std::remove("test_VectorRoundTrip_Double.h5");
}

TEST(Eigen3HDF5_VectorRoundTrip, Int) {
  Eigen::VectorXi mat(12), mat2;
  mat << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12;
#ifdef LOGGING
  std::cout << mat << std::endl;
#endif
  {
    H5::H5File file("test_VectorRoundTrip_Int.h5", H5F_ACC_TRUNC);
    EigenHDF5::save(file, "int_vector", mat);
  }
  {
    H5::H5File file("test_VectorRoundTrip_Int.h5", H5F_ACC_RDONLY);
    EigenHDF5::load(file, "int_vector", mat2);
  }
  ASSERT_EQ(mat, mat2);
  std::remove("test_VectorRoundTrip_Int.h5");
}

} /* namespace Serenity */
