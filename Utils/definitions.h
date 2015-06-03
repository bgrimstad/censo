/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef GENERALDEFINITIONS_H
#define GENERALDEFINITIONS_H

#include <iostream>
#include <iomanip>

#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace CENSO
{

// Eigen vectors
typedef Eigen::VectorXd DenseVector;
typedef Eigen::SparseVector<double> SparseVector;

// Eigen matrices
typedef Eigen::MatrixXd DenseMatrix;
typedef Eigen::SparseMatrix<double> SparseMatrix; // declares a column-major sparse matrix type of double

// Infinity constant
const double INF = std::numeric_limits<double>::infinity();

// Print functions
template<typename T>
std::ostream& operator<<(std::ostream &os, const std::vector<T> &vec)
{
    for (auto val : vec)
        os << std::setprecision(6) << std::setw(8) << std::left << val << " ";
    os << std::endl;
    return os;
}

template<typename T>
std::ostream& operator<<(std::ostream &os, const std::vector< std::vector<T> > &vec)
{
    for (auto val : vec)
        os << val;
    return os;
}

// Compare two numbers
template<typename T>
bool assertNear(T x, T y, double tolAbs = 1e-8, double tolRel = 1e-8)
{
    double dx = std::abs(x - y);
    double xAbs = 0.5*(std::abs(x) + std::abs(y));
    double err = std::max(tolAbs, tolRel*xAbs);

    /*
     * NOTE: Must be strictly less than (<)
     * when comparing infinite numbers
     */
    return dx < err;
}

// Integer test
bool isInteger(double value);

// Randomizer
static bool randomSeedSet = false;
int randomInteger(int min, int max);

// Linspace
std::vector<double> linspace(double start, double stop, unsigned int points);

// Converters
std::vector<double> eigenToStdVec(const DenseVector vec);
DenseVector stdToEigenVec(const std::vector<double> vec);

} // namespace CENSO

#endif // GENERALDEFINITIONS_H
