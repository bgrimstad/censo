/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "eigen_utils.h"


namespace CENSO {

std::vector<double> eigenToStdVec(const DenseVector vec)
{
    std::vector<double> out(vec.size(), 0);
    for (int i = 0; i < vec.size(); i++)
        out.push_back(vec(i));
    return out;
}

DenseVector stdToEigenVec(const std::vector<double> vec)
{
    DenseVector out = DenseVector::Zero(vec.size());
    for (unsigned int i = 0; i < vec.size(); i++)
        out(i) = vec.at(i);
    return out;
}

std::vector<std::vector<double>> eigMatToStdVecVec(const DenseMatrix &mat)
{
    std::vector<std::vector<double>> vec(mat.rows());

    for(size_t i = 0; i < (size_t) mat.rows(); ++i)
    {
        for(size_t j = 0; j < (size_t) mat.cols(); ++j)
        {
            vec.at(i).push_back(mat(i, j));
        }
    }

    return vec;
}

DenseMatrix stdVecVecToEigMat(const std::vector<std::vector<double>> &vec)
{
    size_t numRows = vec.size();
    size_t numCols = numRows > 0 ? vec.at(0).size() : 0;

    DenseMatrix mat(numRows, numCols);

    for(size_t i = 0; i < numRows; ++i)
    {
        for(size_t j = 0; j < numCols; ++j)
        {
            mat(i, j) = vec.at(i).at(j);
        }
    }

    return mat;
}

} // namespace CENSO