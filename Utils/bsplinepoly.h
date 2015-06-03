/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef BSPLINEPOLY_H
#define BSPLINEPOLY_H

#include "Utils/definitions.h"

namespace CENSO
{

DenseMatrix getBSplineBasisCoefficients(DenseVector c, DenseMatrix E, std::vector<double> lb, std::vector<double> ub);

DenseMatrix getPowerBasisCoefficients(DenseVector c, std::vector<unsigned int> degrees, std::vector<double> lb, std::vector<double> ub);

DenseMatrix getTransformationMatrix(std::vector<unsigned int> powers, std::vector<double> lb, std::vector<double> ub);

DenseMatrix getReparameterizationMatrixND(std::vector<unsigned int> degrees, std::vector<double> lb, std::vector<double> ub);

DenseMatrix getReparameterizationMatrix1D(int p, double a, double b);

DenseMatrix getBSplineToPowerBasisMatrix1D(unsigned int p);

double binomialCoeff(int n, int k);

std::vector< std::vector<double> > getRegularKnotVectors(std::vector<unsigned int> deg, std::vector<double> lb, std::vector<double> ub);

std::vector<std::vector<double> > getUniqueKnots(std::vector<std::vector<double> > knots);

DenseMatrix getPowersMatrix(std::vector<unsigned int> degrees);

DenseMatrix getPermutationMatrix(std::vector<unsigned int> num);

} // namespace CENSO

#endif // BSPLINEPOLY_H
