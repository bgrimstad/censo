/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef PWLMODELS_H
#define PWLMODELS_H

#include "Utils/definitions.h"
#include "datatable.h"
#include "OptimizationProblem/constraint.h"

using SPLINTER::DataTable;
using namespace CENSO;

/*
 * A piecewise linear model
 * has a linear constraint set
 * in the variables
 * x in R(dimension),
 * z in R(1),
 * lambda in R(numAuxiliary),
 * y in Z(numBinary).
 */
struct PiecewiseLinearModel
{
    PiecewiseLinearModel(ConstraintPtr constraints, unsigned int dimension, unsigned int numAuxiliary, unsigned int numBinary)
        : constraints(constraints),
          dimension(dimension),
          numAuxiliary(numAuxiliary),
          numBinary(numBinary)
    {
        assert(constraints->getNumVariables() == dimension+1+numAuxiliary+numBinary);
    }

    ConstraintPtr constraints;
    unsigned int dimension, numAuxiliary, numBinary;
};

constexpr double pi();
double sinewave1D(double x);
double sinewave2D(double x);
DataTable sampleFunction1D();
DataTable sampleFunction2D();
DenseMatrix makePointMatrix(const DataTable &points);

void run_mip_models();

// Monovariate problems
unsigned int P00a_MIP(DataTable samples);
unsigned int P00a_global(DataTable samples);

// Bivariate problems
unsigned int P00b_MIP(DataTable samples);
unsigned int P00b_global(DataTable samples);

// Constrained problem
unsigned int P01_MIP(int numSamples);
unsigned int P01_global(int numSamples);

/*
 * Multiple Choice (MC) model
 * for the multivariate case
 */
template <unsigned int N>
PiecewiseLinearModel constraintPiecewiseLinearMC(Variables vars, const DataTable &samples);

/*
 * Convex Combination (CC) model
 * for the multivariate case
 */
template <unsigned int N>
PiecewiseLinearModel constraintPiecewiseLinearCC(Variables vars, const DataTable &samples);

/*
 * Disaggregated Convex Combination (DCC) model
 * for the multivariate case
 */
template <unsigned int N>
PiecewiseLinearModel constraintPiecewiseLinearDCC(Variables vars, const DataTable &samples);

/*
 * Logarithmic Disaggregated Convex Combination (DLog) model
 * for the multivariate case
 */
template <unsigned int N>
PiecewiseLinearModel constraintPiecewiseLinearDLog(Variables vars, const DataTable &samples);

#endif // PWLMODELS_H
