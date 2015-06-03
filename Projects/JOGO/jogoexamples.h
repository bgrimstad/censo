/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef JOGOEXAMPLES_H
#define JOGOEXAMPLES_H

#include "OptimizationProblem/constraintset.h"
#include "OptimizationProblem/constraintlinear.h"
#include "OptimizationProblem/constraintquadratic.h"

#include "datatable.h"

namespace CENSO
{

/*
 * Examples used in the JOGO paper.
 */

void saveDataTable(Splinter::DataTable &data, std::string filename);
void sampleMichalewicz();

void samplePump();
void pumpSynthesis(unsigned int grid);

void optControl1();
void optControl2();
void optControl3();

void Ptest();

double sixHumpCamelFunction(DenseVector x);

void polynomialOptimization();
void sixHumpCamelBackPoly();
void sixHumpCamelBackPoly2();
//void bilinearConstraintBound(std::vector<double> lb, std::vector<double> ub, double &clb, double &cub);
//ConstraintLinear* bilinearConstraint(std::vector<double> lb, std::vector<double> ub);
//ConstraintQuadratic* quadraticIneqConstraint();
//ConstraintLinear* lidConstraint(double lb, double ub);
//void sixHumpCamelBackPoly3();
//void sixHumpCamelBackPoly4();

void subdivision_example();

void cubicSpline();

} // namespace CENSO

#endif // JOGOEXAMPLES_H
