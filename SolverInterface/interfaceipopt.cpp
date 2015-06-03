/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "interfaceipopt.h"
#include "solveripopt.h"
#include <cmath>
#include <iostream>

using std::cout;
using std::endl;

const double IPOPT_UNBOUNDED = 2e19;

InterfaceIpopt::InterfaceIpopt(CENSO::ConstraintPtr constraints)
    : constraints(constraints),
      objectiveValue(CENSO::INF)
{
    getJacobianInfo();
    getHessianInfo();
}

bool InterfaceIpopt::get_nlp_info(Index &n, Index &m, Index &nnz_jac_g, Index &nnz_h_lag, IndexStyleEnum &index_style)
{
    n = constraints->getNumVariables();
    m = constraints->getNumConstraints();
    nnz_jac_g = nnzConstraintJacobian;
    nnz_h_lag = nnzHessian;
    index_style = C_STYLE;

    return true;
}

bool InterfaceIpopt::get_bounds_info(Index n, Number *x_l, Number *x_u, Index m, Number *g_l, Number *g_u)
{
    auto vars = constraints->getVariables();

    for (int i = 0; i < n; i++)
    {
        if (-IPOPT_UNBOUNDED < vars.at(i)->getLowerBound())
        {
            x_l[i] = vars.at(i)->getLowerBound();
        }
        else
        {
            x_l[i] = -IPOPT_UNBOUNDED;
        }

        if (vars.at(i)->getUpperBound() < IPOPT_UNBOUNDED)
        {
            x_u[i] = vars.at(i)->getUpperBound();
        }
        else
        {
            x_u[i] = IPOPT_UNBOUNDED;
        }
    }

    // Constraints bounds
    std::vector<double> constraintLowerBound, constraintUpperBound;
    constraints->getConstraintBounds(constraintLowerBound, constraintUpperBound);

    for (int i = 0; i < m; i++)
    {
        if (-IPOPT_UNBOUNDED < constraintLowerBound.at(i))
        {
            g_l[i] = constraintLowerBound.at(i);
        }
        else
        {
            g_l[i] = -IPOPT_UNBOUNDED;
        }

        if (constraintUpperBound.at(i) < IPOPT_UNBOUNDED)
        {
            g_u[i] = constraintUpperBound.at(i);
        }
        else
        {
            g_u[i] = IPOPT_UNBOUNDED;
        }
    }

    return true;
}

bool InterfaceIpopt::get_starting_point(Index n, bool init_x, Number *x, bool init_z, Number *z_L, Number *z_U, Index m, bool init_lambda, Number *lambda)
{
    //    cout<<"C"<<endl;
    assert(init_x == true);
    assert(init_z == false);
    assert(init_lambda == false);

    auto vars = constraints->getVariables();

    for (int i = 0; i < n; i++)
    {
        x[i] = vars.at(i)->getValue();
    }

    return true;
}

bool InterfaceIpopt::eval_f(Index n, const Number *x, bool new_x, Number &obj_value)
{
    CENSO::DenseVector xv;
    numToVec(x,n,xv);

    double value = 0;

    auto vars = constraints->getVariables();

    for (int i = 0; i < n; i++)
    {
        value += vars.at(i)->getCost()*xv(i);
    }

    obj_value = value;

    return true;
}

bool InterfaceIpopt::eval_grad_f(Index n, const Number *x, bool new_x, Number *grad_f)
{
    CENSO::DenseVector xv,dx;
    numToVec(x,n,xv);
    //objective->evalGradient(xv,dx);

    dx = CENSO::DenseVector::Zero(n);

    auto vars = constraints->getVariables();

    for (int i = 0; i < n; i++)
    {
        dx(i) = vars.at(i)->getCost();
    }

    vecToNum(dx,n,grad_f);

    return true;
}

bool InterfaceIpopt::eval_g(Index n, const Number *x, bool new_x, Index m, Number *g)
{
    CENSO::DenseVector xv,gv;
    numToVec(x,n,xv);

    gv = constraints->eval(xv);

    assert(gv.size() == m);
    vecToNum(gv,m,g);

    return true;
}

bool InterfaceIpopt::eval_jac_g(Index n, const Number *x, bool new_x, Index m, Index nele_jac, Index *iRow, Index *jCol, Number *values)
{
    if (values == NULL)
    {
        vecToIndex(iRowGradient,nele_jac,iRow);
        vecToIndex(jColGradient,nele_jac,jCol);
    }
    else
    {
        // Return the values of the jacobian of the constraints
        CENSO::DenseVector xv,valuesv;
        numToVec(x,n,xv);
        valuesv = constraints->evalJacobian(xv);
        vecToNum(valuesv,nele_jac,values);
    }

    return true;
}

bool InterfaceIpopt::eval_h(Index n, const Number *x, bool new_x, Number obj_factor, Index m, const Number *lambda, bool new_lambda,
                            Index nele_hess, Index *iRow, Index *jCol, Number *values)
{
    if (values == NULL)
    {
        vecToIndex(iRowvHessian,nele_hess,iRow);
        vecToIndex(jColvHessian,nele_hess,jCol);
    }
    else
    {
        // Return the values of the Lagrangian Hessian matrix
        CENSO::DenseVector xv, val;
        numToVec(x,n,xv);

        val = constraints->evalHessian(xv);
        assert(val.size() == valuePositionHessian.size());

        for (int i = 0; i < nele_hess; i++)
        {
            values[i] = 0;
        }

        for (unsigned int i = 0; i < valuePositionHessian.size(); i++)
        {
            values[valuePositionHessian.at(i)] += lambda[eqnrHessian.at(i)]*val[i];
        }
    }

    return true;
}

void InterfaceIpopt::finalize_solution(SolverReturn status,
                                       Index n, const Number* x, const Number* z_L, const Number* z_U,
                                       Index m, const Number* g, const Number* lambda,
                                       Number obj_value,
                                       const IpoptData* ip_data,
                                       IpoptCalculatedQuantities* ip_cq)
{
    //    iterations = ip_data->iter_count();

    //    cout << "Ipopt return status = " << status <<endl;

    // Clear old values
    optimalSolution.clear();
    lowerBoundDuals.clear();
    upperBoundDuals.clear();
    constraintDuals.clear();

    for (Index i = 0; i < n; i++)
    {
        //double xi = 0;
        //if (abs(x[i]) > 1e-8) xi = x[i];
        optimalSolution.push_back(x[i]);
        lowerBoundDuals.push_back(z_L[i]);
        upperBoundDuals.push_back(z_U[i]);
    }

    for (Index i = 0; i < m; i++)
    {
        constraintDuals.push_back(lambda[i]);
    }

    objectiveValue = obj_value;

    bool print_solution = false;

    if (print_solution)
    {
        cout << endl << endl << "Solution of the primal variables, x" << endl;
        for (Index i=0; i<n; i++) {

            std::cout << "x[" << i << "] = " << optimalSolution[i] << endl;
        }

        cout << endl << endl << "Solution of the bound multipliers, z_L and z_U" << endl;
        for (Index i=0; i<n; i++) {
            double li = 0;
            if (std::abs(z_L[i]) > 1e-8) li = z_L[i];
            cout << "z_L[" << i << "] = " << li << "\t  ("<<z_L[i]<<")"<< endl;
        }
        for (Index i=0; i<n; i++) {
            double li = 0;
            if (std::abs(z_U[i]) > 1e-8) li = z_U[i];
            cout << "z_U[" << i << "] = " << li << "\t  ("<<z_U[i]<<")"<< endl;
        }

        cout << std::endl << std::endl << "Solution of the lagrange multipliers, lambda" << endl;
        for (Index i=0; i<m; i++) {
            double li = 0;
            if (std::abs(lambda[i]) > 1e-8) li = lambda[i];
            cout << "lambda[" << i << "] = " << li << "\t  ("<<lambda[i]<<")"<< endl;
        }

        cout << endl << endl << "Objective value" << endl;
        double fx = 0;
        if (std::abs(obj_value) > 1e-8) fx = obj_value;
        cout << "f(x*) = " << fx << "\t  ("<<obj_value<<")"<< endl;

        //    std::cout << std::endl << "Final value of the constraints:" << std::endl;
        //    for (Index i=0; i<m ;i++) {
        //        double gi = 0;
        //        if (std::abs(g[i]) > 1e-8) gi = g[i];

        //        cout << "g[" << i << "] = " << gi << "\t  ("<<g[i]<<")"<< endl;

        //        //        cout << "g(" << i << ") = " << g[i] << endl;
        //    }

    }
}

bool InterfaceIpopt::get_variables_linearity(Index n, TNLP::LinearityType *var_types)
{
    if (constraints->isHessianCalculated())
    {
        for (int i = 0; i < n; i++)
        {
            var_types[i] = TNLP::LINEAR;
        }
        std::vector<int> eqnr;
        std::vector<int> iRow;
        std::vector<int> jCol;
        constraints->structureHessian(eqnr, iRow, jCol);
        for (unsigned int i = 0; i < iRow.size(); i++)
        {
            int dx1 = iRow.at(i);
            int dx2 = jCol.at(i);
            var_types[dx1] = TNLP::NON_LINEAR;
            var_types[dx2] = TNLP::NON_LINEAR;
        }
    }
    else
    {
        for (int i = 0; i < n; i++)
        {
            var_types[i] = TNLP::NON_LINEAR;
        }
    }
    return true;
}

bool InterfaceIpopt::get_constraints_linearity(Index m, TNLP::LinearityType *const_types)
{
    std::vector<CENSO::ConstraintType> constraintTypes = constraints->getConstraintTypes();
    assert((int)constraintTypes.size() == m);
    for (int i = 0; i < m; i++)
    {
        if (constraintTypes.at(i) == CENSO::ConstraintType::LINEAR)
        {
            const_types[i] = TNLP::LINEAR;
        }
        else if (constraintTypes.at(i) == CENSO::ConstraintType::NONLINEAR_CONVEX
                 || constraintTypes.at(i) == CENSO::ConstraintType::NONLINEAR_NONCONVEX)
        {
            const_types[i] = TNLP::NON_LINEAR;
        }
        else
        {
            cout << "ERROR: Invalid constraint property." << endl;
            const_types[i] = TNLP::NON_LINEAR;
        }
    }
    return true;
}

Index InterfaceIpopt::get_number_of_nonlinear_variables()
{
    int count_nonlinear_variables = 0;
    int numVariables = constraints->getNumVariables();
    TNLP::LinearityType* ltv = new TNLP::LinearityType[numVariables];
    get_variables_linearity(numVariables,ltv);
    for (int i = 0; i < numVariables; i++)
    {
        if (ltv[i] == TNLP::NON_LINEAR)
        {
            count_nonlinear_variables++;
        }
    }
    return count_nonlinear_variables;
}

bool InterfaceIpopt::get_list_of_nonlinear_variables(Index num_nonlin_vars, Index *pos_nonlin_vars)
{
    int count_nonlinear_variables = 0;
    int numVariables = constraints->getNumVariables();
    TNLP::LinearityType* ltv = new TNLP::LinearityType[numVariables];
    get_variables_linearity(numVariables,ltv);
    for (int i = 0; i < numVariables; i++)
    {
        if (ltv[i] == TNLP::NON_LINEAR)
        {
            pos_nonlin_vars[count_nonlinear_variables++] = i;
        }
    }
    return true;
}

void InterfaceIpopt::getJacobianInfo()
{
    constraints->structureJacobian(iRowGradient, jColGradient);
    nnzConstraintJacobian = constraints->getNumNonZerosJacobian();
}

void InterfaceIpopt::getHessianInfo()
{
    std::vector<int> iv, jv;
    int count = 0;
    constraints->structureHessian(eqnrHessian, iv, jv);

    for (unsigned int i = 0; i < eqnrHessian.size(); i++)
    {

        bool existingPosition = false;

        for (unsigned int j = 0; j < iRowvHessian.size(); j++)
        {
            if (iRowvHessian.at(j) == iv.at(i) && jColvHessian.at(j) == jv.at(i))
            {
                valuePositionHessian.push_back(j);
                existingPosition = true;
            }
        }

        if (!existingPosition)
        {
            iRowvHessian.push_back(iv.at(i));
            jColvHessian.push_back(jv.at(i));
            valuePositionHessian.push_back(count);
            count++;
        }

    }

    nnzHessian = iRowvHessian.size();
}

void InterfaceIpopt::vecToNum(CENSO::DenseVector &v, int sizen, Number *n)
{
    for (int i = 0; i < sizen; i++)
    {
        n[i] = v(i);
    }
}

void InterfaceIpopt::numToVec(const Number *nr, int n, CENSO::DenseVector &v)
{
    v.resize(n);
    for (int i = 0; i < n; i++)
    {
        v(i) = nr[i];
    }
}

void InterfaceIpopt::vecToNum(std::vector<double> &v, int sizen, Number* n)
{
    for (int i = 0; i < sizen; i++)
    {
        n[i] = v.at(i);
    }
}

void InterfaceIpopt::vecToIndex(std::vector<int> &v,int sizein, Index *in)
{
    for (int i = 0; i < sizein; i++)
    {
        in[i] = v.at(i);
    }
}

void InterfaceIpopt::numToVec(const Number *nr, int n, std::vector<double> &v)
{
    for (int i = 0; i < n; i++)
    {
        v.push_back(nr[i]);
    }
}

void ipoptInterfaceTesting(CENSO::ConstraintPtr constraints, std::vector<double> lb, std::vector<double> ub)
{
    InterfaceIpopt ip(constraints);

    Index m = constraints->getNumConstraints();
    TNLP::LinearityType* lt = new TNLP::LinearityType[m];
    bool returnvalue1 = ip.get_constraints_linearity(m,lt);
    for (int i = 0; i < m; i++)
    {
        if (lt[i] == TNLP::LINEAR)
        {
            cout << "Constraint nr " << i << " is linear" << endl;
        }
        else if (lt[i] == TNLP::NON_LINEAR)
        {
            cout << "Constraint nr " << i << " is nonlinear" << endl;
        }
    }

    Index n = constraints->getNumVariables();
    TNLP::LinearityType* ltv = new TNLP::LinearityType[n];
    bool returnvalue2 = ip.get_variables_linearity(n,ltv);
    for (int i = 0; i < n; i++)
    {
        if (ltv[i] == TNLP::LINEAR)
        {
            cout << "Variable nr " << i << " is linear" << endl;
        }
        else if (ltv[i] == TNLP::NON_LINEAR)
        {
            cout << "Variable nr " << i << " is nonlinear" << endl;
        }
    }
}



