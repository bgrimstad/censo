/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "interfacebonmin.h"
#include "iostream"
using std::cout;
using std::endl;

InterfaceBonmin::InterfaceBonmin(CENSO::ConstraintPtr constraints)
{
    ipoptInterface = new InterfaceIpopt(constraints);
}

/*
 * NOTE: Bonmin has also defined VariableType
 */
bool InterfaceBonmin::get_variables_types(Index n, VariableType* var_types)
{
    auto vars = ipoptInterface->getConstraints()->getVariables();

    for (int i = 0; i < n; i++)
    {
        if (vars.at(i)->getType() == CENSO::VariableType::BINARY)
        {
            var_types[i] = Bonmin::TMINLP::BINARY;
        }
        else if (vars.at(i)->getType() == CENSO::VariableType::INTEGER)
        {
            var_types[i] = Bonmin::TMINLP::INTEGER;
        }
        else if (vars.at(i)->getType() == CENSO::VariableType::CONTINUOUS)
        {
            var_types[i] = Bonmin::TMINLP::CONTINUOUS;
        }
        else
        {
            cout << "Unknown variable type!" << endl;
            exit(1);
        }
    }

    return true;
}

bool InterfaceBonmin::get_variables_linearity(Index n, Ipopt::TNLP::LinearityType* var_types)
{
    ipoptInterface->get_variables_linearity(n,var_types);
    return true;
}

bool InterfaceBonmin::get_constraints_linearity(Index m, Ipopt::TNLP::LinearityType* const_types)
{
    ipoptInterface->get_constraints_linearity(m,const_types);
    return true;
}

bool InterfaceBonmin::get_nlp_info(Index& n, Index&m, Index& nnz_jac_g,
                                   Index& nnz_h_lag, TNLP::IndexStyleEnum& index_style)
{
    ipoptInterface->get_nlp_info(n, m, nnz_jac_g, nnz_h_lag, index_style);
    return true;
}

bool InterfaceBonmin::get_bounds_info(Index n, Number* x_l, Number* x_u,
                                 Index m, Number* g_l, Number* g_u)
{
    ipoptInterface->get_bounds_info(n, x_l, x_u, m, g_l, g_u);
    return true;
}

bool InterfaceBonmin::get_starting_point(Index n, bool init_x, Number* x,
                                    bool init_z, Number* z_L, Number* z_U,
                                    Index m, bool init_lambda,
                                    Number* lambda)
{
    ipoptInterface->get_starting_point(n, init_x, x, init_z, z_L, z_U, m,init_lambda, lambda);
    return true;
}

bool InterfaceBonmin::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
    ipoptInterface->eval_f(n, x, new_x, obj_value);
    return true;
}

bool InterfaceBonmin::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
    ipoptInterface->eval_grad_f(n, x, new_x, grad_f);
    return true;
}

bool InterfaceBonmin::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
    ipoptInterface->eval_g(n, x, new_x, m, g);
    return true;
}

bool InterfaceBonmin::eval_jac_g(Index n, const Number* x, bool new_x,
                            Index m, Index nnz_jac, Index* iRow, Index *jCol,
                            Number* values)
{
    ipoptInterface->eval_jac_g(n, x, new_x, m, nnz_jac, iRow, jCol, values);
    return true;
}

bool InterfaceBonmin::eval_h(Index n, const Number* x, bool new_x,
                        Number obj_factor, Index m, const Number* lambda,
                        bool new_lambda, Index nele_hess, Index* iRow,
                        Index* jCol, Number* values)
{
    ipoptInterface->eval_h(n, x, new_x, obj_factor, m, lambda, new_lambda, nele_hess, iRow, jCol, values);
    return true;
}

void InterfaceBonmin::finalize_solution(TMINLP::SolverReturn status,
                                   Index n, const Number* x, Number obj_value)
{
    objectiveValue = obj_value;
    optimalSolution.clear();

    cout << "Problem status: " << status << endl; // 0 = Success
    cout << "Objective value: " << obj_value << endl;
    //cout << "Number n: " << n << endl;
    if (status == TMINLP::SolverReturn::SUCCESS)
    {
        for (int i = 0 ; i < n ; i++)
            optimalSolution.push_back(x[i]);
    }
    else
    {
        optimalSolution = std::vector<double>(n,0);
    }
}
