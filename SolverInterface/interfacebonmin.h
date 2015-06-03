/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef BONMININTERFACE_H
#define BONMININTERFACE_H

#include "solver.h"
#include "Utils/definitions.h"

#include "interfaceipopt.h"
#include "solveripopt.h"

#include "BonTMINLP.hpp"

// NOTE: namespaces should not overlap
using namespace  Ipopt;
//using namespace Bonmin;

class InterfaceBonmin : public Bonmin::TMINLP
{
public:

    InterfaceBonmin(CENSO::ConstraintPtr constraints);

    InterfaceBonmin(const InterfaceBonmin &copy) = delete;
    InterfaceBonmin& operator=(const InterfaceBonmin &assign) = delete;

    // Virtual destructor
    ~InterfaceBonmin() override {}

    // Copy constructor
    // InterfaceBonmin(const InterfaceBonmin &other) {}
    // Assignment operator. no data = nothing to assign
    //BonminInterface& operator=(const BonminInterface&) {}

    /** Pass the type of the variables (INTEGER, BINARY, CONTINUOUS) to the optimizer.
        \param n size of var_types (has to be equal to the number of variables in the problem)
        \param var_types types of the variables (has to be filled by function).
    */
    bool get_variables_types(Index n, VariableType* var_types) override;

    /** Pass info about linear and nonlinear variables.*/
    bool get_variables_linearity(Index n, Ipopt::TNLP::LinearityType* var_types) override;

    /** Pass the type of the constraints (LINEAR, NON_LINEAR) to the optimizer.
        \param m size of const_types (has to be equal to the number of constraints in the problem)
        \param const_types types of the constraints (has to be filled by function).
    */
    bool get_constraints_linearity(Index m, Ipopt::TNLP::LinearityType* const_types) override;

    /** \name Overloaded functions defining a TNLP.
        This group of function implement the various elements needed to define and solve a TNLP.
        They are the same as those in a standard Ipopt NLP problem
    */

    /** Method to pass the main dimensions of the problem to Ipopt.
        \param n number of variables in problem.
        \param m number of constraints.
        \param nnz_jac_g number of nonzeroes in Jacobian of constraints system.
        \param nnz_h_lag number of nonzeroes in Hessian of the Lagrangean.
        \param index_style indicate wether arrays are numbered from 0 (C-style) or
        from 1 (Fortran).
        \return true in case of success.
    */
     bool get_nlp_info(Index& n, Index&m, Index& nnz_jac_g, Index& nnz_h_lag, TNLP::IndexStyleEnum& index_style) override;

    /** Method to pass the bounds on variables and constraints to Ipopt.
        \param n size of x_l and x_u (has to be equal to the number of variables in the problem)
        \param x_l lower bounds on variables (function should fill it).
        \param x_u upper bounds on the variables (function should fill it).
        \param m size of g_l and g_u (has to be equal to the number of constraints in the problem).
        \param g_l lower bounds of the constraints (function should fill it).
        \param g_u upper bounds of the constraints (function should fill it).
    \return true in case of success.*/
    bool get_bounds_info(Index n, Number* x_l, Number* x_u, Index m, Number* g_l, Number* g_u) override;

    /** Method to to pass the starting point for optimization to Ipopt.
        \param init_x do we initialize primals?
        \param x pass starting primal points (function should fill it if init_x is 1).
        \param m size of lambda (has to be equal to the number of constraints in the problem).
        \param init_lambda do we initialize duals of constraints?
        \param lambda lower bounds of the constraints (function should fill it).
        \return true in case of success.*/
    bool get_starting_point(Index n, bool init_x, Number* x, bool init_z, Number* z_L, Number* z_U, Index m, bool init_lambda, Number* lambda) override;

    /** Method which compute the value of the objective function at point x.
        \param n size of array x (has to be the number of variables in the problem).
        \param x point where to evaluate.
        \param new_x Is this the first time we evaluate functions at this point?
        (in the present context we don't care).
        \param obj_value value of objective in x (has to be computed by the function).
        \return true in case of success.*/
    bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value) override;

    /** Method which compute the gradient of the objective at a point x.
        \param n size of array x (has to be the number of variables in the problem).
        \param x point where to evaluate.
        \param new_x Is this the first time we evaluate functions at this point?
        (in the present context we don't care).
        \param grad_f gradient of objective taken in x (function has to fill it).
        \return true in case of success.*/
    bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f) override;

    /** Method which compute the value of the functions defining the constraints at a point
        x.
        \param n size of array x (has to be the number of variables in the problem).
        \param x point where to evaluate.
        \param new_x Is this the first time we evaluate functions at this point?
        (in the present context we don't care).
        \param m size of array g (has to be equal to the number of constraints in the problem)
        \param grad_f values of the constraints (function has to fill it).
        \return true in case of success.*/
    bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g) override;

    /** Method to compute the Jacobian of the functions defining the constraints.
        If the parameter values==NULL fill the arrays iCol and jRow which store the position of
        the non-zero element of the Jacobian.
        If the paramenter values!=NULL fill values with the non-zero elements of the Jacobian.
        \param n size of array x (has to be the number of variables in the problem).
        \param x point where to evaluate.
        \param new_x Is this the first time we evaluate functions at this point?
        (in the present context we don't care).
        \param m size of array g (has to be equal to the number of constraints in the problem)
        \param grad_f values of the constraints (function has to fill it).
        \return true in case of success.*/
    bool eval_jac_g(Index n, const Number* x, bool new_x, Index m, Index nele_jac, Index* iRow, Index *jCol, Number* values) override;

    /** Method to compute the Jacobian of the functions defining the constraints.
        If the parameter values==NULL fill the arrays iCol and jRow which store the position of
        the non-zero element of the Jacobian.
        If the paramenter values!=NULL fill values with the non-zero elements of the Jacobian.
        \param n size of array x (has to be the number of variables in the problem).
        \param x point where to evaluate.
        \param new_x Is this the first time we evaluate functions at this point?
        (in the present context we don't care).
        \param m size of array g (has to be equal to the number of constraints in the problem)
        \param grad_f values of the constraints (function has to fill it).
        \return true in case of success.*/
    bool eval_h(Index n, const Number* x, bool new_x,
                Number obj_factor, Index m, const Number* lambda,
                bool new_lambda, Index nele_hess, Index* iRow,
                Index* jCol, Number* values) override;

    /** Method called by Ipopt at the end of optimization.*/
    void finalize_solution(TMINLP::SolverReturn status,
                           Index n, const Number* x, Number obj_value) override;


    /** Not implemented yet */
    const SosInfo * sosConstraints() const override {return NULL;}
    const BranchingInfo* branchingInfo() const override {return NULL;}

    // Get functions
    std::vector<double> getOptimalSolution() const { return optimalSolution; }
    double getObjectiveValue() const { return objectiveValue; }

private:

    SmartPtr<InterfaceIpopt> ipoptInterface;

    std::vector<double> optimalSolution;
    double objectiveValue;

};

#endif // BONMININTERFACE_H
