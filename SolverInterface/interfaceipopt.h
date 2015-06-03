/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef IPOPTINTERFACE_H
#define IPOPTINTERFACE_H

#include "Utils/definitions.h"
#include "solver.h"
#include "IpTNLP.hpp"

// NOTE: namespaces should not overlap
using namespace Ipopt;

class InterfaceIpopt : public TNLP
{
public:

    InterfaceIpopt(CENSO::ConstraintPtr constraints);
    InterfaceIpopt(const InterfaceIpopt &copy) = delete;
    InterfaceIpopt& operator=(const InterfaceIpopt &assign) = delete;

    ~InterfaceIpopt() override {}

    bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g, Index& nnz_h_lag, IndexStyleEnum& index_style) override;

    bool get_bounds_info(Index n, Number* x_l, Number* x_u, Index m, Number* g_l, Number* g_u) override;

    bool get_starting_point(Index n, bool init_x, Number* x,bool init_z, Number* z_L, Number* z_U, Index m, bool init_lambda, Number* lambda) override;

    bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value) override;

    bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f) override;

    bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g) override;

    /** Method to return:
       *   1) The structure of the jacobian (if "values" is NULL)
       *   2) The values of the jacobian (if "values" is not NULL)
       */
    bool eval_jac_g(Index n, const Number* x, bool new_x,Index m, Index nele_jac, Index* iRow, Index *jCol,Number* values) override;

    /** Method to return:
       *   1) The structure of the hessian of the lagrangian (if "values" is NULL)
       *   2) The values of the hessian of the lagrangian (if "values" is not NULL)
       */
    bool eval_h(Index n, const Number* x, bool new_x,Number obj_factor, Index m, const Number* lambda,bool new_lambda, Index nele_hess, Index* iRow,Index* jCol, Number* values) override;

    /** This method is called when the algorithm is complete so the TNLP can store/write the solution */
    void finalize_solution(SolverReturn status, Index n, const Number* x, const Number* z_L, const Number* z_U,Index m, const Number* g, const Number* lambda, Number obj_value, const IpoptData* ip_data, IpoptCalculatedQuantities* ip_cq) override;

    /** overload this method to return the variables linearity
     * (TNLP::Linear or TNLP::NonLinear). The var_types
     *  array should be allocated with length at least n. (default implementation
     *  just return false and does not fill the array).*/
    bool get_variables_linearity(Index n, LinearityType* var_types) override;

    /** overload this method to return the constraint linearity.
     * array should be alocated with length at least n. (default implementation
     *  just return false and does not fill the array).*/
    bool get_constraints_linearity(Index m, LinearityType* const_types) override;


    /** @name Methods for quasi-Newton approximation.  If the second
     *  derivatives are approximated by Ipopt, it is better to do this
     *  only in the space of nonlinear variables.  The following
     *  methods are call by Ipopt if the quasi-Newton approximation is
     *  selected.  If -1 is returned as number of nonlinear variables,
     *  Ipopt assumes that all variables are nonlinear.  Otherwise, it
     *  calls get_list_of_nonlinear_variables with an array into which
     *  the indices of the nonlinear variables should be written - the
     *  array has the lengths num_nonlin_vars, which is identical with
     *  the return value of get_number_of_nonlinear_variables().  It
     *  is assumed that the indices are counted starting with 1 in the
     *  FORTRAN_STYLE, and 0 for the C_STYLE. */
    Index get_number_of_nonlinear_variables() override;
    bool get_list_of_nonlinear_variables(Index num_nonlin_vars, Index* pos_nonlin_vars) override;

    /*
     * Get functions
     */
    CENSO::ConstraintPtr getConstraints() const { return constraints; }
    std::vector<double> getOptimalSolution() const { return optimalSolution; }
    double getObjectiveValue() const { return objectiveValue; }

    std::vector<double> getLowerBoundDuals() const { return lowerBoundDuals; }
    std::vector<double> getUpperBoundDuals() const { return upperBoundDuals; }
    std::vector<double> getConstraintDuals() const { return constraintDuals; }

//    DoubleVec optimalSolution;
//    int iterations;

//    virtual bool intermediate_callback(AlgorithmMode mode,
//                                       Index iter, Number obj_value,
//                                       Number inf_pr, Number inf_du,
//                                       Number mu, Number d_norm,
//                                       Number regularization_size,
//                                       Number alpha_du, Number alpha_pr,
//                                       Index ls_trials,
//                                       const IpoptData* ip_data,
//                                       IpoptCalculatedQuantities* ip_cq)
//    {
//      iterations = iter;
//      return true;
//    }

private:

    CENSO::ConstraintPtr constraints;

    int nnzConstraintJacobian;
    std::vector<int> iRowGradient;
    std::vector<int> jColGradient;

    int nnzHessian;
    std::vector<int> eqnrHessian;
    std::vector<int> iRowvHessian;
    std::vector<int> jColvHessian;
    std::vector<int> valuePositionHessian;

    std::vector<double> optimalSolution;
    double objectiveValue;
    std::vector<double> lowerBoundDuals;
    std::vector<double> upperBoundDuals;
    std::vector<double> constraintDuals;

    void vecToNum(std::vector<double> &v, int sizen, Number *n);
    void vecToIndex(std::vector<int> &v, int sizen, Index *in);
    void numToVec(const Number *nr, int n, std::vector<double> &v);
    void getJacobianInfo();
    void getHessianInfo();

    void vecToNum(CENSO::DenseVector &v, int sizen, Number *n);
    void numToVec(const Number *nr, int n, CENSO::DenseVector &v);

};

void ipoptInterfaceTesting(CENSO::ConstraintPtr constraints, std::vector<double> lb, std::vector<double> ub);

#endif // IPOPTINTERFACE_H
