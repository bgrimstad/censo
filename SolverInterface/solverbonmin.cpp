/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "solverbonmin.h"

using std::cout;
using std::endl;

namespace CENSO
{

void SolverBonmin::initialize()
{
    bonminInterface = new InterfaceBonmin(constraints);

    bonmin.initializeOptionsAndJournalist();
    bonmin.readOptionsFile();// This reads the default file "bonmin.opt"

//    bonmin.options()->SetIntegerValue("bb_log_interval",0); // 0 - inf
//    bonmin.options()->SetIntegerValue("bb_log_level",0);    // 0 - 5
//    bonmin.options()->SetIntegerValue("nlp_log_level",0);   // 0 - 2
//    bonmin.options()->SetIntegerValue("lp_log_level",0);    // 0 - 4
//    bonmin.options()->SetIntegerValue("nlp_log_at_root",0); // 0 - 12

    bonmin.options()->SetStringValue("algorithm","B-BB"); // B-BB (default), B-iFP
    bonmin.options()->SetIntegerValue("num_resolve_at_root", 2);
    bonmin.options()->SetIntegerValue("num_resolve_at_node", 1);
    //bonmin.options()->SetNumericValue("allowable_fraction_gap", -1.0);

    // Primal heuristics
    bonmin.options()->SetStringValue("heuristic_feasibility_pump", "yes");
    bonmin.options()->SetStringValue("pump_for_minlp", "yes");

    // Calculation of Jacobian
    if (constraints->isJacobianCalculated())
    {
        if (constraints->isConstraintLinear())
        {
            // Linear constraints with constant Jacobian
            bonmin.options()->SetStringValue("jac_c_constant", "yes"); // Constant Jacobian equality constraints
            bonmin.options()->SetStringValue("jac_d_constant", "yes"); // Constant Jacobian inequality constraints

            //cout << "SolverIpopt: Solving problem with exact, constant constraint Jacobians." << endl;
        }
    }
    else
    {
        cout << "SolverBonmin: Warning! Solving problem with approximated constraint Jacobians." << endl;

        bonmin.options()->SetStringValue("jacobian_approximation", "finite-difference-values");
        bonmin.options()->SetNumericValue("findiff_perturbation", 1e-8);

        // TEMPORARY SETTINGS WHILE TESTING WITH NO JACOBIANS
//        bonmin.options()->SetIntegerValue("max_iter", 10000);
//        bonmin.options()->SetIntegerValue("max_consecutive_infeasible", 5);
//        bonmin.options()->SetIntegerValue("num_resolve_at_infeasibles", 2);
    }

    /*
     * Variable selection strategies:
     * most-fractional: Choose most fractional variable
     * strong-branching: Perform strong branching
     * reliability-branching: Use reliability branching
     * qp-strong-branching: Perform strong branching with QP approximation
     * lp-strong-branching: Perform strong branching with LP approximation
     * nlp-strong-branching: Perform strong branching with NLP approximation
     * osi-simple: Osi method to do simple branching
     * osi-strong: Osi method to do strong branching
     * random: Method to choose branching variable randomly
     */
    //bonmin.options()->SetStringValue("variable_selection","most-fractional");
    //bonmin.options()->SetIntegerValue("number_strong_branch", 0);

//    std::string algo("");
//    bonmin.options()->GetStringValue("algorithm", algo, "bonmin");
//    std::cout << algo << std::endl;

    // Set Ipopt options
    bonmin.options()->SetNumericValue("constr_viol_tol", 1e-6); // Default is 1e-4
    bonmin.options()->SetNumericValue("acceptable_constr_viol_tol", 1e-6); // Default is 1e-4
    bonmin.options()->SetNumericValue("bound_relax_factor", 0); // Important when using B-splines
//    double value = 0;
//    bonmin.options()->GetNumericValue("bound_relax_factor", value, "ipopt");
//    std::cout << "Bonmin: Ipopt bound relax factor: " << value << std::endl;

    // Now initialize from tminlp
    bonmin.initialize(GetRawPtr(bonminInterface));
    setInitialized(true);
}

CENSO::SolverResult SolverBonmin::runOptimizer()
{
    if (!getInitialized())
        initialize();

    Bonmin::Bab bb;
    //    bb(bonmin);
    bb.branchAndBound(bonmin);

    /*
     * MipStatuses
     *
     * FeasibleOptimal: Optimum solution has been found and its optimality proved.
     * ProvenInfeasible: Problem has been proven to be infeasible.
     * Feasible: An integer solution to the problem has been found.
     * UnboundedOrInfeasible: Coninuous relaxation is unbounded.
     * NoSolutionKnown: No feasible solution to the problem is known.
     */

    CENSO::SolverStatus status = CENSO::SolverStatus::ERROR;

    double optimalValue = bonminInterface->getObjectiveValue();
    std::vector<double> optimalPoint = bonminInterface->getOptimalSolution();

    if (bb.mipStatus() == Bonmin::Bab::FeasibleOptimal)
    {
        status = CENSO::SolverStatus::OPTIMAL;
    }
    else if (bb.mipStatus() == Bonmin::Bab::Feasible)
    {
        status = CENSO::SolverStatus::FEASIBLE;
    }
    else if (bb.mipStatus() == Bonmin::Bab::ProvenInfeasible
             || bb.mipStatus() == Bonmin::Bab::NoSolutionKnown)
    {
        status = CENSO::SolverStatus::INFEASIBLE;
    }
    else if (bb.mipStatus() == Bonmin::Bab::UnboundedOrInfeasible)
    {
        status = CENSO::SolverStatus::INFEASIBLE;
    }
    else
    {
        status = CENSO::SolverStatus::ERROR;
    }

    CENSO::SolverResult result(status, optimalValue, optimalPoint);
//    result.lowerBoundDualVariables = bonminInterface->getLowerBoundDuals();
//    result.upperBoundDualVariables = bonminInterface->getUpperBoundDuals();
//    result.constraintDualVariables = bonminInterface->getConstraintDuals();

    return result;
}

} // namespace CENSO
