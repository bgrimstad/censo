/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "solveripopt.h"
#include "interfaceipopt.h"

using std::cout;
using std::endl;

namespace CENSO
{

SolverIpopt::SolverIpopt(ConstraintPtr constraints)
    : Solver(constraints),
      boundRelaxFactor(0)
{
    isDualAdequate = true;
}

// Initialize the optimizer
void SolverIpopt::initialize()
{
    // Setting up the Ipopt interface
    ipoptInterface = new InterfaceIpopt(constraints);

    // Create an instance of the IpoptApplication
    // We are using the factory, since this allows us to compile this
    // example with an Ipopt Windows DLL
    ipoptApplication = IpoptApplicationFactory();

    // Set Ipopt options
    ipoptApplication->Options()->SetNumericValue("tol", 1e-8); // Default is 1e-8
    ipoptApplication->Options()->SetNumericValue("constr_viol_tol", 1e-6); // Default is 1e-4
    ipoptApplication->Options()->SetNumericValue("acceptable_constr_viol_tol", 1e-6); // Default is 1e-4
    ipoptApplication->Options()->SetStringValue("mu_strategy", "adaptive"); // monotone (default), adaptive
    //ipoptApplication->Options()->SetStringValue("output_file", "ipopt.out");
    ipoptApplication->Options()->SetIntegerValue("print_level", 0);
    ipoptApplication->Options()->SetIntegerValue("max_iter", 500); // Default is 3000
    //ipoptApplication->Options()->SetStringValue("nlp_scaling_method", "none"); // none, gradient-based
    ipoptApplication->Options()->SetNumericValue("bound_relax_factor", boundRelaxFactor); // 1e-8 (Ipopt default). Important when using B-splines with closed domain

    // Calculation of Jacobian
    if (constraints->isJacobianCalculated())
    {
        if (constraints->isConstraintLinear())
        {
            // Linear constraints with constant Jacobian
            ipoptApplication->Options()->SetStringValue("jac_c_constant", "yes"); // Constant Jacobian equality constraints
            ipoptApplication->Options()->SetStringValue("jac_d_constant", "yes"); // Constant Jacobian inequality constraints

            //cout << "SolverIpopt: Solving problem with exact, constant constraint Jacobians." << endl;
        }
    }
    else
    {
        // NOTE: finite difference is implemented as a constraint decorator and could be applied here
        //ipoptApplication->Options()->SetStringValue("jacobian_approximation", "finite-difference-values");

        cout << "SolverIpopt: Warning! Solving problem with approximated constraint Jacobians." << endl;
    }

    // Calculation of Hessian
    if (constraints->isHessianCalculated())
    {
        // Exact Hessian
        ipoptApplication->Options()->SetStringValue("hessian_approximation", "exact"); // exact (default, no approx) or limited-memory (quasi-Newton)

        if (constraints->isConstraintLinear())
        {
            // Linear constraints with constant Hessian
            ipoptApplication->Options()->SetStringValue("hessian_constant", "yes");

            //cout << "SolverIpopt: Solving problem with exact, constant constraint Hessian." << endl;
        }
    }
    else
    {
        ipoptApplication->Options()->SetStringValue("hessian_approximation", "limited-memory"); // exact (default, no approx) or limited-memory (quasi-Newton)
        //ipoptApplication->Options()->SetStringValue("limited_memory_update_type", "bfgs"); // BFGS (default) or SR1 (not working well)

        //cout << "SolverIpopt: Solving problem with Hessian approximation (BFGS)." << endl;
    }

    // Derivative checks
    ipoptApplication->Options()->SetStringValue("check_derivatives_for_naninf", "no"); // no (default) or yes (may produce a lot of output)
    ipoptApplication->Options()->SetStringValue("derivative_test", "none"); // none (default), first-order, second-order

    // Branch-and-Bound specific (same as Bonmin uses)
    // NOTE: these options may increase solution time if on. Especially expect_infeasible_problem when set to "yes".
    //ipoptApplication->Options()->SetStringValue("expect_infeasible_problem", "yes"); // Detect infeasibility faster in BB
    //ipoptApplication->Options()->SetNumericValue("required_infeasibility_reduction", 0.1); // Speed up restoration phase

    // Initialize the IpoptApplication and process the options
    ApplicationReturnStatus status;
    status = ipoptApplication->Initialize();
    if (status != Solve_Succeeded)
    {
        cout << endl << endl << "*** Error during initialization!" << endl;
        //return (int) status;
        exit(1);
    }
    else
    {
        // Changing the status to initialized
        setInitialized(true);
    }
}

// Start the optimization
SolverResult SolverIpopt::runOptimizer()
{
    // Checking if everything is initialized
    if(!getInitialized())
        initialize();

    // Set up done, now let's run Ipopt
    ApplicationReturnStatus ipoptStatus;
    ipoptStatus = ipoptApplication->OptimizeTNLP(GetRawPtr(ipoptInterface));

//    Solve_Succeeded=0,
//    Solved_To_Acceptable_Level=1,
//    Infeasible_Problem_Detected=2,
//    Search_Direction_Becomes_Too_Small=3,
//    Diverging_Iterates=4,
//    User_Requested_Stop=5,
//    Feasible_Point_Found=6,

//    Maximum_Iterations_Exceeded=-1,
//    Restoration_Failed=-2,
//    Error_In_Step_Computation=-3,
//    Maximum_CpuTime_Exceeded=-4,
//    Not_Enough_Degrees_Of_Freedom=-10,
//    Invalid_Problem_Definition=-11,
//    Invalid_Option=-12,
//    Invalid_Number_Detected=-13,

//    Unrecoverable_Exception=-100,
//    NonIpopt_Exception_Thrown=-101,
//    Insufficient_Memory=-102,
//    Internal_Error=-199


    SolverStatus status = SolverStatus::ERROR;
    double optimalValue = ipoptInterface->getObjectiveValue();
    std::vector<double> optimalPoint = ipoptInterface->getOptimalSolution();

    if (ipoptStatus == Solve_Succeeded || ipoptStatus == Solved_To_Acceptable_Level)
    {
        status = SolverStatus::OPTIMAL;
    }
    else if (ipoptStatus == Feasible_Point_Found)
    {
        status = SolverStatus::FEASIBLE;
    }
    else if (ipoptStatus == Infeasible_Problem_Detected)
    {
        status = SolverStatus::INFEASIBLE;
        optimalValue = INF;
    }
    else if (ipoptStatus == Diverging_Iterates)
    {
        status = SolverStatus::UNBOUNDED;
        optimalValue = -INF;
    }
    else if(ipoptStatus == Maximum_Iterations_Exceeded)
    {
        status = SolverStatus::LIMIT_EXCEEDED;
    }
    else
    {
        // Error solving problem
        status = SolverStatus::ERROR;
    }

    SolverResult result(status, optimalValue, optimalPoint);
    result.lowerBoundDualVariables = ipoptInterface->getLowerBoundDuals();
    result.upperBoundDualVariables = ipoptInterface->getUpperBoundDuals();
    result.constraintDualVariables = ipoptInterface->getConstraintDuals();

    return result;
}

} // namespace CENSO
