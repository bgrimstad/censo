/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "solvergurobi.h"

using std::cout;
using std::endl;

namespace CENSO
{

SolverGurobi::SolverGurobi(ConstraintPtr constraints)
    : Solver(constraints)
{
    isDualAdequate = true;
}

// Initialize the optimizer
void SolverGurobi::initialize()
{
    ProblemClass probClass = constraints->assessClass();
    if (probClass != ProblemClass::LP
            && probClass != ProblemClass::MILP)
    {
        cout << "Problem class not LP or MILP: cannot use Gurobi!" << endl;
        exit(1);
    }

    setInitialized(true);
}

// Start the optimization
SolverResult SolverGurobi::runOptimizer()
{
    if (!getInitialized())
        initialize();

    try
    {
        // Create Gurobi environment and set parameters
        GRBEnv env = GRBEnv();
        env.set(GRB_IntParam_OutputFlag, 0);

        GRBModel model = GRBModel(env);

        // Get problem info
        int numVars = constraints->getNumVariables();
        int numConstraints = constraints->getNumConstraints();

        // Get variables
        auto variables = constraints->getVariables();

        // Create array of model variables
        GRBVar vars[numVars];
        for (int i = 0; i < numVars; i++)
        {
            //vars[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);

            // Set variable type
            char type = GRB_CONTINUOUS;
            if (variables.at(i)->getType() == VariableType::BINARY)
            {
                type = GRB_BINARY;
            }
            else if (variables.at(i)->getType() == VariableType::INTEGER)
            {
                type = GRB_INTEGER;
            }
            vars[i] = model.addVar(variables.at(i)->getLowerBound(),
                                   variables.at(i)->getUpperBound(),
                                   variables.at(i)->getCost(),
                                   type);
        }

        // Integrate variables into model
        model.update();

        // Set starting points (does not help much...)
        for (int i = 0; i < numVars; i++)
            vars[i].set(GRB_DoubleAttr_Start, variables.at(i)->getValue());

        /*
         * Add constraints Ax <= b (or Ax = b)
         * by evaluating gradient and build A matrix
         */
        DenseVector x = DenseVector::Zero(numVars);
        DenseVector dx = constraints->evalJacobian(x);

        // Get constraint bounds
        std::vector<double> clb;
        std::vector<double> cub;
        constraints->getConstraintBounds(clb,cub);

        std::vector<int> rowGradient, colGradient;
        constraints->structureJacobian(rowGradient, colGradient);
        int nnzJacobian = constraints->getNumNonZerosJacobian();

        // Add constraints one row at the time
        for (int row = 0; row < numConstraints; row++)
        {
            // Build constraint
            GRBLinExpr expr = 0;

            // Loop through all non-zeros (inefficient)
            for (int i = 0; i < nnzJacobian; i++)
            {
                if (rowGradient.at(i) == row)
                {
                    int j = colGradient.at(i);
                    expr += dx(i)*vars[j];
                }
            }

            // Add constraint to model
            if (clb.at(row) == cub.at(row))
            {
                model.addConstr(expr, GRB_EQUAL, cub.at(row));
            }
            else
            {
                model.addConstr(expr, GRB_LESS_EQUAL, cub.at(row));
            }
        }

        // More efficient method - avoids dense matrix
//        std::vector<int> rows = {1,1,1,2,2,3,4,4,4,4,4,5};
//        std::vector<int>::iterator start,stop;
//        start = rows.begin();
//        stop = start;
//        while (start != rows.end())
//        {
//            while (stop != rows.end())
//            {
//                if (*stop == *start)
//                    ++stop;
//                else
//                    break;
//            }
//            for (std::vector<int>::iterator it = start; it != stop; ++it)
//                cout << *it << endl;

//            start = stop;
//        }

        model.update();

        assert(numVars == model.get(GRB_IntAttr_NumVars));
        assert(numConstraints == model.get(GRB_IntAttr_NumConstrs));

        // Optimize model
        model.optimize();

        // Check status
        int optimstatus = model.get(GRB_IntAttr_Status);

        if (optimstatus == GRB_INF_OR_UNBD)
        {
            model.getEnv().set(GRB_IntParam_Presolve, 0);
            model.optimize();
            optimstatus = model.get(GRB_IntAttr_Status);
        }

        // Create result object
        SolverResult result(SolverStatus::ERROR, INF, std::vector<double>(numVars,0));

        // Check Gurobi status
        if (optimstatus == GRB_OPTIMAL)
        {
            result.status = SolverStatus::OPTIMAL;

            // Get solution info
            result.objectiveValue = model.get(GRB_DoubleAttr_ObjVal);

            std::vector<double> optimalSolution;
            for (int i = 0; i < numVars; i++)
            {
                optimalSolution.push_back(vars[i].get(GRB_DoubleAttr_X));
            }

            result.primalVariables = optimalSolution;

            /*
             * Reduced costs and constraint duals are
             * only available for continuous models
             */
            std::vector<double> reducedCosts;
            std::vector<double> constraintDuals;
            if (!model.get(GRB_IntAttr_IsMIP))
            {
                for (int i = 0; i < numVars; i++)
                {
                    // Get reduced costs (related to range constraint duals)
                    reducedCosts.push_back(vars[i].get(GRB_DoubleAttr_RC));
                }

                for (int i = 0; i < numConstraints; i++)
                {
                    GRBConstr c = model.getConstr(i);
                    double pi = c.get(GRB_DoubleAttr_Pi);
                    constraintDuals.push_back(pi);
                }
            }

            result.lowerBoundDualVariables = reducedCosts;
            result.upperBoundDualVariables = reducedCosts;
            result.constraintDualVariables = constraintDuals;

            return result;
        }
        else if (optimstatus == GRB_INFEASIBLE)
        {
            result.status = SolverStatus::INFEASIBLE;
            result.objectiveValue = INF;
            // compute and write out IIS
            // model.computeIIS();
            // model.write("problem.lp");
            return result;
        }
        else if (optimstatus == GRB_UNBOUNDED)
        {
            result.status = SolverStatus::UNBOUNDED;
            result.objectiveValue = -INF;
            return result;
        }
        else
        {
            result.status = SolverStatus::ERROR;
            result.objectiveValue = INF;
            return result;
        }
    }
    catch(GRBException e)
    {
        cout << "SolverGurobi: Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
        return SolverResult(SolverStatus::ERROR, INF, std::vector<double>(constraints->getNumVariables(),0));
    }
    catch (...)
    {
        cout << "SolverGurobi: Error during optimization!" << endl;
        return SolverResult(SolverStatus::ERROR, INF, std::vector<double>(constraints->getNumVariables(),0));
    }
}

} // namespace CENSO
