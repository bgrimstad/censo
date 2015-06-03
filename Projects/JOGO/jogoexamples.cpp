/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "jogoexamples.h"
#include "Utils/definitions.h"
#include "OptimizationProblem/constraintlinear.h"
#include "OptimizationProblem/constraintbspline.h"
#include "OptimizationProblem/constraintquadratic.h"
#include "OptimizationProblem/constraintpolynomial.h"
#include "Utils/bsplinepoly.h"
#include "bspline.h"
#include "datatable.h"
#include "BranchAndBound/branchandbound.h"
#include "SolverInterface/solvergurobi.h"
#include "SolverInterface/solveripopt.h"

using std::cout;
using std::endl;
using Splinter::BSpline;
using Splinter::BSplineType;
using Splinter::DataTable;

namespace CENSO
{

void saveDataTable(DataTable &data, std::string filename)
{
    std::vector<std::set<double>> grid = data.getGrid();

    std::ofstream myfile(filename);

    // Write to file
    if (myfile.is_open())
    {
        for (unsigned int i=0; i<grid.size(); i++)
        {
            std::set<double> x = grid.at(i);
            int j = 0;
            for (auto& xj : x)
            {
                if (j>0)
                    myfile << ",";
                myfile << xj;
                j++;
            }
            myfile << "\n";
        }

        std::vector<double> y = data.getVectorY();

        for (unsigned int i=0; i<y.size(); i++)
        {
            myfile << y.at(i) << "\n";
        }

        myfile.flush();
    }
    else cout << "Unable to open file";
}

/*
 * Used in black-box example
 */
void sampleMichalewicz()
{
    double pi = atan(1)*4;

    auto micha = [pi](double x1, double x2)
    {
        return -std::sin(x1)*std::pow(std::sin(x1*x1/pi), 20) - std::sin(x2)*std::pow(std::sin(2*x2*x2/pi), 20);
    };

    std::vector<double> costs = {0, 0, 1};
    std::vector<double> lb {0, 0, -INF};
    std::vector<double> ub {pi, pi, INF};

    std::vector<VariablePtr> vars;

    for(int i = 0; i < 3; i++)
    {
        auto var = std::make_shared<Variable>(costs.at(i), lb.at(i), ub.at(i));
        vars.push_back(var);
    }

    DataTable data;

    auto x1v = linspace(lb.at(0), ub.at(0), 50);
    auto x2v = linspace(lb.at(1), ub.at(1), 50);
    for (auto x1 : x1v)
    {
        for (auto x2 : x2v)
        {
            std::vector<double> xvec = {x1, x2};
            double y = micha(x1, x2);

            data.addSample(xvec, y);
        }
    }

    // Create B-spline
    BSpline bs(data, BSplineType::CUBIC_FREE);

    // Error
    DataTable error;

    // Test accuracy
    auto x1v2 = linspace(lb.at(0), ub.at(0), 500);
    auto x2v2 = linspace(lb.at(1), ub.at(1), 500);
    double eabs_max = 0;
    double eabs_avg = 0;
    double erel_max = 0;
    double erel_avg = 0;
    double numPoints = 0;
    for (auto x1 : x1v2)
    {
        for (auto x2 : x2v2)
        {
            DenseVector xvec(2);
            xvec(0) = x1;
            xvec(1) = x2;

            double ybs = bs.eval(xvec);

            double y = micha(x1, x2);

            error.addSample(xvec,y-ybs);

            double eabs = std::abs(y-ybs);
            eabs_avg += eabs;
            if (eabs > eabs_max)
                eabs_max = eabs;

            double erel = 0;
            if (std::abs(y) > 0 || std::abs(ybs) > 0)
            {
                erel = std::abs(y-ybs)/std::max(std::abs(y), std::abs(ybs));
                if (erel > erel_max)
                    erel_max = erel;
                erel_avg += erel;
            }
            numPoints++;
        }
    }

    eabs_avg = eabs_avg/numPoints;
    erel_avg /= numPoints;

    if (eabs_max < 1e-3)
    {
        cout << "Test Passed!" << endl;
    }
    else
    {
        cout << "Test Failed!" << endl;
    }
    cout << "Max absolute error: " << eabs_max << endl;
    cout << "Avg absolute error: " << eabs_avg << endl;
    cout << "Max relative error: " << erel_max << endl;
    cout << "Avg relative error: " << erel_avg << endl;

    saveDataTable(error, "micha.txt");

    // Constraints
    ConstraintSetPtr constraints = std::make_shared<ConstraintSet>();

    ConstraintPtr cbspline = std::make_shared<ConstraintBSpline>(vars, bs, false);

    constraints->add(cbspline);

    // Attempt to solve to global optimality
    BB::BranchAndBound bnb(constraints);
    SolverResult result = bnb.optimize();
    cout << result.objectiveValue << endl;
    cout << result.primalVariables;
}

/*
 * Used in black-box example
 */
void samplePump()
{
    std::vector<double> alpha = {19.9, 1.21, 6.52};
    std::vector<double> beta = {0.161, 0.0644, 0.102};
    std::vector<double> gamma = {-0.000561, -0.000564, -0.000232};

    std::vector<double> a = {629.0, 215.0, 361.0};
    std::vector<double> b = {0.696, 2.95, 0.53};
    std::vector<double> c = {-0.0116, -0.115, -0.00946};

    double vtot = 350;
    double wmax = 2950;

    auto power = [wmax](double v, double w, double alpha, double beta, double gamma)
    {
        return alpha*std::pow(w/wmax, 3) + beta*std::pow(w/wmax, 2)*v + gamma*(w/wmax)*v*v;
    };

    auto pressure = [wmax](double v, double w, double a, double b, double c)
    {
        return a*std::pow(w/wmax, 2) + b*(w/wmax)*v + c*v*v;
    };

    DataTable data_pow;
    DataTable data_pres;

    int pump = 1; // 1,2,3

    auto v = linspace(0, vtot, 30);
    auto w = linspace(0, wmax, 30);
    for (auto vi : v)
    {
        for (auto wi : w)
        {
            std::vector<double> x = {vi, wi};
            double pow = power(vi, wi, alpha.at(pump), beta.at(pump), gamma.at(pump));
            double dp = pressure(vi, wi, a.at(pump), b.at(pump), c.at(pump));

            data_pow.addSample(x, pow);
            data_pres.addSample(x, dp);
        }
    }

    // Create B-spline
    BSpline bspline_pow(data_pow, BSplineType::QUADRATIC_FREE);
    BSpline bspline_dp(data_pres, BSplineType::QUADRATIC_FREE);

    // Error
    DataTable error_pow;
    DataTable error_dp;

    // Test accuracy
    double eabs_max_pow = 0;
    double eabs_avg_pow = 0;
    double erel_max_pow = 0;
    double erel_avg_pow = 0;
    double eabs_max_dp = 0;
    double eabs_avg_dp = 0;
    double erel_max_dp = 0;
    double erel_avg_dp = 0;

//    double erel_max = 0;
//    double erel_avg = 0;

    double num_points = 0;

    auto v2 = linspace(0, vtot, 500);
    auto w2 = linspace(0, wmax, 500);
    for (auto vi : v2)
    {
        for (auto wi : w2)
        {
            DenseVector x(2);
            x(0) = vi;
            x(1) = wi;

            double bs_pow = bspline_pow.eval(x);
            double bs_dp = bspline_dp.eval(x);

            double pow = power(vi, wi, alpha.at(pump), beta.at(pump), gamma.at(pump));
            double dp = pressure(vi, wi, a.at(pump), b.at(pump), c.at(pump));

            error_pow.addSample(x,pow-bs_pow);
            error_dp.addSample(x,dp-bs_dp);

            double eabs_pow = std::abs(pow-bs_pow);
            double eabs_dp = std::abs(dp-bs_dp);
            eabs_avg_pow += eabs_pow;
            eabs_avg_dp += eabs_dp;

            if (eabs_pow > eabs_max_pow)
                eabs_max_pow = eabs_pow;

            if (eabs_dp > eabs_max_dp)
                eabs_max_dp = eabs_dp;

            double erel_pow = 0;
            if (std::abs(pow) > 0 && std::abs(bs_pow) > 0)
                erel_pow = std::abs(pow-bs_pow)/std::max(std::abs(pow),std::abs(bs_pow));
            if (erel_pow > erel_max_pow)
                erel_max_pow = erel_pow;
            erel_avg_pow += erel_pow;

            double erel_dp = 0;
            if (std::abs(dp) > 0 && std::abs(bs_dp) > 0)
                erel_dp = std::abs(dp-bs_dp)/std::max(std::abs(dp),std::abs(bs_dp));
            if (erel_dp > erel_max_dp)
                erel_max_dp = erel_dp;
            erel_avg_dp += erel_dp;

            num_points++;
        }
    }

    eabs_avg_pow /= num_points;
    eabs_avg_dp /= num_points;
    erel_avg_pow /= num_points;
    erel_avg_dp /= num_points;
    //erel_avg_dp = std::sqrt(erel_avg_dp);

    cout << "Avg absolute error pow: " << eabs_avg_pow << endl;
    cout << "Max absolute error pow: " << eabs_max_pow << endl;
    cout << "Avg relative error pow: " << erel_avg_pow << endl;
    cout << "Max relative error pow: " << erel_max_pow << endl;
    cout << endl;
    cout << "Avg absolute error pres: " << eabs_avg_dp << endl;
    cout << "Max absolute error pres: " << eabs_max_dp << endl;
    cout << "Avg relative error pres: " << erel_avg_dp << endl;
    cout << "Max relative error pres: " << erel_max_dp << endl;

//    saveDataTable(error, "micha.txt");
}

void pumpSynthesis(unsigned int grid)
{
    const int N = 3;
    const int numCont = N*5;
    const int numInt = N*3;
    const int numVars = numCont + numInt;
    const int numAux = 2*N+N+N+N+N;
    const int totVars = numVars + numAux;

    std::vector<double> cost = {38900, 15300, 20700}; // 20100 in Adjiman
    std::vector<double> ci = {6329.03, 2489.31, 3270.27};
    std::vector<double> ciprime = {1800, 1800, 1800};

    std::vector<double> alpha = {19.9, 1.21, 6.52};
    std::vector<double> beta = {0.161, 0.0644, 0.102};
    std::vector<double> gamma = {-0.000561, -0.000564, -0.000232};

    std::vector<double> ap = {629.0, 215.0, 361.0};
    std::vector<double> bp = {0.696, 2.95, 0.53};
    std::vector<double> cp = {-0.0116, -0.115, -0.00946};

    std::vector<double> pmax = {80, 25, 45}; // 35 in Adjiman

    double vtot = 350;
    double dptot = 400;
    double wmax = 2950;

    /*
     * i = {1,...,N}
     *
     * Cont. variables:
     * x_i
     * v_i
     * w_i
     * p_i
     * dp_i
     *
     * Int. variables:
     * Np_i
     * Ns_i
     * z_i
     *
     * In total: N*(5+3) variables
     */
    int ix = 0;
    int iv = N;
    int iw = 2*N;
    int ip = 3*N;
    int idp = 4*N;
    int inp = 5*N;
    int ins = 6*N;
    int iz = 7*N;
    int il1 = 8*N; // l1 = np*ns*z
    int il2 = 9*N; // l2 = l1*p
    int il3 = 10*N; // l3 = alpha*(w/wmax)^3 + beta*(w/wmax)^2*v+gamma*(w/wmax)*v^2
    int il4 = 11*N; // l4 = a*(w/wmax)^2 + b*(w/wmax)v + c*v^2
    int il5 = 12*N; // l5 = v*np
    int il6 = 13*N; // l6 = dp*ns

    std::vector<double> costs(totVars, 0);
    std::vector<double> lb(totVars, 0);
    std::vector<double> ub(totVars, INF);
    std::vector<VariableType> vt(totVars, VariableType::CONTINUOUS);

    for (int i = 0; i < N; i++)
    {
        ub.at(ix+i) = 1;
        ub.at(iv+i) = vtot;
        ub.at(iw+i) = wmax;
        ub.at(ip+i) = pmax.at(i);
        ub.at(idp+i) = dptot;
        ub.at(inp+i) = 3;
        ub.at(ins+i) = 3;
        ub.at(iz+i) = 1;
        ub.at(il1+i) = 3*3*1;
    }

    for (int i = 0; i < N; i++)
    {
        costs.at(il1+i) = ci.at(i);
        costs.at(il2+i) = ciprime.at(i);
    }

    for (int i = 0; i < N; i++)
    {
        vt.at(inp+i) = VariableType::INTEGER;
        vt.at(ins+i) = VariableType::INTEGER;
        vt.at(iz+i) = VariableType::BINARY;
    }

    std::vector<VariablePtr> vars;
    for (int i = 0; i < totVars; i++)
    {
        auto var = std::make_shared<Variable>(costs.at(i), lb.at(i), ub.at(i), vt.at(i));
        vars.push_back(var);
    }

    // Constraint set
    ConstraintSetPtr cs = std::make_shared<ConstraintSet>();
    ConstraintSetPtr cs_approx = std::make_shared<ConstraintSet>();

    {
        // sum x_i = 1
        std::vector<VariablePtr> cvars;
        for (int i = 0; i < N; i++)
            cvars.push_back(vars.at(ix+i)); // x(i)

        DenseMatrix A = DenseMatrix::Ones(1,N);
        DenseVector b = DenseVector::Ones(1);

        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, true);

        cs->add(lincon);
        cs_approx->add(lincon);
    }

    {
        // p - l3 = 0
        for (int i = 0; i < N; i++)
        {
            std::vector<VariablePtr> cvars;
            cvars.push_back(vars.at(ip+i)); // p(i)
            cvars.push_back(vars.at(il3+i)); // l3(i)

            DenseMatrix A = DenseMatrix::Ones(1,2);
            A(0,1) = -1;
            DenseVector b = DenseVector::Zero(1);

            ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, true);

            cs->add(lincon);
            cs_approx->add(lincon);
        }
    }

    {
        // dp - l4 = 0
        for (int i = 0; i < N; i++)
        {
            std::vector<VariablePtr> cvars;
            cvars.push_back(vars.at(idp+i)); // dp(i)
            cvars.push_back(vars.at(il4+i)); // l4(i)

            DenseMatrix A = DenseMatrix::Ones(1,2);
            A(0,1) = -1;
            DenseVector b = DenseVector::Zero(1);

            ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, true);

            cs->add(lincon);
            cs_approx->add(lincon);
        }
    }

    {
        // l5 - vtot*x = 0
        for (int i = 0; i < N; i++)
        {
            std::vector<VariablePtr> cvars;
            cvars.push_back(vars.at(il5+i)); // l5(i)
            cvars.push_back(vars.at(ix+i)); // x(i)

            DenseMatrix A = DenseMatrix::Ones(1,2);
            A(0,1) = -vtot;
            DenseVector b = DenseVector::Zero(1);

            ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, true);

            cs->add(lincon);
            cs_approx->add(lincon);
        }
    }

    {
        // dptot*z - l6 = 0
        for (int i = 0; i < N; i++)
        {
            std::vector<VariablePtr> cvars;
            cvars.push_back(vars.at(iz+i)); // z(i)
            cvars.push_back(vars.at(il6+i)); // l6(i)

            DenseMatrix A = DenseMatrix::Ones(1,2);
            A(0,0) = dptot; A(0,1) = -1;
            DenseVector b = DenseVector::Zero(1);

            ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, true);

            cs->add(lincon);
            cs_approx->add(lincon);
        }
    }

    {
        // p - pmax*z <= 0
        for (int i = 0; i < N; i++)
        {
            std::vector<VariablePtr> cvars;
            cvars.push_back(vars.at(ip+i)); // p(i)
            cvars.push_back(vars.at(iz+i)); // z(i)

            DenseMatrix A = DenseMatrix::Ones(1,2);
            A(0,1) = -pmax.at(i);
            DenseVector b = DenseVector::Zero(1);

            ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, false);

            cs->add(lincon);
            cs_approx->add(lincon);
        }
    }

    {
        // dp - dptot*z <= 0
        for (int i = 0; i < N; i++)
        {
            std::vector<VariablePtr> cvars;
            cvars.push_back(vars.at(idp+i)); // dp(i)
            cvars.push_back(vars.at(iz+i)); // z(i)

            DenseMatrix A = DenseMatrix::Ones(1,2);
            A(0,1) = -dptot;
            DenseVector b = DenseVector::Zero(1);

            ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, false);

            cs->add(lincon);
            cs_approx->add(lincon);
        }
    }

    {
        // v - vtot*z <= 0
        for (int i = 0; i < N; i++)
        {
            std::vector<VariablePtr> cvars;
            cvars.push_back(vars.at(iv+i)); // v(i)
            cvars.push_back(vars.at(iz+i)); // z(i)

            DenseMatrix A = DenseMatrix::Ones(1,2);
            A(0,1) = -vtot;
            DenseVector b = DenseVector::Zero(1);

            ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, false);

            cs->add(lincon);
            cs_approx->add(lincon);
        }
    }

    {
        // x - z <= 0
        for (int i = 0; i < N; i++)
        {
            std::vector<VariablePtr> cvars;
            cvars.push_back(vars.at(ix+i)); // x(i)
            cvars.push_back(vars.at(iz+i)); // z(i)

            DenseMatrix A = DenseMatrix::Ones(1,2);
            A(0,1) = -1;
            DenseVector b = DenseVector::Zero(1);

            ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, false);

            cs->add(lincon);
            cs_approx->add(lincon);
        }
    }

    {
        // w - wmax*z <= 0
        for (int i = 0; i < N; i++)
        {
            std::vector<VariablePtr> cvars;
            cvars.push_back(vars.at(iw+i)); // w(i)
            cvars.push_back(vars.at(iz+i)); // z(i)

            DenseMatrix A = DenseMatrix::Ones(1,2);
            A(0,1) = -wmax;
            DenseVector b = DenseVector::Zero(1);

            ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, false);

            cs->add(lincon);
            cs_approx->add(lincon);
        }
    }

    {
        // np - 3*z <= 0
        for (int i = 0; i < N; i++)
        {
            std::vector<VariablePtr> cvars;
            cvars.push_back(vars.at(inp+i)); // np(i)
            cvars.push_back(vars.at(iz+i)); // z(i)

            DenseMatrix A = DenseMatrix::Ones(1,2);
            A(0,1) = -3;
            DenseVector b = DenseVector::Zero(1);

            ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, false);

            cs->add(lincon);
            cs_approx->add(lincon);
        }
    }

    {
        // ns - 3*z <= 0
        for (int i = 0; i < N; i++)
        {
            std::vector<VariablePtr> cvars;
            cvars.push_back(vars.at(ins+i)); // ns(i)
            cvars.push_back(vars.at(iz+i)); // z(i)

            DenseMatrix A = DenseMatrix::Ones(1,2);
            A(0,1) = -3;
            DenseVector b = DenseVector::Zero(1);

            ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, false);

            cs->add(lincon);
            cs_approx->add(lincon);
        }
    }

    {
        // l1 = np*ns*z
        for (int i = 0; i < N; i++)
        {
            std::vector<VariablePtr> cvars = {
                vars.at(inp+i),
                vars.at(ins+i),
                vars.at(iz+i),
                vars.at(il1+i)
            };

            std::vector<double> clb, cub;
            for (unsigned int j = 0; j < cvars.size()-1; j++)
            {
                clb.push_back(cvars.at(j)->getLowerBound());
                cub.push_back(cvars.at(j)->getUpperBound());
            }

            std::vector<unsigned int> deg = {1,1,1};

            DenseVector c = DenseVector::Zero(8);
            c(7) = 1; // trilinear term

            DenseMatrix T = getTransformationMatrix(deg, clb, cub);
            DenseMatrix coeffs = T*c;

            std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, clb, cub);

            BSpline bs(coeffs.transpose(), knots, deg);

            ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

            cs->add(cbs);
            cs_approx->add(cbs);
        }
    }

    {
        // l2 = p*l1
        for (int i = 0; i < N; i++)
        {
            std::vector<VariablePtr> cvars = {
                vars.at(ip+i),
                vars.at(il1+i),
                vars.at(il2+i)
            };

            std::vector<double> clb, cub;
            for (unsigned int j = 0; j < cvars.size()-1; j++)
            {
                clb.push_back(cvars.at(j)->getLowerBound());
                cub.push_back(cvars.at(j)->getUpperBound());
            }

            std::vector<unsigned int> deg = {1,1};

            DenseVector c = DenseVector::Zero(4);
            c(3) = 1; // bilinear term

            DenseMatrix T = getTransformationMatrix(deg, clb, cub);
            DenseMatrix coeffs = T*c;

            std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, clb, cub);

            BSpline bs(coeffs.transpose(), knots, deg);

            ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

            cs->add(cbs);
            cs_approx->add(cbs);
        }
    }

    {
        // l3 = alpha*(w/wmax)^3 + beta*(w/wmax)^2*v + gamma*(w/wmax)*v^2
        for (int i = 0; i < N; i++)
        {
            std::vector<VariablePtr> cvars = {
                vars.at(iw+i),
                vars.at(iv+i),
                vars.at(il3+i)
            };

            std::vector<double> clb, cub;
            for (unsigned int j = 0; j < cvars.size()-1; j++)
            {
                clb.push_back(cvars.at(j)->getLowerBound());
                cub.push_back(cvars.at(j)->getUpperBound());
            }

            std::vector<unsigned int> deg = {3,2};

            DenseVector c = DenseVector::Zero(12);
            c(5) = gamma.at(i)/wmax; // w*v^2
            c(7) = beta.at(i)/(wmax*wmax); // w^2*v
            c(9) = alpha.at(i)/(wmax*wmax*wmax); // w^3

            DenseMatrix T = getTransformationMatrix(deg, clb, cub);
            DenseMatrix coeffs = T*c;

            std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, clb, cub);

            BSpline bs(coeffs.transpose(), knots, deg);

            ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

            cs->add(cbs);
        }

        // Black-box approach
        for (int i = 0; i < N; i++)
        {
            std::vector<VariablePtr> cvars = {
                vars.at(iw+i),
                vars.at(iv+i),
                vars.at(il3+i)
            };

            std::vector<double> clb, cub;
            for (unsigned int j = 0; j < cvars.size()-1; j++)
            {
                clb.push_back(cvars.at(j)->getLowerBound());
                cub.push_back(cvars.at(j)->getUpperBound());
            }

            auto f_pow = [wmax](double v, double w, double alpha, double beta, double gamma)
            {
                return alpha*std::pow(w/wmax, 3) + beta*std::pow(w/wmax, 2)*v + gamma*(w/wmax)*v*v;
            };

            DataTable samples;

//            auto w = linspace(clb.at(0), cub.at(0), 10);
//            auto v = linspace(clb.at(1), cub.at(1), 10);
            auto w = linspace(clb.at(0), cub.at(0), grid);
            auto v = linspace(clb.at(1), cub.at(1), grid);
            for (auto wi : w)
            {
                for (auto vi : v)
                {
                    std::vector<double> x = {wi,vi};
                    double y = f_pow(vi, wi, alpha.at(i), beta.at(i), gamma.at(i));
                    //double dp = f_dp(vi, wi, a.at(i), b.at(i), c.at(i));

                    samples.addSample(x, y);
                }
            }

            // Create B-spline
            BSpline bs(samples, BSplineType::CUBIC_FREE);

            ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

            cs_approx->add(cbs);
        }
    }

    {
        // l4 = a*(w/wmax)^2 + b*(w/wmax)*v + c*v^2
        for (int i = 0; i < N; i++)
        {
            std::vector<VariablePtr> cvars = {
                vars.at(iw+i),
                vars.at(iv+i),
                vars.at(il4+i)
            };

            std::vector<double> clb, cub;
            for (unsigned int j = 0; j < cvars.size()-1; j++)
            {
                clb.push_back(cvars.at(j)->getLowerBound());
                cub.push_back(cvars.at(j)->getUpperBound());
            }

            std::vector<unsigned int> deg = {2,2};

            DenseVector c = DenseVector::Zero(9);
            c(2) = cp.at(i); // v^2
            c(4) = bp.at(i)/wmax; // w*v
            c(6) = ap.at(i)/(wmax*wmax); // w^2

            DenseMatrix T = getTransformationMatrix(deg, clb, cub);
            DenseMatrix coeffs = T*c;

            std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, clb, cub);

            BSpline bs(coeffs.transpose(), knots, deg);

            ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

            cs->add(cbs);
        }

        // Black-box approach
        for (int i = 0; i < N; i++)
        {
            std::vector<VariablePtr> cvars = {
                vars.at(iw+i),
                vars.at(iv+i),
                vars.at(il4+i)
            };

            std::vector<double> clb, cub;
            for (unsigned int j = 0; j < cvars.size()-1; j++)
            {
                clb.push_back(cvars.at(j)->getLowerBound());
                cub.push_back(cvars.at(j)->getUpperBound());
            }

            auto f_dp = [wmax](double v, double w, double a, double b, double c)
            {
                return a*std::pow(w/wmax, 2) + b*(w/wmax)*v + c*v*v;
            };

            DataTable samples;

//            auto w = linspace(clb.at(0), cub.at(0), 10);
//            auto v = linspace(clb.at(1), cub.at(1), 10);
            auto w = linspace(clb.at(0), cub.at(0), grid);
            auto v = linspace(clb.at(1), cub.at(1), grid);
            for (auto wi : w)
            {
                for (auto vi : v)
                {
                    std::vector<double> x = {wi,vi};
                    //double y = f_pow(vi, wi, alpha.at(i), beta.at(i), gamma.at(i));
                    double y = f_dp(vi, wi, ap.at(i), bp.at(i), cp.at(i));

                    samples.addSample(x, y);
                }
            }

            // Create B-spline
            BSpline bs(samples, BSplineType::CUBIC_FREE);

            ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

            cs_approx->add(cbs);
        }
    }

    {
        // l5 = v*np
        for (int i = 0; i < N; i++)
        {
            std::vector<VariablePtr> cvars = {
                vars.at(iv+i),
                vars.at(inp+i),
                vars.at(il5+i)
            };

            std::vector<double> clb, cub;
            for (unsigned int j = 0; j < cvars.size()-1; j++)
            {
                clb.push_back(cvars.at(j)->getLowerBound());
                cub.push_back(cvars.at(j)->getUpperBound());
            }

            std::vector<unsigned int> deg = {1,1};

            DenseVector c = DenseVector::Zero(4);
            c(3) = 1; // v*np

            DenseMatrix T = getTransformationMatrix(deg, clb, cub);
            DenseMatrix coeffs = T*c;

            std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, clb, cub);

            BSpline bs(coeffs.transpose(), knots, deg);

            ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

            cs->add(cbs);
            cs_approx->add(cbs);
        }
    }

    {
        // l6 = dp*ns
        for (int i = 0; i < N; i++)
        {
            std::vector<VariablePtr> cvars = {
                vars.at(idp+i),
                vars.at(ins+i),
                vars.at(il6+i)
            };

            std::vector<double> clb, cub;
            for (unsigned int j = 0; j < cvars.size()-1; j++)
            {
                clb.push_back(cvars.at(j)->getLowerBound());
                cub.push_back(cvars.at(j)->getUpperBound());
            }

            std::vector<unsigned int> deg = {1,1};

            DenseVector c = DenseVector::Zero(4);
            c(3) = 1; // dp*ns

            DenseMatrix T = getTransformationMatrix(deg, clb, cub);
            DenseMatrix coeffs = T*c;

            std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, clb, cub);

            BSpline bs(coeffs.transpose(), knots, deg);

            ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

            cs->add(cbs);
            cs_approx->add(cbs);
        }
    }

    //BB::BranchAndBound solver(cs_approx);
    BB::BranchAndBound solver(cs);
    auto res = solver.optimize();
    cout << res << endl;

    auto x = res.primalVariables;
    int count = 0;
    for (auto xi : x)
    {
        cout << std::setprecision(10) << count++ << ": " << xi << endl;
    }

    {
        // Get variables (but may be in different order than vars)
        auto cs_approx_vars = cs_approx->getVariables();

        int counter = 0;
        for (auto var : cs_approx_vars)
        {
            // Fix integer variables
            if (var->getType() == VariableType::BINARY
                    || var->getType() == VariableType::INTEGER)
            {
                auto xval = x.at(counter);
                var->setValue(xval);
                var->setLowerBound(xval);
                var->setUpperBound(xval);
            }
            counter++;
        }
    }

    // Set initial values and type to continuous
    {
        for (unsigned int i = 0; i < vars.size(); i++)
        {
            vars.at(i)->setValue(x.at(i));
            vars.at(i)->setType(VariableType::CONTINUOUS);
        }
    }

    // Do local search
    SolverIpopt nlp(cs);
    auto res2 = nlp.optimize();

    cout << res2 << endl;
    auto x2 = res2.primalVariables;
    int count2 = 0;
    for (auto xi : x2)
    {
        cout << std::setprecision(10) << count2++ << ": " << xi << endl;
    }
}

/*
 * Optimal control problem 1 from Liao and Schoemaker (1992)
 */
void optControl1()
{
    const int N = 10;   // horizon
    const int n = 2;  // states
    const int m = 1;   // inputs
    const double mu = 1.0/2.0; // nonlinearity parameter

    DenseMatrix A = DenseMatrix::Zero(n,n);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i == j)
                A(i,j) = 0.5;
            else if (j == i+1)
                A(i,j) = 0.25;
            else if (j == i-1)
                A(i,j) = -0.25;
        }
    }

    DenseMatrix B = DenseMatrix::Zero(n,m);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            double num = (double)(i-j);
            double den = (double)(n+m);
            B(i,j) = num/den;
        }
    }

    DenseMatrix C = DenseMatrix::Zero(n,m);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            double num = (double)(i+1+j+1); // (i+1) and (j+1) since 0-indexed
            double den = (double)(n+m);
            C(i,j) = mu*num/den;
        }
    }

    /*
     * Optimization variables:
     * N*n states
     * (N-1)*m inputs
     * N*n + (N-1)*m auxiliary variables
     */
    int numVars = N*n + (N-1)*m;
    int numObjAuxVars = N*n + (N-1)*m;
    int numConAuxVars = (N-1)*n*m;
    int totVars = numVars + numObjAuxVars + numConAuxVars;

    // Variables
    std::vector<double> costs(totVars, 0);
    for (int i = 0; i < numObjAuxVars; i++)
        costs.at(numVars+i) = 1;

    std::vector<VariablePtr> vars;
    for (int i = 0; i < totVars; i++)
    {
        auto var = std::make_shared<Variable>(costs.at(i));
        vars.push_back(var);
    }

    // Constraint set
    ConstraintSetPtr cs = std::make_shared<ConstraintSet>();

    {
        // x(k+1) = A x(k) + B u(k) + (x(k)^T C u(k))e
        for (int k = 0; k < N-1; k++)
        {
            std::vector<VariablePtr> cvars;
            for (int i = 0; i < n; i++)
                cvars.push_back(vars.at(n*(k+1)+i)); // x(k+1)
            for (int i = 0; i < n; i++)
                cvars.push_back(vars.at(n*k+i)); // x(k)
            for (int i = 0; i < m; i++)
                cvars.push_back(vars.at(N*n+k*m+i)); // u(k)
            for (int i = 0; i < n*m; i++)
                cvars.push_back(vars.at(numVars+numObjAuxVars+k*(n*m)+i)); // lambda(k)

            DenseMatrix I = DenseMatrix::Identity(n,n);
            DenseMatrix ones = DenseMatrix::Ones(n,n*m);

            DenseMatrix IABC = DenseMatrix::Zero(n,n+n+m+n*m);

            IABC.block(0,0,n,n) = I;
            IABC.block(0,n,n,n) = -A;
            IABC.block(0,n+n,n,m) = -B;
            IABC.block(0,n+n+m,n,n*m) = -ones;

            DenseVector b = DenseVector::Zero(n);

            ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, IABC, b, true);

            cs->add(lincon);
        }
    }

    {
        // x(1) = 0
        std::vector<VariablePtr> cvars;
        for (int i = 0; i < n; i++)
            cvars.push_back(vars.at(i)); // x(1)

        DenseMatrix I = DenseMatrix::Identity(n,n);
        DenseVector b = DenseVector::Zero(n);

        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, I, b, true);

        cs->add(lincon);
    }

    {
        // Nonlinear objective terms
        for (int i = 0; i < numObjAuxVars; i++)
        {
            std::vector<VariablePtr> cvars = {vars.at(i), vars.at(numVars+i)};

            std::vector<double> lb = {-10};
            std::vector<double> ub = {10};

            std::vector<unsigned int> deg = {4};

            double a = 0.25;
            if (i >= N*n)
                a = 0.5; // u terms
            DenseVector c = DenseVector::Zero(5);
            c(0) = a*a*a*a;
            c(1) = 4*a*a*a;
            c(2) = 6*a*a;
            c(3) = 4*a;
            c(4) = 1;

            DenseMatrix T = getTransformationMatrix(deg, lb, ub);
            DenseMatrix coeffs = T*c;

            std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, lb, ub);

            BSpline bs(coeffs.transpose(), knots, deg);

            ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

            cs->add(cbs);

            for (int j = 0; j < 100; j++)
            {
                double frac = j/100.0;
                double x = lb.at(0) + frac*(ub.at(0)-lb.at(0));
                DenseVector xv(1); xv(0) = x;

                auto y = bs.eval(xv);
                auto yreal = std::pow(x+a,4);

                assert(assertNear(y, yreal, 1e-03));
             }
        }

    }

    {
        // Nonlinear constraint terms
        int counter = 0;

        for (int k = 0; k < N-1; k++)
        {
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < m; j++)
                {
                    // lambda = x(k,i)*C(i,j)*u(k,j)

                    std::vector<VariablePtr> cvars = {
                        vars.at(k*n+i),
                        vars.at(N*n+k*m+j),
                        vars.at(numVars+numObjAuxVars+counter)
                    };

                    std::vector<double> lb = {-10,-10};
                    std::vector<double> ub = {10,10};

                    std::vector<unsigned int> deg = {1,1};

                    DenseVector c = DenseVector::Zero(4);
                    c(3) = C(i,j);

                    DenseMatrix T = getTransformationMatrix(deg, lb, ub);
                    DenseMatrix coeffs = T*c;

                    std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, lb, ub);

                    BSpline bs(coeffs.transpose(), knots, deg);

                    ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

                    cs->add(cbs);

                    counter++;
                }
            }
        }
        assert(counter == numConAuxVars);
    }

//    std::vector<double> lb,ub;
//    cs->getDomainBounds(lb,ub);
//    for (int i = 0; i < numVars; i++)
//    {
//        lb.at(i) = -1000;
//        ub.at(i) = 1000;
//    }
//    cs->setDomainBounds(lb,ub);

    BB::BranchAndBound solver(cs);
    auto res = solver.optimize();

    cout << res << endl;

}

/*
 * Optimal control problem 2 from Liao and Schoemaker (1992)
 */
void optControl2()
{
    const int n = 1;    // states
    const int m = 1;    // inputs
    const int N = 10;   // horizon

    /*
     * Optimization variables:
     * N*n states
     * (N-1)*m inputs
     * N*n + (N-1)*m auxiliary variables
     */
    int numStates = n*N;
    int numInputs = m*(N-1);
    int numAux = N;
    int numVars = numStates + numInputs + numAux;

    // Variables
    std::vector<double> costs(numVars, 0);
    costs.at(numStates+numInputs) = 1;

    std::vector<VariablePtr> vars;
    for (int i = 0; i < numVars; i++)
    {
        auto var = std::make_shared<Variable>(costs.at(i));
        vars.push_back(var);
    }

    // Initial condition x(1) = 1
    vars.at(0)->setLowerBound(0);
    vars.at(0)->setUpperBound(0);

    // Constraint set
    ConstraintSetPtr cs = std::make_shared<ConstraintSet>();

    // (1/N) sum_{i = 1,...,N-1} (x(k)^2 + u(k)^2)
    {
        std::vector<VariablePtr> cvars;
        for (int i = 0; i < numStates-1; i++)
            cvars.push_back(vars.at(i)); // x(k)
        for (int i = numStates; i < numStates+numInputs; i++)
            cvars.push_back(vars.at(i)); // u(k)
        cvars.push_back(vars.at(numStates+numInputs)); // lambda(1)

        int nv = 2*(N-1)+1;
        DenseMatrix A = (1.0/N)*DenseMatrix::Identity(nv,nv);
        A(nv-1,nv-1) = 0;

        DenseMatrix b = DenseMatrix::Zero(nv,1);
        b(nv-1,0) = -1;

        ConstraintPtr quadcon = std::make_shared<ConstraintQuadratic>(cvars, A, b, 0, -INF, 0);

        assert(quadcon->isConstraintConvex());

        cs->add(quadcon);
    }

    // x(k+1) = x(k) + (1/N)(lambda(k) - u(k)), for k = 1,...,N-1
    {
        for (int k = 0; k < N-1; k++)
        {
            std::vector<VariablePtr> cvars;
            cvars.push_back(vars.at(k+1)); // x(k+1)
            cvars.push_back(vars.at(k)); // x(k)
            cvars.push_back(vars.at(numStates+numInputs+1+k)); // lambda(k)
            cvars.push_back(vars.at(numStates+k)); // u(k)

            double invN = (1.0/N);
            DenseMatrix A = DenseMatrix::Zero(1,4);
            A(0,0) = 1;     // x(k+1)
            A(0,1) = -1;    // x(k)
            A(0,2) = -invN; // lambda(k)
            A(0,3) = invN;  // u(k)

            DenseVector b = DenseVector::Zero(1);

            ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, true);

            cs->add(lincon);
        }
    }

    // lambda(k) = x(k)^2, for k = 1,...,N-1
    {
        // Nonlinear terms
        for (int i = 0; i < N-1; i++)
        {
            std::vector<VariablePtr> cvars = {
                vars.at(i),
                vars.at(numStates+numInputs+1+i)
            };

            std::vector<double> lb = {-100};
            std::vector<double> ub = {100};

            std::vector<unsigned int> deg = {2};

            DenseVector c = DenseVector::Zero(3);
            c(2) = 1;

            DenseMatrix T = getTransformationMatrix(deg, lb, ub);
            DenseMatrix coeffs = T*c;

            std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, lb, ub);

            BSpline bs(coeffs.transpose(), knots, deg);

            ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

            cs->add(cbs);

            for (int j = 0; j < 100; j++)
            {
                double frac = j/100.0;
                double x = lb.at(0) + frac*(ub.at(0)-lb.at(0));
                DenseVector xv(1); xv(0) = x;

                auto y = bs.eval(xv);
                auto yreal = x*x;

                if (!(std::abs(y-yreal) < 1e-10))
                {
                    cout << std::setprecision(10) << yreal << endl;
                    cout << std::setprecision(10) << y << endl;
                }
                assert(std::abs(y-yreal) < 1e-10);
                //assert(isApprox(y, yreal, 1e-03));
             }
        }
    }

    BB::BranchAndBound solver(cs);
    auto res = solver.optimize();

    cout << res << endl;
}

/*
 * Optimal control problem 2 from Liao and Schoemaker (1992)
 */
void optControl3()
{
    const int n = 1;    // states
    const int m = 1;    // inputs
    const int N = 20;   // horizon

    /*
     * Optimization variables:
     * N*n states
     * (N-1)*m inputs
     * N*n + (N-1)*m auxiliary variables
     */
    int numStates = n*N;
    int numInputs = m*(N-1);
    int numAux = 1 + 2*(N-1);
    int numVars = numStates + numInputs + numAux;

    // Variables
    std::vector<double> costs(numVars, 0);
    costs.at(numStates+numInputs) = 1;

    std::vector<VariablePtr> vars;
    for (int i = 0; i < numVars; i++)
    {
        auto var = std::make_shared<Variable>(costs.at(i));
        vars.push_back(var);
    }

    // Initial condition x(1) = 1
    vars.at(0)->setLowerBound(1);
    vars.at(0)->setUpperBound(1);

    // Constraint set
    ConstraintSetPtr cs = std::make_shared<ConstraintSet>();

    // (1/N) sum_{i = 1,...,N-1} (lambda1(k) + lambda2(k)) = lambda0
    {
        std::vector<VariablePtr> cvars;
        for (int i = 0; i < numAux; i++)
            cvars.push_back(vars.at(numStates+numInputs+i)); // lambda0, lambda1(k) and lambda2(k)

        int nv = 2*(N-1)+1;
        DenseMatrix A = (1.0/N)*DenseMatrix::Ones(1,nv);
        A(0,0) = -1; // lambda0

        DenseVector b = DenseVector::Zero(1,1);

        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, false);

        //assert(quadcon->isConstraintConvex());

        cs->add(lincon);
    }

    // x(k+1) = x(k) + (1/N)(lambda1(k) - u(k)), for k = 1,...,N-1
    {
        for (int k = 0; k < N-1; k++)
        {
            std::vector<VariablePtr> cvars;
            cvars.push_back(vars.at(k+1)); // x(k+1)
            cvars.push_back(vars.at(k)); // x(k)
            cvars.push_back(vars.at(numStates+numInputs+1+k)); // lambda(k)
            cvars.push_back(vars.at(numStates+k)); // u(k)

            double invN = (1.0/N);
            DenseMatrix A = DenseMatrix::Zero(1,4);
            A(0,0) = 1;     // x(k+1)
            A(0,1) = -1;    // x(k)
            A(0,2) = -invN; // lambda1(k)
            A(0,3) = invN;  // u(k)

            DenseVector b = DenseVector::Zero(1);

            ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, true);

            cs->add(lincon);
        }
    }

    // lambda1(k) = x(k)^2, for k = 1,...,N-1
    {
        // Nonlinear terms
        for (int i = 0; i < N-1; i++)
        {
            std::vector<VariablePtr> cvars = {vars.at(i), vars.at(numStates+numInputs+1+i)};

            std::vector<double> lb = {-100};
            std::vector<double> ub = {100};

            std::vector<unsigned int> deg = {2};

            DenseVector c = DenseVector::Zero(3);
            c(2) = 1;

            DenseMatrix T = getTransformationMatrix(deg, lb, ub);
            DenseMatrix coeffs = T*c;

            std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, lb, ub);

            BSpline bs(coeffs.transpose(), knots, deg);

            ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

            cs->add(cbs);

            for (int j = 0; j < 100; j++)
            {
                double frac = j/100.0;
                double x = lb.at(0) + frac*(ub.at(0)-lb.at(0));
                DenseVector xv(1); xv(0) = x;

                auto y = bs.eval(xv);
                auto yreal = x*x;

                if (!(std::abs(y-yreal) < 1e-10))
                {
                    cout << std::setprecision(10) << yreal << endl;
                    cout << std::setprecision(10) << y << endl;
                }
                assert(std::abs(y-yreal) < 1e-10);
                //assert(isApprox(y, yreal, 1e-03));
             }
        }
    }

    // lambda2(k) = u(k)^2, for k = 1,...,N-1
    {
        // Nonlinear terms
        for (int i = 0; i < N-1; i++)
        {
            std::vector<VariablePtr> cvars = {
                vars.at(numStates+i),
                vars.at(numStates+numInputs+N+i)
            };

            std::vector<double> lb = {-100};
            std::vector<double> ub = {100};

            std::vector<unsigned int> deg = {2};

            DenseVector c = DenseVector::Zero(3);
            c(2) = 1;

            DenseMatrix T = getTransformationMatrix(deg, lb, ub);
            DenseMatrix coeffs = T*c;

            std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, lb, ub);

            BSpline bs(coeffs.transpose(), knots, deg);

            ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

            cs->add(cbs);

            for (int j = 0; j < 100; j++)
            {
                double frac = j/100.0;
                double x = lb.at(0) + frac*(ub.at(0)-lb.at(0));
                DenseVector xv(1); xv(0) = x;

                auto y = bs.eval(xv);
                auto yreal = x*x;

                if (!(std::abs(y-yreal) < 1e-10))
                {
                    cout << std::setprecision(10) << yreal << endl;
                    cout << std::setprecision(10) << y << endl;
                }
                assert(std::abs(y-yreal) < 1e-10);
                //assert(isApprox(y, yreal, 1e-03));
             }
        }
    }

    //SolverIpopt solver(cs);
    BB::BranchAndBound solver(cs);
    auto res = solver.optimize();

    cout << res << endl;
}

void Ptest()
{
    // Test problem
    cout << "\n\nSolving problem Ptest..." << endl;

    int dim = 2+3;

    // x1,x2,l1,l2
    std::vector<double> lb = {-1,-1,-INF,-INF,-INF};
    std::vector<double> ub = {2,2,INF,INF,INF};
    std::vector<double> costs = {0, 0, -1, 0, 0};

    std::vector<VariablePtr> vars;
    for (int i = 0; i < dim; i++)
    {
        auto var = std::make_shared<Variable>(costs.at(i), lb.at(i), ub.at(i));
        vars.push_back(var);
    }

    ConstraintSetPtr cs = std::make_shared<ConstraintSet>();

    { // w2 = x1^2

        std::vector<VariablePtr> cvars = {vars.at(0), vars.at(3)};

        std::vector<double> thislb = {lb.at(0)};
        std::vector<double> thisub = {ub.at(0)};

        std::vector<unsigned int> deg = {2};

        DenseVector c(3);
        c.setZero();
        c(2) = 1;

        DenseMatrix T = getTransformationMatrix(deg, thislb, thisub);
        DenseMatrix coeffs = T*c;

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

        BSpline bs(coeffs.transpose(), knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    { // w3 = x2^2

        std::vector<VariablePtr> cvars = {vars.at(1), vars.at(4)};

        std::vector<double> thislb = {lb.at(1)};
        std::vector<double> thisub = {ub.at(1)};

        std::vector<unsigned int> deg = {2};

        DenseVector c(3);
        c.setZero();
        c(2) = 1;

        DenseMatrix T = getTransformationMatrix(deg, thislb, thisub);
        DenseMatrix coeffs = T*c;

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

        BSpline bs(coeffs.transpose(), knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    { // w1 = w2*w3

        std::vector<VariablePtr> cvars = {vars.at(3), vars.at(4), vars.at(2)};

        std::vector<double> thislb = {lb.at(3), lb.at(4)};
        std::vector<double> thisub = {ub.at(3), ub.at(4)};

        std::vector<unsigned int> deg = {1,1};

        DenseVector c(4);
        c.setZero();
        c(3) = 1;

        DenseMatrix T = getTransformationMatrix(deg, thislb, thisub);
        DenseMatrix coeffs = T*c;

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

        BSpline bs(coeffs.transpose(), knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }


    //std::vector<int> vt= {INTEGER,INTEGER,CONTINUOUS,CONTINUOUS}; // Integer problem

    BB::BranchAndBound bnb(cs);
    Timer timer;
    timer.start();
    bnb.optimize();
    timer.stop();
    cout << "Time: " << timer.getMilliSeconds() << " (ms)" << endl;
}



double sixHumpCamelFunction(DenseVector x)
{
    assert(x.rows() == 2);
    return (4 - 2.1*x(0)*x(0) + (1/3.)*x(0)*x(0)*x(0)*x(0))*x(0)*x(0) + x(0)*x(1) + (-4 + 4*x(1)*x(1))*x(1)*x(1);
}

void polynomialOptimization()
{
    std::vector<VariablePtr> vars =
    {
        std::make_shared<Variable>(0, -1, 1),
        std::make_shared<Variable>(0, -0.5, 2),
        std::make_shared<Variable>(1, -INF, INF)
    };

    std::vector<unsigned int> degrees = {3,2};
    unsigned int dim = 2;
    int terms = 2;

    std::vector<double> lb = {-1,-0.5};
    std::vector<double> ub = {1,2};

    DenseVector c(terms); c.setZero();
    c(0) = -1;
    c(1) = 1;

    DenseMatrix E(terms,dim); E.setZero();
    E(0,0) = 1; E(0,1) = 2;
    E(1,0) = 3; E(1,1) = 2;

    cout << E << endl;

    // Coefficients
    DenseMatrix coeffs = getBSplineBasisCoefficients(c, E, lb, ub);

    cout << "Coefficients: " << endl;
    cout << coeffs << endl;

    // Knots
    std::vector< std::vector<double> > knots;
    for (unsigned int j = 0; j < dim; j++)
    {
        std::vector<double> ks;
        for (unsigned int i = 0; i < degrees.at(j)+1; i++)
        {
            ks.push_back(lb.at(j));
        }
        for (unsigned int i = 0; i < degrees.at(j)+1; i++)
        {
            ks.push_back(ub.at(j));
        }

        knots.push_back(ks);
        cout << ks;
    }

    BSpline bs(coeffs.transpose(), knots, degrees);

    ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(vars, bs, true);

    ConstraintSetPtr cs = std::make_shared<ConstraintSet>();

    cs->add(cbs);

    BB::BranchAndBound bnb(cs);
    Timer timer;
    timer.start();
    bnb.optimize();
    timer.stop();
    cout << "Time: " << timer.getMilliSeconds() << endl;

//    cout << "Starting test..." << endl;
//    for (double x = a; x <= b; x=x+0.001)
//    {
//        DenseVector xv(1);
//        xv(0) = x;
//        DenseVector yv = bs.evaluate(xv);

//        double yv1 = yv(0);
//        double yv2 = 1 - x + x*x - x*x*x + 2*x*x*x*x;

//        if (std::abs(yv1-yv2) > 1e-16)
//            cout << yv1 - yv2 << endl;

//    }

}

void sixHumpCamelBackPoly()
{
    // Six hump camel

    double d = 3; // Hyperrectangle "radius"
    std::vector<double> lb = {-d,-d};
    std::vector<double> ub = {d,d};

    // Variables
    std::vector<VariablePtr> vars =
    {
        std::make_shared<Variable>(0, -d, d),
        std::make_shared<Variable>(0, -d, d),
        std::make_shared<Variable>(1, -INF, INF)
    };

    // Degrees and coefficients
    std::vector<unsigned int> deg = {6,4};

    DenseVector c(35);
    c.setZero();
    c(10) = 4;
    c(20) = -2.1;
    c(30) = 1/3.;
    c(6) = 1;
    c(2) = -4;
    c(4) = 4;

    DenseMatrix T = getTransformationMatrix(deg,lb,ub);

    DenseMatrix coeffs = T*c;

    std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, lb, ub);

//    printVector(knots);

    BSpline bs(coeffs.transpose(), knots, deg);

    ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(vars, bs, true);

    ConstraintSetPtr cs = std::make_shared<ConstraintSet>();
    cs->add(cbs);

    cout << "Num control points: ";
    DenseMatrix cpts = bs.getControlPoints();
//    cout << cpts << endl;
    cout << cpts.cols() << endl;
    //exit(1);

    // Relaxed problem
    ConstraintPtr csr = cs->getConvexRelaxation();

    SolverGurobi solver(csr);
    SolverResult result = solver.optimize();
    cout << "frelaxed below: " << endl;
    cout << result << endl;
    cout << "fopt: " << -1.0316 << endl;
    cout << "gap: " << -1.0316 - result.objectiveValue << endl;

//    // Solve problem
//    std::vector<int> vt(dim+1,CONTINUOUS);
//    std::vector<int> bv = {0,1};
//    std::vector<double> z0 = {0,0,0};

//    BB::BranchAndBound bnb(obj, cs, z0, vt, bv);
//    Timer timer;
//    timer.start();
//    bnb.optimize();
//    timer.stop();
//    cout << "Time: " << timer.getMilliSeconds() << endl;
}

void sixHumpCamelBackPoly2()
{
    // Six hump camelback
    // Only the nonconvex terms are relaxed with a B-spline
    // The convex terms are kept
    int dim = 2;
    int aux = 4; // num aux vars
    std::vector<unsigned int> deg = {4,2};

    double d = 10;

    std::vector<double> costs = {0, 0, 1, 4, 1./3., 4};
    std::vector<double> clb = {-d, -d, -INF, -INF, -INF, -INF};
    std::vector<double> cub = {d, d, INF, INF, INF, INF};

    std::vector<VariablePtr> vars;
    for (int i = 0; i < dim+aux; i++)
    {
        auto var = std::make_shared<Variable>(costs.at(i), clb.at(i), cub.at(i));
        vars.push_back(var);
    }

    ConstraintSetPtr cs = std::make_shared<ConstraintSet>();

    {
        std::vector<VariablePtr> cvars = {
            vars.at(0),
            vars.at(1),
            vars.at(2)
        };

        std::vector<double> lb = {-d,-d};
        std::vector<double> ub = {d,d};

        DenseVector c(15);
        c.setZero();
        c(2) = -4;
        c(4) = 1;
        c(12) = -2.1;

        DenseMatrix T = getTransformationMatrix(deg,lb,ub);

        DenseMatrix coeffs = T*c;

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, lb, ub);

        BSpline bs(coeffs.transpose(), knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        // Add B-spline
        cs->add(cbs);
    }

    // Add convex terms
    {
        std::vector<VariablePtr> cvars = {
            vars.at(0),
            vars.at(3)
        };

        DenseMatrix E;
        E.setZero(2,2);
        E(0,0) = 2;
        E(1,1) = 1;

        DenseVector c(2);
        c(0) = 1;
        c(1) = -1;

        ConstraintPtr poly = std::make_shared<ConstraintPolynomial>(cvars, c, E, false);

        if (poly->isConstraintConvex())
            cout << "Convex!" << endl;
        else
            cout << "Nonconvex!" << endl;

        cs->add(poly);
    }

    {
        std::vector<VariablePtr> cvars = {vars.at(0), vars.at(4)};

        DenseMatrix E;
        E.setZero(2,2);
        E(0,0) = 6;
        E(1,1) = 1;

        DenseVector c(2);
        c(0) = 1;
        c(1) = -1;

        ConstraintPtr poly = std::make_shared<ConstraintPolynomial>(cvars, c, E, false);

        if (poly->isConstraintConvex())
            cout << "Convex!" << endl;
        else
            cout << "Nonconvex!" << endl;

        cs->add(poly);
    }

    {
        std::vector<VariablePtr> cvars = {vars.at(1), vars.at(5)};
        DenseMatrix E;
        E.setZero(2,2);
        E(0,0) = 4;
        E(1,1) = 1;

        DenseVector c(2);
        c(0) = 1;
        c(1) = -1;

        ConstraintPtr poly = std::make_shared<ConstraintPolynomial>(cvars, c, E, false);

        if (poly->isConstraintConvex())
            cout << "Convex!" << endl;
        else
            cout << "Nonconvex!" << endl;

        cs->add(poly);
    }

    // Relaxed problem
    ConstraintPtr csr = cs->getConvexRelaxation();

    // Solve relaxed problem
    SolverIpopt solver(csr);
    SolverResult result = solver.optimize();
    cout << "frelaxed below:" << endl;
    cout << result << endl;
    cout << "fopt: " << -1.0316 << endl;
    cout << "gap: " << -1.0316 - result.objectiveValue << endl;
}

//void bilinearConstraintBound(std::vector<double> lb, std::vector<double> ub, double &clb, double &cub)
//{
//    assert(lb.size() == 2);
//    assert(lb.size() == ub.size());

//    double l1 = lb.at(0);
//    double l2 = lb.at(1);
//    double u1 = ub.at(0);
//    double u2 = ub.at(1);

//    clb = l1*l2;
//    clb = std::min(clb, l1*u2);
//    clb = std::min(clb, u1*l2);
//    clb = std::min(clb, u1*u2);

//    cub = l1*l2;
//    cub = std::max(cub, l1*u2);
//    cub = std::max(cub, u1*l2);
//    cub = std::max(cub, u1*u2);
//}

//ConstraintPtr bilinearConstraint(std::vector<double> lb, std::vector<double> ub)
//{
//    assert(lb.size() == 2);
//    assert(lb.size() == ub.size());

//    DenseMatrix Ac(4,3); Ac.setZero();
//    Ac(0,0) = lb.at(1);     Ac(0,1) = lb.at(0);     Ac(0,2) = -1;
//    Ac(1,0) = ub.at(1);     Ac(1,1) = ub.at(0);     Ac(1,2) = -1;
//    Ac(2,0) = -lb.at(1);    Ac(2,1) = -ub.at(0);    Ac(2,2) = 1;
//    Ac(3,0) = -ub.at(1);    Ac(3,1) = -lb.at(0);    Ac(3,2) = 1;

//    DenseVector bc(4); bc.setZero();
//    bc(0) = lb.at(0)*lb.at(1);
//    bc(1) = ub.at(0)*ub.at(1);
//    bc(2) = -ub.at(0)*lb.at(1);
//    bc(3) = -lb.at(0)*ub.at(1);

//    return std::make_shared<ConstraintLinear>(Ac, bc, false);
//}

//// x^2 <= y
//ConstraintQuadratic* quadraticIneqConstraint()
//{
//    DenseMatrix A;
//    A.setZero(2,2);
//    A(0,0) = 1;

//    DenseMatrix b;
//    b.setZero(2,1);
//    b(1,0) = -1;

//    return new ConstraintQuadratic(A,b,0,-INF,0);
//}

//ConstraintLinear* lidConstraint(double lb, double ub)
//{
//    double a = (ub*ub - lb*lb)/(ub - lb);

//    DenseMatrix A;
//    A.setZero(1,2);
//    A(0,0) = -a;
//    A(0,1) = 1;

//    DenseVector b;
//    b.setZero(1);
//    b(0) = -a*lb + lb*lb;

//    return new ConstraintLinear(A,b,false);
//}

//void sixHumpCamelBackPoly3()
//{
//    // Six-hump camelback relaxation using binary tree for function expansion (standard form) and McCormick bilinear relaxation
//    int dim = 8;
//    double d = 10; // hyperrectangle "radius"
//    std::vector<double> lb = {-d,-d};
//    std::vector<double> ub = {d,d};

//    // Domains of relaxation variables (McCormick)
//    std::vector<double> clb = {lb.at(0), lb.at(1), -INF, -INF, -INF, -INF, -INF, -INF};
//    std::vector<double> cub = {ub.at(0), ub.at(1), INF, INF, INF, INF, INF, INF};

//    // Constraint composite
//    ConstraintCompositePtr cs(new ConstraintComposite(dim, clb, cub));

//    // Add constraint w2 = w0^2
//    {
//        std::vector<int> v = {0,0,2};
//        std::vector<double> thislb = {clb.at(v.at(0)), clb.at(v.at(1))};
//        std::vector<double> thisub = {cub.at(v.at(0)), cub.at(v.at(1))};

//        // Update bounds
//        double auxlb, auxub;
//        bilinearConstraintBound(thislb, thisub, auxlb, auxub);
//        clb.at(v.at(2)) = auxlb;
//        cub.at(v.at(2)) = auxub;

//        // Add linearized bilinear constraint
//        ConstraintPtr lincon = bilinearConstraint(thislb, thisub);
//        cs->add(lincon, v);
//    }

//    // Add constraint w3 = w1^2
//    {
//        std::vector<int> v = {1,1,3};
//        std::vector<double> thislb = {clb.at(v.at(0)), clb.at(v.at(1))};
//        std::vector<double> thisub = {cub.at(v.at(0)), cub.at(v.at(1))};

//        // Update bounds
//        double auxlb, auxub;
//        bilinearConstraintBound(thislb, thisub, auxlb, auxub);
//        clb.at(v.at(2)) = auxlb;
//        cub.at(v.at(2)) = auxub;

//        // Add linear constraint
//        ConstraintPtr lincon = bilinearConstraint(thislb, thisub);
//        cs->add(lincon, v);
//    }

//    // Add constraint w4 = w2*w3
//    {
//        std::vector<int> v = {2,3,4};
//        std::vector<double> thislb = {clb.at(v.at(0)), clb.at(v.at(1))};
//        std::vector<double> thisub = {cub.at(v.at(0)), cub.at(v.at(1))};

//        // Update bounds
//        double auxlb, auxub;
//        bilinearConstraintBound(thislb, thisub, auxlb, auxub);
//        clb.at(v.at(2)) = auxlb;
//        cub.at(v.at(2)) = auxub;

//        // Add linear constraint
//        ConstraintPtr lincon = bilinearConstraint(thislb, thisub);
//        cs->add(lincon, v);
//    }

//    // Add constraint w5 = w0*w1
//    {
//        std::vector<int> v = {0,1,5};
//        std::vector<double> thislb = {clb.at(v.at(0)), clb.at(v.at(1))};
//        std::vector<double> thisub = {cub.at(v.at(0)), cub.at(v.at(1))};

//        // Update bounds
//        double auxlb, auxub;
//        bilinearConstraintBound(thislb, thisub, auxlb, auxub);
//        clb.at(v.at(2)) = auxlb;
//        cub.at(v.at(2)) = auxub;

//        // Add linear constraint
//        ConstraintPtr lincon = bilinearConstraint(thislb, thisub);
//        cs->add(lincon, v);
//    }

//    // Add constraint w6 = w1^2
//    {
//        std::vector<int> v = {1,1,6};
//        std::vector<double> thislb = {clb.at(v.at(0)), clb.at(v.at(1))};
//        std::vector<double> thisub = {cub.at(v.at(0)), cub.at(v.at(1))};

//        // Update bounds
//        double auxlb, auxub;
//        bilinearConstraintBound(thislb, thisub, auxlb, auxub);
//        clb.at(v.at(2)) = auxlb;
//        cub.at(v.at(2)) = auxub;

//        // Add linear constraint
//        ConstraintPtr lincon = bilinearConstraint(thislb, thisub);
//        cs->add(lincon, v);
//    }

//    // Add constraint w7 = w6^2
//    {
//        std::vector<int> v = {6,6,7};
//        std::vector<double> thislb = {clb.at(v.at(0)), clb.at(v.at(1))};
//        std::vector<double> thisub = {cub.at(v.at(0)), cub.at(v.at(1))};

//        // Update bounds
//        double auxlb, auxub;
//        bilinearConstraintBound(thislb, thisub, auxlb, auxub);
//        clb.at(v.at(2)) = auxlb;
//        cub.at(v.at(2)) = auxub;

//        // Add linear constraint
//        ConstraintPtr lincon = bilinearConstraint(thislb, thisub);
//        cs->add(lincon, v);
//    }

//    // Update domain bounds
//    cs->setDomainBounds(clb, cub);

//    DenseMatrix cobj(1,dim);
//    cobj.setZero();
//    cobj(0,2) = 4;
//    cobj(0,3) = -2.1;
//    cobj(0,4) = 1.0/3.0;
//    cobj(0,5) = 1;
//    cobj(0,6) = -4;
//    cobj(0,7) = 4;

//    ObjectivePtr obj(new ObjectiveLinear(cobj));

//    std::vector<double> z0(dim, 0);
//    ProblemPtr prob(new Problem(obj, cs));
//    SolverGurobi solver(prob, z0);
//    SolverResult result = solver.optimize();

//    cout << "frelaxed below:" << endl;
//    cout << result << endl;
//    cout << "fopt: " << -1.0316 << endl;
//    cout << "gap: " << -1.0316 - result.objectiveValue << endl;

//    std::vector<double> x = result.primalVariables;
//    //printVec(x);
//}

//void sixHumpCamelBackPoly4()
//{
//    // Six-hump camelback relaxation using binary tree for function expansion (standard form)
//    // The problem is relaxed using McCormick bilinear relaxation and the convex envelope of quadratic terms
//    int dim = 8;
//    double d = 3; // hyperrectangle "radius"
//    std::vector<double> lb = {-d,-d};
//    std::vector<double> ub = {d,d};

//    // Domains of relaxation variables (McCormick)
//    std::vector<double> clb = {lb.at(0), lb.at(1), -INF, -INF, -INF, -INF, -INF, -INF};
//    std::vector<double> cub = {ub.at(0), ub.at(1), INF, INF, INF, INF, INF, INF};

//    // Constraint composite
//    ConstraintCompositePtr cs(new ConstraintComposite(dim, clb, cub));

//    // Add constraint w2 = w0^2
//    {
//        std::vector<int> v = {0,2};

//        double lb = clb.at(0);
//        double ub = cub.at(0);

//        // Update bounds on w2
//        clb.at(2) = 0; // since lb < 0 < ub
//        cub.at(2) = std::max(ub*ub,lb*lb);

//        // Add linearized bilinear constraint
//        ConstraintPtr quadcon(quadraticIneqConstraint());
//        cs->add(quadcon, v);

//        ConstraintPtr lidcon(lidConstraint(lb,ub));
//        cs->add(lidcon, v);
//    }

//    // Add constraint w3 = w1^2
//    {
//        std::vector<int> v = {1,3};

//        double lb = clb.at(1);
//        double ub = cub.at(1);

//        // Update bounds on w2
//        clb.at(3) = 0; // since lb < 0 < ub
//        cub.at(3) = std::max(ub*ub,lb*lb);

//        // Add linearized bilinear constraint
//        ConstraintPtr quadcon(quadraticIneqConstraint());
//        cs->add(quadcon, v);

//        ConstraintPtr lidcon(lidConstraint(lb,ub));
//        cs->add(lidcon, v);
//    }

//    // Add constraint w4 = w2*w3
//    {
//        std::vector<int> v = {2,3,4};
//        std::vector<double> thislb = {clb.at(v.at(0)), clb.at(v.at(1))};
//        std::vector<double> thisub = {cub.at(v.at(0)), cub.at(v.at(1))};

//        // Update bounds
//        double auxlb, auxub;
//        bilinearConstraintBound(thislb, thisub, auxlb, auxub);
//        clb.at(v.at(2)) = auxlb;
//        cub.at(v.at(2)) = auxub;

//        // Add linear constraint
//        ConstraintPtr lincon = bilinearConstraint(thislb, thisub);
//        cs->add(lincon, v);
//    }

//    // Add constraint w5 = w0*w1
//    {
//        std::vector<int> v = {0,1,5};
//        std::vector<double> thislb = {clb.at(v.at(0)), clb.at(v.at(1))};
//        std::vector<double> thisub = {cub.at(v.at(0)), cub.at(v.at(1))};

//        // Update bounds
//        double auxlb, auxub;
//        bilinearConstraintBound(thislb, thisub, auxlb, auxub);
//        clb.at(v.at(2)) = auxlb;
//        cub.at(v.at(2)) = auxub;

//        // Add linear constraint
//        ConstraintPtr lincon = bilinearConstraint(thislb, thisub);
//        cs->add(lincon, v);
//    }

//    // Add constraint w6 = w1^2
//    {
//        std::vector<int> v = {1,6};

//        double lb = clb.at(1);
//        double ub = cub.at(1);

//        // Update bounds on w2
//        clb.at(6) = 0; // since lb < 0 < ub
//        cub.at(6) = std::max(ub*ub,lb*lb);

//        // Add linearized bilinear constraint
//        ConstraintPtr quadcon(quadraticIneqConstraint());
//        cs->add(quadcon, v);

//        ConstraintPtr lidcon(lidConstraint(lb,ub));
//        cs->add(lidcon, v);
//    }

//    // Add constraint w7 = w6^2
//    {
//        std::vector<int> v = {6,7};

//        double lb = clb.at(6);
//        double ub = cub.at(6);

//        // Update bounds on w2
//        clb.at(7) = 0; // since lb < 0 < ub
//        cub.at(7) = std::max(ub*ub,lb*lb);

//        // Add linearized bilinear constraint
//        ConstraintPtr quadcon(quadraticIneqConstraint());
//        cs->add(quadcon, v);

//        ConstraintPtr lidcon(lidConstraint(lb,ub));
//        cs->add(lidcon, v);
//    }

//    // Update domain bounds
//    cs->setDomainBounds(clb, cub);

//    cout << "Num variables: " << cs->getDimensionDomain() << endl;
//    cout << "Num constraints: " << cs->getDimensionCodomain() << endl;

//    DenseMatrix cobj(1,dim);
//    cobj.setZero();
//    cobj(0,2) = 4;
//    cobj(0,3) = -2.1;
//    cobj(0,4) = 1.0/3.0;
//    cobj(0,5) = 1;
//    cobj(0,6) = -4;
//    cobj(0,7) = 4;

//    ObjectivePtr obj(new ObjectiveLinear(cobj));

//    std::vector<double> z0(dim, 0);
//    ProblemPtr prob(new Problem(obj, cs));
//    SolverIpopt solver(prob, z0);
//    SolverResult result = solver.optimize();

//    cout << "frelaxed below:" << endl;
//    cout << result << endl;
//    cout << "fopt: " << -1.0316 << endl;
//    cout << "gap: " << -1.0316 - result.objectiveValue << endl;

//    std::vector<double> x = result.primalVariables;
//    //printVec(x);
//}

void subdivision_example()
{
    DenseMatrix coeffs(1,3);
    coeffs << 1, -1, 1;
    std::vector<double> knots = {0,0,0,1,1,1};
    std::vector<unsigned int> degrees = {2};

    std::vector< std::vector<double> > knots2;
    knots2.push_back(knots);

    BSpline bs(coeffs, knots2, degrees);

    // Subdivision
    bs.insertKnots(0.5,0);
    bs.insertKnots(0.5,0);
    bs.insertKnots(0.5,0);

    // Refinement
    bs.insertKnots(0.25,0);
    bs.insertKnots(0.75,0);

    DenseVector x(1); x(0) = 0.5;
    auto y = bs.eval(x);
    cout << y << endl;

    DenseMatrix C = bs.getControlPoints();

    cout << C << endl;
}

void cubicSpline()
{
    //std::vector<double> x = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21};
    //std::vector<double> y = {0, 0, 0, 0, 0.86, 0.98, 1, 1, 1, 1, 1, 0.65, 0.52, 0.5, 0.5, 0.5, 0.48, 0.35, 0, 0, 0, 0};

    std::vector<double> x = {0,     1,      2,      3,      4,      6,      8,          10,     12,     14,     16,     18,     19,     20};
    std::vector<double> y = {1,     1.02,   0.98,    1.0,   0.4,    0,   0.025,      0.05,      0.45,    0.50,   0.75,   1.02,   0.98,    1};

    assert(x.size() == y.size());

    DataTable table;
    for (unsigned int i = 0; i < x.size(); i++)
    {
        table.addSample(x.at(i), y.at(i));
    }

    //table.printTable();

    BSpline bs(table, BSplineType::CUBIC_FREE);

    // Add knots
    bs.insertKnots(1,0);
    bs.insertKnots(2,0);
    //bs.insertKnot(2,0);
    bs.insertKnots(2.7,0);
    bs.insertKnots(3.5,0);
    bs.insertKnots(5,0);
    bs.insertKnots(9,0);
    bs.insertKnots(11,0);
    bs.insertKnots(13,0);
    bs.insertKnots(17,0);
    bs.insertKnots(18.8,0);

    DenseMatrix cp = bs.getControlPoints();

    cout << "Control points:" << endl;
    cout << cp << endl;

}

} // namespace CENSO
