/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "pwlmodels.h"

#include "OptimizationProblem/constraintset.h"
#include "OptimizationProblem/constraintbspline.h"
#include "OptimizationProblem/constraintlinear.h"
#include "bspline.h"
#include "BranchAndBound/branchandbound.h"
#include "SolverInterface/solvergurobi.h"
#include "polytope.h"
#include "triangulation.h"

using std::cout;
using std::endl;

using SPLINTER::DataTable;
using SPLINTER::BSpline;
using SPLINTER::BSplineType;

constexpr double pi()
{
    return std::atan(1)*4;
}

// f* = -1 at x = 5
double sinewave1D(double x)
{
    return std::sin(pi()*(x+0.5)) + 0.1*std::pow(x-5,2);
}

// f* = -1 at x = [5, 0]
double sinewave2D(double x1, double x2)
{
    return std::sin(pi()*(x1+0.5)) + 0.1*std::pow(x1-5,2) + x2*x2;
}

DataTable sampleFunction1D()
{
    std::vector<double> xs = linspace(0, 4*pi(), 1000);

    DataTable samples;
    for (const auto &x : xs)
    {
        double y = sinewave1D(x);
        samples.addSample(x,y);
    }
    return samples;
}

DataTable sampleFunction2D()
{
    std::vector<double> x1s = linspace(0, 4*pi(), 160);
    std::vector<double> x2s = linspace(-2, 2, 20);

    DataTable samples;
    for (const auto &x1 : x1s)
    {
        for (const auto &x2 : x2s)
        {
            std::vector<double> xv;
            xv.push_back(x1);
            xv.push_back(x2);
            double y = sinewave2D(x1, x2);
            samples.addSample(xv,y);
        }
    }
    return samples;
}

DenseMatrix makePointMatrix(const DataTable &points)
{
    DenseMatrix A = DenseMatrix::Zero(points.getNumSamples(), points.getNumVariables()+1);

    int i = 0;
    for (auto it = points.cbegin(); it != points.cend(); ++it)
    {
        auto x = it->getX();
        auto y = it->getY();

        int j = 0;
        for (const auto &xi : x)
        {
            A(i,j) = xi;
            j++;
        }
        A(i,j) = y;
        i++;
    }
    return A;
}

void run_mip_models()
{
//    {
//        cout << "*********************************" << endl;
//        cout << "SOLVING MONOVARIATE CASE" << endl;
//        cout << "*********************************" << endl;

//        DataTable samples = sampleFunction1D();
//        cout << "N: " << samples.getNumSamples() << endl;

//        const int loops = 3;

//        int sum = 0;
//        for (int i = 0; i < loops; i++)
//            sum += P00a_MIP(samples);
//        cout << "Avg. time MIP: " << (double)sum/loops << endl;

//        sum = 0;
//        for (int i = 0; i < loops; i++)
//            sum += P00a_global(samples);
//        cout << "Avg. time Glo: " << (double)sum/loops << endl;
//    }

//    {
//        cout << "*********************************" << endl;
//        cout << "SOLVING BIVARIATE CASE" << endl;
//        cout << "*********************************" << endl;

//        DataTable samples = sampleFunction2D();
//        cout << "N: " << samples.getNumSamples() << endl;

//        P00b_MIP(samples);

//        P00b_global(samples);
//    }
//    return;

    {
        // Number of samples
        std::vector<int> nums = {25, 50, 100, 200, 400, 800, 1600};
        nums.clear(); nums = {1600};

        int loops = 2;
        for (const auto n : nums)
        {
            int sum = 0;
            for (int i = 0; i < loops; i++)
                sum += P01_MIP(n);
            cout << "Avg. time P01_mip(" << n << "): " << (double)sum/loops << endl;
        }

        for (const auto n : nums)
        {
            int sum = 0;
            for (int i = 0; i < loops; i++)
                sum += P01_global(n);
            cout << "Avg. time P01_global(" << n << "): " << (double)sum/loops << endl;
        }

        return;
    }

    // Test
    std::vector<Point2D> points;
    points.push_back(Point2D({0,0},1));
    points.push_back(Point2D({0,1}));
    points.push_back(Point2D({0,2}));
    points.push_back(Point2D({1,0}));
    points.push_back(Point2D({1,1}));
    points.push_back(Point2D({1,2}));
    points.push_back(Point2D({2,0}));
    points.push_back(Point2D({2,1}));
    points.push_back(Point2D({2,2}));
    auto polytopes = triangulate<2>(points);

    auto poly = polytopes.at(0);
    cout << poly << endl;
    DenseMatrix mc = poly.computeInterpolationCoefficients();
    cout << mc << endl;
}

/*
 * Global optimization model
 * for the monovariate case
 */
unsigned int P00a_global(DataTable samples)
{
    //cout << "Running global model..." << endl;

    std::vector<VariablePtr> vars = {
        std::make_shared<Variable>(0),
        std::make_shared<Variable>(1)
    };

    ConstraintSetPtr cs = std::make_shared<ConstraintSet>();

    BSpline bspline(samples, BSplineType::LINEAR);
    ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(vars, bspline, true);
    cs->add(cbs);

    BB::BranchAndBound bnb(cs);
    Timer timer;
    timer.start();
    auto res = bnb.optimize();
    timer.stop();
    cout << res << endl;
    //cout << "Time: " << timer.getMilliSeconds() << " (ms)" << endl;

    return timer.getMilliSeconds();
}

unsigned int P00a_MIP(DataTable samples)
{
    // Create some variables
    Variables vars;
    for (unsigned int i = 0; i < samples.getNumVariables()+1; ++i)
    {
        auto var = std::make_shared<Variable>(0);
        vars.push_back(var);
    }

    // Build MIP model
    PiecewiseLinearModel pwl = constraintPiecewiseLinearDLog<1>(vars, samples);

    assert(pwl.dimension == samples.getNumVariables());
    unsigned int dim = pwl.dimension;
    unsigned int nl = pwl.numAuxiliary; // number of lambdas (cont)
    unsigned int ny = pwl.numBinary; // log number of binary variables y

    // Objective cost
    VariablePtr var = pwl.constraints->getVariableAt(dim);
    var->setCost(1);

    // Constraints
    ConstraintPtr cs = pwl.constraints;
    assert(cs->getNumVariables() == dim+1+nl+ny);

    SolverGurobi solver(cs);
    Timer timer;
    timer.start();
    SolverResult res = solver.optimize();
    timer.stop();
    assert(res.status == SolverStatus::OPTIMAL);

    //cout << "Time: " << timer.getMilliSeconds() << " (ms)" << endl;
    cout << res << endl;
    auto sol = res.primalVariables;
    cout << "Sol: x = " << sol.at(0) << ", z = " << sol.at(1) << endl;

    return timer.getMilliSeconds();
}

/*
 * Global optimization model
 * for the bivariate case
 */
unsigned int P00b_global(DataTable samples)
{
    //cout << "Running global model..." << endl;

    unsigned int dim = samples.getNumVariables();

    std::vector<VariablePtr> vars = {
        std::make_shared<Variable>(0),
        std::make_shared<Variable>(0),
        std::make_shared<Variable>(1)
    };

    ConstraintSetPtr cs = std::make_shared<ConstraintSet>();

    BSpline bspline(samples, BSplineType::LINEAR);
    ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(vars, bspline, true);
    cs->add(cbs);

    BB::BranchAndBound bnb(cs);
    Timer timer;
    timer.start();
    bnb.optimize();
    timer.stop();
    //cout << "Time: " << timer.getMilliSeconds() << " (ms)" << endl;

    return timer.getMilliSeconds();
}

unsigned int P00b_MIP(DataTable samples)
{
    // Create some variables
    Variables vars;
    for (unsigned int i = 0; i < samples.getNumVariables()+1; ++i)
    {
        auto var = std::make_shared<Variable>(0);
        vars.push_back(var);
    }

    // Build MIP model
    PiecewiseLinearModel pwl = constraintPiecewiseLinearDLog<2>(vars, samples);

    ConstraintPtr cs = pwl.constraints;
    unsigned int dim = pwl.dimension; // dimension
    unsigned int nl = pwl.numAuxiliary; // number of lambdas (cont)
    unsigned int ny = pwl.numBinary; // number of ys (binary)

    assert(cs->getNumVariables() == dim+1+nl+ny);
    assert(pwl.dimension == samples.getNumVariables());

    // Objective costs
    VariablePtr var = cs->getVariableAt(dim);
    var->setCost(1);

    SolverGurobi solver(cs);
//    BB::BranchAndBound solver(cs);
    Timer timer;
    timer.start();
    SolverResult result = solver.optimize();
    timer.stop();
    assert(result.status == SolverStatus::OPTIMAL);

    cout << "Time: " << timer.getMilliSeconds() << " (ms)" << endl;
    cout << result << endl;
    auto sol = result.primalVariables;
    cout << "Sol: x1 = " << sol.at(0) << ", x2 = " << sol.at(1) << ", z = " << sol.at(2) << endl;

    return timer.getMilliSeconds();
}

unsigned int P01_MIP(int numSamples)
{
    assert(numSamples >= 2);
    std::vector<double> xs = linspace(0, 3, numSamples);

    auto con1 = [](double x) { return -2 -8*x*x +8*x*x*x -2*x*x*x*x; };
    auto con2 = [](double x) { return -36 +96*x -88*x*x +32*x*x*x -4*x*x*x*x; };

    DataTable samples1, samples2;
    for (const auto &xv : xs)
    {
        double y1 = con1(xv);
        double y2 = con2(xv);

        samples1.addSample(xv, y1);
        samples2.addSample(xv, y2);
    }

    // Create some variables
    Variables vars = {
        std::make_shared<Variable>(-1, 0, 3),
        std::make_shared<Variable>(1, -4, 0)
    };

    PiecewiseLinearModel pwl1 = constraintPiecewiseLinearDCC<1>(vars, samples1);
    PiecewiseLinearModel pwl2 = constraintPiecewiseLinearDCC<1>(vars, samples2);

    ConstraintPtr cs1 = pwl1.constraints;
    ConstraintPtr cs2 = pwl2.constraints;

    ConstraintSetPtr cs = std::make_shared<ConstraintSet>();

    cs->add(cs1);
    cs->add(cs2);

    cs->writeToGAMS("test.gms");

    SolverGurobi solver(cs);
    Timer timer;
    timer.start();
    SolverResult result = solver.optimize();
    timer.stop();
    assert(result.status == SolverStatus::OPTIMAL);

    cout << "Time: " << timer.getMilliSeconds() << " (ms)" << endl;
    cout << result << endl;
    auto sol = result.primalVariables;
    cout << "Sol: x1 = " << sol.at(0) << ", x2 = " << sol.at(1) << endl;
    //printVector(sol);

    DenseVector xsol(sol.size());
    for (unsigned int i = 0; i < sol.size(); i++)
        xsol(i) = sol.at(i);
    if (!cs->checkFeasibility(xsol))
        cout << "Optimal point not feasible!" << endl;

    return timer.getMilliSeconds();
}

unsigned int P01_global(int numSamples)
{
    assert(numSamples >= 2);
    std::vector<double> xs = linspace(0, 3, numSamples);

    auto con1 = [](double x) { return -2 -8*x*x +8*x*x*x -2*x*x*x*x; };
    auto con2 = [](double x) { return -36 +96*x -88*x*x +32*x*x*x -4*x*x*x*x; };

    DataTable samples1, samples2;
    for (const auto &xv : xs)
    {
        double y1 = con1(xv);
        double y2 = con2(xv);

        samples1.addSample(xv, y1);
        samples2.addSample(xv, y2);
    }

    int dim = 2;
    int aux = 2;

    Variables vars = {
        std::make_shared<Variable>(-1, 0, 3),
        std::make_shared<Variable>(1, -4, 0),
        std::make_shared<Variable>(0),
        std::make_shared<Variable>(0)
    };

    ConstraintSetPtr cs = std::make_shared<ConstraintSet>();

    // Build B-splines
    BSpline bspline1(samples1, BSplineType::LINEAR);
    BSpline bspline2(samples2, BSplineType::LINEAR);

    std::vector<VariablePtr> vars1 = {vars.at(0), vars.at(2)};
    std::vector<VariablePtr> vars2 = {vars.at(0), vars.at(3)};

    ConstraintPtr cbs1 = std::make_shared<ConstraintBSpline>(vars1, bspline1, true);
    ConstraintPtr cbs2 = std::make_shared<ConstraintBSpline>(vars2, bspline2, true);

    cs->add(cbs1);
    cs->add(cbs2);

    DenseMatrix A = DenseMatrix::Zero(2,4);
    A(0,1) = -1; A(0,2) = 1;
    A(1,1) = -1; A(1,3) = 1;
    DenseVector b; b.setZero(2);
    ConstraintPtr lincon = std::make_shared<ConstraintLinear>(vars, A, b, false);
    cs->add(lincon);

//    SolverGurobi solver(cs);
    BB::BranchAndBound solver(cs);
    Timer timer;
    timer.start();
    SolverResult result = solver.optimize();
    timer.stop();
    assert(result.status == SolverStatus::OPTIMAL);

    cout << "Time: " << timer.getMilliSeconds() << " (ms)" << endl;
    cout << result << endl;
    auto sol = result.primalVariables;
    cout << "Sol: x1 = " << sol.at(0) << ", x2 = " << sol.at(1) << endl;
    //printVector(sol);

    // Check feasibility (note that BB may use a different threshold than constraint)
    DenseVector xsol(sol.size());
    for (unsigned int i = 0; i < sol.size(); i++)
        xsol(i) = sol.at(i);
    if (!cs->checkFeasibility(xsol))
        cout << "Optimal point not feasible!" << endl;

    return timer.getMilliSeconds();
}

/*
 * Returns constraint set with
 * N+1+N*|P|+|P| variables [x z lambda y],
 * where x in R(N), z in R(1), lambda in R(N*|P|),
 * and y in Z(|P|).
 */
template <unsigned int N>
PiecewiseLinearModel constraintPiecewiseLinearMC(Variables vars, const DataTable &samples)
{
    assert(samples.getNumVariables() == N);
    assert(vars.size() == N+1);

    DenseMatrix A = makePointMatrix(samples);

    auto points = makePointVector<N>(A);

    auto polytopes = triangulate<N>(points);

    int dim = N; // Dimension
    int nl = dim*polytopes.size(); // Number of xs (cont) is N*|P|
    int ny = polytopes.size(); // Number of ys (binary) is |P|

    /*
     * Constraint set for the MC model
     * The variables are ordered as follows:
     * [x, z, x_p, y_p], where
     * x in R(N), z in R(1), x_p in R(N*|P|), y_p in R(|P|)
     */
    for (int i = 0; i < nl; ++i)
    {
        auto var = std::make_shared<Variable>(0);
        vars.push_back(var);
    }
    for (int i = 0; i < ny; ++i)
    {
        auto var = std::make_shared<Variable>(0, 0, 1, VariableType::BINARY);
        vars.push_back(var);
    }

    // Create constraint set with variables (variable order maintained)
    ConstraintSetPtr cs = std::make_shared<ConstraintSet>(vars);

    {
        // Add constraint sum x_p = x
        DenseMatrix A = DenseMatrix::Ones(dim,dim+nl);
        DenseMatrix I = DenseMatrix::Identity(dim,dim);
        A.block(0,0,dim,dim) = -I;
        for (unsigned int i = 1; i <= polytopes.size(); i++)
            A.block(0,i*dim,dim,dim) = I;

        DenseVector b = DenseMatrix::Zero(dim,1);

        std::vector<VariablePtr> cvars;
        for (int i = 0; i < dim; i++)
            cvars.push_back(vars.at(i)); // Add x
        for (int i = 0; i < nl; i++)
            cvars.push_back(vars.at(dim+1+i)); // Add x_p

        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, true);
        cs->add(lincon);
    }

    {
        // Add constraint sum_p (m_p x_p + c_p y_p) <= z
        DenseMatrix A = DenseMatrix::Zero(1, nl+ny+1);
        A(0,0) = -1;
        DenseVector b = DenseMatrix::Zero(1,1);

        std::vector<VariablePtr> cvars;
        cvars.push_back(vars.at(dim)); // z

        int i = 0;
        for (const auto &poly : polytopes)
        {
            DenseMatrix mc = poly.computeInterpolationCoefficients();
            mc.transposeInPlace(); // Transpose to row vector
            assert(mc.rows() == 1);
            assert(mc.cols() == dim+1);
            A.block(0, 1+(dim+1)*i, 1, dim+1) = mc;
            for (int j = 0; j < dim; j++)
                cvars.push_back(vars.at(dim+1+dim*i+j)); // x_p
            cvars.push_back(vars.at(dim+1+nl+i)); // y_p

            i++;
        }

        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, false);
        cs->add(lincon);
    }

    {
        // Add constraints A_p x_p <= b_p y_p
        int i = 0;
        for (auto &poly : polytopes)
        {
            DenseMatrix A;
            DenseVector b;
            poly.computeFacets(A, b);
            int npf = A.rows(); // Facets of polytope
            assert(npf == b.rows());
            assert(A.cols() == dim);
            assert(b.cols() == 1);

            DenseMatrix A2 = DenseMatrix::Zero(npf, dim+1);
            A2.block(0, 0, npf, dim) = A;
            A2.block(0, dim, npf, 1) = -b;
            DenseVector b2 = DenseMatrix::Zero(npf,1);

            std::vector<VariablePtr> cvars;
            for (int j = 0; j < dim; j++)
                cvars.push_back(vars.at(dim+1+dim*i+j)); // x_p
            cvars.push_back(vars.at(dim+1+nl+i)); // y_p

            ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A2, b2, false);
            cs->add(lincon);

            i++;
        }
    }

    {
        // Add constraint sum_p y_p = 1
        DenseMatrix A = DenseMatrix::Ones(1, ny);
        DenseVector b(1); b(0) = 1;

        std::vector<VariablePtr> cvars;
        for (int i = 0; i < ny; i++)
            cvars.push_back(vars.at(dim+1+nl+i)); // y_p

        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, true);
        cs->add(lincon);
    }

    PiecewiseLinearModel pwl(cs, dim, nl, ny);

    return pwl;
}

/*
 * Returns constraint set with
 * N+1+N*|P|+|P| variables [x z lambda y],
 * where x in R(N), z in R(1), lambda in R(V(P)),
 * and y in Z(|P|).
 */
template <unsigned int N>
PiecewiseLinearModel constraintPiecewiseLinearCC(Variables vars, const DataTable &samples)
{
    assert(samples.getNumVariables() == N);
    assert(vars.size() == N+1);

    DenseMatrix A = makePointMatrix(samples);

    auto points = makePointVector<N>(A);

    auto polys = triangulate<N>(points);

    int dim = N; // dimension
    int nl = points.size(); // number of lambdas (cont)
    int ny = polys.size(); // number of ys (binary)

    // Create variables
    for (int i = 0; i < nl; ++i)
    {
        auto var = std::make_shared<Variable>(0);
        vars.push_back(var);
    }
    for (int i = 0; i < ny; ++i)
    {
        auto var = std::make_shared<Variable>(0, 0, 1, VariableType::BINARY);
        vars.push_back(var);
    }

    // Create map (point, index) for lambda variables
    std::map<Point<N>, VariablePtr> var_map;
    for (unsigned int i = 0; i < points.size(); i++)
        var_map.emplace(points.at(i), vars.at(dim+1+i));

    // Create constraint set with variables (variable order maintained)
    ConstraintSetPtr cs = std::make_shared<ConstraintSet>(vars);

    {
        // Add constraint sum_v (lambda_v v) = x
        DenseMatrix A = DenseMatrix::Ones(dim, dim+nl);
        A.block(0,0,dim,dim) = -1*DenseMatrix::Identity(dim,dim);
        DenseVector b = DenseMatrix::Zero(dim,1);

        std::vector<VariablePtr> cvars;
        for (int i = 0; i < dim; i++)
            cvars.push_back(vars.at(i)); // x
        for (int i = 0; i < nl; i++)
            cvars.push_back(vars.at(dim+1+i)); // lambda

        for (int i = 0; i < nl; i++)
        {
            auto x = points.at(i).x;
            for (int j = 0; j < dim; j++)
            {
                A(j,i+dim) = x.at(j);
            }
        }

        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, true);
        cs->add(lincon);
    }

    {
        // Add constraint sum_v lambda_v (m_p v + c_p) <= z
        DenseMatrix A = DenseMatrix::Zero(1, nl+1);
        A(0,0) = -1;
        DenseVector b = DenseMatrix::Zero(1,1);

        std::vector<VariablePtr> cvars;
        cvars.push_back(vars.at(dim)); // z

        for (unsigned int i = 0; i < points.size(); i++)
        {
            auto p = points.at(i);
            A(0,1+i) = p.y; // y = m_p v + c_p

            auto varit = var_map.find(p);
            assert(varit != var_map.end());
            auto var = varit->second;
            cvars.push_back(var); // lambda_v
        }

        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, false);
        cs->add(lincon);
    }

    {
        // Add constraint sum_v lambda_v = 1
        DenseMatrix A = DenseMatrix::Ones(1, points.size());
        DenseVector b(1); b(0) = 1;

        std::vector<VariablePtr> cvars;
        for (unsigned int i = 0; i < points.size(); i++)
        {
            auto p = points.at(i);
            auto varit = var_map.find(p);
            assert(varit != var_map.end());
            auto var = varit->second;
            cvars.push_back(var); // lambda_v
        }

        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, true);
        cs->add(lincon);
    }

    {
        // Add constraints lambda_v <= sum_p y_p, for all v

        for (unsigned int i = 0; i < points.size(); i++)
        {
            auto p = points.at(i);
            auto varit = var_map.find(p);
            assert(varit != var_map.end());
            VariablePtr var = varit->second;

            // Find polytopes that p is a vertex in
            std::vector<int> polys_p;
            int j = 0;
            for (const auto &poly : polys)
            {
                if (poly.isVertex(p))
                    polys_p.push_back(j);
                j++;
            }

            DenseMatrix A = DenseMatrix::Zero(1,polys_p.size()+1);
            A(0,0) = 1;
            DenseVector b = DenseMatrix::Zero(1,1);

            std::vector<VariablePtr> cvars;
            cvars.push_back(var); // x_p

            for (unsigned int k = 0; k < polys_p.size(); k++)
            {
                A(0,1+k) = -1;
                int yp = polys_p.at(k);
                cvars.push_back(vars.at(dim+1+nl+yp)); // y_p
            }

            ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, false);
            cs->add(lincon);
        }
    }

    {
        // Add constraint sum_p y_p = 1
        DenseMatrix A = DenseMatrix::Ones(1, ny);
        DenseVector b(1); b(0) = 1;

        std::vector<VariablePtr> cvars;
        for (int i = 0; i < ny; i++)
            cvars.push_back(vars.at(dim+1+nl+i)); // y_p

        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, true);
        cs->add(lincon);
    }

    // lambda_v >= 0
    for (int i = 0; i < nl; i++)
        vars.at(dim+1+i)->setLowerBound(0);

    PiecewiseLinearModel pwl(cs, dim, nl, ny);

    return pwl;
}

/*
 * Disaggregated Convex Combination (DCC) model
 * for the multivariate case
 */
template <unsigned int N>
PiecewiseLinearModel constraintPiecewiseLinearDCC(Variables vars, const DataTable &samples)
{
    assert(samples.getNumVariables() == N);
    assert(vars.size() == N+1);

    DenseMatrix A = makePointMatrix(samples);

    auto points = makePointVector<N>(A);

    auto polys = triangulate<N>(points);

    int dim = N; // dimension
    int nl = 0; // number of lambdas (cont)
    int ny = polys.size(); // number of ys (binary)

    // One auxilliary variable for each vertex in each polytope
    for (const auto &p : polys)
        nl += p.vertices.size();

    // Create variables
    for (int i = 0; i < nl; ++i)
    {
        auto var = std::make_shared<Variable>(0);
        vars.push_back(var);
    }
    for (int i = 0; i < ny; ++i)
    {
        auto var = std::make_shared<Variable>(0, 0, 1, VariableType::BINARY);
        vars.push_back(var);
    }

    // Create constraint set with variables (variable order maintained)
    ConstraintSetPtr cs = std::make_shared<ConstraintSet>(vars);

    {
        // Add constraint sum_p sum_v lambda_{p,v} v = x
        DenseMatrix A = DenseMatrix::Ones(dim, dim+nl);
        DenseMatrix I = DenseMatrix::Identity(dim,dim);
        A.block(0,0,dim,dim) = -I;
        DenseVector b = DenseMatrix::Zero(dim,1);

        std::vector<VariablePtr> cvars;
        for (int i = 0; i < dim; i++)
            cvars.push_back(vars.at(i)); // x
        for (int i = 0; i < nl; i++)
            cvars.push_back(vars.at(dim+1+i)); // lambdas

        int k = 0;
        for (const auto &p : polys)
        {
            for (const auto &v : p.vertices)
            {
                // Vertex
                auto x = v.x;
                for (int l = 0; l < dim; l++)
                    A(l,dim+k) = x.at(l);
                k++;
            }
        }

        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, true);
        cs->add(lincon);
    }

    {
        // Add constraint sum_p sum_v lambda_v (m_p v + c_p) <= z
        DenseMatrix A = DenseMatrix::Zero(1, 1+nl);
        A(0,0) = -1;
        DenseVector b = DenseMatrix::Zero(1,1);

        std::vector<VariablePtr> cvars;
        cvars.push_back(vars.at(dim)); // z
        for (int i = 0; i < nl; i++)
            cvars.push_back(vars.at(dim+1+i)); // lambdas

        int k = 0;
        for (const auto &p : polys)
        {
            for (const auto &v : p.vertices)
            {
                // Vertex
                A(0,k+1) = v.y; // y = m_p v + c_p
                k++;
            }
        }

        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, false);
        cs->add(lincon);
    }

    {
        // Add constraint sum_v lambda_{P,v} = y_p
        int ip = 0; // polytope counter
        int il = 0; // lambda counter

        // Assuming most polytopes have an equal number of vertices
        DenseMatrix A0 = DenseMatrix::Ones(1,1+polys.front().vertices.size());
        A0(0,0) = -1;
        DenseVector b = DenseMatrix::Zero(1,1);

        for (const auto &p : polys)
        {
            std::vector<VariablePtr> cvars;
            cvars.push_back(vars.at(dim+1+nl+ip)); // y

            for (const auto &v : p.vertices)
            {
                cvars.push_back(vars.at(dim+1+il));
                il++;
            }

            // Check if preallocated matrix can be used
            if (A0.cols() != 1+p.vertices.size())
            {
                DenseMatrix A = DenseMatrix::Ones(1,1+p.vertices.size());
                A(0,0) = -1;
                ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, true);
                cs->add(lincon);
            }
            else
            {
                ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A0, b, true);
                cs->add(lincon);
            }

            ip++;
        }
    }

    {
        // Add constraint sum_p y_p = 1
        DenseMatrix A = DenseMatrix::Ones(1, ny);
        DenseVector b(1); b(0) = 1;

        std::vector<VariablePtr> cvars;
        for (int i = 0; i < ny; i++)
            cvars.push_back(vars.at(dim+1+nl+i)); // y_p

        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, true);
        cs->add(lincon);
    }

    // lambda_v >= 0
    for (int i = 0; i < nl; i++)
        vars.at(dim+1+i)->setLowerBound(0);

    PiecewiseLinearModel pwl(cs, dim, nl, ny);

    return pwl;
}

/*
 * Disaggregated convex combination Logarithmic (DLog) model
 * for the multivariate case
 *
 * For a problem with |P| polygons
 * nlog = ceil(log2(|P|)) binary variables are required
 *
 * An injective mapping B : P -> {0,1}^nlog is needed
 * to model the constraints. Here each P is given a number
 * l in {1,...,nlog}. Below, B is an injective mapping that
 * takes this number to its binary form.
 */
template <unsigned int N>
PiecewiseLinearModel constraintPiecewiseLinearDLog(Variables vars, const DataTable &samples)
{
    assert(samples.getNumVariables() == N);
    assert(vars.size() == N+1);

    DenseMatrix A = makePointMatrix(samples);

    auto points = makePointVector<N>(A);

    auto polys = triangulate<N>(points);

    int dim = N; // dimension
    int nl = 0; // number of lambdas (cont)
    int ny = std::ceil(std::log2(polys.size())); // log number of binary variables y

    // One auxilliary variable for each vertex of each polytope
    for (const auto &p : polys)
        nl += p.vertices.size();

    // Create variables
    for (int i = 0; i < nl; ++i)
    {
        auto var = std::make_shared<Variable>(0);
        vars.push_back(var);
    }
    for (int i = 0; i < ny; ++i)
    {
        auto var = std::make_shared<Variable>(0, 0, 1, VariableType::BINARY);
        vars.push_back(var);
    }

    /*
     * Check if digit l of integer b is 1,
     * where l is nonnegative (zero-indexed).
     */
    auto isDigitOne = [](int b, int l)
    {
        assert(l >= 0);
        int mask = 1 << l;
        int res = (b & mask) >> l;
        return (res == 1);
    };

    // Create constraint set with variables (variable order maintained)
    ConstraintSetPtr cs = std::make_shared<ConstraintSet>(vars);

    {
        // Add constraint sum_p sum_v lambda_{p,v} v = x
        DenseMatrix A = DenseMatrix::Ones(dim, dim+nl);
        DenseMatrix I = DenseMatrix::Identity(dim,dim);
        A.block(0,0,dim,dim) = -I;
        DenseVector b = DenseMatrix::Zero(dim,1);

        std::vector<VariablePtr> cvars;
        for (int i = 0; i < dim; i++)
            cvars.push_back(vars.at(i)); // x
        for (int i = 0; i < nl; i++)
            cvars.push_back(vars.at(dim+1+i)); // lambdas

        int k = 0;
        for (const auto &p : polys)
        {
            for (const auto &v : p.vertices)
            {
                // Vertex
                auto x = v.x;
                for (int i = 0; i < dim; i++)
                    A(i,dim+k) = x.at(i);
                k++;
            }
        }

        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, true);
        cs->add(lincon);
    }

    {
        // Add constraint sum_p sum_v lambda_{p,v} (m_p v + c_p) <= z
        DenseMatrix A = DenseMatrix::Zero(1, nl+1);
        A(0,0) = -1;
        DenseVector b = DenseMatrix::Zero(1,1);

        std::vector<VariablePtr> cvars;
        cvars.push_back(vars.at(dim)); // z
        for (int i = 0; i < nl; i++)
            cvars.push_back(vars.at(dim+1+i)); // lambdas

        int k = 0;
        for (const auto &p : polys)
        {
            for (const auto &v : p.vertices)
            {
                // Vertex
                A(0,1+k) = v.y; // y = m_p v + c_p
                k++;
            }
        }

        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, false);
        cs->add(lincon);
    }

    {
        // Add constraint sum_p sum_v lambda_{p,v} = 1
        DenseMatrix A = DenseMatrix::Ones(1, nl);
        DenseVector b(1); b(0) = 1;

        std::vector<VariablePtr> cvars;
        for (int i = 0; i < nl; i++)
            cvars.push_back(vars.at(dim+1+i)); // lambda_{p,v}

        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, true);
        cs->add(lincon);
    }

    {
        // For all l in L(P) = {1,...,ceil(log2(|P|))}
        // add constraint sum_p+ sum_v lambda_{p,v} <= y_l
        for (int l = 0; l < ny; l++)
        {
            std::vector<VariablePtr> cvars;
            cvars.push_back(vars.at(dim+1+nl+l)); // y_l

            int numl = 0;
            for (unsigned int i = 0; i < polys.size(); i++)
            {
                auto p = polys.at(i);
                if (!isDigitOne(i,l))
                {
                    numl += p.vertices.size();
                    continue;
                }

                for (unsigned int j = 0; j < p.vertices.size(); j++)
                {
                    cvars.push_back(vars.at(dim+1+numl)); // lambda_{p,v}
                    numl++;
                }
            }
            assert(numl == nl);

            DenseMatrix A = DenseMatrix::Ones(1,cvars.size());
            A(0,0) = -1;
            DenseVector b = DenseMatrix::Zero(1,1);

            ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, false);
            cs->add(lincon);
        }
    }

    {
        // For all l in L(P) = {1,...,ceil(log2(|P|))}
        // add constraint sum_p0 sum_v lambda_{p,v} <= (1 - y_l)
        for (int l = 0; l < ny; l++)
        {
            std::vector<VariablePtr> cvars;
            cvars.push_back(vars.at(dim+1+nl+l)); // y_l

            int numl = 0;
            for (unsigned int i = 0; i < polys.size(); i++)
            {
                auto p = polys.at(i);
                if (isDigitOne(i,l))
                {
                    numl += p.vertices.size();
                    continue;
                }

                for (unsigned int j = 0; j < p.vertices.size(); j++)
                {
                    cvars.push_back(vars.at(dim+1+numl)); // lambda_{p,v}
                    numl++;
                }
            }
            assert(numl == nl);

            DenseMatrix A = DenseMatrix::Ones(1,cvars.size());
            A(0,0) = 1;
            DenseVector b = DenseMatrix::Ones(1,1);

            ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, false);
            cs->add(lincon);
        }
    }

    // lambda_v >= 0
    for (int i = 0; i < nl; i++)
        vars.at(dim+1+i)->setLowerBound(0);

    PiecewiseLinearModel pwl(cs, dim, nl, ny);

    return pwl;
}
