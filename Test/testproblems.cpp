/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <iostream>
#include <fstream>
#include <string>

#include "Utils/definitions.h"
#include "OptimizationProblem/constraint.h"
#include "Utils/definitions.h"

#include "OptimizationProblem/constraintlinear.h"
#include "OptimizationProblem/constraintbilinear.h"
#include "OptimizationProblem/constraintset.h"

//void testFeasibilityCheck()
//{
//    // Problem 4 (Nataraj)

//    int dim = 2;

//    // x1,x2
//    std::vector<double> lb = {-1,-1};
//    std::vector<double> ub = {10,10};

//    ConstraintCompositePtr cs(new ConstraintComposite(dim, lb, ub));

//    {
//        DenseMatrix A(2,2);
//        A.setZero();
//        A(0,0) = 1; A(0,1) = 1;
//        A(1,0) = -1; A(1,1) = -1;

//        DenseVector b; b.setZero(2);
//        b(0) = 10; b(1) = 10;

//        ConstraintPtr c(new ConstraintLinear(A,b,false));

//        std::vector<int> vars = {0,1};
//        cs->add(c,vars);
//    }

//    DenseVector x(2);
//    x(0) = -0.00000000;
//    x(1) = -10;

//    if (cs->checkFeasibility(x))
//    {
//        cout << "Feasible!" << endl;
//    }
//    else
//    {
//        cout << "Infeasible!" << endl;
//    }

//    std::vector<double> x0 = {10,10};

//    DenseMatrix c(1,2);
//    c.setOnes();
//    c(0) = -2;
//    ObjectivePtr obj(new ObjectiveLinear(c));

//    //SolverIpopt solver(obj, cs, x0);
//    SolverGurobi solver(obj, cs, x0);
//    solver.optimize();
//}

//void diamondExample()
//{
//    int d = 2; // dim
//    int n = 4; // number of points

//    // A holds vertices of convex hull
//    // Diamond shape example
//    DenseMatrix A; A.setZero(d,n);
//    A(0,0) = -1;
//    A(1,0) = 0;

//    A(0,1) = 1;
//    A(1,1) = 0;

//    A(0,2) = 0;
//    A(1,2) = -1;

//    A(0,3) = 0;
//    A(1,3) = 1;

//    DenseMatrix ones;
//    ones.setOnes(1,n);

//    DenseMatrix A2(d+1,n);
//    A2.block(0, 0,  d,  n) = A;
//    A2.block(d, 0,  1,  n) = ones;

//    // Objective
//    DenseMatrix c; c.setOnes(1,n);
//    ObjectivePtr obj(new ObjectiveLinear(c));

//    // [0,1] constraints on variables
//    DoubleVec lb(n,0);
//    DoubleVec ub(n,1);
//    DoubleVec z0(n,0); //z0.at(0) = 1;

//    int points_inside = 0;
//    int points_total = 1000;

//    for (int k = 0; k < points_total; k++)
//    {
//        // Ax = b, where b holds a points (last element is for the convexity requirement)
//        VecD b2(d+1);
//        b2(0) = randomInteger(-1,1); // Random integer point
//        b2(1) = randomInteger(-1,1);
//        b2(2) = 1;

//        ConstraintPtr cHull(new ConstraintLinear(A2, b2, true));
//        cHull->setDomainBounds(lb, ub);

//        SolverIpopt ipopt(obj, cHull, z0);

//        int status = ipopt.optimize();
//        if (status == 1)
//        {
//            //cout << "Point in convex hull!" << endl;
//            points_inside++;
//        }
//        else
//        {
//            //cout << "Point not in convex hull!" << endl;
//        }

//    }

//    // Expect that on average 5/9 of the points are inside the convex hull
//    // This is due to the diamond shape and the random integers
//    // (9 possible points, 5 are inside (on the boundary) of the convex hull)
//    cout << "Points inside: " << points_inside << " / " << points_total << endl;
//}

// Test problems from paper
//void P01()
//{
//    // Problem 1 (Nataraj)
//    cout << "\n\nSolving problem P01..." << endl;

//    int dim = 4;

//    // x1,x2,l1,l2
//    std::vector<double> lb = {0,0,-INF,-INF};
//    std::vector<double> ub = {3,4,INF,INF};

//    ConstraintCompositePtr cs(new ConstraintComposite(dim, lb, ub));

//    { // 2*x1^4 - 8*x1^3 + 8*x1^2 + 2
//        int var = 0;
//        std::vector<double> thislb = {lb.at(var)};
//        std::vector<double> thisub = {ub.at(var)};

//        std::vector<int> deg = {4};

//        DenseVector c(5); c.setZero();
//        c(0) = 2;
//        c(1) = 0;
//        c(2) = 8;
//        c(3) = -8;
//        c(4) = 2;

//        DenseMatrix T = getTransformationMatrix(deg, thislb, thisub);
//        DenseMatrix coeffs = T*c;

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {var,2};

//        cs->add(cbs,vars);
//    }

//    { // 4*x1^4 - 32*x1^3 + 88*x1^2 - 96*x1 + 36
//        int var = 0;
//        std::vector<double> thislb = {lb.at(var)};
//        std::vector<double> thisub = {ub.at(var)};

//        std::vector<int> deg = {4};

//        DenseVector c(5); c.setZero();
//        c(0) = 36;
//        c(1) = -96;
//        c(2) = 88;
//        c(3) = -32;
//        c(4) = 4;

//        DenseMatrix T = getTransformationMatrix(deg, thislb, thisub);
//        DenseMatrix coeffs = T*c;

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {var,3};

//        cs->add(cbs,vars);
//    }

//    { // x2 - l1 <= 0, x2 - l2 <= 0

//        DenseMatrix A(2,3);
//        A.setZero();
//        A(0,0) = 1; A(0,1) = -1;
//        A(1,0) = 1; A(1,2) = -1;

//        DenseVector b; b.setZero(2);

//        ConstraintPtr c(new ConstraintLinear(A,b,false));

//        std::vector<int> vars = {1,2,3};
//        cs->add(c,vars);
//    }

//    DenseMatrix cobj(1,dim); cobj.setZero();
//    cobj(0,0) = -1;
//    cobj(0,1) = -1;
//    ObjectivePtr obj(new ObjectiveLinear(cobj));

//    std::vector<int> vt(dim,CONTINUOUS);
//    //std::vector<int> vt= {INTEGER,INTEGER,CONTINUOUS,CONTINUOUS}; // Integer problem
//    std::vector<int> bv = {0,1};
//    std::vector<double> z0(dim,0);

//    BB::BranchAndBound bnb(obj, cs, z0, vt, bv);
//    Timer timer;
//    timer.start();
//    bnb.optimize();
//    timer.stop();
//    cout << "Time: " << timer.getMilliSeconds() << " (ms)" << endl;
//}

//void P02()
//{
//    // Problem 2 (Nataraj)
//    cout << "\n\nSolving problem P02..." << endl;

//    int dim = 7;

//    // x1,x2,l1,...,l5
//    std::vector<double> lb = {13,0,-INF,-INF,0,0,0};
//    std::vector<double> ub = {100,100,INF,INF,INF,INF,INF};

//    ConstraintCompositePtr cs(new ConstraintComposite(dim, lb, ub));

//    { // (x1-10)^3
//        int var = 0;
//        std::vector<double> thislb = {lb.at(var)};
//        std::vector<double> thisub = {ub.at(var)};

//        std::vector<int> deg = {3};

//        DenseVector c(4); c.setZero();
//        c(0) = -1000;
//        c(1) = 300;
//        c(2) = -30;
//        c(3) = 1;

//        DenseMatrix T = getTransformationMatrix(deg, thislb, thisub);
//        DenseMatrix coeffs = T*c;

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {var,2};

//        cs->add(cbs,vars);
//    }

//    { // (x2-20)^3
//        int var = 1;
//        std::vector<double> thislb = {lb.at(var)};
//        std::vector<double> thisub = {ub.at(var)};

//        std::vector<int> deg = {3};

//        DenseVector c(4); c.setZero();
//        c(0) = -8000;
//        c(1) = 1200;
//        c(2) = -60;
//        c(3) = 1;

//        DenseMatrix T = getTransformationMatrix(deg, thislb, thisub);
//        DenseMatrix coeffs = T*c;

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {var,3};

//        cs->add(cbs,vars);
//    }

//    { // (x1-5)^2
//        int var = 0;
//        std::vector<double> thislb = {lb.at(var)};
//        std::vector<double> thisub = {ub.at(var)};

//        std::vector<int> deg = {2};

//        DenseVector c(3); c.setZero();
//        c(0) = 25;
//        c(1) = -10;
//        c(2) = 1;

//        DenseMatrix T = getTransformationMatrix(deg,thislb,thisub);
//        DenseMatrix coeffs = T*c;

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {var,4};

//        cs->add(cbs,vars);
//    }

//    { // (x2-5)^2
//        int var = 1;
//        std::vector<double> thislb = {lb.at(var)};
//        std::vector<double> thisub = {ub.at(var)};

//        std::vector<int> deg = {2};

//        DenseVector c(3); c.setZero();
//        c(0) = 25;
//        c(1) = -10;
//        c(2) = 1;

//        DenseMatrix T = getTransformationMatrix(deg,thislb,thisub);
//        DenseMatrix coeffs = T*c;

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {var,5};

//        cs->add(cbs,vars);
//    }

//    { // (x1-6)^2
//        int var = 0;
//        std::vector<double> thislb = {lb.at(var)};
//        std::vector<double> thisub = {ub.at(var)};

//        std::vector<int> deg = {2};

//        DenseVector c(3); c.setZero();
//        c(0) = 36;
//        c(1) = -12;
//        c(2) = 1;

//        DenseMatrix T = getTransformationMatrix(deg,thislb,thisub);
//        DenseMatrix coeffs = T*c;

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {var,6};

//        cs->add(cbs,vars);
//    }

//    { // -l3 - l4 <= -100

//        DenseMatrix A(1,2); A(0,0) = -1; A(0,1) = -1;

//        DenseVector b(1); b(0) = -100;

//        ConstraintPtr c(new ConstraintLinear(A,b,false));

//        std::vector<int> vars = {4,5};
//        cs->add(c,vars);
//    }

//    { // l5 + l4 <= 82.81

//        DenseMatrix A(1,2); A(0,0) = 1; A(0,1) = 1;

//        DenseVector b(1); b(0) = 82.81;

//        ConstraintPtr c(new ConstraintLinear(A,b,false));

//        std::vector<int> vars = {6,5};
//        cs->add(c,vars);
//    }

//    DenseMatrix cobj(1,dim); cobj.setZero();
//    cobj(0,2) = 1;
//    cobj(0,3) = 1;
//    ObjectivePtr obj(new ObjectiveLinear(cobj));

//    std::vector<int> vt(dim,CONTINUOUS);
//    //std::vector<int> vt = {INTEGER,INTEGER,CONTINUOUS,CONTINUOUS,CONTINUOUS,CONTINUOUS,CONTINUOUS};
//    std::vector<int> bv = {0,1};
//    std::vector<double> z0(dim,0);

//    BB::BranchAndBound bnb(obj, cs, z0, vt, bv);
//    Timer timer;
//    timer.start();
//    bnb.optimize();
//    timer.stop();
//    cout << "Time: " << timer.getMilliSeconds() << " (ms)" << endl;
//}

//void P03()
//{
//    // Problem 3 (Nataraj)
//    cout << "\n\nSolving problem P03..." << endl;

//    int dim = 4;

//    // x1,x2,l1,l2
//    std::vector<double> lb = {-10,-10,0,-1000};
//    std::vector<double> ub = {10,10,100,1000};

//    ConstraintCompositePtr cs(new ConstraintComposite(dim, lb, ub));

//    { // x1^2 = l1
//        int var = 0;
//        std::vector<double> thislb = {lb.at(var)};
//        std::vector<double> thisub = {ub.at(var)};

//        std::vector<int> deg = {2};

//        DenseVector c(3);
//        c.setZero();
//        c(0) = 0;
//        c(1) = 0;
//        c(2) = 1;

//        DenseMatrix T = getTransformationMatrix(deg, thislb, thisub);
//        DenseMatrix coeffs = T*c;

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {var,2};

//        cs->add(cbs,vars);
//    }

//    { // x1^3 = l1
//        int var = 0;
//        std::vector<double> thislb = {lb.at(var)};
//        std::vector<double> thisub = {ub.at(var)};

//        std::vector<int> deg = {3};

//        DenseVector c(4); c.setZero();
//        c(0) = 0;
//        c(1) = 0;
//        c(2) = 0;
//        c(3) = 1;

//        DenseMatrix T = getTransformationMatrix(deg, thislb, thisub);
//        DenseMatrix coeffs = T*c;

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {var,3};

//        cs->add(cbs,vars);
//    }

//    { // x2 - l1 <= 0, x2 - l2 <= 0

//        DenseMatrix A(2,3);
//        A.setZero();
//        A(0,0) = -1; A(0,1) = 1;
//        A(1,0) = 1; A(1,1) = 2; A(1,2) = -1;

//        DenseVector b; b.setZero(2); b(1) = -1e-5;

//        ConstraintPtr c(new ConstraintLinear(A,b,false));

//        std::vector<int> vars = {1,2,3};
//        cs->add(c,vars);
//    }

//    DenseMatrix cobj(1,dim); cobj.setZero();
//    cobj(0,0) = 1;
//    ObjectivePtr obj(new ObjectiveLinear(cobj));

//    std::vector<int> vt(dim,CONTINUOUS);
//    std::vector<int> bv = {0,1};
//    std::vector<double> z0(dim,0);

//    BB::BranchAndBound bnb(obj, cs, z0, vt, bv);
//    Timer timer;
//    timer.start();
//    bnb.optimize();
//    timer.stop();
//    cout << "Time: " << timer.getMilliSeconds() << " (ms)" << endl;
//}

//void P04()
//{
//    // Problem 4 (Nataraj)
//    cout << "\n\nSolving problem P04..." << endl;

//    int dim = 3+1;

//    // x1,x2,x3,l1, l1 >= 0
//    std::vector<double> lb = {0,0,0,0};
//    std::vector<double> ub = {2,10,3,INF};

//    ConstraintCompositePtr cs(new ConstraintComposite(dim, lb, ub));

//    { // Quadratic constraint = l1
//        std::vector<double> thislb = {lb.at(0), lb.at(1), lb.at(2)};
//        std::vector<double> thisub = {ub.at(0), ub.at(1), ub.at(2)};

//        std::vector<int> deg = {2,2,2};

//        DenseVector c(27);
//        c.setZero();
//        c(0) = 24; // Constant
//        c(1) = -13; // x3
//        c(2) = 2; // x3^2
//        c(3) = 9; // x2
//        c(4) = -2; // x2*x3
//        c(6) = 2; // x2^2
//        c(9) = -20; // x1
//        c(10) = 4; // x1*x3
//        c(12) = -4; // x1*x2
//        c(18) = 4; // x1^2

//        DenseMatrix T = getTransformationMatrix(deg, thislb, thisub);
//        DenseMatrix coeffs = T*c;

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {0,1,2,3};

//        cs->add(cbs,vars);
//    }

//    { // Linear constraints

//        DenseMatrix A(2,3);
//        A.setZero();
//        A(0,0) = 1; A(0,1) = 1; A(0,2) = 1;
//        A(1,1) = 3; A(1,2) = 1;

//        DenseVector b;
//        b.setZero(2);
//        b(0) = 4;
//        b(1) = 6;

//        ConstraintPtr c(new ConstraintLinear(A,b,false));

//        std::vector<int> vars = {0,1,2};
//        cs->add(c,vars);
//    }

//    DenseMatrix cobj(1,dim); cobj.setZero();
//    cobj(0,0) = -2;
//    cobj(0,1) = 1;
//    cobj(0,2) = -1;
//    ObjectivePtr obj(new ObjectiveLinear(cobj));

//    std::vector<int> vt(dim,CONTINUOUS);
//    std::vector<int> bv = {0,1,2};
//    std::vector<double> z0(dim,0);

//    BB::BranchAndBound bnb(obj, cs, z0, vt, bv);
//    Timer timer;
//    timer.start();
//    bnb.optimize();
//    timer.stop();
//    cout << "Time: " << timer.getMilliSeconds() << " (ms)" << endl;
//}

//void P05()
//{
//    // Problem 5 (Nataraj), Himmelblau root finding problem
//    cout << "\n\nSolving problem P05..." << endl;

//    int dim = 3+2;

//    // x1,x2,x3,l1,l2
//    std::vector<double> lb = {-5,-5,-5,-INF,-INF};
//    std::vector<double> ub = {5,5,5,INF,INF};

//    ConstraintCompositePtr cs(new ConstraintComposite(dim, lb, ub));

//    { // f1(x1,x2) = 2*x1^2 + 4*x1*x2 - 42*x1 +4*x1^3 = l1
//        std::vector<double> thislb = {lb.at(0), lb.at(1)};
//        std::vector<double> thisub = {ub.at(0), ub.at(1)};

//        std::vector<int> deg = {3,1};

//        DenseVector c(8);
//        c.setZero();
//        c(2) = -42; // x1
//        c(3) = 4; // x1*x2
//        c(4) = 2; // x1^2
//        c(6) = 4; // x1^3

//        DenseMatrix T = getTransformationMatrix(deg, thislb, thisub);
//        DenseMatrix coeffs = T*c;

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {0,1,3};

//        cs->add(cbs,vars);
//    }

//    { // f2(x1,x2) = 2*x1^2 + 4*x1*x2 - 26*x2 +4*x2^3 = l2
//        std::vector<double> thislb = {lb.at(0), lb.at(1)};
//        std::vector<double> thisub = {ub.at(0), ub.at(1)};

//        std::vector<int> deg = {2,3};

//        DenseVector c(12);
//        c.setZero();
//        c(1) = -26; // x2
//        c(3) = 4; // x2^3
//        c(5) = 4; // x1*x2
//        c(8) = 2; // x1^2

//        DenseMatrix T = getTransformationMatrix(deg, thislb, thisub);
//        DenseMatrix coeffs = T*c;

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {0,1,4};

//        cs->add(cbs,vars);
//    }

//    { // Linear constraints

//        DenseMatrix A(4,5);
//        A.setZero();
//        A(0,2) = -1; A(0,3) = 1;
//        A(1,2) = -1; A(1,3) = -1;
//        A(2,2) = -1; A(2,4) = 1;
//        A(3,2) = -1; A(3,4) = -1;

//        DenseVector b;
//        b.setZero(4);
//        b(0) = 14;
//        b(1) = -14;
//        b(2) = 22;
//        b(3) = -22;

//        ConstraintPtr c(new ConstraintLinear(A,b,false));
//        cs->add(c);
//    }

//    DenseMatrix cobj(1,dim);
//    cobj.setZero();
//    cobj(0,2) = 1;
//    ObjectivePtr obj(new ObjectiveLinear(cobj));

//    std::vector<int> vt(dim,CONTINUOUS);
//    std::vector<int> bv = {0,1};
//    std::vector<double> z0(dim,0);

//    BB::BranchAndBound bnb(obj, cs, z0, vt, bv);
//    Timer timer;
//    timer.start();
//    bnb.optimize();
//    timer.stop();
//    cout << "Time: " << timer.getMilliSeconds() << " (ms)" << endl;
//}

//void P06()
//{
//    // Problem 6 (Nataraj)
//    // 5 2-D B-splines
//    int dim = 4+2;

//    // x1,x2,x3,x4,l1,l2
//    std::vector<double> lb = {1,0.625,47.5,90,-INF,-INF};
//    std::vector<double> ub = {1.1375,1,52.5,112,INF,INF};

//    // x1,x2,x3,x4,l1,l2,l3,l4,l5
////    std::vector<double> lb = {1,0.625,47.5,90,-INF,-INF,-INF,-INF,-INF};
////    std::vector<double> ub = {1.1375,1,52.5,112,INF,INF,INF,INF,INF};

//    ConstraintCompositePtr cs(new ConstraintComposite(dim, lb, ub));

//    { // obj = l1
//        std::vector<double> thislb = {lb.at(0), lb.at(1), lb.at(2), lb.at(3)};
//        std::vector<double> thisub = {ub.at(0), ub.at(1), ub.at(2), ub.at(3)};

//        std::vector<int> deg = {2,1,2,1};

//        // Poly coeffs
//        DenseVector c(4);
//        c.setZero();
//        c(0) = 0.6224;
//        c(1) = 1.7781;
//        c(2) = 3.1661;
//        c(3) = 19.84;

//        // Poly exponents
//        DenseMatrix E(4,4);
//        E.setZero();
//        E(0,2) = 1; E(0,3) = 1;
//        E(1,1) = 1; E(1,2) = 2;
//        E(2,0) = 2; E(2,3) = 1;
//        E(3,0) = 2; E(3,2) = 1;

//        DenseMatrix coeffs = getBsplineCoefficients(c, E, thislb, thisub);

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {0,1,2,3,4};

//        cs->add(cbs,vars);

//        DenseMatrix cpoints = bs.getControlPoints();
//        cout << cpoints << endl;
//    }

//    { // noncon = l2
//        std::vector<double> thislb = {lb.at(2), lb.at(3)};
//        std::vector<double> thisub = {ub.at(2), ub.at(3)};

//        std::vector<int> deg = {3,1};

//        // Poly coeffs
//        DenseVector c(2);
//        c.setZero();
//        c(0) = -1;
//        c(1) = -(4.0/3.0);

//        // Poly exponents
//        DenseMatrix E(2,2);
//        E.setZero();
//        E(0,0) = 2; E(0,1) = 1;
//        E(1,0) = 3;

//        DenseMatrix coeffs = getBsplineCoefficients(c, E, thislb, thisub);

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {2,3,5};

//        cs->add(cbs,vars);
//    }

//    { // Linear constraints of auxiliary variables
//        DenseMatrix A(4,6);
//        A.setZero();
//        A(0,0) = -1; A(0,2) = 0.0193;
//        A(1,1) = -1; A(1,2) = 0.00954;
//        A(2,5) = 1;
//        A(3,3) = 1;

//        DenseVector b;
//        b.setZero(4);
//        b(2) = -750.1728/3.14159265359;
//        b(3) = 240;

//        ConstraintPtr c(new ConstraintLinear(A,b,false));
//        cs->add(c);
//    }

//    DenseMatrix cobj(1,dim);
//    cobj.setZero();
//    cobj(0,4) = 1;
//    ObjectivePtr obj(new ObjectiveLinear(cobj));

//    std::vector<int> vt(dim,CONTINUOUS);
//    std::vector<int> bv = {0,1,2,3};
//    std::vector<double> z0(dim,0);

//    BB::BranchAndBound bnb(obj, cs, z0, vt, bv);
//    Timer timer;
//    timer.start();
//    bnb.optimize();
//    timer.stop();
//    cout << "Time: " << timer.getMilliSeconds() << " (ms)" << endl;
//    cout << "Time: " << timer.getMicroSeconds() << " (ms)" << endl;

//    //1 0.625 47.5 90 6395.50783 -345958.333
//}

//void P07()
//{
//    // Problem 7 (Nataraj)
//    // Two 2-D B-splines
//    cout << "\n\nSolving problem P07..." << endl;

//    int dim = 4+2;

//    // x1,x2,x3,x4,l1,l2
//    std::vector<double> lb = {0,0,0,0,-INF,-INF};
//    std::vector<double> ub = {5,5,5,5,INF,INF};

//    ConstraintCompositePtr cs(new ConstraintComposite(dim, lb, ub));

//    { // x1^4*x2^4 - x1^4 = l1
//        std::vector<double> thislb = {lb.at(0), lb.at(1)};
//        std::vector<double> thisub = {ub.at(0), ub.at(1)};

//        std::vector<int> deg = {4,4};

//        DenseVector c(5*5); c.setZero();
//        c(20) = -1;
//        c(24) = 1;

//        DenseMatrix T = getTransformationMatrix(deg, thislb, thisub);
//        DenseMatrix coeffs = T*c;

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {0,1,4};

//        cs->add(cbs,vars);
//    }

//    { // x2^4*x3 = l1
//        std::vector<double> thislb = {lb.at(1), lb.at(2)};
//        std::vector<double> thisub = {ub.at(1), ub.at(2)};

//        std::vector<int> deg = {4,1};

//        DenseVector c(10); c.setZero();
//        c(9) = 1;

//        DenseMatrix T = getTransformationMatrix(deg, thislb, thisub);
//        DenseMatrix coeffs = T*c;

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {1,2,5};

//        cs->add(cbs,vars);
//    }

//    { // Linear constraints of auxiliary variables
//        DenseMatrix A(1,2);
//        A(0,0) = 1; A(0,1) = -1;

//        DenseVector b; b.setZero(1);

//        ConstraintPtr c(new ConstraintLinear(A,b,true));

//        std::vector<int> vars = {4,5};
//        cs->add(c,vars);
//    }

//    { // Linear constraints

//        DenseMatrix A(6,4);
//        A.setZero();
//        A(0,0) = -1; A(0,3) = -0.25;
//        A(1,0) = 1; A(1,3) = -0.25;
//        A(2,1) = -1; A(2,3) = -0.2;
//        A(3,1) = 1; A(3,3) = -0.2;
//        A(4,2) = -1; A(4,3) = -0.2;
//        A(5,2) = 1; A(5,3) = -0.2;

//        DenseVector b; b.setZero(6);
//        b(0) = -1.4;
//        b(1) = 1.4;
//        b(2) = -1.5;
//        b(3) = 1.5;
//        b(4) = -0.8;
//        b(5) = 0.8;

//        ConstraintPtr c(new ConstraintLinear(A,b,false));

//        std::vector<int> vars = {0,1,2,3};
//        cs->add(c,vars);
//    }

//    DenseMatrix cobj(1,dim); cobj.setZero();
//    cobj(0,3) = 1;
//    ObjectivePtr obj(new ObjectiveLinear(cobj));

//    std::vector<int> vt(dim,CONTINUOUS);
//    std::vector<int> bv = {0,1,2};
//    std::vector<double> z0(dim,0);

//    BB::BranchAndBound bnb(obj, cs, z0, vt, bv);
//    Timer timer;
//    timer.start();
//    bnb.optimize();
//    timer.stop();
//    cout << "Time: " << timer.getMilliSeconds() << " (ms)" << endl;
//}

//void P07_2()
//{
//    // Problem 7 (Nataraj)
//    // One 3-D B-spline

//    int dim = 4+1;

//    // x1,x2,x3,x4,l1,l2
//    std::vector<double> lb = {0,0,0,0,-INF};
//    std::vector<double> ub = {5,5,5,5,INF};

//    ConstraintCompositePtr cs(new ConstraintComposite(dim, lb, ub));

//    { // x1^4*x2^4 - x1^4 = l1
//        std::vector<double> thislb = {lb.at(0), lb.at(1), lb.at(2)};
//        std::vector<double> thisub = {ub.at(0), ub.at(1), ub.at(2)};

//        std::vector<int> deg = {4,4,1};

//        // Poly coeffs
//        DenseVector c(3); c.setZero();
//        c(0) = 1;
//        c(1) = -1;
//        c(2) = -1;

//        // Poly exponents
//        DenseMatrix E(3,3); E.setZero();
//        E(0,0) = 4; E(0,1) = 4; E(0,2) = 0;
//        E(1,0) = 4; E(1,1) = 0; E(1,2) = 0;
//        E(2,0) = 0; E(2,1) = 4; E(2,2) = 1;

//        DenseMatrix coeffs = getBsplineCoefficients(c, E, thislb, thisub);

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {0,1,2,4};

//        cs->add(cbs,vars);
//    }

//    { // Linear constraints of auxiliary variables
//        DenseMatrix A(1,1);
//        A(0,0) = 1;

//        DenseVector b; b.setZero(1);

//        ConstraintPtr c(new ConstraintLinear(A,b,true));

//        std::vector<int> vars = {4};
//        cs->add(c,vars);
//    }

//    { // Linear constraints

//        DenseMatrix A(6,4);
//        A.setZero();
//        A(0,0) = -1; A(0,3) = -0.25;
//        A(1,0) = 1; A(1,3) = -0.25;
//        A(2,1) = -1; A(2,3) = -0.2;
//        A(3,1) = 1; A(3,3) = -0.2;
//        A(4,2) = -1; A(4,3) = -0.2;
//        A(5,2) = 1; A(5,3) = -0.2;

//        DenseVector b; b.setZero(6);
//        b(0) = -1.4;
//        b(1) = 1.4;
//        b(2) = -1.5;
//        b(3) = 1.5;
//        b(4) = -0.8;
//        b(5) = 0.8;

//        ConstraintPtr c(new ConstraintLinear(A,b,false));

//        std::vector<int> vars = {0,1,2,3};
//        cs->add(c,vars);
//    }

//    DenseMatrix cobj(1,dim); cobj.setZero();
//    cobj(0,3) = 1;
//    ObjectivePtr obj(new ObjectiveLinear(cobj));

//    std::vector<int> vt(dim,CONTINUOUS);
//    std::vector<int> bv = {0,1,2};
//    std::vector<double> z0(dim,0);

//    BB::BranchAndBound bnb(obj, cs, z0, vt, bv);
//    Timer timer;
//    timer.start();
//    bnb.optimize();
//    timer.stop();
//    cout << "Time: " << timer.getMilliSeconds() << " (ms)" << endl;
//}

//void P08()
//{
//    // Problem 8 (Nataraj)
//    // Three 4-D B-splines
//    // Runs fast when refinement is turned off!
//    // A lot of iterations when using bounding box relaxation (as expected)
//    cout << "\n\nSolving problem P08..." << endl;

//    int dim = 4+3;

//    // x1,x2,x3,x4,l1,l2,l3
//    std::vector<double> lb = {3,2,0.125,0.25,-INF,-INF,-INF};
//    std::vector<double> ub = {20,15,0.75,1.25,INF,INF,INF};

//    ConstraintCompositePtr cs(new ConstraintComposite(dim, lb, ub));

//    { // f(x) = l3
//        std::vector<double> thislb = {lb.at(0), lb.at(1), lb.at(2), lb.at(3)};
//        std::vector<double> thisub = {ub.at(0), ub.at(1), ub.at(2), ub.at(3)};

//        std::vector<int> deg = {1,1,1,1};

//        // Poly coeffs
//        double a = 27.264;
//        DenseVector c(3); c.setZero();
//        c(0) = 2*a;
//        c(1) = a;
//        c(2) = -2*a;

//        // Poly exponents
//        DenseMatrix E(3,4); E.setZero();
//        E(0,0) = 0; E(0,1) = 1; E(0,2) = 0; E(0,3) = 1;
//        E(1,0) = 1; E(1,1) = 0; E(1,2) = 1; E(1,3) = 0;
//        E(2,0) = 0; E(2,1) = 0; E(2,2) = 1; E(2,3) = 1;

//        DenseMatrix coeffs = getBsplineCoefficients(c, E, thislb, thisub);

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {0,1,2,3,6};

//        cs->add(cbs,vars);
//    }

//    { // I1(x) = l1
//        std::vector<double> thislb = {lb.at(0), lb.at(1), lb.at(2), lb.at(3)};
//        std::vector<double> thisub = {ub.at(0), ub.at(1), ub.at(2), ub.at(3)};

//        std::vector<int> deg = {3,1,3,1};

//        // Poly coeffs
//        DenseVector c(7); c.setZero();
//        c(0) = 6;
//        c(1) = -12;
//        c(2) = 8;
//        c(3) = 1;
//        c(4) = -6;
//        c(5) = 12;
//        c(6) = -8;

//        // Poly exponents
//        DenseMatrix E(7,4); E.setZero();
//        E(0,0) = 2; E(0,1) = 1; E(0,2) = 1; E(0,3) = 0;
//        E(1,0) = 1; E(1,1) = 1; E(1,2) = 2; E(1,3) = 0;
//        E(2,0) = 0; E(2,1) = 1; E(2,2) = 3; E(2,3) = 0;
//        E(3,0) = 3; E(3,1) = 0; E(3,2) = 0; E(3,3) = 1;
//        E(4,0) = 2; E(4,1) = 0; E(4,2) = 1; E(4,3) = 1;
//        E(5,0) = 1; E(5,1) = 0; E(5,2) = 2; E(5,3) = 1;
//        E(6,0) = 0; E(6,1) = 0; E(6,2) = 3; E(6,3) = 1;

//        DenseMatrix coeffs = getBsplineCoefficients(c, E, thislb, thisub);

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {0,1,2,3,4};

//        cs->add(cbs,vars);
//    }

//    { // I2(x) = l2
//        std::vector<double> thislb = {lb.at(0), lb.at(1), lb.at(2), lb.at(3)};
//        std::vector<double> thisub = {ub.at(0), ub.at(1), ub.at(2), ub.at(3)};

//        std::vector<int> deg = {3,1,4,2};

//        // Poly coeffs
//        double a = -3.5;
//        DenseVector c(12); c.setZero();
//        c(0) = 6*a;
//        c(1) = -12*a;
//        c(2) = 8*a;
//        c(3) = 1*a;
//        c(4) = -6*a;
//        c(5) = 12*a;
//        c(6) = -8*a;

//        c(7) = 1;
//        c(8) = -1;
//        c(9) = 1;
//        c(10) = 1;
//        c(11) = -2;

//        // Poly exponents
//        DenseMatrix E(12,4); E.setZero();
//        E(0,0) = 2; E(0,1) = 1; E(0,2) = 2; E(0,3) = 0;
//        E(1,0) = 1; E(1,1) = 1; E(1,2) = 3; E(1,3) = 0;
//        E(2,0) = 0; E(2,1) = 1; E(2,2) = 4; E(2,3) = 0;
//        E(3,0) = 3; E(3,1) = 0; E(3,2) = 1; E(3,3) = 1;
//        E(4,0) = 2; E(4,1) = 0; E(4,2) = 2; E(4,3) = 1;
//        E(5,0) = 1; E(5,1) = 0; E(5,2) = 3; E(5,3) = 1;
//        E(6,0) = 0; E(6,1) = 0; E(6,2) = 4; E(6,3) = 1;

//        E(7,0) = 1;  E(7,1) = 1;  E(7,2) = 0;  E(7,3) = 1;
//        E(8,0) = 0;  E(8,1) = 1;  E(8,2) = 0;  E(8,3) = 2;
//        E(9,0) = 2;  E(9,1) = 0;  E(9,2) = 1;  E(9,3) = 0;
//        E(10,0) = 0; E(10,1) = 0; E(10,2) = 1; E(10,3) = 2;
//        E(11,0) = 1; E(11,1) = 0; E(11,2) = 1; E(11,3) = 1;

//        DenseMatrix coeffs = getBsplineCoefficients(c, E, thislb, thisub);

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {0,1,2,3,5};

//        cs->add(cbs,vars);
//    }

//    { // Linear constraints of auxiliary variables
//        DenseMatrix A = DenseMatrix::Zero(7,6);
//        A(0,4) = -1;
//        A(1,0) = 8; A(1,4) = -1;
//        A(2,5) = 1;
//        A(3,0) = 1; A(3,1) = -3;
//        A(4,1) = 2; A(4,0) = -1;
//        A(5,2) = 1; A(5,3) = -1.5;
//        A(6,3) = 0.5; A(6,2) = -1;

//        DenseVector b; b.setZero(7);
//        b(0) = -61.01627586;

//        ConstraintPtr c(new ConstraintLinear(A,b,false));

//        std::vector<int> vars = {0,1,2,3,4,5};
//        cs->add(c,vars);
//    }

//    DenseMatrix cobj(1,dim); cobj.setZero();
//    cobj(0,6) = 1;
//    ObjectivePtr obj(new ObjectiveLinear(cobj));

//    std::vector<int> vt(dim,CONTINUOUS);
//    std::vector<int> bv = {0,1,2,3};
//    std::vector<double> z0(dim,0);


////    DenseVector x(7);
////    x(0) = 3.76117027;
////    x(1) = 2;
////    x(2) = 0.125;
////    x(3) = 0.25;
////    x(4) = 30.662118;
////    x(5) = -10.1180522;
////    x(6) = 38.3780683;

////    bool feas = cs->checkFeasibility(x, 1e-5);
////    if (feas) cout << "YAAAI!" << endl;
////    else cout << "NEEEY!" << endl;
////    exit(1);

//    BB::BranchAndBound bnb(obj, cs, z0, vt, bv);
//    Timer timer;
//    timer.start();
//    bnb.optimize();
//    timer.stop();
//    cout << "Time: " << timer.getMilliSeconds() << " (ms)" << endl;
//}

//void P09()
//{
//    // Problem 9 (Nataraj), QP
//    // Five 1-D B-splines
//    cout << "\n\nSolving problem P09..." << endl;

//    int dim = 5+5;

//    // x1,x2,x3,x4,l1,l2,l3
//    std::vector<double> lb = {0,0,0,0,0,-INF,-INF,-INF,-INF,-INF};
//    std::vector<double> ub = {1,1,1,1,1,INF,INF,INF,INF,INF};

//    ConstraintCompositePtr cs(new ConstraintComposite(dim, lb, ub));

//    std::vector<double> a = {42,44,45,47,47.5};

//    // Add one B-spline for each variable
//    for (int i = 0; i < 5; i++)
//    {
//        std::vector<double> thislb = {lb.at(i)};
//        std::vector<double> thisub = {ub.at(i)};

//        std::vector<int> deg = {2};

//        // Poly coeffs
//        DenseVector c(3);
//        c.setZero();
//        c(0) = 0;
//        c(1) = a.at(i);
//        c(2) = -50; // -50 or -0.5 (Floudas' problem has -50)

//        DenseMatrix T = getTransformationMatrix(deg, thislb, thisub);
//        DenseMatrix coeffs = T*c;

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {i,i+5};

//        cs->add(cbs,vars);
//    }

//    { // Linear constraints
//        DenseMatrix A = DenseMatrix::Zero(1,5);
//        A(0,0) = 20;
//        A(0,1) = 12;
//        A(0,2) = 11;
//        A(0,3) = 7;
//        A(0,4) = 4;

//        DenseVector b;
//        b.setZero(1);
//        b(0) = 40;

//        ConstraintPtr c(new ConstraintLinear(A,b,false));

//        std::vector<int> vars = {0,1,2,3,4};
//        cs->add(c,vars);
//    }

//    DenseMatrix cobj(1,dim);
//    cobj.setZero();
//    cobj(0,5) = 1;
//    cobj(0,6) = 1;
//    cobj(0,7) = 1;
//    cobj(0,8) = 1;
//    cobj(0,9) = 1;
//    ObjectivePtr obj(new ObjectiveLinear(cobj));

//    std::vector<int> vt(dim,CONTINUOUS);
//    std::vector<int> bv = {0,1,2,3,4};
//    std::vector<double> z0(dim,0);

//    BB::BranchAndBound bnb(obj, cs, z0, vt, bv);
//    Timer timer;
//    timer.start();
//    bnb.optimize();
//    timer.stop();
//    cout << "Time: " << timer.getMilliSeconds() << " (ms)" << endl;
//}

//void P10()
//{
//    // Problem 10 (Nataraj), QP
//    // Five 1-D B-splines
//    cout << "\n\nSolving problem P10..." << endl;

//    int dim = 6+5; // x1,..,x5,y,l1,...,l5

//    // x1,x2,x3,x4,l1,l2,l3
//    std::vector<double> lb = {0,0,0,0,0,0,-INF,-INF,-INF,-INF,-INF};
//    std::vector<double> ub = {1,1,1,1,1,INF,INF,INF,INF,INF,INF};

//    ConstraintCompositePtr cs(new ConstraintComposite(dim, lb, ub));

//    std::vector<double> a = {-10.5,-7.5,-3.5,-2.5,-1.5};

//    // Add one B-spline for each variable
//    for (int i = 0; i < 5; i++)
//    {
//        std::vector<double> thislb = {lb.at(i)};
//        std::vector<double> thisub = {ub.at(i)};

//        std::vector<int> deg = {2};

//        // Poly coeffs
//        DenseVector c(3);
//        c.setZero();
//        c(0) = 0;
//        c(1) = a.at(i);
//        c(2) = -0.5; // -50 or -0.5 (Floudas' problem has -50)

//        DenseMatrix T = getTransformationMatrix(deg, thislb, thisub);
//        DenseMatrix coeffs = T*c;

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {i,i+6};

//        cs->add(cbs,vars);
//    }

//    { // Linear constraints
//        DenseMatrix A = DenseMatrix::Zero(2,6);
//        A(0,0) = 6; A(0,1) = 3; A(0,2) = 3; A(0,3) = 2; A(0,4) = 1;
//        A(1,0) = 10; A(1,2) = 10; A(1,5) = 1;

//        DenseVector b;
//        b.setZero(2);
//        b(0) = 6.5;
//        b(1) = 20;

//        ConstraintPtr c(new ConstraintLinear(A,b,false));

//        std::vector<int> vars = {0,1,2,3,4,5};
//        cs->add(c,vars);
//    }

//    DenseMatrix cobj(1,dim);
//    cobj.setZero();
//    cobj(0,5) = -10;
//    cobj(0,6) = 1;
//    cobj(0,7) = 1;
//    cobj(0,8) = 1;
//    cobj(0,9) = 1;
//    cobj(0,10) = 1;
//    ObjectivePtr obj(new ObjectiveLinear(cobj));

//    std::vector<int> vt(dim,CONTINUOUS);
//    std::vector<int> bv = {0,1,2,3,4};
//    std::vector<double> z0(dim,0);

//    BB::BranchAndBound bnb(obj, cs, z0, vt, bv);
//    Timer timer;
//    timer.start();
//    bnb.optimize();
//    timer.stop();
//    cout << "Time: " << timer.getMilliSeconds() << " (ms)" << endl;
//    cout << "Time: " << timer.getMicroSeconds() << " (us)" << endl;
//}

//void P11()
//{
//    // Problem 11 (Nataraj) - Test problem 3, Ch. 3.3 (Floudas)
//    cout << "\n\nSolving problem P11..." << endl;

//    int dim = 14;

//    // x1,...,x6,l1,...,l8
//    std::vector<double> lb = {0,0,1,0,1,0,0,0,0,0,0,0,0,0};
//    std::vector<double> ub = {5,5,5,6,5,10,INF,INF,INF,INF,INF,INF,INF,INF};

//    ConstraintCompositePtr cs(new ConstraintComposite(dim, lb, ub));

//    { // (x1-2)^2
//        int var = 0;
//        std::vector<double> thislb = {lb.at(var)};
//        std::vector<double> thisub = {ub.at(var)};

//        std::vector<int> deg = {2};

//        DenseVector c(3); c.setZero();
//        c(0) = 4;
//        c(1) = -4;
//        c(2) = 1;

//        DenseMatrix T = getTransformationMatrix(deg,thislb,thisub);
//        DenseMatrix coeffs = T*c;

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {var,6};

//        cs->add(cbs,vars);
//    }

//    { // (x2-2)^2
//        int var = 1;
//        std::vector<double> thislb = {lb.at(var)};
//        std::vector<double> thisub = {ub.at(var)};

//        std::vector<int> deg = {2};

//        DenseVector c(3); c.setZero();
//        c(0) = 4;
//        c(1) = -4;
//        c(2) = 1;

//        DenseMatrix T = getTransformationMatrix(deg,thislb,thisub);
//        DenseMatrix coeffs = T*c;

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {var,7};

//        cs->add(cbs,vars);
//    }

//    { // (x3-1)^2
//        int var = 2;
//        std::vector<double> thislb = {lb.at(var)};
//        std::vector<double> thisub = {ub.at(var)};

//        std::vector<int> deg = {2};

//        DenseVector c(3); c.setZero();
//        c(0) = 1;
//        c(1) = -2;
//        c(2) = 1;

//        DenseMatrix T = getTransformationMatrix(deg,thislb,thisub);
//        DenseMatrix coeffs = T*c;

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {var,8};

//        cs->add(cbs,vars);
//    }

//    { // (x4-4)^2
//        int var = 3;
//        std::vector<double> thislb = {lb.at(var)};
//        std::vector<double> thisub = {ub.at(var)};

//        std::vector<int> deg = {2};

//        DenseVector c(3); c.setZero();
//        c(0) = 16;
//        c(1) = -8;
//        c(2) = 1;

//        DenseMatrix T = getTransformationMatrix(deg,thislb,thisub);
//        DenseMatrix coeffs = T*c;

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {var,9};

//        cs->add(cbs,vars);
//    }

//    { // (x5-1)^2
//        int var = 4;
//        std::vector<double> thislb = {lb.at(var)};
//        std::vector<double> thisub = {ub.at(var)};

//        std::vector<int> deg = {2};

//        DenseVector c(3); c.setZero();
//        c(0) = 1;
//        c(1) = -2;
//        c(2) = 1;

//        DenseMatrix T = getTransformationMatrix(deg,thislb,thisub);
//        DenseMatrix coeffs = T*c;

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {var,10};

//        cs->add(cbs,vars);
//    }

//    { // (x6-4)^2
//        int var = 5;
//        std::vector<double> thislb = {lb.at(var)};
//        std::vector<double> thisub = {ub.at(var)};

//        std::vector<int> deg = {2};

//        DenseVector c(3); c.setZero();
//        c(0) = 16;
//        c(1) = -8;
//        c(2) = 1;

//        DenseMatrix T = getTransformationMatrix(deg,thislb,thisub);
//        DenseMatrix coeffs = T*c;

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {var,11};

//        cs->add(cbs,vars);
//    }

//    { // (x3-3)^2
//        int var = 2;
//        std::vector<double> thislb = {lb.at(var)};
//        std::vector<double> thisub = {ub.at(var)};

//        std::vector<int> deg = {2};

//        DenseVector c(3); c.setZero();
//        c(0) = 9;
//        c(1) = -6;
//        c(2) = 1;

//        DenseMatrix T = getTransformationMatrix(deg,thislb,thisub);
//        DenseMatrix coeffs = T*c;

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {var,12};

//        cs->add(cbs,vars);
//    }

//    { // (x5-3)^2
//        int var = 4;
//        std::vector<double> thislb = {lb.at(var)};
//        std::vector<double> thisub = {ub.at(var)};

//        std::vector<int> deg = {2};

//        DenseVector c(3); c.setZero();
//        c(0) = 9;
//        c(1) = -6;
//        c(2) = 1;

//        DenseMatrix T = getTransformationMatrix(deg,thislb,thisub);
//        DenseMatrix coeffs = T*c;

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {var,13};

//        cs->add(cbs,vars);
//    }

//    { // l7 + x4 >= 4

//        DenseMatrix A(1,2); A(0,0) = -1; A(0,1) = -1;

//        DenseVector b(1); b(0) = -4;

//        ConstraintPtr c(new ConstraintLinear(A,b,false));

//        std::vector<int> vars = {12,3};
//        cs->add(c,vars);
//    }

//    { // l8 + x6 >= 4

//        DenseMatrix A(1,2); A(0,0) = -1; A(0,1) = -1;

//        DenseVector b(1); b(0) = -4;

//        ConstraintPtr c(new ConstraintLinear(A,b,false));

//        std::vector<int> vars = {13,5};
//        cs->add(c,vars);
//    }

//    { // x1 - 3x2 <= 2, -x1 + x2 <= 2, x1 + x2 <= 6, x1 + x2 >= 2

//        DenseMatrix A(4,2);
//        A(0,0) = 1; A(0,1) = -3;
//        A(1,0) = -1; A(1,1) = 1;
//        A(2,0) = 1; A(2,1) = 1;
//        A(3,0) = -1; A(3,1) = -1;

//        DenseVector b(4);
//        b(0) = 2;
//        b(1) = 2;
//        b(2) = 6;
//        b(3) = -2;

//        ConstraintPtr c(new ConstraintLinear(A,b,false));

//        std::vector<int> vars = {0,1};
//        cs->add(c,vars);
//    }

//    DenseMatrix cobj(1,dim); cobj.setZero();
//    cobj(0,6) = -25;
//    cobj(0,7) = -1;
//    cobj(0,8) = -1;
//    cobj(0,9) = -1;
//    cobj(0,10) = -1;
//    cobj(0,11) = -1;
//    ObjectivePtr obj(new ObjectiveLinear(cobj));

//    std::vector<int> vt(dim,CONTINUOUS);
//    std::vector<int> bv = {0,1,2,3,4,5};
//    std::vector<double> z0(dim,0);

//    BB::BranchAndBound bnb(obj, cs, z0, vt, bv);
//    Timer timer;
//    timer.start();
//    bnb.optimize();
//    timer.stop();
//    cout << "Time: " << timer.getMilliSeconds() << " (ms)" << endl;
//}

//void P12()
//{
//    // Problem 12 (Nataraj), Bilinear problem
//    // Six 2-D B-splines
//    cout << "\n\nSolving problem P12..." << endl;

//    int dim = 7+6; // x1,..,x7,l1,...,l6

//    std::vector<double> lb = {0,0,0,0,0,0,0,            0,0,0,0,0,0};
//    std::vector<double> ub = {1,1,1,100,200,100,200,    100,100,100,200,200,200};

//    ConstraintCompositePtr cs(new ConstraintComposite(dim, lb, ub));

//    // x_i * x_4 = li
//    for (int i = 0; i < 3; i++)
//    {
//        std::vector<double> thislb = {lb.at(i), lb.at(3)};
//        std::vector<double> thisub = {ub.at(i), ub.at(3)};

//        std::vector<int> deg = {1,1};

//        // Poly coeffs
//        DenseVector c(4);
//        c.setZero();
//        c(3) = 1;

//        DenseMatrix T = getTransformationMatrix(deg, thislb, thisub);
//        DenseMatrix coeffs = T*c;

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {i,3,i+7};

//        cs->add(cbs,vars);
//    }

//    // x_i * x_5 = li
//    for (int i = 0; i < 3; i++)
//    {
//        std::vector<double> thislb = {lb.at(i), lb.at(4)};
//        std::vector<double> thisub = {ub.at(i), ub.at(4)};

//        std::vector<int> deg = {1,1};

//        // Poly coeffs
//        DenseVector c(4);
//        c.setZero();
//        c(3) = 1;

//        DenseMatrix T = getTransformationMatrix(deg, thislb, thisub);
//        DenseMatrix coeffs = T*c;

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {i,4,i+7+3};

//        cs->add(cbs,vars);
//    }

//    { // Linear inequality constraints
//        DenseMatrix A = DenseMatrix::Zero(5,dim);
//        A(0,9) = 1; A(0,12) = 1;
//        A(1,3) = 1; A(1,5) = 1;
//        A(2,4) = 1; A(2,6) = 1;
//        A(3,7) = 3; A(3,8) = 1; A(3,9) = 1; A(3,3) = -2.5; A(3,5) = -0.5;
//        A(4,10) = 3; A(4,11) = 1; A(4,12) = 1; A(4,4) = -1.5; A(4,6) = 0.5;


//        DenseVector b;
//        b.setZero(5);
//        b(0) = 50;
//        b(1) = 100;
//        b(2) = 200;

//        ConstraintPtr c(new ConstraintLinear(A,b,false));

//        //std::vector<int> vars = {0,1,2,3,4,5};
//        cs->add(c); // All variables
//    }

//    { // Linear equality constraints
//        DenseMatrix A = DenseMatrix::Zero(1,3);
//        A(0,0) = 1; A(0,1) = 1; A(0,2) = 1;

//        DenseVector b(1);
//        b(0) = 1;

//        ConstraintPtr c(new ConstraintLinear(A,b,true));

//        std::vector<int> vars = {0,1,2};
//        cs->add(c,vars);
//    }

//    DenseMatrix cobj(1,dim);
//    cobj.setZero();
//    cobj(0,3) = -9;
//    cobj(0,7) = 6;
//    cobj(0,8) = 16;
//    cobj(0,9) = 15;
//    cobj(0,4) = -15;
//    cobj(0,10) = 6;
//    cobj(0,11) = 16;
//    cobj(0,12) = 15;
//    cobj(0,5) = 1;
//    cobj(0,6) = -5;

//    ObjectivePtr obj(new ObjectiveLinear(cobj));

//    std::vector<int> vt(dim,CONTINUOUS);
//    std::vector<int> bv = {0,1,2,3,4};
//    std::vector<double> z0(dim,0);

//    BB::BranchAndBound bnb(obj, cs, z0, vt, bv);
//    Timer timer;
//    timer.start();
//    bnb.optimize();
//    timer.stop();
//    cout << "Time: " << timer.getMilliSeconds() << " (ms)" << endl;
//}

//void P12_2()
//{
//    // Problem 12 (Nataraj), Bilinear problem
//    cout << "\n\nSolving problem P12..." << endl;

//    int dim = 7+5; // x1,..,x7,l1,...,l5

//    std::vector<double> lb = {0,0,0,0,0,0,0,            -INF,-INF,-INF,-INF,-INF};
//    std::vector<double> ub = {1,1,1,100,200,100,200,    INF,INF,INF,INF,INF};

//    ConstraintCompositePtr cs(new ConstraintComposite(dim, lb, ub));

//    // x4*(9 - 6*x1 - 16*x2 - 15*x3) = l1
//    {
//        std::vector<double> thislb = {lb.at(0), lb.at(1), lb.at(2), lb.at(3)};
//        std::vector<double> thisub = {ub.at(0), ub.at(1), ub.at(2), ub.at(3)};

//        std::vector<int> deg = {1,1,1,1};

//        // Poly coeffs
//        DenseVector c(4);
//        c.setZero();
//        c(0) = 9;
//        c(1) = -6;
//        c(2) = -16;
//        c(3) = -15;

//        // Poly exponents
//        DenseMatrix E(4,4); E.setZero();
//        E(0,3) = 1;
//        E(1,0) = 1; E(1,3) = 1;
//        E(2,1) = 1; E(2,3) = 1;
//        E(3,2) = 1; E(3,3) = 1;

//        DenseMatrix coeffs = getBsplineCoefficients(c, E, thislb, thisub);

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {0,1,2,3,7};

//        cs->add(cbs,vars);
//    }

//    // x5*(15 - 6*x1 - 16*x2 - 15*x3) = l2
//    {
//        std::vector<double> thislb = {lb.at(0), lb.at(1), lb.at(2), lb.at(4)};
//        std::vector<double> thisub = {ub.at(0), ub.at(1), ub.at(2), ub.at(4)};

//        std::vector<int> deg = {1,1,1,1};

//        // Poly coeffs
//        DenseVector c(4); c.setZero();
//        c(0) = 15;
//        c(1) = -6;
//        c(2) = -16;
//        c(3) = -15;

//        // Poly exponents
//        DenseMatrix E(4,4); E.setZero();
//        E(0,3) = 1;
//        E(1,0) = 1; E(1,3) = 1;
//        E(2,1) = 1; E(2,3) = 1;
//        E(3,2) = 1; E(3,3) = 1;

//        DenseMatrix coeffs = getBsplineCoefficients(c, E, thislb, thisub);

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {0,1,2,4,8};

//        cs->add(cbs,vars);
//    }

//    // x3*(x4 + x5) = l3
//    {
//        std::vector<double> thislb = {lb.at(2), lb.at(3), lb.at(4)};
//        std::vector<double> thisub = {ub.at(2), ub.at(3), ub.at(4)};

//        std::vector<int> deg = {1,1,1};

//        // Poly coeffs
//        DenseVector c(2); c.setZero();
//        c(0) = 1;
//        c(1) = 1;

//        // Poly exponents
//        DenseMatrix E(2,3); E.setZero();
//        E(0,0) = 1; E(0,1) = 1;
//        E(1,0) = 1; E(1,2) = 1;

//        DenseMatrix coeffs = getBsplineCoefficients(c, E, thislb, thisub);

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {2,3,4,9};

//        cs->add(cbs,vars);
//    }

//    // x4*(3*x1 + x2 + x3 - 2.5) = l4
//    {
//        std::vector<double> thislb = {lb.at(0), lb.at(1), lb.at(2), lb.at(3)};
//        std::vector<double> thisub = {ub.at(0), ub.at(1), ub.at(2), ub.at(3)};

//        std::vector<int> deg = {1,1,1,1};

//        // Poly coeffs
//        DenseVector c(4); c.setZero();
//        c(0) = 3;
//        c(1) = 1;
//        c(2) = 1;
//        c(3) = -2.5;

//        // Poly exponents
//        DenseMatrix E(4,4); E.setZero();
//        E(0,0) = 1; E(0,3) = 1;
//        E(1,1) = 1; E(1,3) = 1;
//        E(2,2) = 1; E(2,3) = 1;
//        E(3,3) = 1;

//        DenseMatrix coeffs = getBsplineCoefficients(c, E, thislb, thisub);

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {0,1,2,3,10};

//        cs->add(cbs,vars);
//    }

//    // x5*(3*x1 + x2 + x3 - 1.5) = l4
//    {
//        std::vector<double> thislb = {lb.at(0), lb.at(1), lb.at(2), lb.at(4)};
//        std::vector<double> thisub = {ub.at(0), ub.at(1), ub.at(2), ub.at(4)};

//        std::vector<int> deg = {1,1,1,1};

//        // Poly coeffs
//        DenseVector c(4); c.setZero();
//        c(0) = 3;
//        c(1) = 1;
//        c(2) = 1;
//        c(3) = -1.5;

//        // Poly exponents
//        DenseMatrix E(4,4); E.setZero();
//        E(0,0) = 1; E(0,3) = 1;
//        E(1,1) = 1; E(1,3) = 1;
//        E(2,2) = 1; E(2,3) = 1;
//        E(3,3) = 1;

//        DenseMatrix coeffs = getBsplineCoefficients(c, E, thislb, thisub);

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {0,1,2,4,11};

//        cs->add(cbs,vars);
//    }

//    { // Linear inequality constraints
//        DenseMatrix A = DenseMatrix::Zero(5,dim);
//        A(0,9) = 1;
//        A(1,3) = 1; A(1,5) = 1;
//        A(2,4) = 1; A(2,6) = 1;
//        A(3,10) = 1; A(3,5) = -0.5;
//        A(4,11) = 1; A(4,6) = 0.5;

//        DenseVector b;
//        b.setZero(5);
//        b(0) = 50;
//        b(1) = 100;
//        b(2) = 200;

//        ConstraintPtr c(new ConstraintLinear(A,b,false));

//        //std::vector<int> vars = {0,1,2,3,4,5};
//        cs->add(c); // All variables
//    }

//    { // Linear equality constraints
//        DenseMatrix A = DenseMatrix::Zero(1,3);
//        A(0,0) = 1; A(0,1) = 1; A(0,2) = 1;

//        DenseVector b(1);
//        b(0) = 1;

//        ConstraintPtr c(new ConstraintLinear(A,b,true));

//        std::vector<int> vars = {0,1,2};
//        cs->add(c,vars);
//    }

//    DenseMatrix cobj(1,dim);
//    cobj.setZero();
//    cobj(0,5) = 1;
//    cobj(0,6) = -5;
//    cobj(0,7) = -1;
//    cobj(0,8) = -1;

//    ObjectivePtr obj(new ObjectiveLinear(cobj));

//    std::vector<int> vt(dim,CONTINUOUS);
//    std::vector<int> bv = {0,1,2,3,4};
//    std::vector<double> z0(dim,0);

//    BB::BranchAndBound bnb(obj, cs, z0, vt, bv);
//    Timer timer;
//    timer.start();
//    bnb.optimize();
//    timer.stop();
//    cout << "Time: " << timer.getMilliSeconds() << " (ms)" << endl;
//}

//void P12_3()
//{
//    // Problem 12 (Nataraj), Bilinear problem
//    // Six 2-D B-splines
//    cout << "\n\nSolving problem P12..." << endl;

//    int dim = 7+6; // x1,..,x7,l1,...,l6

//    std::vector<double> lb = {0,0,0,0,0,0,0,            0,0,0,0,0,0};
//    std::vector<double> ub = {1,1,1,1,1,1,1,    1,1,1,1,1,1};

//    ConstraintCompositePtr cs(new ConstraintComposite(dim, lb, ub));

//    // x_i * x_4 = li
//    for (int i = 0; i < 3; i++)
//    {
//        std::vector<double> thislb = {lb.at(i), lb.at(3)};
//        std::vector<double> thisub = {ub.at(i), ub.at(3)};

//        std::vector<int> deg = {1,1};

//        // Poly coeffs
//        DenseVector c(4);
//        c.setZero();
//        c(3) = 1;

//        DenseMatrix T = getTransformationMatrix(deg, thislb, thisub);
//        DenseMatrix coeffs = T*c;

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {i,3,i+7};

//        cs->add(cbs,vars);
//    }

//    // x_i * x_5 = li
//    for (int i = 0; i < 3; i++)
//    {
//        std::vector<double> thislb = {lb.at(i), lb.at(4)};
//        std::vector<double> thisub = {ub.at(i), ub.at(4)};

//        std::vector<int> deg = {1,1};

//        // Poly coeffs
//        DenseVector c(4);
//        c.setZero();
//        c(3) = 1;

//        DenseMatrix T = getTransformationMatrix(deg, thislb, thisub);
//        DenseMatrix coeffs = T*c;

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {i,4,i+7+3};

//        cs->add(cbs,vars);
//    }

//    { // Linear inequality constraints
//        DenseMatrix A = DenseMatrix::Zero(5,dim);
//        A(0,9) = 0.5; A(0,12) = 1;
//        A(1,3) = 1; A(1,5) = 1;
//        A(2,4) = 1; A(2,6) = 1;
//        A(3,7) = 3; A(3,8) = 1; A(3,9) = 1; A(3,3) = -2.5; A(3,5) = -0.5;
//        A(4,10) = 3; A(4,11) = 1; A(4,12) = 1; A(4,4) = -1.5; A(4,6) = 0.5;

//        DenseVector b;
//        b.setZero(5);
//        b(0) = 50/200.0;
//        b(1) = 1;
//        b(2) = 1;

//        ConstraintPtr c(new ConstraintLinear(A,b,false));

//        //std::vector<int> vars = {0,1,2,3,4,5};
//        cs->add(c); // All variables
//    }

//    { // Linear equality constraints
//        DenseMatrix A = DenseMatrix::Zero(1,3);
//        A(0,0) = 1; A(0,1) = 1; A(0,2) = 1;

//        DenseVector b(1);
//        b(0) = 1;

//        ConstraintPtr c(new ConstraintLinear(A,b,true));

//        std::vector<int> vars = {0,1,2};
//        cs->add(c,vars);
//    }

//    DenseMatrix cobj(1,dim);
//    cobj.setZero();
//    cobj(0,3) = -9;
//    cobj(0,7) = 6;
//    cobj(0,8) = 16;
//    cobj(0,9) = 15;
//    cobj(0,4) = -15*2;
//    cobj(0,10) = 6*2;
//    cobj(0,11) = 16*2;
//    cobj(0,12) = 15*2;
//    cobj(0,5) = 1;
//    cobj(0,6) = -5*2;

//    ObjectivePtr obj(new ObjectiveLinear(cobj));

//    std::vector<int> vt(dim,CONTINUOUS);
//    std::vector<int> bv = {0,1,2,3,4};
//    std::vector<double> z0(dim,0);

//    BB::BranchAndBound bnb(obj, cs, z0, vt, bv);
//    Timer timer;
//    timer.start();
//    bnb.optimize();
//    timer.stop();
//    cout << "Time: " << timer.getMilliSeconds() << " (ms)" << endl;
//}

//void P13()
//{
//    // Problem 13 (Nataraj)
//    cout << "\n\nSolving problem P13..." << endl;

//    int dim = 7+8; // x1,..,x7,l1,...,l5

//    std::vector<double> lb = {2.6,0.7,17,7.3,7.3,2.9,5,-INF,-INF,-INF,-INF,-INF,-INF,-INF,-INF};
//    std::vector<double> ub = {3.6,0.8,28,8.3,8.3,3.9,5.5,INF,INF,INF,INF,INF,INF,INF,INF};

//    // Testing manual bounds tightening
//    lb.at(0) = 3.5;
//    ub.at(1) = 0.72;
//    lb.at(4) = 7.4;

//    ConstraintCompositePtr cs(new ConstraintComposite(dim, lb, ub));

//    // obj = l1
//    {
//        std::vector<double> thislb = {lb.at(0), lb.at(1), lb.at(2), lb.at(3), lb.at(4), lb.at(5), lb.at(6)};
//        std::vector<double> thisub = {ub.at(0), ub.at(1), ub.at(2), ub.at(3), ub.at(4), ub.at(5), ub.at(6)};

//        std::vector<int> deg = {1,2,2,1,1,3,3};

//        double a1 = 0.7854;
//        double a2 = 3.3333;
//        double a3 = 14.9334;
//        double a4 = -43.0934;
//        double a5 = -1.508;
//        double a6 = 7.477;
//        double a7 = 0.7854;

//        // Poly coeffs
//        DenseVector c(9);
//        c.setZero();
//        c(0) = a1*a2;
//        c(1) = a1*a3;
//        c(2) = a1*a4;
//        c(3) = a5;
//        c(4) = a5;
//        c(5) = a6;
//        c(6) = a6;
//        c(7) = a7;
//        c(8) = a7;

//        // Poly exponents
//        DenseMatrix E(9,7);
//        E.setZero();
//        E(0,0) = 1; E(0,1) = 2; E(0,2) = 2;
//        E(1,0) = 1; E(1,1) = 2; E(1,2) = 1;
//        E(2,0) = 1; E(2,1) = 2;
//        E(3,0) = 1; E(3,5) = 2;
//        E(4,0) = 1; E(4,6) = 2;
//        E(5,5) = 3;
//        E(6,6) = 3;
//        E(7,3) = 1; E(7,5) = 2;
//        E(8,4) = 1; E(8,6) = 2;

//        DenseMatrix coeffs = getBsplineCoefficients(c, E, thislb, thisub);

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {0,1,2,3,4,5,6,7};

//        cs->add(cbs,vars);
//    }

//    // x1*x2*x2*x3 = l2
//    {
//        std::vector<double> thislb = {lb.at(0), lb.at(1), lb.at(2)};
//        std::vector<double> thisub = {ub.at(0), ub.at(1), ub.at(2)};

//        std::vector<int> deg = {1,2,1};

//        // Poly coeffs
//        DenseVector c(1);
//        c.setZero();
//        c(0) = 1;

//        // Poly exponents
//        DenseMatrix E(1,3);
//        E.setZero();
//        E(0,0) = 1; E(0,1) = 2; E(0,2) = 1;

//        DenseMatrix coeffs = getBsplineCoefficients(c, E, thislb, thisub);

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {0,1,2,8};

//        cs->add(cbs,vars);
//    }

//    // x1*x2*x2*x3*x3 = l3
//    {
//        std::vector<double> thislb = {lb.at(0), lb.at(1), lb.at(2)};
//        std::vector<double> thisub = {ub.at(0), ub.at(1), ub.at(2)};

//        std::vector<int> deg = {1,2,2};

//        // Poly coeffs
//        DenseVector c(1);
//        c.setZero();
//        c(0) = 1;

//        // Poly exponents
//        DenseMatrix E(1,3);
//        E.setZero();
//        E(0,0) = 1; E(0,1) = 2; E(0,2) = 2;

//        DenseMatrix coeffs = getBsplineCoefficients(c, E, thislb, thisub);

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {0,1,2,9};

//        cs->add(cbs,vars);
//    }

//    // x2*x6^4*x3 - 1.93*x4^3 = l4
//    {
//        std::vector<double> thislb = {lb.at(1), lb.at(2), lb.at(3), lb.at(5)};
//        std::vector<double> thisub = {ub.at(1), ub.at(2), ub.at(3), ub.at(5)};

//        std::vector<int> deg = {1,1,3,4};

//        // Poly coeffs
//        DenseVector c(2);
//        c.setZero();
//        c(0) = 1;
//        c(1) = -1.93;

//        // Poly exponents
//        DenseMatrix E(2,4);
//        E.setZero();
//        E(0,0) = 1; E(0,1) = 1; E(0,3) = 4;
//        E(1,2) = 3;

//        DenseMatrix coeffs = getBsplineCoefficients(c, E, thislb, thisub);

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {1,2,3,5,10};

//        cs->add(cbs,vars);
//    }

//    // x2*x7^4*x3 - 1.93*x5^3 = l5
//    {
//        std::vector<double> thislb = {lb.at(1), lb.at(2), lb.at(4), lb.at(6)};
//        std::vector<double> thisub = {ub.at(1), ub.at(2), ub.at(4), ub.at(6)};

//        std::vector<int> deg = {1,1,3,4};

//        // Poly coeffs
//        DenseVector c(2);
//        c.setZero();
//        c(0) = 1;
//        c(1) = -1.93;

//        // Poly exponents
//        DenseMatrix E(2,4);
//        E.setZero();
//        E(0,0) = 1; E(0,1) = 1; E(0,3) = 4;
//        E(1,2) = 3;

//        DenseMatrix coeffs = getBsplineCoefficients(c, E, thislb, thisub);

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {1,2,4,6,11};

//        cs->add(cbs,vars);
//    }

//    // g(x2,x3,x4,x6) = l6
//    {
//        std::vector<double> thislb = {lb.at(1), lb.at(2), lb.at(3), lb.at(5)};
//        std::vector<double> thisub = {ub.at(1), ub.at(2), ub.at(3), ub.at(5)};

//        std::vector<int> deg = {2,2,2,6};

//        // Poly coeffs
//        DenseVector c(3);
//        c.setZero();
//        c(0) = 745*745/16.911e6;
//        c(1) = 1;
//        c(2) = -110*110/16.911e6;

//        // Poly exponents
//        DenseMatrix E(3,4);
//        E.setZero();
//        E(0,2) = 2;
//        E(1,0) = 2; E(1,1) = 2;
//        E(2,0) = 2; E(2,1) = 2; E(2,3) = 6;

//        DenseMatrix coeffs = getBsplineCoefficients(c, E, thislb, thisub);

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {1,2,3,5,12};

//        cs->add(cbs,vars);
//    }

//    // g(x2,x3,x5,x7) = l7
//    {
//        std::vector<double> thislb = {lb.at(1), lb.at(2), lb.at(4), lb.at(6)};
//        std::vector<double> thisub = {ub.at(1), ub.at(2), ub.at(4), ub.at(6)};

//        std::vector<int> deg = {2,2,2,6};

//        // Poly coeffs
//        DenseVector c(3);
//        c.setZero();
//        c(0) = 745*745/157.50e6;
//        c(1) = 1;
//        c(2) = -85*85/157.50e6;

//        // Poly exponents
//        DenseMatrix E(3,4);
//        E.setZero();
//        E(0,2) = 2;
//        E(1,0) = 2; E(1,1) = 2;
//        E(2,0) = 2; E(2,1) = 2; E(2,3) = 6;

//        DenseMatrix coeffs = getBsplineCoefficients(c, E, thislb, thisub);

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {1,2,4,6,13};

//        cs->add(cbs,vars);
//    }

//    // x2*x3 = l8
//    {
//        std::vector<double> thislb = {lb.at(1), lb.at(2)};
//        std::vector<double> thisub = {ub.at(1), ub.at(2)};

//        std::vector<int> deg = {1,1};

//        // Poly coeffs
//        DenseVector c(4);
//        c.setZero();
//        c(3) = 1;

//        DenseMatrix T = getTransformationMatrix(deg, thislb, thisub);
//        DenseMatrix coeffs = T*c;

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {1,2,14};

//        cs->add(cbs,vars);
//    }

//    { // Linear inequality constraints
//        DenseMatrix A = DenseMatrix::Zero(11,dim);
//        A(0,8) = -1;
//        A(1,9) = -1;
//        A(2,10) = -1;
//        A(3,11) = -1;
//        A(4,12) = 1;
//        A(5,13) = 1;
//        A(6,14) = 1;
//        A(7,0) = -1; A(7,1) = 5;
//        A(8,0) = 1; A(8,1) = -12;
//        A(9,3) = -1; A(9,5) = 1.5;
//        A(10,4) = -1; A(10,6) = 1.1;

//        DenseVector b;
//        b.setZero(11);
//        b(0) = -27;
//        b(1) = -397.5;
//        b(6) = 40;
//        b(9) = -1.9;
//        b(10) = -1.9;

//        ConstraintPtr c(new ConstraintLinear(A,b,false));

//        //std::vector<int> vars = {0,1,2,3,4,5};
//        cs->add(c); // All variables
//    }

//    DenseMatrix cobj(1,dim);
//    cobj.setZero();
//    cobj(0,7) = 1;

//    ObjectivePtr obj(new ObjectiveLinear(cobj));

//    std::vector<int> vt(dim,CONTINUOUS);
//    std::vector<int> bv = {0,1,2,3,4,5,6};
//    //std::vector<int> bv = {0,2,3,4}; // With BT
//    std::vector<double> z0(dim,0);

//    BB::BranchAndBound bnb(obj, cs, z0, vt, bv);
//    Timer timer;
//    timer.start();
//    bnb.optimize();
//    timer.stop();
//    cout << "Time: " << timer.getMilliSeconds() << " (ms)" << endl;
//}

//void P13_2()
//{
//    // Problem 13 (Nataraj)
//    cout << "\n\nSolving problem P13..." << endl;

//    int dim = 7+15; // x1,..,x7,l1,...,l15

//    std::vector<double> lb = {2.6,0.7,17,7.3,7.3,2.9,5,-INF,-INF,-INF,-INF,-INF,-INF,-INF,-INF,-INF,-INF,-INF,-INF,-INF,-INF,-INF};
//    std::vector<double> ub = {3.6,0.8,28,8.3,8.3,3.9,5.5,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF};

//    // Testing manual bounds tightening
//    lb.at(0) = 3.5;
//    ub.at(1) = 0.72;
//    lb.at(4) = 7.4;

//    ConstraintCompositePtr cs(new ConstraintComposite(dim, lb, ub));

//    // obj = l1
//    {
//        std::vector<double> thislb = {lb.at(0), lb.at(1), lb.at(2)};
//        std::vector<double> thisub = {ub.at(0), ub.at(1), ub.at(2)};

//        std::vector<int> deg = {1,2,2};

//        double a1 = 0.7854;
//        double a2 = 3.3333;
//        double a3 = 14.9334;
//        double a4 = -43.0934;

//        // Poly coeffs
//        DenseVector c(3);
//        c.setZero();
//        c(0) = a1*a2;
//        c(1) = a1*a3;
//        c(2) = a1*a4;

//        // Poly exponents
//        DenseMatrix E(3,3);
//        E.setZero();
//        E(0,0) = 1; E(0,1) = 2; E(0,2) = 2;
//        E(1,0) = 1; E(1,1) = 2; E(1,2) = 1;
//        E(2,0) = 1; E(2,1) = 2;

//        DenseMatrix coeffs = getBsplineCoefficients(c, E, thislb, thisub);

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {0,1,2,7};

//        cs->add(cbs,vars);
//    }

//    // obj = l2
//    {
//        std::vector<double> thislb = {lb.at(0), lb.at(5), lb.at(6)};
//        std::vector<double> thisub = {ub.at(0), ub.at(5), ub.at(6)};

//        std::vector<int> deg = {1,3,3};

//        double a1 = -1.508;
//        double a2 = 7.477;

//        // Poly coeffs
//        DenseVector c(4);
//        c.setZero();
//        c(0) = a1;
//        c(1) = a1;
//        c(2) = a2;
//        c(3) = a2;

//        // Poly exponents
//        DenseMatrix E(4,3);
//        E.setZero();
//        E(0,0) = 1; E(0,1) = 2;
//        E(1,0) = 1; E(1,2) = 2;
//        E(2,1) = 3;
//        E(3,2) = 3;

//        DenseMatrix coeffs = getBsplineCoefficients(c, E, thislb, thisub);

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {0,5,6,8};

//        cs->add(cbs,vars);
//    }

//    // obj = l3
//    {
//        std::vector<double> thislb = {lb.at(3), lb.at(5)};
//        std::vector<double> thisub = {ub.at(3), ub.at(5)};

//        std::vector<int> deg = {1,2};

//        double a1 = 0.7854;

//        // Poly coeffs
//        DenseVector c(1);
//        c.setZero();
//        c(0) = a1;

//        // Poly exponents
//        DenseMatrix E(1,2);
//        E.setZero();
//        E(0,0) = 1; E(0,1) = 2;

//        DenseMatrix coeffs = getBsplineCoefficients(c, E, thislb, thisub);

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {3,5,9};

//        cs->add(cbs,vars);
//    }

//    // obj = l4
//    {
//        std::vector<double> thislb = {lb.at(4), lb.at(6)};
//        std::vector<double> thisub = {ub.at(4), ub.at(6)};

//        std::vector<int> deg = {1,2};

//        double a1 = 0.7854;

//        // Poly coeffs
//        DenseVector c(1);
//        c.setZero();
//        c(0) = a1;

//        // Poly exponents
//        DenseMatrix E(1,2);
//        E.setZero();
//        E(0,0) = 1; E(0,1) = 2;

//        DenseMatrix coeffs = getBsplineCoefficients(c, E, thislb, thisub);

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {4,6,10};

//        cs->add(cbs,vars);
//    }

//    // x1*x2*x2*x3 = l5
//    {
//        std::vector<double> thislb = {lb.at(0), lb.at(1), lb.at(2)};
//        std::vector<double> thisub = {ub.at(0), ub.at(1), ub.at(2)};

//        std::vector<int> deg = {1,2,1};

//        // Poly coeffs
//        DenseVector c(1);
//        c.setZero();
//        c(0) = 1;

//        // Poly exponents
//        DenseMatrix E(1,3);
//        E.setZero();
//        E(0,0) = 1; E(0,1) = 2; E(0,2) = 1;

//        DenseMatrix coeffs = getBsplineCoefficients(c, E, thislb, thisub);

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {0,1,2,11};

//        cs->add(cbs,vars);
//    }

//    // x1*x2*x2*x3*x3 = l6
//    {
//        std::vector<double> thislb = {lb.at(0), lb.at(1), lb.at(2)};
//        std::vector<double> thisub = {ub.at(0), ub.at(1), ub.at(2)};

//        std::vector<int> deg = {1,2,2};

//        // Poly coeffs
//        DenseVector c(1);
//        c.setZero();
//        c(0) = 1;

//        // Poly exponents
//        DenseMatrix E(1,3);
//        E.setZero();
//        E(0,0) = 1; E(0,1) = 2; E(0,2) = 2;

//        DenseMatrix coeffs = getBsplineCoefficients(c, E, thislb, thisub);

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {0,1,2,12};

//        cs->add(cbs,vars);
//    }

//    // x2*x3*x6^4 = l7
//    {
//        std::vector<double> thislb = {lb.at(1), lb.at(2), lb.at(5)};
//        std::vector<double> thisub = {ub.at(1), ub.at(2), ub.at(5)};

//        std::vector<int> deg = {1,1,4};

//        // Poly coeffs
//        DenseVector c(1);
//        c.setZero();
//        c(0) = 1;

//        // Poly exponents
//        DenseMatrix E(1,3);
//        E.setZero();
//        E(0,0) = 1; E(0,1) = 1; E(0,2) = 4;

//        DenseMatrix coeffs = getBsplineCoefficients(c, E, thislb, thisub);

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {1,2,5,13};

//        cs->add(cbs,vars);
//    }

//    // x4^3 = l8
//    {
//        std::vector<double> thislb = {lb.at(3)};
//        std::vector<double> thisub = {ub.at(3)};

//        std::vector<int> deg = {3};

//        // Poly coeffs
//        DenseVector c(1);
//        c.setZero();
//        c(0) = 1;

//        // Poly exponents
//        DenseMatrix E(1,1);
//        E.setZero();
//        E(0,0) = 3;

//        DenseMatrix coeffs = getBsplineCoefficients(c, E, thislb, thisub);

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {3,14};

//        cs->add(cbs,vars);
//    }

//    // x2*x3*x7^4 = l9
//    {
//        std::vector<double> thislb = {lb.at(1), lb.at(2), lb.at(6)};
//        std::vector<double> thisub = {ub.at(1), ub.at(2), ub.at(6)};

//        std::vector<int> deg = {1,1,4};

//        // Poly coeffs
//        DenseVector c(1);
//        c.setZero();
//        c(0) = 1;

//        // Poly exponents
//        DenseMatrix E(1,3);
//        E.setZero();
//        E(0,0) = 1; E(0,1) = 1; E(0,2) = 4;

//        DenseMatrix coeffs = getBsplineCoefficients(c, E, thislb, thisub);

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {1,2,6,15};

//        cs->add(cbs,vars);
//    }

//    // x5^3 = l10
//    {
//        std::vector<double> thislb = {lb.at(4)};
//        std::vector<double> thisub = {ub.at(4)};

//        std::vector<int> deg = {3};

//        // Poly coeffs
//        DenseVector c(1);
//        c.setZero();
//        c(0) = 1;

//        // Poly exponents
//        DenseMatrix E(1,1);
//        E.setZero();
//        E(0,0) = 3;

//        DenseMatrix coeffs = getBsplineCoefficients(c, E, thislb, thisub);

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {4,16};

//        cs->add(cbs,vars);
//    }

//    // x4^2 = l11
//    {
//        std::vector<double> thislb = {lb.at(3)};
//        std::vector<double> thisub = {ub.at(3)};

//        std::vector<int> deg = {2};

//        // Poly coeffs
//        DenseVector c(1);
//        c.setZero();
//        c(0) = 1;

//        // Poly exponents
//        DenseMatrix E(1,1);
//        E.setZero();
//        E(0,0) = 2;

//        DenseMatrix coeffs = getBsplineCoefficients(c, E, thislb, thisub);

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {3,17};

//        cs->add(cbs,vars);
//    }

//    // l12
//    {
//        std::vector<double> thislb = {lb.at(1), lb.at(2), lb.at(5)};
//        std::vector<double> thisub = {ub.at(1), ub.at(2), ub.at(5)};

//        std::vector<int> deg = {2,2,6};

//        double a1 = 12100/(745.0*745.0);
//        double a2 = -16.91e6/(745.0*745.0);

//        // Poly coeffs
//        DenseVector c(2);
//        c.setZero();
//        c(0) = a1;
//        c(1) = a2;

//        // Poly exponents
//        DenseMatrix E(2,3);
//        E.setZero();
//        E(0,0) = 2; E(0,1) = 2; E(0,2) = 6;
//        E(1,0) = 2; E(1,1) = 2;

//        DenseMatrix coeffs = getBsplineCoefficients(c, E, thislb, thisub);

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {1,2,5,18};

//        cs->add(cbs,vars);
//    }

//    // x5^2 = l13
//    {
//        std::vector<double> thislb = {lb.at(4)};
//        std::vector<double> thisub = {ub.at(4)};

//        std::vector<int> deg = {2};

//        // Poly coeffs
//        DenseVector c(1);
//        c.setZero();
//        c(0) = 1;

//        // Poly exponents
//        DenseMatrix E(1,1);
//        E.setZero();
//        E(0,0) = 2;

//        DenseMatrix coeffs = getBsplineCoefficients(c, E, thislb, thisub);

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {4,19};

//        cs->add(cbs,vars);
//    }

//    // l14
//    {
//        std::vector<double> thislb = {lb.at(1), lb.at(2), lb.at(6)};
//        std::vector<double> thisub = {ub.at(1), ub.at(2), ub.at(6)};

//        std::vector<int> deg = {2,2,6};

//        double a1 = 7225/(745.0*745.0);
//        double a2 = -157.5e6/(745.0*745.0);

//        // Poly coeffs
//        DenseVector c(2);
//        c.setZero();
//        c(0) = a1;
//        c(1) = a2;

//        // Poly exponents
//        DenseMatrix E(2,3);
//        E.setZero();
//        E(0,0) = 2; E(0,1) = 2; E(0,2) = 6;
//        E(1,0) = 2; E(1,1) = 2;

//        DenseMatrix coeffs = getBsplineCoefficients(c, E, thislb, thisub);

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {1,2,6,20};

//        cs->add(cbs,vars);
//    }

//    // x2*x3 = l15
//    {
//        std::vector<double> thislb = {lb.at(1), lb.at(2)};
//        std::vector<double> thisub = {ub.at(1), ub.at(2)};

//        std::vector<int> deg = {1,1};

//        // Poly coeffs
//        DenseVector c(4);
//        c.setZero();
//        c(3) = 1;

//        DenseMatrix T = getTransformationMatrix(deg, thislb, thisub);
//        DenseMatrix coeffs = T*c;

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {1,2,21};

//        cs->add(cbs,vars);
//    }

//    { // Linear inequality constraints
//        DenseMatrix A = DenseMatrix::Zero(11,dim);
//        A(0,11) = -1;
//        A(1,12) = -1;
//        A(2,13) = -1; A(2,14) = 1.93;
//        A(3,15) = -1; A(3,16) = 1.93;
//        A(4,17) = 1; A(4,18) = -1;
//        A(5,19) = 1; A(5,20) = -1;
//        A(6,21) = 1;
//        A(7,0) = -1; A(7,1) = 5;
//        A(8,0) = 1; A(8,1) = -12;
//        A(9,3) = -1; A(9,5) = 1.5;
//        A(10,4) = -1; A(10,6) = 1.1;

//        DenseVector b;
//        b.setZero(11);
//        b(0) = -27;
//        b(1) = -397.5;
//        b(6) = 40;
//        b(9) = -1.9;
//        b(10) = -1.9;

//        ConstraintPtr c(new ConstraintLinear(A,b,false));

//        //std::vector<int> vars = {0,1,2,3,4,5};
//        cs->add(c); // All variables
//    }

//    DenseMatrix cobj(1,dim);
//    cobj.setZero();
//    cobj(0,7) = 1;
//    cobj(0,8) = 1;
//    cobj(0,9) = 1;
//    cobj(0,10) = 1;

//    ObjectivePtr obj(new ObjectiveLinear(cobj));

//    std::vector<int> vt(dim,CONTINUOUS);
//    std::vector<int> bv = {0,1,2,3,4,5,6};
//    //std::vector<int> bv = {0,2,3,4}; // With BT
//    std::vector<double> z0(dim,0);

//    BB::BranchAndBound bnb(obj, cs, z0, vt, bv);
//    Timer timer;
//    timer.start();
//    bnb.optimize();
//    timer.stop();
//    cout << "Time: " << timer.getMilliSeconds() << " (ms)" << endl;
//}

//void P13_3()
//{
//    // Problem 13 (Nataraj)
//    cout << "\n\nSolving problem P13..." << endl;

//    int dim = 7+15; // x1,..,x7,l1,...,l15

//    double x3s = 28.0;
//    double x4s = 8.3;
//    double x5s = 8.3;
//    std::vector<double> lb = {2.6,0.7,17/x3s,7.3/x4s,7.3/x5s,2.9,5,-INF,-INF,-INF,-INF,-INF,-INF,-INF,-INF,-INF,-INF,-INF,-INF,-INF,-INF,-INF};
//    std::vector<double> ub = {3.6,0.8,28/x3s,8.3/x4s,8.3/x5s,3.9,5.5,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF};

//    // Testing manual bounds tightening
//    lb.at(0) = 3.5;
//    ub.at(1) = 0.72;
//    lb.at(4) = 7.4/x5s;

//    ConstraintCompositePtr cs(new ConstraintComposite(dim, lb, ub));

//    // obj = l1
//    {
//        std::vector<double> thislb = {lb.at(0), lb.at(1), lb.at(2)};
//        std::vector<double> thisub = {ub.at(0), ub.at(1), ub.at(2)};

//        std::vector<int> deg = {1,2,2};

//        double a1 = 0.7854;
//        double a2 = 3.3333;
//        double a3 = 14.9334;
//        double a4 = -43.0934;

//        // Poly coeffs
//        DenseVector c(3);
//        c.setZero();
//        c(0) = a1*a2*x3s*x3s;
//        c(1) = a1*a3*x3s;
//        c(2) = a1*a4;

//        // Poly exponents
//        DenseMatrix E(3,3);
//        E.setZero();
//        E(0,0) = 1; E(0,1) = 2; E(0,2) = 2;
//        E(1,0) = 1; E(1,1) = 2; E(1,2) = 1;
//        E(2,0) = 1; E(2,1) = 2;

//        DenseMatrix coeffs = getBsplineCoefficients(c, E, thislb, thisub);

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {0,1,2,7};

//        cs->add(cbs,vars);
//    }

//    // obj = l2
//    {
//        std::vector<double> thislb = {lb.at(0), lb.at(5), lb.at(6)};
//        std::vector<double> thisub = {ub.at(0), ub.at(5), ub.at(6)};

//        std::vector<int> deg = {1,3,3};

//        double a1 = -1.508;
//        double a2 = 7.477;

//        // Poly coeffs
//        DenseVector c(4);
//        c.setZero();
//        c(0) = a1;
//        c(1) = a1;
//        c(2) = a2;
//        c(3) = a2;

//        // Poly exponents
//        DenseMatrix E(4,3);
//        E.setZero();
//        E(0,0) = 1; E(0,1) = 2;
//        E(1,0) = 1; E(1,2) = 2;
//        E(2,1) = 3;
//        E(3,2) = 3;

//        DenseMatrix coeffs = getBsplineCoefficients(c, E, thislb, thisub);

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {0,5,6,8};

//        cs->add(cbs,vars);
//    }

//    // obj = l3
//    {
//        std::vector<double> thislb = {lb.at(3), lb.at(5)};
//        std::vector<double> thisub = {ub.at(3), ub.at(5)};

//        std::vector<int> deg = {1,2};

//        double a1 = 0.7854;

//        // Poly coeffs
//        DenseVector c(1);
//        c.setZero();
//        c(0) = a1*x4s;

//        // Poly exponents
//        DenseMatrix E(1,2);
//        E.setZero();
//        E(0,0) = 1; E(0,1) = 2;

//        DenseMatrix coeffs = getBsplineCoefficients(c, E, thislb, thisub);

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {3,5,9};

//        cs->add(cbs,vars);
//    }

//    // obj = l4
//    {
//        std::vector<double> thislb = {lb.at(4), lb.at(6)};
//        std::vector<double> thisub = {ub.at(4), ub.at(6)};

//        std::vector<int> deg = {1,2};

//        double a1 = 0.7854*x5s;

//        // Poly coeffs
//        DenseVector c(1);
//        c.setZero();
//        c(0) = a1;

//        // Poly exponents
//        DenseMatrix E(1,2);
//        E.setZero();
//        E(0,0) = 1; E(0,1) = 2;

//        DenseMatrix coeffs = getBsplineCoefficients(c, E, thislb, thisub);

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {4,6,10};

//        cs->add(cbs,vars);
//    }

//    // x1*x2*x2*x3 = l5
//    {
//        std::vector<double> thislb = {lb.at(0), lb.at(1), lb.at(2)};
//        std::vector<double> thisub = {ub.at(0), ub.at(1), ub.at(2)};

//        std::vector<int> deg = {1,2,1};

//        // Poly coeffs
//        DenseVector c(1);
//        c.setZero();
//        c(0) = 1;

//        // Poly exponents
//        DenseMatrix E(1,3);
//        E.setZero();
//        E(0,0) = 1; E(0,1) = 2; E(0,2) = 1;

//        DenseMatrix coeffs = getBsplineCoefficients(c, E, thislb, thisub);

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {0,1,2,11};

//        cs->add(cbs,vars);
//    }

//    // x1*x2*x2*x3*x3 = l6
//    {
//        std::vector<double> thislb = {lb.at(0), lb.at(1), lb.at(2)};
//        std::vector<double> thisub = {ub.at(0), ub.at(1), ub.at(2)};

//        std::vector<int> deg = {1,2,2};

//        // Poly coeffs
//        DenseVector c(1);
//        c.setZero();
//        c(0) = 1;

//        // Poly exponents
//        DenseMatrix E(1,3);
//        E.setZero();
//        E(0,0) = 1; E(0,1) = 2; E(0,2) = 2;

//        DenseMatrix coeffs = getBsplineCoefficients(c, E, thislb, thisub);

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {0,1,2,12};

//        cs->add(cbs,vars);
//    }

//    // x2*x3*x6^4 = l7
//    {
//        std::vector<double> thislb = {lb.at(1), lb.at(2), lb.at(5)};
//        std::vector<double> thisub = {ub.at(1), ub.at(2), ub.at(5)};

//        std::vector<int> deg = {1,1,4};

//        // Poly coeffs
//        DenseVector c(1);
//        c.setZero();
//        c(0) = 1;

//        // Poly exponents
//        DenseMatrix E(1,3);
//        E.setZero();
//        E(0,0) = 1; E(0,1) = 1; E(0,2) = 4;

//        DenseMatrix coeffs = getBsplineCoefficients(c, E, thislb, thisub);

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {1,2,5,13};

//        cs->add(cbs,vars);
//    }

//    // x4^3 = l8
//    {
//        std::vector<double> thislb = {lb.at(3)};
//        std::vector<double> thisub = {ub.at(3)};

//        std::vector<int> deg = {3};

//        // Poly coeffs
//        DenseVector c(1);
//        c.setZero();
//        c(0) = 1*x4s*x4s*x4s;

//        // Poly exponents
//        DenseMatrix E(1,1);
//        E.setZero();
//        E(0,0) = 3;

//        DenseMatrix coeffs = getBsplineCoefficients(c, E, thislb, thisub);

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {3,14};

//        cs->add(cbs,vars);
//    }

//    // x2*x3*x7^4 = l9
//    {
//        std::vector<double> thislb = {lb.at(1), lb.at(2), lb.at(6)};
//        std::vector<double> thisub = {ub.at(1), ub.at(2), ub.at(6)};

//        std::vector<int> deg = {1,1,4};

//        // Poly coeffs
//        DenseVector c(1);
//        c.setZero();
//        c(0) = 1;

//        // Poly exponents
//        DenseMatrix E(1,3);
//        E.setZero();
//        E(0,0) = 1; E(0,1) = 1; E(0,2) = 4;

//        DenseMatrix coeffs = getBsplineCoefficients(c, E, thislb, thisub);

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {1,2,6,15};

//        cs->add(cbs,vars);
//    }

//    // x5^3 = l10
//    {
//        std::vector<double> thislb = {lb.at(4)};
//        std::vector<double> thisub = {ub.at(4)};

//        std::vector<int> deg = {3};

//        // Poly coeffs
//        DenseVector c(1);
//        c.setZero();
//        c(0) = 1*x5s*x5s;

//        // Poly exponents
//        DenseMatrix E(1,1);
//        E.setZero();
//        E(0,0) = 3;

//        DenseMatrix coeffs = getBsplineCoefficients(c, E, thislb, thisub);

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {4,16};

//        cs->add(cbs,vars);
//    }

//    // x4^2 = l11
//    {
//        std::vector<double> thislb = {lb.at(3)};
//        std::vector<double> thisub = {ub.at(3)};

//        std::vector<int> deg = {2};

//        // Poly coeffs
//        DenseVector c(1);
//        c.setZero();
//        c(0) = 1*x4s*x4s;

//        // Poly exponents
//        DenseMatrix E(1,1);
//        E.setZero();
//        E(0,0) = 2;

//        DenseMatrix coeffs = getBsplineCoefficients(c, E, thislb, thisub);

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {3,17};

//        cs->add(cbs,vars);
//    }

//    // l12
//    {
//        std::vector<double> thislb = {lb.at(1), lb.at(2), lb.at(5)};
//        std::vector<double> thisub = {ub.at(1), ub.at(2), ub.at(5)};

//        std::vector<int> deg = {2,2,6};

//        double a1 = 12100/(745.0*745.0)*x3s*x3s;
//        double a2 = -16.91e6/(745.0*745.0)*x3s*x3s;

//        // Poly coeffs
//        DenseVector c(2);
//        c.setZero();
//        c(0) = a1;
//        c(1) = a2;

//        // Poly exponents
//        DenseMatrix E(2,3);
//        E.setZero();
//        E(0,0) = 2; E(0,1) = 2; E(0,2) = 6;
//        E(1,0) = 2; E(1,1) = 2;

//        DenseMatrix coeffs = getBsplineCoefficients(c, E, thislb, thisub);

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {1,2,5,18};

//        cs->add(cbs,vars);
//    }

//    // x5^2 = l13
//    {
//        std::vector<double> thislb = {lb.at(4)};
//        std::vector<double> thisub = {ub.at(4)};

//        std::vector<int> deg = {2};

//        // Poly coeffs
//        DenseVector c(1);
//        c.setZero();
//        c(0) = 1*x5s*x5s;

//        // Poly exponents
//        DenseMatrix E(1,1);
//        E.setZero();
//        E(0,0) = 2;

//        DenseMatrix coeffs = getBsplineCoefficients(c, E, thislb, thisub);

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {4,19};

//        cs->add(cbs,vars);
//    }

//    // l14
//    {
//        std::vector<double> thislb = {lb.at(1), lb.at(2), lb.at(6)};
//        std::vector<double> thisub = {ub.at(1), ub.at(2), ub.at(6)};

//        std::vector<int> deg = {2,2,6};

//        double a1 = 7225/(745.0*745.0)*x3s*x3s;
//        double a2 = -157.5e6/(745.0*745.0)*x3s*x3s;

//        // Poly coeffs
//        DenseVector c(2);
//        c.setZero();
//        c(0) = a1;
//        c(1) = a2;

//        // Poly exponents
//        DenseMatrix E(2,3);
//        E.setZero();
//        E(0,0) = 2; E(0,1) = 2; E(0,2) = 6;
//        E(1,0) = 2; E(1,1) = 2;

//        DenseMatrix coeffs = getBsplineCoefficients(c, E, thislb, thisub);

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {1,2,6,20};

//        cs->add(cbs,vars);
//    }

//    // x2*x3 = l15
//    {
//        std::vector<double> thislb = {lb.at(1), lb.at(2)};
//        std::vector<double> thisub = {ub.at(1), ub.at(2)};

//        std::vector<int> deg = {1,1};

//        // Poly coeffs
//        DenseVector c(4);
//        c.setZero();
//        c(3) = 1;

//        DenseMatrix T = getTransformationMatrix(deg, thislb, thisub);
//        DenseMatrix coeffs = T*c;

//        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        Bspline bs(coeffs.transpose(), knots, deg);

//        ConstraintPtr cbs(new ConstraintBspline(bs,true));

//        std::vector<int> vars = {1,2,21};

//        cs->add(cbs,vars);
//    }

//    { // Linear inequality constraints
//        DenseMatrix A = DenseMatrix::Zero(11,dim);
//        A(0,11) = -1;
//        A(1,12) = -1;
//        A(2,13) = -1*x3s; A(2,14) = 1.93;
//        A(3,15) = -1*x3s; A(3,16) = 1.93;
//        A(4,17) = 1; A(4,18) = -1;
//        A(5,19) = 1; A(5,20) = -1;
//        A(6,21) = 1;
//        A(7,0) = -1; A(7,1) = 5;
//        A(8,0) = 1; A(8,1) = -12;
//        A(9,3) = -1*x4s; A(9,5) = 1.5;
//        A(10,4) = -1*x5s; A(10,6) = 1.1;

//        DenseVector b;
//        b.setZero(11);
//        b(0) = -27/x3s;
//        b(1) = -397.5/(x3s*x3s);
//        b(6) = 40/x3s;
//        b(9) = -1.9;
//        b(10) = -1.9;

//        ConstraintPtr c(new ConstraintLinear(A,b,false));

//        //std::vector<int> vars = {0,1,2,3,4,5};
//        cs->add(c); // All variables
//    }

//    DenseMatrix cobj(1,dim);
//    cobj.setZero();
//    cobj(0,7) = 1;
//    cobj(0,8) = 1;
//    cobj(0,9) = 1;
//    cobj(0,10) = 1;

//    ObjectivePtr obj(new ObjectiveLinear(cobj));

//    std::vector<int> vt(dim,CONTINUOUS);
//    std::vector<int> bv = {0,1,2,3,4,5,6};
//    //std::vector<int> bv = {0,2,3,4}; // With BT
//    std::vector<double> z0(dim,0);

//    BB::BranchAndBound bnb(obj, cs, z0, vt, bv);
//    Timer timer;
//    timer.start();
//    bnb.optimize();
//    timer.stop();
//    cout << "Time: " << timer.getMilliSeconds() << " (ms)" << endl;
//}


/*** START code from Anders S. ***/
//template<class T>
//void recSampling(std::vector<double> lb, std::vector<double> ub, std::vector<int> index, std::vector<int> numSamples, DataTable& table, T f, std::vector<double> x={}, int iter = 0)
//{
//    if( iter < index.size())
//    {
//        int i = index.at(iter);
//        for(const auto xji : linspace(lb.at(i), ub.at(i), numSamples.at(iter)))
//        {
//            auto xi = x;
//            xi.push_back(xji);
//            recSampling(lb, ub, index, numSamples, table, f, xi, iter+1);
//        }
//    }
//    else
//    {
//        table.addSample(x, f(x));
//    }
//}

//void P_14_1_4()
//{
//    double e = 2.718281828459045;
//    double pi = 3.141592653589793;

//    std::vector<std::vector<double>> xstar =
//    {{0.29945, 2.83693},
//    { 0.5,     3.14159}};

//    cout << "Test problem P_14_1_4" << endl;

//    int numVars = 2;
//    int numSlacks =  1;

//    int dim = numVars + numSlacks;

//    // x1, x2, s
//    std::vector<double> lb = {0.25, 1.5, 0};
//    std::vector<double> ub = {1, 2*pi, INF};

//    ConstraintCompositePtr cc(new ConstraintComposite(dim, lb, ub));
//    int ns = 10;
//    {
//        std::vector<int> indexTable = {0, 1};
//        std::vector<int> numSamples = {ns, ns};
//        auto fun = [=](std::vector<double> x)
//        {
//            auto x1 = x.at(0);
//            auto x2 = x.at(1);
//            return 0.5*std::sin(x1*x2) - 0.25*x2/pi - 0.5*x1;
//        };

//        DataTable table;
//        recSampling(lb, ub, indexTable, numSamples, table, fun);

//        std::vector<int> indexConstraint = indexTable;
//        indexConstraint.push_back(numVars);

//        for(const auto ic : indexConstraint) cout << ic << ",  "; cout << endl;

//        Bspline bs(table, BsplineType::CUBIC_FREE);
//        for(unsigned int solnr = 0; solnr < xstar.size(); ++solnr)
//        {
//            std::vector<double> xtest = xstar.at(solnr);
//            DenseVector xdv(xtest.size()); for(unsigned int i = 0; i < xtest.size(); ++i) xdv[i] = xtest.at(i);

//            auto y = bs.eval(xdv);
//            cout << "Equation nr  1" << ", solution nr " << solnr << "\t fun(xtest) = " << fun(xtest) << ", bs(xtest) = " << y << endl;
//        }
//        cc->add({new ConstraintBspline(bs, true)}, indexConstraint);
//    }
//    {
//        std::vector<int> indexTable = {0, 1};
//        std::vector<int> numSamples = {ns, ns};
//        auto fun = [=](std::vector<double> x)
//        {
//            auto x1 = x.at(0);
//            auto x2 = x.at(1);
//            return (1 - 0.25/pi)*(std::exp(2*x1) - e) + e*x2/pi - 2*e*x1;
//        };

//        DataTable table;
//        recSampling(lb, ub, indexTable, numSamples, table, fun);

//        std::vector<int> indexConstraint = indexTable;
//        indexConstraint.push_back(numVars);

//        for(const auto ic : indexConstraint) cout << ic << ",  "; cout << endl;

//        Bspline bs(table, BsplineType::CUBIC_FREE);
//        for(unsigned int solnr = 0; solnr < xstar.size(); ++solnr)
//        {
//            std::vector<double> xtest = xstar.at(solnr);
//            DenseVector xdv(xtest.size()); for(unsigned int i = 0; i < xtest.size(); ++i) xdv[i] = xtest.at(i);

//            auto y = bs.eval(xdv);
//            cout << "Equation nr  2" << ", solution nr " << solnr << "\t fun(xtest) = " << fun(xtest) << ", bs(xtest) = " << y << endl;
//        }
//        cc->add({new ConstraintBspline(bs, true)}, indexConstraint);
//    }
//    DenseMatrix cobj(1,dim);
//    cobj.setZero();
//    cobj(0, numVars) = 1;
//    ObjectivePtr obj(new ObjectiveLinear(cobj));
//    cout << "obj = " << cobj << endl;

//    std::vector<int> vt(dim,CONTINUOUS);
//    std::vector<int> bv(numVars); std::iota(begin(bv), end(bv), 0);
//    std::vector<double> z0(dim,0);

//    BB::BranchAndBound bnb(obj, cc, z0, vt, bv);
//    Timer timer;
//    timer.start();
//    bnb.optimize();
//    timer.stop();
//    cout << "Time: " << timer.getMilliSeconds() << " (ms)" << endl;

//    auto sol = bnb.getOptimalSolution();
//    auto solKnown = xstar.at(1);
//    cout << "Solution: " << endl;
//    for(int i = 0; i < numVars; ++i) cout << "x" << i << " = " << sol.at(i) << ", \t Error = " << sol.at(i) - solKnown.at(i) << endl;
//    cout << "Slack variables: " << endl;
//    for(int i = 0; i < numSlacks; ++i) cout << "s" << i << " = " << sol.at(numVars + i) << endl;

//}

//void P_14_1_7()
//{
//    // Test Problem 7, Chapter 14, Floudas 1999 - Handbook of test problems ...
//    std::vector<double> xstar = {0.89999, 0.44999,1.00001,2.00007,7.99997, 7.99969,5.00003,0.99999,2.00005, 0, 0, 0, 0, 0, 0, 0, 0, 0};

//    std::vector<std::vector<double>> g =
//    {{0.4850, 0.7520, 0.8690, 0.9820},
//    {0.3690, 1.2540, 0.7030, 1.4550},
//    {5.2095, 10.0677, 22.9274, 20.2153},
//    {23.3037, 101.7790, 111.4610, 191.2670},
//    {28.5132, 111.8467, 134.3884, 211.4823}};

//    cout << "Test problem P_14_1_7" << endl;

//    int numX = 9;
//    int numK = 4;
//    int numS = 2*numK + 1;

//    int dim = numX + numS;

//    // x1, ..., x9, s1, ..., s9
//    double xlb = 0;
//    double xub = 10;
//    std::vector<double> lb(dim, xlb);
//    std::vector<double> ub(dim, xub);
////    for(unsigned int i = 0; i < numX; ++i)
////    {
////        lb.at(i) = xstar.at(i) - 0.5;
////        ub.at(i) = xstar.at(i) + 0.5;
////    }

//    ConstraintCompositePtr cc(new ConstraintComposite(dim, lb, ub));
//    int ns = 6;
//    {
//        std::vector<int> indexTable = {0, 1, 2, 4, 6, 7};
//        std::vector<int> numSamples = {4, 4, 4, ns, ns, ns};



//        for(unsigned int k = 0; k < numK; ++k)
//        {
//            auto fun = [g, k](std::vector<double> x)
//            {
//                int j = -1;
//                double x1 = x.at(++j);
//                double x2 = x.at(++j);
//                double x3 = x.at(++j);
//                double x5 = x.at(++j);
//                double x7 = x.at(++j);
//                double x8 = x.at(++j);
//                double g1 = g.at(0).at(k);
//                double g3 = g.at(2).at(k);
//                double g4 = g.at(3).at(k);
//                double g5 = g.at(4).at(k);

//                auto curly = exp(x5*(g1 - g3*x7/1000- g5*x8/1000)) - 1;
//                return (1 - x1*x2)*x3 * curly  - g5 + g4*x2;
//            };


//            DataTable table;
//            recSampling(lb, ub, indexTable, numSamples, table, fun);

//            std::vector<int> indexConstraint = indexTable;
//            indexConstraint.push_back(numX + k);

//            //            for(const auto ic : indexConstraint) cout << ic << ",  "; cout << endl;


//            std::vector<double> xtest = {0.89999,0.44999,1.00001,7.99997,5.00003,0.99999};
//            //            for(const auto it : indexTable) xtest.push_back( xstar.at(it) );
//            cout << "C1: k = " << k << "\t fun(xtest) = " << fun(xtest) << endl;
//            Bspline bs(table, BsplineType::CUBIC_FREE);
//            cc->add({new ConstraintBspline(bs, true)}, indexConstraint);
//        }
//    }
//    {
//        std::vector<int> indexTable = {0, 1, 3, 5, 6, 8};
//        std::vector<int> numSamples = {4, 4, 4, ns, ns, ns};
//        for(int k = 0; k < numK; ++k)
//        {
//            auto fun = [g, k](std::vector<double> x)
//            {
//                int j = -1;
//                double x1 = x.at(++j);
//                double x2 = x.at(++j);
//                double x4 = x.at(++j);
//                double x6 = x.at(++j);
//                double x7 = x.at(++j);
//                double x9 = x.at(++j);
//                double g1 = g.at(0).at(k);
//                double g2 = g.at(1).at(k);
//                double g3 = g.at(2).at(k);
//                double g4 = g.at(3).at(k);
//                double g5 = g.at(4).at(k);

//                auto curly = exp(x6*(g1 -g2 -g3*x7*1e-3 + g4*x9*1e-3)) - 1;
//                return (1 - x1*x2)*x4 * curly  - g5*x1 + g4;
//            };
//            DataTable table;
//            recSampling(lb, ub, indexTable, numSamples, table, fun);

//            std::vector<int> indexConstraint = indexTable;
//            indexConstraint.push_back(numX + numK + k);


//            std::vector<double> xtest;
//            for(const auto it : indexTable) xtest.push_back( xstar.at(it) );
//            //            for(const auto ic : indexConstraint) cout << ic << ",  "; cout << endl;
//            cout << "C2: k = " << k << "\t fun(xtest) = " << fun(xtest) << endl;
//            Bspline bs(table, BsplineType::CUBIC_FREE);
//            cc->add({new ConstraintBspline(bs, true)}, indexConstraint);
//        }
//    }

//    {
//        std::vector<int> indexTable = {0, 1, 2, 3};
//        std::vector<int> numSamples = {2, 2, 2, 2};
//        auto fun = [](std::vector<double> x)
//        {
//            int j = -1;
//            double x1 = x.at(++j);
//            double x2 = x.at(++j);
//            double x3 = x.at(++j);
//            double x4 = x.at(++j);
//            return x1*x3 - x2*x4;
//        };
//        DataTable table;
//        recSampling(lb, ub, indexTable, numSamples, table, fun);

//        std::vector<int> indexConstraint = indexTable;
//        indexConstraint.push_back(numX + 2*numK );
//        for(const auto ic : indexConstraint) cout << ic << ",  "; cout << endl;


//        std::vector<double> xtest;
//        for(const auto it : indexTable) xtest.push_back( xstar.at(it) );
//        cout << "C3: " <<"\t fun(xtest) = " << fun(xtest) << endl;
//        Bspline bs(table, BsplineType::LINEAR);
//        cc->add({new ConstraintBspline(bs, true)}, indexConstraint);
//    }

//    DenseMatrix cobj(1,dim); cobj.setZero();
//    for(int i = 0; i < numS; ++i) cobj(0, numX + i) = 1;
//    ObjectivePtr obj(new ObjectiveLinear(cobj));

//    cout << "obj = " << cobj << endl;

//    std::vector<int> vt(dim,CONTINUOUS);
//    std::vector<int> bv(numX);
//    std::iota(begin(bv), end(bv), 0);
//    std::vector<double> z0(dim,0);
//    z0 = {0.89999, 0.44999,1.00001,2.00007,7.99997, 7.99969,5.00003,0.99999,2.00005, 0, 0, 0, 0, 0, 0, 0, 0, 0};

//    BB::BranchAndBound bnb(obj, cc, z0, vt, bv);
//    Timer timer;
//    timer.start();
//    bnb.optimize();
//    timer.stop();
//    cout << "Time: " << timer.getMilliSeconds() << " (ms)" << endl;

//    auto sol = bnb.getOptimalSolution();
//    std::vector<double> solKnown = {0.89999, 0.44999,1.00001,2.00007,7.99997, 7.99969,5.00003,0.99999,2.00005};

//    cout << "Solution: " << endl;
//    for(int i = 0; i < numX; ++i) cout << "x" << i << " = " << sol.at(i) << ", \t Error = " << sol.at(i) - solKnown.at(i) << endl;
//    cout << "Slack variables: " << endl;
//    for(int i = 0; i < numS; ++i) cout << "s" << i << " = " << sol.at(numX + i) << endl;

//}

//void P_14_2_4()
//{
//    double P = 760;
//    double R = 1.98721; //Unused?
//    std::vector<double> a = {16.388,    16.268,     18.607};
//    std::vector<double> b = {2787.50,   2665.54,    3643.31};
//    std::vector<double> c = {229.66,    219.73,     239.73};
//    std::vector<std::vector<double>> Lambda =
//    {{1.0,       0.48,       0.768},
//    {1.55,      1.0,        0.544},
//    {0.566,     0.65,       1.0}};

//    std::vector<std::vector<double>> xstar =
//    { {0.532,     0.468,  0.000},
//    {0.747,     0.000,  0.253},
//    {0.000,     0.677,  0.323},
//    {0.272,     0.465,  0.253}};
//    std::vector<double> Tstar = {55.675,    54.505,     54.356,     54.254};

//    cout << "Test problem P_14_2_4" << endl;

//    int N = 3;
//    int numVar = N + 1;
//    int numS =  N + 1;

//    int dim = numVar + numS;

//    // x1, x2, x3, T
//    double xlb = 0;
//    double xub = 1;
//    std::vector<double> lb(dim, xlb);
//    std::vector<double> ub(dim, xub);
//    lb.at(N) = 50;
//    ub.at(N) = 60;//no bounds on T in book formulation?

//    //    for(unsigned int i = 0; i < numX; ++i)
//    //    {
//    //        lb.at(i) = xstar.at(i) - 0.5;
//    //        ub.at(i) = xstar.at(i) + 0.5;
//    //    }

//    ConstraintCompositePtr cc(new ConstraintComposite(dim, lb, ub));
//    int ns = 6;
//    {
//        std::vector<int> indexTable = {0, 1, 2, 3};
//        std::vector<int> numSamples = {ns, ns, ns, ns};
//        for(int i = 0; i < N; ++i)
//        {
//            auto fun = [=](std::vector<double> x)
//            {
//                //quickfix div by zero
//                auto maxval = *(std::max_element(begin(x), begin(x)+N));
//                double eps = 0.01;
//                if (maxval < eps) for(unsigned int i = 0; i < N; ++i) x.at(i) = eps;

//                //Temperature eq
//                auto T = x.at(N);
//                auto lnGammaI_TempEq = std::log(P) - (a.at(i) - b.at(i)/(T + c.at(i)));

//                //Wilson:
//                // ln sum
//                auto Li = Lambda.at(i);
//                auto sumprod_Lij_xj = std::inner_product(begin(Li), end(Li), begin(x), 0.0);
//                auto lnsum = std::log(std::max(1e-5, sumprod_Lij_xj));
//                //sum of fraction
//                std::vector<double> xjLji(N);
//                std::transform(begin(Lambda), end(Lambda), begin(x), begin(xjLji), [&](std::vector<double> Lj, double xj){ return Lj.at(i)*xj;});
//                std::vector<double> inv_sum_xkLjk(N);
//                std::transform(begin(Lambda), end(Lambda), begin(inv_sum_xkLjk), [&](std::vector<double> Lj){return 1.0 / std::inner_product(begin(Lj), end(Lj), begin(x), 0.0);});
//                auto sumfrac = std::inner_product(begin(xjLji), end(xjLji), begin(inv_sum_xkLjk), 0.0);
//                auto lnGammaI_WilsonEq = 1 - lnsum - sumfrac;



//                return lnGammaI_TempEq - lnGammaI_WilsonEq;
//            };

//            DataTable table;
//            recSampling(lb, ub, indexTable, numSamples, table, fun);

//            std::vector<int> indexConstraint = indexTable;
//            indexConstraint.push_back(numVar + i);

//            for(const auto ic : indexConstraint) cout << ic << ",  "; cout << endl;

//            Bspline bs(table, BsplineType::CUBIC_FREE);
//            for(unsigned int solnr = 0; solnr < Tstar.size(); ++solnr)
//            {
//                std::vector<double> xtest = xstar.at(solnr);
//                xtest.push_back(Tstar.at(solnr));

//                DenseVector xdv(xtest.size());
//                for(unsigned int i = 0; i < xtest.size(); ++i) xdv[i] = xtest.at(i);

//                auto y = bs.eval(xdv);
//                cout << "Equation nr  " << i << ", solution nr " << solnr << "\t fun(xtest) = " << fun(xtest) << ", bs(xtest) = " << y << endl;

//            }
//            cc->add({new ConstraintBspline(bs, true)}, indexConstraint);
//        }
//    }
//    {

//        DenseMatrix A; A.setOnes(1,N);
//        DenseVector b; b.setOnes(1);
//        cc->add({new ConstraintLinear(A, b, true)}, {0, 1, 2});
//    }

//    DenseMatrix cobj(1,dim); cobj.setZero();
//    for(int i = 0; i < numS; ++i) cobj(0, numVar + i) = 1;
//    ObjectivePtr obj(new ObjectiveLinear(cobj));
//    cout << "obj = " << cobj << endl;

//    std::vector<int> vt(dim,CONTINUOUS);
//    std::vector<int> bv(numVar);
//    std::iota(begin(bv), end(bv), 0);
//    std::vector<double> z0(dim,0);

//    BB::BranchAndBound bnb(obj, cc, z0, vt, bv);
//    Timer timer;
//    timer.start();
//    bnb.optimize();
//    timer.stop();
//    cout << "Time: " << timer.getMilliSeconds() << " (ms)" << endl;

//    //    auto sol = bnb.getOptimalSolution();

//    //    cout << "Solution: " << endl;
//    //    for(unsigned int i = 0; i < numVar; ++i) cout << "x" << i << " = " << sol.at(i) << ", \t Error = " << sol.at(i) - solKnown.at(i) << endl;
//    //    cout << "Slack variables: " << endl;
//    //    for(unsigned int i = 0; i < numS; ++i) cout << "s" << i << " = " << sol.at(numX + i) << endl;

//}
/*** END code from Anders S. ***/
