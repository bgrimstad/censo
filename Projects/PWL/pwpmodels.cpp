/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "pwpmodels.h"

#include "OptimizationProblem/constraintset.h"
#include "OptimizationProblem/constraintlinear.h"

#include "OptimizationProblem/constraintbspline.h"
#include "BranchAndBound/branchandbound.h"

SparseMatrix getSelectionMatrix(unsigned int k, unsigned int degree, unsigned int knotIntervals)
{
    assert(k < knotIntervals);
    SparseMatrix S((degree+1)*knotIntervals,degree+1);
    S.reserve(degree+1);

    for (unsigned int i = 0; i < degree + 1; ++i)
    {
        S.insert((degree+1)*k+i,i) = 1;
    }

    return S;
}

template <unsigned int N>
ConstraintPtr constraintPiecewisePolynomialDCC(Variables vars, SPLINTER::BSpline bspline)
{
    unsigned int dim = bspline.getNumVariables();
    assert(dim == N);
    assert(vars.size() == dim + 1);

    // Decompose B-spline
    PiecewisePolynomial<N> pwp = decomposeBSplinePower<N>(bspline);
    auto polynomials = pwp.polynomials;
    auto polytopes = pwp.polytopes;

    // Compute number of variables
    // Order of variables: [x z x_P z_P lambda_{P,v} y_P]
    int npol = polynomials.size(); // Number of polytopes
    int nl = 0; // Number of auxiliary (cont)
    int ny = npol; // Number of binary (y_P)

    // N+1 auxiliary variables for each polytope (x_P and z_P)
    nl += (N+1)*npol;

    // One auxilliary variable for each vertex in each polytope (lambda_{P,v})
    unsigned int nlambda = std::pow(2,N)*npol; // 2^N vertices for hyperrectangles
    nl += nlambda;

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

    auto cs = std::make_shared<ConstraintSet>(vars);

    // Update polynomial variables
    for (unsigned int i = 0; i < npol; ++i)
    {
        Variables pVars;
        for (unsigned int j = 0; j < dim; ++j)
            pVars.push_back(vars.at(dim+1+i*dim+j)); // x_P
        pVars.push_back(vars.at(dim+1+npol*dim+i)); // z_P

        polynomials[i]->setVariables(pVars);
    }

    {
        // Add constraint sum_P x_P = x
        DenseMatrix A = DenseMatrix::Zero(dim, dim+dim*npol);
        DenseMatrix I = DenseMatrix::Identity(dim,dim);
        A.block(0,0,dim,dim) = -I;

        for (int k = 1; k <= npol; ++k)
            A.block(0,k*dim,dim,dim) = I;

        DenseVector b = DenseMatrix::Zero(dim,1);

        std::vector<VariablePtr> cvars;
        for (int i = 0; i < dim; ++i)
            cvars.push_back(vars.at(i)); // x
        for (int i = 0; i < dim*npol; ++i)
            cvars.push_back(vars.at(dim+1+i)); // x_P for all P

        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, true);
        cs->add(lincon);
    }

    {
        // Add constraints sum_v lambda_{P,v} v = x_P for all P
        int ip = 0; // Polytope counter
        int il = 0; // Lambda counter
        for (const auto &p : polytopes)
        {
            int nv = p.vertices.size();

            DenseMatrix A = DenseMatrix::Ones(dim, dim+nv);
            DenseMatrix I = DenseMatrix::Identity(dim,dim);
            A.block(0,0,dim,dim) = -I;
            DenseVector b = DenseMatrix::Zero(dim,1);

            std::vector<VariablePtr> cvars;
            for (int k = 0; k < dim; ++k)
                cvars.push_back(vars.at(dim+1+dim*ip+k)); // x_P
            for (int k = 0; k < nv; ++k)
                cvars.push_back(vars.at((dim+1)*(npol+1)+il+k)); // lambda_{P,v}

            int k = 0;
            for (const auto &v : p.vertices)
            {
                // Vertex
                auto x = v.x;
                for (int l = 0; l < dim; l++)
                    A(l,dim+k) = x.at(l);
                ++k;
            }

            ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, true);
            cs->add(lincon);

            il += nv;
            ++ip;
        }

        assert(il == nlambda);
    }

    {
        // Add constraint sum_P y_P*z_P <= z
        DenseMatrix E = DenseMatrix::Zero(npol+1, 2*npol+1);
        for (int k = 0; k < npol; ++k)
        {
            E(k, k) = 1; // z_P
            E(k, npol+k) = 1; // y_p
        }
        E(npol,2*npol) = 1; // z

        DenseVector c = DenseVector::Ones(npol+1);
        c(npol) = -1; // -z

        std::vector<VariablePtr> cvars;
        for (int k = 0; k < npol; ++k)
            cvars.push_back(vars.at(dim+1+dim*npol+k)); // z_P
        for (int k = 0; k < npol; ++k)
            cvars.push_back(vars.at(dim+1+nl+k)); // y_P
        cvars.push_back(vars.at(dim)); // z

        ConstraintPtr lincon = std::make_shared<ConstraintPolynomial>(cvars, c, E, false);
        cs->add(lincon);
    }

    {
        // Add polynomial constraints f_P(x_P) <= z_P for all P
        for (auto &pol : polynomials)
            cs->add(pol);
    }

    {
        // lambda_{P,v} >= 0
        for (int k = 0; k < nlambda; ++k)
            vars.at((dim+1)*(npol+1)+k)->setLowerBound(0);
    }

    {
        // Add constraint sum_v lambda_{P,v} = y_P for all P

        int ip = 0; // polytope counter
        int il = 0; // lambda counter

        // All polytopes have an equal number of vertices (2^N)
        assert(polytopes.front().vertices.size() == std::pow(2,N));
        DenseMatrix A = DenseMatrix::Ones(1,1+polytopes.front().vertices.size());
        A(0,0) = -1;
        DenseVector b = DenseMatrix::Zero(1,1);

        for (const auto &p : polytopes)
        {
            std::vector<VariablePtr> cvars;
            cvars.push_back(vars.at(dim+1+nl+ip)); // y_P

            for (const auto &v : p.vertices)
            {
                cvars.push_back(vars.at((dim+1)*(npol+1)+il)); // lambda_{P,v}
                il++;
            }

            // Check if preallocated matrix can be used
            assert(A.cols() == 1+p.vertices.size());

            ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, true);
            cs->add(lincon);

            ip++;
        }

        assert(nlambda == il);
    }

    {
        // Add constraint sum_P y_P = 1
        DenseMatrix A = DenseMatrix::Ones(1, ny);
        DenseVector b(1); b(0) = 1;

        std::vector<VariablePtr> cvars;
        for (int k = 0; k < ny; ++k)
            cvars.push_back(vars.at(dim+1+nl+k)); // y_p

        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, true);
        cs->add(lincon);
    }

    return cs;
}

template <unsigned int N>
ConstraintPtr constraintPiecewisePolynomialDCC2(Variables vars, SPLINTER::BSpline bspline)
{
    unsigned int dim = bspline.getNumVariables();
    assert(dim == N);
    assert(vars.size() == dim + 1);

    // Decompose B-spline
    PiecewisePolynomial<N> pwp = decomposeBSplinePower<N>(bspline);
    auto polynomials = pwp.polynomials;
    auto polytopes = pwp.polytopes;

    // Compute number of variables
    // Order of variables: [x z x_P z_P lambda_{P,v} y_P]
    int npol = polynomials.size(); // Number of polytopes
    int nl = 0; // Number of auxiliary (cont)
    int ny = npol; // Number of binary (y_P)

    // 1 auxiliary variable (z_P) for each polytope
    nl += npol;

    // One auxilliary variable for each vertex in each polytope (lambda_{P,v})
    unsigned int nlambda = std::pow(2,N)*npol; // 2^N vertices for hyperrectangles
    nl += nlambda;

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

    auto cs = std::make_shared<ConstraintSet>(vars);

    // Update polynomial variables
    for (unsigned int i = 0; i < npol; ++i)
    {
        Variables pVars;
        for (unsigned int j = 0; j < dim; ++j)
            pVars.push_back(vars.at(j)); // x
        pVars.push_back(vars.at(dim+1+i)); // z_P

        polynomials[i]->setVariables(pVars);
    }

    {
        // Add constraint sum_P sum_v lambda_{P,v} v = x
        DenseMatrix A = DenseMatrix::Ones(dim, dim+nlambda);
        DenseMatrix I = DenseMatrix::Identity(dim,dim);
        A.block(0,0,dim,dim) = -I;
        DenseVector b = DenseMatrix::Zero(dim,1);

        std::vector<VariablePtr> cvars;
        for (int i = 0; i < dim; i++)
            cvars.push_back(vars.at(i)); // x
        for (int i = 0; i < nlambda; i++)
            cvars.push_back(vars.at(dim+1+npol+i)); // lambdas

        int k = 0;
        for (const auto &p : polytopes)
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
        // Add constraint sum_P y_P*z_P <= z
        DenseMatrix E = DenseMatrix::Zero(npol+1, 2*npol+1);
        for (int k = 0; k < npol; ++k)
        {
            E(k, k) = 1; // z_P
            E(k, npol+k) = 1; // y_p
        }
        E(npol,2*npol) = 1; // z

        DenseVector c = DenseVector::Ones(npol+1);
        c(npol) = -1; // -z

        std::vector<VariablePtr> cvars;
        for (int k = 0; k < npol; ++k)
            cvars.push_back(vars.at(dim+1+k)); // z_P
        for (int k = 0; k < npol; ++k)
            cvars.push_back(vars.at(dim+1+nl+k)); // y_P
        cvars.push_back(vars.at(dim)); // z

        ConstraintPtr lincon = std::make_shared<ConstraintPolynomial>(cvars, c, E, false);
        cs->add(lincon);
    }

    {
        // Add polynomial constraints f_P(x_P) <= z_P for all P
        for (auto &pol : polynomials)
            cs->add(pol);
    }

    {
        // lambda_{P,v} >= 0
        for (int k = 0; k < nlambda; ++k)
            vars.at(dim+1+npol+k)->setLowerBound(0);
    }

    {
        // Add constraint sum_v lambda_{P,v} = y_P for all P

        int ip = 0; // polytope counter
        int il = 0; // lambda counter

        // All polytopes have an equal number of vertices (2^N)
        assert(polytopes.front().vertices.size() == std::pow(2,N));
        DenseMatrix A = DenseMatrix::Ones(1,1+polytopes.front().vertices.size());
        A(0,0) = -1;
        DenseVector b = DenseMatrix::Zero(1,1);

        for (const auto &p : polytopes)
        {
            std::vector<VariablePtr> cvars;
            cvars.push_back(vars.at(dim+1+nl+ip)); // y_P

            for (const auto &v : p.vertices)
            {
                cvars.push_back(vars.at(dim+1+npol+il)); // lambda_{P,v}
                il++;
            }

            // Check if preallocated matrix can be used
            assert(A.cols() == 1+p.vertices.size());

            ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, true);
            cs->add(lincon);

            ip++;
        }

        assert(nlambda == il);
    }

    {
        // Add constraint sum_P y_P = 1
        DenseMatrix A = DenseMatrix::Ones(1, ny);
        DenseVector b(1); b(0) = 1;

        std::vector<VariablePtr> cvars;
        for (int k = 0; k < ny; ++k)
            cvars.push_back(vars.at(dim+1+nl+k)); // y_p

        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, true);
        cs->add(lincon);
    }

    return cs;
}

template <unsigned int N>
ConstraintPtr constraintPiecewisePolynomialDLog(Variables vars, SPLINTER::BSpline bspline)
{
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

    unsigned int dim = bspline.getNumVariables();
    assert(dim == N);
    assert(vars.size() == dim + 1);

    // Decompose B-spline
    PiecewisePolynomial<N> pwp = decomposeBSplinePower<N>(bspline);
    auto polynomials = pwp.polynomials;
    auto polytopes = pwp.polytopes;

    // Compute number of variables
    // Order of variables: [x z x_P x_P z_P Lambda_P lambda_{P,v} y_P]
    int npol = polytopes.size(); // Number of polytopes
    int ncaux = 0; // Number of auxiliary (cont)
    int nbaux = std::ceil(std::log2(npol)); // Logarithmic number of binary variables (y_P)

    // 2 auxiliary variables (z_P and Lambda_P) for each polytope
    ncaux += (dim+2)*npol;

    // One auxilliary variable for each vertex in each polytope (lambda_{P,v})
    unsigned int nlambda = std::pow(2,N)*npol; // 2^N vertices for hyperrectangles
    ncaux += nlambda;

    // Create variables
    for (int i = 0; i < ncaux; ++i)
    {
        auto var = std::make_shared<Variable>(0);
        vars.push_back(var);
    }
    for (int i = 0; i < nbaux; ++i)
    {
        auto var = std::make_shared<Variable>(0, 0, 1, VariableType::BINARY);
        vars.push_back(var);
    }

    auto cs = std::make_shared<ConstraintSet>(vars);

    // Update polynomial variables
    for (unsigned int i = 0; i < npol; ++i)
    {
        Variables pVars;
        for (unsigned int j = 0; j < dim; ++j)
            pVars.push_back(vars.at(dim+1+dim*i+j)); // x_P
        pVars.push_back(vars.at(dim+1+dim*npol+i)); // z_P

        polynomials[i]->setVariables(pVars);
    }

    {
        // Add constraint sum_P x_P = x
        DenseMatrix A = DenseMatrix::Zero(dim, dim+dim*npol);
        DenseMatrix I = DenseMatrix::Identity(dim,dim);
        A.block(0,0,dim,dim) = -I;

        for (int k = 1; k <= npol; ++k)
            A.block(0,k*dim,dim,dim) = I;

        DenseVector b = DenseMatrix::Zero(dim,1);

        std::vector<VariablePtr> cvars;
        for (int i = 0; i < dim; ++i)
            cvars.push_back(vars.at(i)); // x
        for (int i = 0; i < dim*npol; ++i)
            cvars.push_back(vars.at(dim+1+i)); // x_P for all P

        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, true);
        cs->add(lincon);
    }

    {
        // Add constraints sum_v lambda_{P,v} v = x_P for all P
        int ip = 0; // Polytope counter
        int il = 0; // Lambda counter
        for (const auto &p : polytopes)
        {
            int nv = p.vertices.size();

            DenseMatrix A = DenseMatrix::Ones(dim, dim+nv);
            DenseMatrix I = DenseMatrix::Identity(dim,dim);
            A.block(0,0,dim,dim) = -I;
            DenseVector b = DenseMatrix::Zero(dim,1);

            std::vector<VariablePtr> cvars;
            for (int k = 0; k < dim; ++k)
                cvars.push_back(vars.at(dim+1+dim*ip+k)); // x_P
            for (int k = 0; k < nv; ++k)
                cvars.push_back(vars.at(dim+1+(dim+2)*npol+il+k)); // lambda_{P,v}

            int k = 0;
            for (const auto &v : p.vertices)
            {
                // Vertex
                auto x = v.x;
                for (int l = 0; l < dim; l++)
                    A(l,dim+k) = x.at(l);
                ++k;
            }

            ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, true);
            cs->add(lincon);

            il += nv;
            ++ip;
        }

        assert(il == nlambda);
    }

    {
        // Add constraint sum_P Lambda_P*z_P <= z
        DenseMatrix E = DenseMatrix::Zero(npol+1, 2*npol+1);
        for (int k = 0; k < npol; ++k)
        {
            E(k, k) = 1; // z_P
            E(k, npol+k) = 1; // L_p
        }
        E(npol,2*npol) = 1; // z

        DenseVector c = DenseVector::Ones(npol+1);
        c(npol) = -1; // -z

        std::vector<VariablePtr> cvars;
        for (int k = 0; k < npol; ++k)
            cvars.push_back(vars.at(dim+1+dim*npol+k)); // z_P
        for (int k = 0; k < npol; ++k)
            cvars.push_back(vars.at(dim+1+(dim+1)*npol+k)); // Lambda_P
        cvars.push_back(vars.at(dim)); // z

        ConstraintPtr lincon = std::make_shared<ConstraintPolynomial>(cvars, c, E, false);
        cs->add(lincon);
    }

    {
        // Add polynomial constraints f_P(x) <= z_P for all P
        for (auto &pol : polynomials)
            cs->add(pol);
    }

    {
        // lambda_{P,v} >= 0
        for (int k = 0; k < nlambda; ++k)
            vars.at(dim+1+(dim+2)*npol+k)->setLowerBound(0);

        // L_P >= 0
        for (int k = 0; k < npol; ++k)
            vars.at(dim+1+(dim+1)*npol+k)->setLowerBound(0);
    }

    {
        // Add constraint sum_v lambda_{P,v} = L_P for all P
        int ip = 0;
        int iv = 0;
        for (const auto &p : polytopes)
        {
            int nvert = p.vertices.size();

            DenseMatrix A = DenseMatrix::Ones(1, 1+nvert);
            A(0,0) = -1;
            DenseVector b(1); b(0) = 0;

            std::vector<VariablePtr> cvars;
            cvars.push_back(vars.at(dim+1+(dim+1)*npol+ip)); // L_P
            for (int i = 0; i < nvert; i++)
                cvars.push_back(vars.at(dim+1+(dim+2)*npol+iv+i)); // lambda_{P,v}

            ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, true);
            cs->add(lincon);

            ip++;
            iv += nvert;
        }
        assert(iv == nlambda);
    }

    {
        // Add constraint sum_P L_{P} = 1
        DenseMatrix A = DenseMatrix::Ones(1, npol);
        DenseVector b(1); b(0) = 1;

        std::vector<VariablePtr> cvars;
        for (int i = 0; i < npol; ++i)
            cvars.push_back(vars.at(dim+1+(dim+1)*npol+i)); // L_P

        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, true);
        cs->add(lincon);
    }

    {
        // For all l in L(P) = {1,...,ceil(log2(|P|))}
        // add constraint sum_P+ sum_v lambda_{P,v} <= y_l
        for (int l = 0; l < nbaux; l++)
        {
            std::vector<VariablePtr> cvars;
            cvars.push_back(vars.at(dim+1+ncaux+l)); // y_l

            for (unsigned int i = 0; i < polytopes.size(); i++)
            {
                if (!isDigitOne(i,l)) continue;

                cvars.push_back(vars.at(dim+1+(dim+1)*npol+i)); // Lambda_P
            }

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
        for (int l = 0; l < nbaux; l++)
        {
            std::vector<VariablePtr> cvars;
            cvars.push_back(vars.at(dim+1+ncaux+l)); // y_l

            for (unsigned int i = 0; i < polytopes.size(); i++)
            {
                if (isDigitOne(i,l)) continue;

                cvars.push_back(vars.at(dim+1+(dim+1)*npol+i)); // Lambda_P
            }

            DenseMatrix A = DenseMatrix::Ones(1,cvars.size());
            A(0,0) = 1;
            DenseVector b = DenseMatrix::Ones(1,1);

            ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, false);
            cs->add(lincon);
        }
    }

    return cs;
}

template <unsigned int N>
ConstraintPtr constraintPiecewisePolynomialDLog2(Variables vars, SPLINTER::BSpline bspline)
{
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

    unsigned int dim = bspline.getNumVariables();
    assert(dim == N);
    assert(vars.size() == dim + 1);

    // Decompose B-spline
    PiecewisePolynomial<N> pwp = decomposeBSplinePower<N>(bspline);
    auto polynomials = pwp.polynomials;
    auto polytopes = pwp.polytopes;

    // Compute number of variables
    // Order of variables: [x z x_P z_P Lambda_P lambda_{P,v} y_P]
    int npol = polytopes.size(); // Number of polytopes
    int ncaux = 0; // Number of auxiliary (cont)
    int nbaux = std::ceil(std::log2(npol)); // Logarithmic number of binary variables (y_P)

    // 2 auxiliary variables (z_P and Lambda_P) for each polytope
    ncaux += 2*npol;

    // One auxilliary variable for each vertex in each polytope (lambda_{P,v})
    unsigned int nlambda = std::pow(2,N)*npol; // 2^N vertices for hyperrectangles
    ncaux += nlambda;

    // Create variables
    for (int i = 0; i < ncaux; ++i)
    {
        auto var = std::make_shared<Variable>(0);
        vars.push_back(var);
    }
    for (int i = 0; i < nbaux; ++i)
    {
        auto var = std::make_shared<Variable>(0, 0, 1, VariableType::BINARY);
        vars.push_back(var);
    }

    auto cs = std::make_shared<ConstraintSet>(vars);

    // Update polynomial variables
    for (unsigned int i = 0; i < npol; ++i)
    {
        Variables pVars;
        for (unsigned int j = 0; j < dim; ++j)
            pVars.push_back(vars.at(j)); // x
        pVars.push_back(vars.at(dim+1+i)); // z_P

        polynomials[i]->setVariables(pVars);
    }

    {
        // Add constraint sum_P sum_v lambda_{P,v} v = x
        DenseMatrix A = DenseMatrix::Ones(dim, dim+nlambda);
        DenseMatrix I = DenseMatrix::Identity(dim,dim);
        A.block(0,0,dim,dim) = -I;
        DenseVector b = DenseMatrix::Zero(dim,1);

        std::vector<VariablePtr> cvars;
        for (int i = 0; i < dim; i++)
            cvars.push_back(vars.at(i)); // x
        for (int i = 0; i < nlambda; i++)
            cvars.push_back(vars.at(dim+1+2*npol+i)); // lambda_{P,v}

        int k = 0;
        for (const auto &p : polytopes)
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
        // Add constraint sum_P Lambda_P*z_P <= z
        DenseMatrix E = DenseMatrix::Zero(npol+1, 2*npol+1);
        for (int k = 0; k < npol; ++k)
        {
            E(k, k) = 1; // z_P
            E(k, npol+k) = 1; // L_p
        }
        E(npol,2*npol) = 1; // z

        DenseVector c = DenseVector::Ones(npol+1);
        c(npol) = -1; // -z

        std::vector<VariablePtr> cvars;
        for (int k = 0; k < npol; ++k)
            cvars.push_back(vars.at(dim+1+k)); // z_P
        for (int k = 0; k < npol; ++k)
            cvars.push_back(vars.at(dim+1+npol+k)); // Lambda_P
        cvars.push_back(vars.at(dim)); // z

        ConstraintPtr lincon = std::make_shared<ConstraintPolynomial>(cvars, c, E, false);
        cs->add(lincon);
    }

    {
        // Add polynomial constraints f_P(x) <= z_P for all P
        for (auto &pol : polynomials)
            cs->add(pol);
    }

    {
        // lambda_{P,v} >= 0
        for (int k = 0; k < nlambda; ++k)
            vars.at(dim+1+2*npol+k)->updateBounds(0,1);

        // L_P >= 0
        for (int k = 0; k < npol; ++k)
            vars.at(dim+1+npol+k)->updateBounds(0,1);
    }

    {
        // Add constraint sum_v lambda_{P,v} = L_P for all P
        int ip = 0;
        int iv = 0;
        for (const auto &p : polytopes)
        {
            int nvert = p.vertices.size();

            DenseMatrix A = DenseMatrix::Ones(1, 1+nvert);
            A(0,0) = -1;
            DenseVector b(1); b(0) = 0;

            std::vector<VariablePtr> cvars;
            cvars.push_back(vars.at(dim+1+npol+ip)); // L_P
            for (int i = 0; i < nvert; i++)
                cvars.push_back(vars.at(dim+1+2*npol+iv+i)); // lambda_{P,v}

            ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, true);
            cs->add(lincon);

            ip++;
            iv += nvert;
        }
        assert(iv == nlambda);
    }

    {
        // Add constraint sum_P L_{P} = 1
        DenseMatrix A = DenseMatrix::Ones(1, npol);
        DenseVector b(1); b(0) = 1;

        std::vector<VariablePtr> cvars;
        for (int i = 0; i < npol; ++i)
            cvars.push_back(vars.at(dim+1+npol+i)); // L_P

        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, true);
        cs->add(lincon);
    }

    {
        // For all l in L(P) = {1,...,ceil(log2(|P|))}
        // add constraint sum_P+ sum_v lambda_{P,v} <= y_l
        for (int l = 0; l < nbaux; l++)
        {
            std::vector<VariablePtr> cvars;
            cvars.push_back(vars.at(dim+1+ncaux+l)); // y_l

            for (unsigned int i = 0; i < polytopes.size(); i++)
            {
                if (!isDigitOne(i,l)) continue;

                cvars.push_back(vars.at(dim+1+npol+i)); // Lambda_P
            }

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
        for (int l = 0; l < nbaux; l++)
        {
            std::vector<VariablePtr> cvars;
            cvars.push_back(vars.at(dim+1+ncaux+l)); // y_l

            for (unsigned int i = 0; i < polytopes.size(); i++)
            {
                if (isDigitOne(i,l)) continue;

                cvars.push_back(vars.at(dim+1+npol+i)); // Lambda_P
            }

            DenseMatrix A = DenseMatrix::Ones(1,cvars.size());
            A(0,0) = 1;
            DenseVector b = DenseMatrix::Ones(1,1);

            ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, false);
            cs->add(lincon);
        }
    }

    return cs;
}

template <unsigned int N>
ConstraintPtr constraintPiecewisePolynomialDLog3(Variables vars, SPLINTER::BSpline bspline)
{
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

    unsigned int dim = bspline.getNumVariables();
    assert(dim == N);
    assert(vars.size() == dim + 1);

    // Decompose B-spline
    PiecewisePolynomial<N> pwp = decomposeBSplinePower<N>(bspline);
    auto polynomials = pwp.polynomials;
    auto polytopes = pwp.polytopes;

    // Get monomial powers
    DenseMatrix M = getPowersMatrix(bspline.getBasisDegrees());

    // Compute monomial bounds xi^L and xi^U
    DenseVector xil = DenseVector::Ones(M.rows());
    DenseVector xiu = DenseVector::Ones(M.rows());

    for (int i = 0; i < M.rows(); ++i)
    {
        bool even = true;

        double ximax = 1;

        for (int j = 0; j < M.cols(); ++j)
        {
            int pow = M(i,j);

            if (pow < 1)
                continue;

            if (pow % 2 != 0)
                even = false;

            double xjl = vars.at(j)->getLowerBound();
            double xju = vars.at(j)->getUpperBound();

            ximax *= std::pow(std::max(xjl, xju), pow);
        }

        xiu(i) = ximax;
        xil(i) = -ximax;

        if (even) xil(i) = 0;
    }

    // Get number of polytopes
    int npol = polytopes.size(); // Number of polytopes

    // Variables: [x z z_P Lambda_P lambda_{P,v} xi phi_P y_P]
    int nz = 1;
    int iz = N;

    int nzp = npol;
    int izp = iz+nz;

    int nLambda = npol;
    int iLambda = izp+nzp;

    int nlambda = std::pow(2,N)*npol; // 2^N vertices for hyperrectangles
    int ilambda = iLambda+nLambda;

    int nxi = M.rows();
    int ixi = ilambda+nlambda;

    int nphi = npol;
    int iphi = ixi+nxi;

    int ny = std::ceil(std::log2(npol)); // Logarithmic number of binary variables (y_P)
    int iy = iphi + nphi;

    int ncaux = nzp + nlambda + nLambda + nxi + nphi; // Number of aux. cont. variables
    int nbaux = ny; // Number of aux. binary variables

    // Create variables
    for (int i = 0; i < ncaux; ++i)
    {
        auto var = std::make_shared<Variable>(0);
        vars.push_back(var);
    }
    for (int i = 0; i < nbaux; ++i)
    {
        auto var = std::make_shared<Variable>(0, 0, 1, VariableType::BINARY);
        vars.push_back(var);
    }

    auto cs = std::make_shared<ConstraintSet>(vars);

    {
        // Update polynomial variables
        Variables pVars;
        for (unsigned int j = 0; j < dim; ++j)
            pVars.push_back(vars.at(j)); // x
        pVars.push_back(vars.at(iz)); // z

        for (auto pol : polynomials)
        {
            pol->setVariables(pVars);
        }
    }

    {
        // xil <= xi <= xiu
        for (int i = 0; i < nxi; ++i)
            vars.at(ixi+i)->updateBounds(xil(i), xiu(i));
    }

    {
        // Add constraint sum_P sum_v lambda_{P,v} v = x
        DenseMatrix A = DenseMatrix::Ones(N, N+nlambda);
        DenseMatrix I = DenseMatrix::Identity(N,N);
        A.block(0,0,N,N) = -I;
        DenseVector b = DenseMatrix::Zero(N,1);

        std::vector<VariablePtr> cvars;
        for (int i = 0; i < N; i++)
            cvars.push_back(vars.at(i)); // x
        for (int i = 0; i < nlambda; i++)
            cvars.push_back(vars.at(ilambda+i)); // lambda_{P,v}

        int k = 0;
        for (const auto &p : polytopes)
        {
            for (const auto &v : p.vertices)
            {
                // Vertex
                auto x = v.x;
                for (int l = 0; l < N; l++)
                    A(l,N+k) = x.at(l);
                k++;
            }
        }

        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, true);
        cs->add(lincon);
    }

    {
        // Add constraint sum_P phi_P <= z
        DenseMatrix A = DenseMatrix::Ones(1,1+nphi);
        A(0,0) = -1;
        DenseVector b(1); b(0) = 0;

        Variables cvars;
        cvars.push_back(vars.at(iz));
        for (int i = 0; i < nphi; ++i)
            cvars.push_back(vars.at(iphi+i));

        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, false);
        cs->add(lincon);
    }

    {
        // Add constraint Lambda_P*z_P = phi_P for all P
        for (int i = 0; i < npol; ++i)
        {
            Variables cvars = {
                vars.at(iphi+i), // phi_P
                vars.at(iLambda+i), // Lambda_P
                vars.at(izp+i) // z_P
            };

            // Create monomial xi = m(x)
            DenseVector ci = DenseVector::Ones(2);
            ci(0) = -1;
            DenseMatrix Ei = DenseMatrix::Zero(2,3);
            Ei(0,0) = 1;
            Ei(1,1) = 1;
            Ei(1,2) = 1;

            ConstraintPtr con = std::make_shared<ConstraintPolynomial>(cvars, ci, Ei, true);
            cs->add(con);
        }
    }

    {
        // Add constraint z_P = a_P*xi for all P
        int ipol = 0;
        for (auto &pol : polynomials)
        {
            DenseMatrix a = pol->getCoefficients();
            assert(a.rows() == nxi+1); // Last coefficient related to z
            assert(a.cols() == 1);
            a.transposeInPlace();

            DenseVector b(1); b(0) = 0;

            Variables cvars;

            for (int i = 0; i < nxi; ++i)
                cvars.push_back(vars.at(ixi+i));
            cvars.push_back(vars.at(izp+ipol));

            ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, a, b, true);
            //cout << "Before: " << vars.at(izp+ipol)->getLowerBound() << " / " << vars.at(izp+ipol)->getUpperBound() << endl;
            assert(lincon->reduceVariableRanges()); // Set bounds on zp
            //cout << "After: " << vars.at(izp+ipol)->getLowerBound() << " / " << vars.at(izp+ipol)->getUpperBound() << endl;
            cs->add(lincon);

            ++ipol;
        }
    }

    {
        // Add monomials
        for (int i = 0; i < nxi; ++i)
        {
            Variables cvars = {
                vars.at(ixi+i) // xi
            };

            std::vector<unsigned int> powers;

            for (int j = 0; j < N; ++j)
            {
                if (M(i,j) > 0)
                {
                    cvars.push_back(vars.at(j)); // x
                    powers.push_back(M(i,j));
                }
            }

            if (powers.size() > 0)
            {
                // Create monomial xi = m(x)
                DenseVector ci = DenseVector::Ones(2);
                ci(0) = -1;
                DenseMatrix Ei = DenseMatrix::Zero(2,1+powers.size());
                Ei(0,0) = 1;

                for (int j = 0; j < powers.size(); ++j)
                {
                    Ei(1,1+j) = powers[j];
                }

                auto mono = std::make_shared<ConstraintPolynomial>(cvars, ci, Ei, true);

                cs->add(mono);
            }
            else
            {
                // Fix variable xi to 1
                cvars.at(0)->setLowerBound(1);
                cvars.at(0)->setUpperBound(1);
            }
        }
    }

    {
        // lambda_{P,v} >= 0
        for (int i = 0; i < nlambda; ++i)
            vars.at(ilambda+i)->updateBounds(0,1);

        // L_P >= 0
        for (int i = 0; i < nLambda; ++i)
            vars.at(iLambda+i)->updateBounds(0,1);

        // xil <= xi <= xiu
        for (int i = 0; i < nxi; ++i)
            vars.at(ixi+i)->updateBounds(xil(i), xiu(i));
    }

    {
        // Add constraint sum_v lambda_{P,v} = L_P for all P
        int ip = 0;
        int iv = 0;
        for (const auto &p : polytopes)
        {
            int nvert = p.vertices.size();

            DenseMatrix A = DenseMatrix::Ones(1, 1+nvert);
            A(0,0) = -1;
            DenseVector b(1); b(0) = 0;

            std::vector<VariablePtr> cvars;
            cvars.push_back(vars.at(iLambda+ip)); // L_P
            for (int i = 0; i < nvert; i++)
                cvars.push_back(vars.at(ilambda+iv+i)); // lambda_{P,v}

            ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, true);
            cs->add(lincon);

            ip++;
            iv += nvert;
        }
        assert(iv == nlambda);
    }

    {
        // Add constraint sum_P L_{P} = 1
        DenseMatrix A = DenseMatrix::Ones(1, nLambda);
        DenseVector b(1); b(0) = 1;

        std::vector<VariablePtr> cvars;
        for (int i = 0; i < nLambda; ++i)
            cvars.push_back(vars.at(iLambda+i)); // L_P

        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, true);
        cs->add(lincon);
    }

    {
        // For all l in L(P) = {1,...,ceil(log2(|P|))}
        // add constraint sum_P+ sum_v lambda_{P,v} <= y_l
        for (int l = 0; l < ny; l++)
        {
            std::vector<VariablePtr> cvars;
            cvars.push_back(vars.at(iy+l)); // y_l

            for (unsigned int i = 0; i < npol; i++)
            {
                if (!isDigitOne(i,l)) continue;

                cvars.push_back(vars.at(iLambda+i)); // Lambda_P
            }

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
            cvars.push_back(vars.at(iy+l)); // y_l

            for (unsigned int i = 0; i < npol; i++)
            {
                if (isDigitOne(i,l)) continue;

                cvars.push_back(vars.at(iLambda+i)); // Lambda_P
            }

            DenseMatrix A = DenseMatrix::Ones(1,cvars.size());
            A(0,0) = 1;
            DenseVector b = DenseMatrix::Ones(1,1);

            ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, false);
            cs->add(lincon);
        }
    }

    return cs;
}

template <unsigned int N>
ConstraintPtr constraintPiecewisePolynomialDLog4(Variables vars, SPLINTER::BSpline bspline)
{
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

    unsigned int dim = bspline.getNumVariables();
    assert(dim == N);
    assert(vars.size() == dim + 1);

    // Decompose B-spline
    PiecewisePolynomial<N> pwp = decomposeBSplinePower<N>(bspline);
    auto polynomials = pwp.polynomials;
    auto polytopes = pwp.polytopes;

    // Get monomial powers
    DenseMatrix M = getPowersMatrix(bspline.getBasisDegrees());

    // Compute monomial bounds xi^L and xi^U
    DenseVector xil = DenseVector::Ones(M.rows());
    DenseVector xiu = DenseVector::Ones(M.rows());

    for (int i = 0; i < M.rows(); ++i)
    {
        bool even = true;

        double ximax = 1;

        for (int j = 0; j < M.cols(); ++j)
        {
            int pow = M(i,j);

            if (pow < 1)
                continue;

            if (pow % 2 != 0)
                even = false;

            double xjl = vars.at(j)->getLowerBound();
            double xju = vars.at(j)->getUpperBound();

            ximax *= std::pow(std::max(xjl, xju), pow);
        }

        xiu(i) = ximax;
        xil(i) = -ximax;

        if (even) xil(i) = 0;
    }

    cout << xil << endl;
    cout << "--" << endl;
    cout << xiu << endl;

    // Get number of polytopes
    int npol = polytopes.size(); // Number of polytopes

    // Variables: [x z Lambda_P lambda_{P,v} xi phi_P y_P]
    int nz = 1;
    int iz = N;

    int nLambda = npol;
    int iLambda = N+1;

    int nlambda = std::pow(2,N)*npol; // 2^N vertices for hyperrectangles
    int ilambda = iLambda+nLambda;

    int nxi = M.rows();
    int ixi = ilambda+nlambda;

    int nphi = nxi*npol;
    int iphi = ixi+nxi;

    int ny = std::ceil(std::log2(npol)); // Logarithmic number of binary variables (y_P)
    int iy = iphi + nphi;

    int ncaux = nlambda + nLambda + nxi + nphi; // Number of aux. cont. variables
    int nbaux = ny; // Number of aux. binary variables

    // Create variables
    for (int i = 0; i < ncaux; ++i)
    {
        auto var = std::make_shared<Variable>(0);
        vars.push_back(var);
    }
    for (int i = 0; i < nbaux; ++i)
    {
        auto var = std::make_shared<Variable>(0, 0, 1, VariableType::BINARY);
        vars.push_back(var);
    }

    auto cs = std::make_shared<ConstraintSet>(vars);

    // Update polynomial variables
    {
        Variables pVars;
        for (unsigned int j = 0; j < dim; ++j)
            pVars.push_back(vars.at(j)); // x
        pVars.push_back(vars.at(iz)); // z

        for (auto pol : polynomials)
        {
            pol->setVariables(pVars);
        }
    }

    {
        // Add constraint sum_P sum_v lambda_{P,v} v = x
        DenseMatrix A = DenseMatrix::Ones(N, N+nlambda);
        DenseMatrix I = DenseMatrix::Identity(N,N);
        A.block(0,0,N,N) = -I;
        DenseVector b = DenseMatrix::Zero(N,1);

        std::vector<VariablePtr> cvars;
        for (int i = 0; i < N; i++)
            cvars.push_back(vars.at(i)); // x
        for (int i = 0; i < nlambda; i++)
            cvars.push_back(vars.at(ilambda+i)); // lambda_{P,v}

        int k = 0;
        for (const auto &p : polytopes)
        {
            for (const auto &v : p.vertices)
            {
                // Vertex
                auto x = v.x;
                for (int l = 0; l < N; l++)
                    A(l,N+k) = x.at(l);
                k++;
            }
        }

        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, true);
        cs->add(lincon);
    }

    {
        // Add constraint sum_P a_P*phi_P <= z
        DenseMatrix A = DenseMatrix::Zero(1, 1+nphi);
        A(0,0) = -1; // -z

        int cphi = 0;
        for (auto &pol : polynomials)
        {
            DenseMatrix ap = pol->getCoefficients();
            DenseMatrix a = ap.block(0,0,ap.rows()-1,ap.cols()); // Last coefficient related to z
            assert(a.rows() == nxi);
            assert(a.cols() == 1);
            a.transposeInPlace();

            A.block(0,1+cphi,1,nxi) = a;

            cphi += nxi;
        }
        assert(cphi == nphi);

        DenseVector b = DenseVector::Zero(1);

        Variables cvars;
        cvars.push_back(vars.at(iz)); // z

        for (int i = 0; i < nphi; ++i)
            cvars.push_back(vars.at(iphi+i)); // phi_p

        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, false);
        cs->add(lincon);
    }

    {
        // Add McCormick constraints A*[phi_P Lambda_P xi] <= b for all P
        DenseMatrix A = DenseMatrix::Zero(4*nxi, nxi+1+nxi);
        DenseVector b = DenseVector::Zero(4*nxi);
        DenseMatrix I = DenseMatrix::Identity(nxi, nxi);

        A.block(0, 0, nxi, nxi) = -I;
        A.block(nxi, 0, nxi, nxi) = -I;
        A.block(2*nxi, 0, nxi, nxi) = I;
        A.block(3*nxi, 0, nxi, nxi) = I;

        A.block(0, nxi, nxi, 1) = xiu;
        A.block(nxi, nxi, nxi, 1) = xil;
        A.block(2*nxi, nxi, nxi, 1) = -xil;
        A.block(3*nxi, nxi, nxi, 1) = -xiu;

        A.block(0, nxi+1, nxi, nxi) = I;
        //A.block(nxi, nxi+1, nxi, nxi) = 0;
        A.block(2*nxi, nxi+1, nxi, nxi) = -I;
        //A.block(3*nxi, nxi+1, nxi, nxi) = 0;

        b.block(0, 0, nxi, 1) = xiu;
        b.block(2*nxi, 0, nxi, 1) = -xil;

        for (int i = 0; i < npol; ++i)
        {
            Variables cvars;
            for (int j = 0; j < nxi; ++j)
                cvars.push_back(vars.at(iphi+i*nxi+j)); // phi

            cvars.push_back(vars.at(iLambda+i)); // Lambda_P

            for (int j = 0; j < nxi; ++j)
                cvars.push_back(vars.at(ixi+j)); // xi

            auto lincon = std::make_shared<ConstraintLinear>(cvars, A, b, false);

            cs->add(lincon);
        }
    }

    {
        // Add monomials
        for (int i = 0; i < nxi; ++i)
        {
            Variables cvars = {
                vars.at(ixi+i) // xi
            };

            std::vector<unsigned int> powers;

            for (int j = 0; j < N; ++j)
            {
                if (M(i,j) > 0)
                {
                    cvars.push_back(vars.at(j)); // x
                    powers.push_back(M(i,j));
                }
            }

            if (powers.size() > 0)
            {
                // Create monomial xi = m(x)
                DenseVector ci = DenseVector::Ones(2);
                ci(0) = -1;
                DenseMatrix Ei = DenseMatrix::Zero(2,1+powers.size());
                Ei(0,0) = 1;

                for (int j = 0; j < powers.size(); ++j)
                {
                    Ei(1,1+j) = powers[j];
                }

                auto mono = std::make_shared<ConstraintPolynomial>(cvars, ci, Ei, true);

                cs->add(mono);
            }
            else
            {
                // Fix variable xi to 1
                cvars.at(0)->setLowerBound(1);
                cvars.at(0)->setUpperBound(1);
            }
        }
    }

    {
        // lambda_{P,v} >= 0
        for (int i = 0; i < nlambda; ++i)
            vars.at(ilambda+i)->updateBounds(0,1);

        // L_P >= 0
        for (int i = 0; i < nLambda; ++i)
            vars.at(iLambda+i)->updateBounds(0,1);

        // xil <= xi <= xiu
        for (int i = 0; i < nxi; ++i)
            vars.at(ixi+i)->updateBounds(xil(i), xiu(i));
    }

    {
        // Add constraint sum_v lambda_{P,v} = L_P for all P
        int ip = 0;
        int iv = 0;
        for (const auto &p : polytopes)
        {
            int nvert = p.vertices.size();

            DenseMatrix A = DenseMatrix::Ones(1, 1+nvert);
            A(0,0) = -1;
            DenseVector b(1); b(0) = 0;

            std::vector<VariablePtr> cvars;
            cvars.push_back(vars.at(iLambda+ip)); // L_P
            for (int i = 0; i < nvert; i++)
                cvars.push_back(vars.at(ilambda+iv+i)); // lambda_{P,v}

            ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, true);
            cs->add(lincon);

            ip++;
            iv += nvert;
        }
        assert(iv == nlambda);
    }

    {
        // Add constraint sum_P L_{P} = 1
        DenseMatrix A = DenseMatrix::Ones(1, nLambda);
        DenseVector b(1); b(0) = 1;

        std::vector<VariablePtr> cvars;
        for (int i = 0; i < nLambda; ++i)
            cvars.push_back(vars.at(iLambda+i)); // L_P

        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, true);
        cs->add(lincon);
    }

    {
        // For all l in L(P) = {1,...,ceil(log2(|P|))}
        // add constraint sum_P+ sum_v lambda_{P,v} <= y_l
        for (int l = 0; l < ny; l++)
        {
            std::vector<VariablePtr> cvars;
            cvars.push_back(vars.at(iy+l)); // y_l

            for (unsigned int i = 0; i < npol; i++)
            {
                if (!isDigitOne(i,l)) continue;

                cvars.push_back(vars.at(iLambda+i)); // Lambda_P
            }

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
            cvars.push_back(vars.at(iy+l)); // y_l

            for (unsigned int i = 0; i < npol; i++)
            {
                if (isDigitOne(i,l)) continue;

                cvars.push_back(vars.at(iLambda+i)); // Lambda_P
            }

            DenseMatrix A = DenseMatrix::Ones(1,cvars.size());
            A(0,0) = 1;
            DenseVector b = DenseMatrix::Ones(1,1);

            ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, false);
            cs->add(lincon);
        }
    }

    return cs;
}

template <unsigned int N>
ConstraintPtr constraintPiecewisePolynomialDLog5(Variables vars, SPLINTER::BSpline bspline, bool equality)
{
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

    unsigned int dim = bspline.getNumVariables();
    assert(dim == N);
    assert(vars.size() == dim + 1);

    // Decompose B-spline
    PiecewisePolynomial<N> pwp = decomposeBSplinePower<N>(bspline);
    auto polynomials = pwp.polynomials;
    auto polytopes = pwp.polytopes;

    // Get monomial powers
    DenseMatrix M = getPowersMatrix(bspline.getBasisDegrees());

    // Compute monomial bounds xi^L and xi^U
    DenseVector xil = DenseVector::Ones(M.rows());
    DenseVector xiu = DenseVector::Ones(M.rows());

    for (int i = 0; i < M.rows(); ++i)
    {
        bool even = true;

        double ximax = 1;

        for (int j = 0; j < M.cols(); ++j)
        {
            int pow = M(i,j);

            if (pow < 1)
                continue;

            if (pow % 2 != 0)
                even = false;

            double xjl = vars.at(j)->getLowerBound();
            double xju = vars.at(j)->getUpperBound();

            ximax *= std::pow(std::max(xjl, xju), pow);
        }

        xiu(i) = ximax;
        xil(i) = -ximax;

        if (even) xil(i) = 0;
    }

    cout << xil << endl;
    cout << "--" << endl;
    cout << xiu << endl;

    // Get number of polytopes
    int npol = polytopes.size(); // Number of polytopes

    // Variables: [x z z_P Lambda_P lambda_{P,v} xi phi_P y_P]
    int nz = 1;
    int iz = N;

    int nzp = npol;
    int izp = iz+nz;

    int nLambda = npol;
    int iLambda = izp+nzp;

    int nlambda = std::pow(2,N)*npol; // 2^N vertices for hyperrectangles
    int ilambda = iLambda+nLambda;

    int nxi = M.rows();
    int ixi = ilambda+nlambda;

    int nphi = npol;
    int iphi = ixi+nxi;

    int ny = std::ceil(std::log2(npol)); // Logarithmic number of binary variables (y_P)
    int iy = iphi + nphi;

    int ncaux = nzp + nlambda + nLambda + nxi + nphi; // Number of aux. cont. variables
    int nbaux = ny; // Number of aux. binary variables

    // Create variables
    for (int i = 0; i < ncaux; ++i)
    {
        auto var = std::make_shared<Variable>(0);
        vars.push_back(var);
    }
    for (int i = 0; i < nbaux; ++i)
    {
        auto var = std::make_shared<Variable>(0, 0, 1, VariableType::BINARY);
        vars.push_back(var);
    }

    auto cs = std::make_shared<ConstraintSet>(vars);

    {
        // Update polynomial variables
        Variables pVars;
        for (unsigned int j = 0; j < dim; ++j)
            pVars.push_back(vars.at(j)); // x
        pVars.push_back(vars.at(iz)); // z

        for (auto pol : polynomials)
        {
            pol->setVariables(pVars);
        }
    }

    {
        // xil <= xi <= xiu
        for (int i = 0; i < nxi; ++i)
            vars.at(ixi+i)->updateBounds(xil(i), xiu(i));
    }

    {
        // Add constraint sum_P sum_v lambda_{P,v} v = x
        DenseMatrix A = DenseMatrix::Ones(N, N+nlambda);
        DenseMatrix I = DenseMatrix::Identity(N,N);
        A.block(0,0,N,N) = -I;
        DenseVector b = DenseMatrix::Zero(N,1);

        std::vector<VariablePtr> cvars;
        for (int i = 0; i < N; i++)
            cvars.push_back(vars.at(i)); // x
        for (int i = 0; i < nlambda; i++)
            cvars.push_back(vars.at(ilambda+i)); // lambda_{P,v}

        int k = 0;
        for (const auto &p : polytopes)
        {
            for (const auto &v : p.vertices)
            {
                // Vertex
                auto x = v.x;
                for (int l = 0; l < N; l++)
                    A(l,N+k) = x.at(l);
                k++;
            }
        }

        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, true);
        cs->add(lincon);
    }

    {
        // Add constraint sum_P phi_P <= z
        DenseMatrix A = DenseMatrix::Ones(1,1+nphi);
        A(0,0) = -1;
        DenseVector b(1); b(0) = 0;

        Variables cvars;
        cvars.push_back(vars.at(iz));
        for (int i = 0; i < nphi; ++i)
            cvars.push_back(vars.at(iphi+i));

        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, equality);
        cs->add(lincon);
    }

//    {
//        // Add constraint sum_P a_P*phi_P <= z
//        DenseMatrix A = DenseMatrix::Zero(1, 1+nphi);
//        A(0,0) = -1; // -z

//        int cphi = 0;
//        for (auto &pol : polynomials)
//        {
//            DenseMatrix ap = pol->getCoefficients();
//            DenseMatrix a = ap.block(0,0,ap.rows()-1,ap.cols()); // Last coefficient related to z
//            assert(a.rows() == nxi);
//            assert(a.cols() == 1);
//            a.transposeInPlace();

//            A.block(0,1+cphi,1,nxi) = a;

//            cphi += nxi;
//        }
//        assert(cphi == nphi);

//        DenseVector b = DenseVector::Zero(1);

//        Variables cvars;
//        cvars.push_back(vars.at(iz)); // z

//        for (int i = 0; i < nphi; ++i)
//            cvars.push_back(vars.at(iphi+i)); // phi_p

//        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, false);
//        cs->add(lincon);
//    }

    {
        // Add constraint z_P = a_P*xi for all P
        int ipol = 0;
        for (auto &pol : polynomials)
        {
            DenseMatrix a = pol->getCoefficients();
            assert(a.rows() == nxi+1); // Last coefficient related to z
            assert(a.cols() == 1);
            a.transposeInPlace();

            DenseVector b(1); b(0) = 0;

            Variables cvars;

            for (int i = 0; i < nxi; ++i)
                cvars.push_back(vars.at(ixi+i));
            cvars.push_back(vars.at(izp+ipol));

            ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, a, b, true);
            //cout << "Before: " << vars.at(izp+ipol)->getLowerBound() << " / " << vars.at(izp+ipol)->getUpperBound() << endl;
            assert(lincon->reduceVariableRanges()); // Set bounds on zp
            //cout << "After: " << vars.at(izp+ipol)->getLowerBound() << " / " << vars.at(izp+ipol)->getUpperBound() << endl;
            cs->add(lincon);

            ++ipol;
        }
    }

    {

        // Add McCormick constraints A*[phi_P Lambda_P z_P] <= b for all P

        for (int i = 0; i < npol; ++i)
        {
            DenseMatrix A = DenseMatrix::Zero(4, 3);
            DenseVector b = DenseVector::Zero(4);
            DenseMatrix I = DenseMatrix::Identity(nxi, nxi);

            double zpl = vars.at(izp+i)->getLowerBound();
            double zpu = vars.at(izp+i)->getUpperBound();

            A(0, 0) = -1;
            A(1, 0) = -1;
            A(2, 0) = 1;
            A(3, 0) = 1;

            A(0, 1) = zpu;
            A(1, 1) = zpl;
            A(2, 1) = -zpl;
            A(3, 1) = -zpu;

            A(0, 2) = 1;
            A(1, 2) = 0;
            A(2, 2) = -1;
            A(3, 2) = 0;

            b(0) = zpu;
            b(2) = -zpl;

            Variables cvars = {
                vars.at(iphi+i), // phi_P
                vars.at(iLambda+i), // Lambda_P
                vars.at(izp+i) // z_P
            };

            auto lincon = std::make_shared<ConstraintLinear>(cvars, A, b, false);

            cs->add(lincon);
        }
    }

    {
        // Add monomials
        for (int i = 0; i < nxi; ++i)
        {
            Variables cvars = {
                vars.at(ixi+i) // xi
            };

            std::vector<unsigned int> powers;

            for (int j = 0; j < N; ++j)
            {
                if (M(i,j) > 0)
                {
                    cvars.push_back(vars.at(j)); // x
                    powers.push_back(M(i,j));
                }
            }

            if (powers.size() > 0)
            {
                // Create monomial xi = m(x)
                DenseVector ci = DenseVector::Ones(2);
                ci(0) = -1;
                DenseMatrix Ei = DenseMatrix::Zero(2,1+powers.size());
                Ei(0,0) = 1;

                for (int j = 0; j < powers.size(); ++j)
                {
                    Ei(1,1+j) = powers[j];
                }

                auto mono = std::make_shared<ConstraintPolynomial>(cvars, ci, Ei, true);

                cs->add(mono);
            }
            else
            {
                // Fix variable xi to 1
                cvars.at(0)->setLowerBound(1);
                cvars.at(0)->setUpperBound(1);
            }
        }
    }

    {
        // lambda_{P,v} >= 0
        for (int i = 0; i < nlambda; ++i)
            vars.at(ilambda+i)->updateBounds(0,1);

        // L_P >= 0
        for (int i = 0; i < nLambda; ++i)
            vars.at(iLambda+i)->updateBounds(0,1);

        // xil <= xi <= xiu
        for (int i = 0; i < nxi; ++i)
            vars.at(ixi+i)->updateBounds(xil(i), xiu(i));
    }

    {
        // Add constraint sum_v lambda_{P,v} = L_P for all P
        int ip = 0;
        int iv = 0;
        for (const auto &p : polytopes)
        {
            int nvert = p.vertices.size();

            DenseMatrix A = DenseMatrix::Ones(1, 1+nvert);
            A(0,0) = -1;
            DenseVector b(1); b(0) = 0;

            std::vector<VariablePtr> cvars;
            cvars.push_back(vars.at(iLambda+ip)); // L_P
            for (int i = 0; i < nvert; i++)
                cvars.push_back(vars.at(ilambda+iv+i)); // lambda_{P,v}

            ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, true);
            cs->add(lincon);

            ip++;
            iv += nvert;
        }
        assert(iv == nlambda);
    }

    {
        // Add constraint sum_P L_{P} = 1
        DenseMatrix A = DenseMatrix::Ones(1, nLambda);
        DenseVector b(1); b(0) = 1;

        std::vector<VariablePtr> cvars;
        for (int i = 0; i < nLambda; ++i)
            cvars.push_back(vars.at(iLambda+i)); // L_P

        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, true);
        cs->add(lincon);
    }

    {
        // For all l in L(P) = {1,...,ceil(log2(|P|))}
        // add constraint sum_P+ sum_v lambda_{P,v} <= y_l
        for (int l = 0; l < ny; l++)
        {
            std::vector<VariablePtr> cvars;
            cvars.push_back(vars.at(iy+l)); // y_l

            for (unsigned int i = 0; i < npol; i++)
            {
                if (!isDigitOne(i,l)) continue;

                cvars.push_back(vars.at(iLambda+i)); // Lambda_P
            }

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
            cvars.push_back(vars.at(iy+l)); // y_l

            for (unsigned int i = 0; i < npol; i++)
            {
                if (isDigitOne(i,l)) continue;

                cvars.push_back(vars.at(iLambda+i)); // Lambda_P
            }

            DenseMatrix A = DenseMatrix::Ones(1,cvars.size());
            A(0,0) = 1;
            DenseVector b = DenseMatrix::Ones(1,1);

            ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, false);
            cs->add(lincon);
        }
    }

    return cs;
}

template <unsigned int N>
ConstraintPtr constraintPiecewisePolynomialDLog6(Variables vars, SPLINTER::BSpline bspline, bool equality)
{
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

    unsigned int dim = bspline.getNumVariables();
    assert(dim == N);
    assert(vars.size() == dim + 1);

    // Decompose B-spline
    PiecewisePolynomial<N> pwp = decomposeBSplinePower<N>(bspline);
    auto polynomials = pwp.polynomials;
    auto polytopes = pwp.polytopes;

    // Get monomial powers
    DenseMatrix M = getPowersMatrix(bspline.getBasisDegrees());

    // Get number of polytopes
    int npol = polytopes.size(); // Number of polytopes

    // Variables: [x z t z_P Lambda_P lambda_{P,v} xi phi_P y_P]
    int nz = 1;
    int iz = N;

    int nt = N; // Scaled x
    int it = iz+nz;

    int nzp = npol;
    int izp = it+nt;

    int nLambda = npol;
    int iLambda = izp+nzp;

    int nlambda = std::pow(2,N)*npol; // 2^N vertices for hyperrectangles
    int ilambda = iLambda+nLambda;

    int nxi = M.rows();
    int ixi = ilambda+nlambda;

    int nphi = npol;
    int iphi = ixi+nxi;

    int ny = std::ceil(std::log2(npol)); // Logarithmic number of binary variables (y_P)
    int iy = iphi + nphi;

    int ncaux = nt + nzp + nlambda + nLambda + nxi + nphi; // Number of aux. cont. variables
    int nbaux = ny; // Number of aux. binary variables

    // Create variables
    for (int i = 0; i < ncaux; ++i)
    {
        auto var = std::make_shared<Variable>(0);
        vars.push_back(var);
    }
    for (int i = 0; i < nbaux; ++i)
    {
        auto var = std::make_shared<Variable>(0, 0, 1, VariableType::BINARY);
        vars.push_back(var);
    }

    auto cs = std::make_shared<ConstraintSet>(vars);

    {
        // Update polynomial variables
        Variables pVars;
        for (unsigned int j = 0; j < dim; ++j)
            pVars.push_back(vars.at(j)); // x
        pVars.push_back(vars.at(iz)); // z

        for (auto pol : polynomials)
        {
            pol->setVariables(pVars);
        }
    }

    {
        // xil <= xi <= xiu
        for (int i = 0; i < nxi; ++i)
            vars.at(ixi+i)->updateBounds(0,1);

        // lambda_{P,v} >= 0
        for (int i = 0; i < nlambda; ++i)
            vars.at(ilambda+i)->updateBounds(0,1);

        // L_P >= 0
        for (int i = 0; i < nLambda; ++i)
            vars.at(iLambda+i)->updateBounds(0,1);

        // 0 <= t_i <= 1
        for (int i = 0; i < nt; ++i)
            vars.at(it+i)->updateBounds(0,1);
    }

    {
        // x_i = (x_i^U-x_i^L)*t_i + x_i^L
        for (int i = 0; i < N; ++i)
        {
            auto xi = vars.at(i);

            Variables cvars = {
                vars.at(i), // x_i
                vars.at(it+i) // t_i
            };

            DenseMatrix A = DenseMatrix::Zero(1,2);
            A(0,0) = 1; A(0,1) = -(xi->getUpperBound() - xi->getLowerBound());

            DenseVector b(1); b(0) = xi->getLowerBound();

            auto con = std::make_shared<ConstraintLinear>(cvars, A, b, true);

            cs->add(con);
        }
    }

    {
        // Add constraint sum_P sum_v lambda_{P,v} v = x
        DenseMatrix A = DenseMatrix::Ones(N, N+nlambda);
        DenseMatrix I = DenseMatrix::Identity(N,N);
        A.block(0,0,N,N) = -I;
        DenseVector b = DenseMatrix::Zero(N,1);

        std::vector<VariablePtr> cvars;
        for (int i = 0; i < N; i++)
            cvars.push_back(vars.at(i)); // x
        for (int i = 0; i < nlambda; i++)
            cvars.push_back(vars.at(ilambda+i)); // lambda_{P,v}

        int k = 0;
        for (const auto &p : polytopes)
        {
            for (const auto &v : p.vertices)
            {
                // Vertex
                auto x = v.x;
                for (int l = 0; l < N; l++)
                    A(l,N+k) = x.at(l);
                k++;
            }
        }

        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, true);
        cs->add(lincon);
    }

    {
        // Add constraint sum_P phi_P <= z
        DenseMatrix A = DenseMatrix::Ones(1,1+nphi);
        A(0,0) = -1;
        DenseVector b(1); b(0) = 0;

        Variables cvars;
        cvars.push_back(vars.at(iz));
        for (int i = 0; i < nphi; ++i)
            cvars.push_back(vars.at(iphi+i));

        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, equality);
        cs->add(lincon);
    }

    {
        // Add constraint z_P = a_P*xi for all P
        int ipol = 0;
        for (auto &pol : polynomials)
        {
            DenseMatrix a = pol->getCoefficients();
            assert(a.rows() == nxi+1); // Last coefficient related to z
            assert(a.cols() == 1);
            a.transposeInPlace();

            DenseVector b(1); b(0) = 0;

            Variables cvars;

            for (int i = 0; i < nxi; ++i)
                cvars.push_back(vars.at(ixi+i));
            cvars.push_back(vars.at(izp+ipol));

            ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, a, b, true);

            assert(lincon->reduceVariableRanges()); // Set bounds on zp

            cs->add(lincon);

            ++ipol;
        }
    }

    {

        // Add McCormick constraints A*[phi_P Lambda_P z_P] <= b for all P

        for (int i = 0; i < npol; ++i)
        {
            DenseMatrix A = DenseMatrix::Zero(4, 3);
            DenseVector b = DenseVector::Zero(4);

            double zpl = vars.at(izp+i)->getLowerBound();
            double zpu = vars.at(izp+i)->getUpperBound();

            A(0, 0) = -1;
            A(1, 0) = -1;
            A(2, 0) = 1;
            A(3, 0) = 1;

            A(0, 1) = zpu;
            A(1, 1) = zpl;
            A(2, 1) = -zpl;
            A(3, 1) = -zpu;

            A(0, 2) = 1;
            A(1, 2) = 0;
            A(2, 2) = -1;
            A(3, 2) = 0;

            b(0) = zpu;
            b(2) = -zpl;

            Variables cvars = {
                vars.at(iphi+i), // phi_P
                vars.at(iLambda+i), // Lambda_P
                vars.at(izp+i) // z_P
            };

            auto lincon = std::make_shared<ConstraintLinear>(cvars, A, b, false);

            cs->add(lincon);
        }
    }

    {
        // Add monomials
        for (int i = 0; i < nxi; ++i)
        {
            Variables cvars = {
                vars.at(ixi+i) // xi
            };

            std::vector<unsigned int> powers;

            for (int j = 0; j < N; ++j)
            {
                if (M(i,j) > 0)
                {
                    cvars.push_back(vars.at(it+j)); // t_i
                    powers.push_back(M(i,j));
                }
            }

            if (powers.size() > 0)
            {
                // Create monomial xi = m(x)
                DenseVector ci = DenseVector::Ones(2);
                ci(0) = -1;
                DenseMatrix Ei = DenseMatrix::Zero(2,1+powers.size());
                Ei(0,0) = 1;

                for (int j = 0; j < powers.size(); ++j)
                {
                    Ei(1,1+j) = powers[j];
                }

                auto mono = std::make_shared<ConstraintPolynomial>(cvars, ci, Ei, true);

                cs->add(mono);
            }
            else
            {
                // Fix variable xi to 1
                cvars.at(0)->setLowerBound(1);
                cvars.at(0)->setUpperBound(1);
            }
        }
    }

    {
        // Add constraint sum_v lambda_{P,v} = L_P for all P
        int ip = 0;
        int iv = 0;
        for (const auto &p : polytopes)
        {
            int nvert = p.vertices.size();

            DenseMatrix A = DenseMatrix::Ones(1, 1+nvert);
            A(0,0) = -1;
            DenseVector b(1); b(0) = 0;

            std::vector<VariablePtr> cvars;
            cvars.push_back(vars.at(iLambda+ip)); // L_P
            for (int i = 0; i < nvert; i++)
                cvars.push_back(vars.at(ilambda+iv+i)); // lambda_{P,v}

            ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, true);
            cs->add(lincon);

            ip++;
            iv += nvert;
        }
        assert(iv == nlambda);
    }

    {
        // Add constraint sum_P L_{P} = 1
        DenseMatrix A = DenseMatrix::Ones(1, nLambda);
        DenseVector b(1); b(0) = 1;

        std::vector<VariablePtr> cvars;
        for (int i = 0; i < nLambda; ++i)
            cvars.push_back(vars.at(iLambda+i)); // L_P

        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, true);
        cs->add(lincon);
    }

    {
        // For all l in L(P) = {1,...,ceil(log2(|P|))}
        // add constraint sum_P+ sum_v lambda_{P,v} <= y_l
        for (int l = 0; l < ny; l++)
        {
            std::vector<VariablePtr> cvars;
            cvars.push_back(vars.at(iy+l)); // y_l

            for (unsigned int i = 0; i < npol; i++)
            {
                if (!isDigitOne(i,l)) continue;

                cvars.push_back(vars.at(iLambda+i)); // Lambda_P
            }

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
            cvars.push_back(vars.at(iy+l)); // y_l

            for (unsigned int i = 0; i < npol; i++)
            {
                if (isDigitOne(i,l)) continue;

                cvars.push_back(vars.at(iLambda+i)); // Lambda_P
            }

            DenseMatrix A = DenseMatrix::Ones(1,cvars.size());
            A(0,0) = 1;
            DenseVector b = DenseMatrix::Ones(1,1);

            ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, false);
            cs->add(lincon);
        }
    }

    return cs;
}

template <unsigned int N>
ConstraintPtr constraintPiecewisePolynomialDLog7(Variables vars, SPLINTER::BSpline bspline, bool equality)
{
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

    unsigned int dim = bspline.getNumVariables();
    assert(dim == N);
    assert(vars.size() == dim + 1);

    // Decompose B-spline
    PiecewisePolynomial<N> pwp = decomposeBSplinePower<N>(bspline);
    auto polynomials = pwp.polynomials;
    auto polytopes = pwp.polytopes;

    // Get monomial powers
    DenseMatrix M = getPowersMatrix(bspline.getBasisDegrees());

    // Get number of polytopes
    int npol = polytopes.size(); // Number of polytopes

    // Variables: [x z t z_P Lambda_P xi phi_P y_P]
    int nz = 1;
    int iz = N;

    int nt = N; // Scaled x
    int it = iz+nz;

    int nzp = npol;
    int izp = it+nt;

    int nLambda = npol;
    int iLambda = izp+nzp;

    int nxi = M.rows();
    int ixi = iLambda+nLambda;

    int nphi = npol;
    int iphi = ixi+nxi;

    int ny = std::ceil(std::log2(npol)); // Logarithmic number of binary variables (y_P)
    int iy = iphi + nphi;

    int ncaux = nt + nzp + nLambda + nxi + nphi; // Number of aux. cont. variables
    int nbaux = ny; // Number of aux. binary variables

    // Create variables
    for (int i = 0; i < ncaux; ++i)
    {
        auto var = std::make_shared<Variable>(0);
        vars.push_back(var);
    }
    for (int i = 0; i < nbaux; ++i)
    {
        auto var = std::make_shared<Variable>(0, 0, 1, VariableType::BINARY);
        vars.push_back(var);
    }

    auto cs = std::make_shared<ConstraintSet>(vars);

    {
        // Update polynomial variables
        Variables pVars;
        for (unsigned int j = 0; j < dim; ++j)
            pVars.push_back(vars.at(j)); // x
        pVars.push_back(vars.at(iz)); // z

        for (auto pol : polynomials)
        {
            pol->setVariables(pVars);
        }
    }

    {
        // xil <= xi <= xiu
        for (int i = 0; i < nxi; ++i)
            vars.at(ixi+i)->updateBounds(0,1);

        // L_P >= 0
        for (int i = 0; i < nLambda; ++i)
            vars.at(iLambda+i)->updateBounds(0,1);

        // 0 <= t_i <= 1
        for (int i = 0; i < nt; ++i)
            vars.at(it+i)->updateBounds(0,1);
    }

    {
        // x_i = (x_i^U-x_i^L)*t_i + x_i^L
        for (int i = 0; i < N; ++i)
        {
            auto xi = vars.at(i);

            Variables cvars = {
                vars.at(i), // x_i
                vars.at(it+i) // t_i
            };

            DenseMatrix A = DenseMatrix::Zero(1,2);
            A(0,0) = 1; A(0,1) = -(xi->getUpperBound() - xi->getLowerBound());

            DenseVector b(1); b(0) = xi->getLowerBound();

            auto con = std::make_shared<ConstraintLinear>(cvars, A, b, true);

            cs->add(con);
        }
    }

    {
        /*
         * x <= Lambda_P x_P^U - (1 - Lambda_P) x^U
         * x >= Lambda_P x_P^L + (1 - Lambda_P) x^L
         */
        int ipoly = 0;
        for (auto poly : polytopes)
        {
            auto verts = poly.vertices;

            for (int i = 0; i < N; ++i)
            {
                auto xi = vars.at(i);

                // Find upper (x_P^U) and lower (x_P^L) bound on rectangle
                double xpi_lb = xi->getUpperBound(); // Start at upper bound
                double xpi_ub = xi->getLowerBound(); // Start at lower bound
                for (auto v : verts)
                {
                    double vxi = v.x.at(i);

                    if (vxi > xpi_ub) xpi_ub = vxi;
                    if (vxi < xpi_lb) xpi_lb = vxi;
                }
                assert(xpi_lb < xpi_ub);

                Variables cvars = {
                    vars.at(i), // x
                    vars.at(iLambda+ipoly) // Lambda_P
                };

                DenseMatrix A_lb = DenseMatrix::Ones(1,2);
                A_lb(0,0) = -1;
                A_lb(0,1) = xpi_lb - xi->getLowerBound();
                DenseVector b_lb(1); b_lb(0) = -xi->getLowerBound();
                auto lincon_lb = std::make_shared<ConstraintLinear>(cvars, A_lb, b_lb, false);
                cs->add(lincon_lb);

                DenseMatrix A_ub = DenseMatrix::Ones(1,2);
                A_ub(0,1) = xi->getUpperBound() - xpi_ub;
                DenseVector b_ub(1); b_ub(0) = xi->getUpperBound();
                auto lincon_ub = std::make_shared<ConstraintLinear>(cvars, A_ub, b_ub, false);
                cs->add(lincon_ub);
            }

            ++ipoly;
        }
    }

    {
        // Add constraint sum_P phi_P <= z
        DenseMatrix A = DenseMatrix::Ones(1,1+nphi);
        A(0,0) = -1;
        DenseVector b(1); b(0) = 0;

        Variables cvars;
        cvars.push_back(vars.at(iz));
        for (int i = 0; i < nphi; ++i)
            cvars.push_back(vars.at(iphi+i));

        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, equality);
        cs->add(lincon);
    }

    {
        // Add constraint z_P = a_P*xi for all P
        int ipol = 0;
        for (auto &pol : polynomials)
        {
            DenseMatrix a = pol->getCoefficients();
            assert(a.rows() == nxi+1); // Last coefficient related to z
            assert(a.cols() == 1);
            a.transposeInPlace();

            DenseVector b(1); b(0) = 0;

            Variables cvars;

            for (int i = 0; i < nxi; ++i)
                cvars.push_back(vars.at(ixi+i));
            cvars.push_back(vars.at(izp+ipol));

            ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, a, b, true);

            assert(lincon->reduceVariableRanges()); // Set bounds on zp

            cs->add(lincon);

            ++ipol;
        }
    }

    {
        // Add McCormick constraints A*[phi_P Lambda_P z_P] <= b for all P
        for (int i = 0; i < npol; ++i)
        {
            DenseMatrix A = DenseMatrix::Zero(4, 3);
            DenseVector b = DenseVector::Zero(4);

            double zpl = vars.at(izp+i)->getLowerBound();
            double zpu = vars.at(izp+i)->getUpperBound();

            A(0, 0) = -1;
            A(1, 0) = -1;
            A(2, 0) = 1;
            A(3, 0) = 1;

            A(0, 1) = zpu;
            A(1, 1) = zpl;
            A(2, 1) = -zpl;
            A(3, 1) = -zpu;

            A(0, 2) = 1;
            A(1, 2) = 0;
            A(2, 2) = -1;
            A(3, 2) = 0;

            b(0) = zpu;
            b(2) = -zpl;

            Variables cvars = {
                vars.at(iphi+i), // phi_P
                vars.at(iLambda+i), // Lambda_P
                vars.at(izp+i) // z_P
            };

            auto lincon = std::make_shared<ConstraintLinear>(cvars, A, b, false);

            cs->add(lincon);
        }
    }

    {
        // Add monomials
        for (int i = 0; i < nxi; ++i)
        {
            Variables cvars = {
                vars.at(ixi+i) // xi
            };

            std::vector<unsigned int> powers;

            for (int j = 0; j < N; ++j)
            {
                if (M(i,j) > 0)
                {
                    cvars.push_back(vars.at(it+j)); // t_i
                    powers.push_back(M(i,j));
                }
            }

            if (powers.size() > 0)
            {
                // Create monomial xi = m(x)
                DenseVector ci = DenseVector::Ones(2);
                ci(0) = -1;
                DenseMatrix Ei = DenseMatrix::Zero(2,1+powers.size());
                Ei(0,0) = 1;

                for (int j = 0; j < powers.size(); ++j)
                {
                    Ei(1,1+j) = powers[j];
                }

                auto mono = std::make_shared<ConstraintPolynomial>(cvars, ci, Ei, true);

                cs->add(mono);
            }
            else
            {
                // Fix variable xi to 1
                cvars.at(0)->setLowerBound(1);
                cvars.at(0)->setUpperBound(1);
            }
        }
    }

    {
        // Add constraint sum_P Lambda_{P} = 1
        DenseMatrix A = DenseMatrix::Ones(1, nLambda);
        DenseVector b(1); b(0) = 1;

        std::vector<VariablePtr> cvars;
        for (int i = 0; i < nLambda; ++i)
            cvars.push_back(vars.at(iLambda+i)); // L_P

        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, true);
        cs->add(lincon);
    }

    {
        // For all l in L(P) = {1,...,ceil(log2(|P|))}
        // add constraint sum_P1 Lambda_P <= y_l
        for (int l = 0; l < ny; l++)
        {
            std::vector<VariablePtr> cvars;
            cvars.push_back(vars.at(iy+l)); // y_l

            for (unsigned int i = 0; i < npol; i++)
            {
                if (!isDigitOne(i,l)) continue;

                cvars.push_back(vars.at(iLambda+i)); // Lambda_P
            }

            DenseMatrix A = DenseMatrix::Ones(1,cvars.size());
            A(0,0) = -1;
            DenseVector b = DenseMatrix::Zero(1,1);

            ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, false);
            cs->add(lincon);
        }
    }

    {
        // For all l in L(P) = {1,...,ceil(log2(|P|))}
        // add constraint sum_P0 Lambda_P <= (1 - y_l)
        for (int l = 0; l < ny; l++)
        {
            std::vector<VariablePtr> cvars;
            cvars.push_back(vars.at(iy+l)); // y_l

            for (unsigned int i = 0; i < npol; i++)
            {
                if (isDigitOne(i,l)) continue;

                cvars.push_back(vars.at(iLambda+i)); // Lambda_P
            }

            DenseMatrix A = DenseMatrix::Ones(1,cvars.size());
            A(0,0) = 1;
            DenseVector b = DenseMatrix::Ones(1,1);

            ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, false);
            cs->add(lincon);
        }
    }

    return cs;
}

void testSplineDecomposition()
{
    // 2*x1^4 - 8*x1^3 + 8*x1^2 + 2
    std::vector<double> xs1 = linspace(0, 3, 10);
    std::vector<double> xs2 = linspace(0, 3, 10);

    auto func = [](double x, double y){ return -8*x*x*x + 8*x*x + 2 - y*y + 10*y + 5*x*y*y; };

    SPLINTER::DataTable data;

    for (auto x1 : xs1)
    {
        for (auto x2 : xs2)
        {
            auto y = func(x1, x2);
            std::vector<double> x = {x1, x2};
            data.addSample(x, y);
        }
    }

    SPLINTER::BSpline bs(data, SPLINTER::BSplineType::CUBIC);

    Variables vars = {
        std::make_shared<Variable>(0),
        std::make_shared<Variable>(0),
        std::make_shared<Variable>(0)
    };

    auto con = constraintPiecewisePolynomialDCC<2>(vars, bs);

}

unsigned int polyTest(int numSamples)
{
    assert(numSamples >= 2);
    std::vector<double> xs = linspace(0, 3, numSamples);

    auto func1 = [](double x) { return -2 -8*x*x +8*x*x*x -2*x*x*x*x; };
    auto func2 = [](double x) { return -36 +96*x -88*x*x +32*x*x*x -4*x*x*x*x; };

    SPLINTER::DataTable samples1, samples2;
    for (const auto &xv : xs)
    {
        double y1 = func1(xv);
        double y2 = func2(xv);

        samples1.addSample(xv, y1);
        samples2.addSample(xv, y2);
    }

    SPLINTER::BSpline bs1(samples1, SPLINTER::BSplineType::CUBIC);
    SPLINTER::BSpline bs2(samples2, SPLINTER::BSplineType::CUBIC);

    // Create some variables
    Variables vars = {
        std::make_shared<Variable>(-1, 0, 3),
        std::make_shared<Variable>(1, -4, 0)
    };

//    ConstraintPtr con1 = std::make_shared<ConstraintBSpline>(vars, bs1, true);
//    ConstraintPtr con2 = std::make_shared<ConstraintBSpline>(vars, bs2, true);

    auto con1 = constraintPiecewisePolynomialDCC<1>(vars, bs1);
    auto con2 = constraintPiecewisePolynomialDCC<1>(vars, bs2);

    ConstraintSetPtr cs = std::make_shared<ConstraintSet>();

    cs->add(con1);
    cs->add(con2);

    cs->writeToGAMS("test.gms");

    return 0;

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

    return timer.getMilliSeconds();
}

void michaTest()
{
    double pi = atan(1)*4;

    auto micha = [pi](double x1, double x2)
    {
        return -std::sin(x1)*std::pow(std::sin(x1*x1/pi), 20) - std::sin(x2)*std::pow(std::sin(2*x2*x2/pi), 20);
    };

    std::vector<double> costs = {0, 0, 1};
    std::vector<double> lb {0, 0, -INF};
    std::vector<double> ub {pi, pi, INF};

    Variables vars;

    for (int i = 0; i < 3; i++)
    {
        auto var = std::make_shared<Variable>(costs.at(i), lb.at(i), ub.at(i));
        vars.push_back(var);
    }

    SPLINTER::DataTable data;

    auto x1v = linspace(lb.at(0), ub.at(0), 30);
    auto x2v = linspace(lb.at(1), ub.at(1), 30);
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
    SPLINTER::BSpline bs(data, SPLINTER::BSplineType::LINEAR);

    //ConstraintPtr con = std::make_shared<ConstraintBSpline>(vars, bs, true);
    auto con = constraintPiecewisePolynomialDLog7<2>(vars, bs);

    ConstraintSetPtr cs = std::make_shared<ConstraintSet>();

    cs->add(con);

    cs->writeToGAMS("test.gms");

    exit(1);

    BB::BranchAndBound solver(cs);

    SolverResult result = solver.optimize();

    assert(result.status == SolverStatus::OPTIMAL);

    cout << result << endl;
    auto sol = result.primalVariables;
    cout << "Sol: x1 = " << sol.at(0) << ", x2 = " << sol.at(1) << endl;

}

void pumpSynthesis2(unsigned int grid)
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

            DenseVector ci = DenseVector::Ones(2);
            ci(1) = -1;
            DenseMatrix Ei = DenseMatrix::Zero(2,4);
            Ei(0,0) = 1; Ei(0,1) = 1; Ei(0,2) = 1;
            Ei(1,3) = 1;

            ConstraintPtr con = std::make_shared<ConstraintPolynomial>(cvars, ci, Ei, true);

            cs->add(con);
            cs_approx->add(con);
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

            DenseVector ci = DenseVector::Ones(2);
            ci(1) = -1;
            DenseMatrix Ei = DenseMatrix::Zero(2,3);
            Ei(0,0) = 1; Ei(0,1) = 1;
            Ei(1,2) = 1;

            ConstraintPtr con = std::make_shared<ConstraintPolynomial>(cvars, ci, Ei, true);

            cs->add(con);
            cs_approx->add(con);
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

            SPLINTER::DataTable samples;

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
            SPLINTER::BSpline bs(samples, SPLINTER::BSplineType::LINEAR);

            ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

            auto cbscon = constraintPiecewisePolynomialDLog7<2>(cvars, bs, true);

            cs_approx->add(cbscon);
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

            SPLINTER::DataTable samples;

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
            SPLINTER::BSpline bs(samples, SPLINTER::BSplineType::LINEAR);

            ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

            auto cbscon = constraintPiecewisePolynomialDLog7<2>(cvars, bs, true);

            cs_approx->add(cbscon);
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

            DenseVector ci = DenseVector::Ones(2);
            ci(1) = -1;
            DenseMatrix Ei = DenseMatrix::Zero(2,3);
            Ei(0,0) = 1; Ei(0,1) = 1;
            Ei(1,2) = 1;

            ConstraintPtr con = std::make_shared<ConstraintPolynomial>(cvars, ci, Ei, true);

            cs->add(con);
            cs_approx->add(con);
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

            DenseVector ci = DenseVector::Ones(2);
            ci(1) = -1;
            DenseMatrix Ei = DenseMatrix::Zero(2,3);
            Ei(0,0) = 1; Ei(0,1) = 1;
            Ei(1,2) = 1;

            ConstraintPtr con = std::make_shared<ConstraintPolynomial>(cvars, ci, Ei, true);

            cs->add(con);
            cs_approx->add(con);
        }
    }

    cs_approx->writeToGAMS("test.gms");
}
