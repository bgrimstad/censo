#include "unittests.h"
#include "Utils/definitions.h"
#include "OptimizationProblem/constraintset.h"
#include "OptimizationProblem/constraintbspline.h"
#include "OptimizationProblem/constraintlinear.h"

using std::cout;
using std::endl;
using namespace CENSO;

/*
 * Create a B-spline and reduce bounds
 */
void testConstraintBSplineRangeReduction1()
{
    DenseMatrix coefficients(1,6);
    coefficients << 1, 1, 0, 0, -1, -1;
    std::vector< std::vector<double> > knots = {{0,0,1,2,3,4,5,5}}; // Should have n+p+1 knots!
    std::vector<unsigned int> degrees = {1};
    Splinter::BSpline bspline(coefficients, knots, degrees);

    std::vector<VariablePtr> vars = {std::make_shared<Variable>(1, 0, 5),
                                     std::make_shared<Variable>(1, 0, 0.99)};

    ConstraintBSpline cbs(vars, bspline, true);

    cbs.reduceVariableRanges();

    for (auto var : vars)
        cout << *var << endl;
}

void testConstraintBSplineRangeReduction2()
{
    DenseMatrix coefficients(1,7);
    coefficients << 1, 1, 1, 0, -1, -1, -1;
    std::vector< std::vector<double> > knots = {{0,0,0,1,2,3,4,5,5,5}}; // Should have n+p+1 knots!
    std::vector<unsigned int> degrees = {2};
    Splinter::BSpline bspline(coefficients, knots, degrees);

    std::vector<VariablePtr> vars = {std::make_shared<Variable>(1, 0, 5),
                                     std::make_shared<Variable>(1, 0, 0.99)};

    ConstraintBSpline cbs(vars, bspline, true);

    cbs.reduceVariableRanges();

    for (auto var : vars)
        cout << *var << endl;
}

void testConstraintBSplineRangeReduction3()
{

    // 2*x1^4 - 8*x1^3 + 8*x1^2 + 2
    std::vector<double> xs = linspace(0, 3, 100);

    auto func = [](double x){ return 2*x*x*x*x - 8*x*x*x + 8*x*x + 2; };

    Splinter::DataTable data;

    for (auto x : xs)
    {
        auto y = func(x);
        data.addSample(x,y);
    }

    Splinter::BSpline bs(data, Splinter::BSplineType::CUBIC_FREE);

    std::vector<VariablePtr> vars = {std::make_shared<Variable>(1, -10, 10),
                                     std::make_shared<Variable>(1, 0, 5)};

    ConstraintBSpline cbs(vars, bs, true);

    cbs.reduceVariableRanges();

    // Should give x in [0, 2.515] which is the best possible reduction
    for (auto var : vars)
        cout << *var << endl;

}

//void test4()
//{
//    // 2*x^4 - 8*x^3 + 8*x^2 + 2 - 10*y^2
//    std::vector<double> xs = linspace(0, 3, 100);
//    std::vector<double> ys = linspace(-1, 1, 20);

//    auto func = [](double x, double y){ return 2*x*x*x*x - 8*x*x*x + 8*x*x + 2 - 10*y*y; };

//    Splinter::DataTable data;

//    for (auto x : xs)
//    {
//        for (auto y : ys)
//        {
//            std::vector<double> p = {x,y};
//            auto z = func(x,y);
//            data.addSample(p,z);
//        }

//    }

//    Splinter::BSpline bs(data, Splinter::BSplineType::CUBIC_FREE);

//    /*
//     * Should tighten bounds to
//     * x in [0, 2.909]
//     * y in [-1, 1]
//     */
//    CENSO::ConstraintBSpline cbs(bs, true);

//    std::vector<double> lb, ub;
//    cbs.getDomainBounds(lb, ub);
//    lb.back() = 0;
//    ub.back() = 5;
//    cbs.setDomainBounds(lb, ub);
//    cout << "Old lb:" << endl;
//    cout << lb << endl;
//    cout << "Old ub:" << endl;
//    cout << ub << endl;

//    cbs.controlPointBoundsDeduction(lb, ub);

//    cout << "New lb:" << endl;
//    cout << lb << endl;
//    cout << "New ub:" << endl;
//    cout << ub << endl;
//}

void testNewConstraintClasses()
{
    // Create some vars
    std::vector<VariablePtr> vars;

    for (int i = 0; i < 3; i++)
    {
        auto var = std::make_shared<Variable>(1, 0, i+1);
        vars.push_back(var);
        cout << *var << endl;
    }

    ConstraintSetPtr cs = std::make_shared<ConstraintSet>();

    DenseMatrix A = DenseMatrix::Zero(3,3);
    DenseVector b = DenseVector::Zero(3);
    b(0) = 0.5;

    A(0,0) = 1;
    ConstraintPtr lincon = std::make_shared<ConstraintLinear>(vars, A, b, true);

    A(1,1) = 2;
    ConstraintPtr lincon2 = std::make_shared<ConstraintLinear>(vars, A, b, true);

    A(2,2) = 3;
    ConstraintPtr lincon3 = std::make_shared<ConstraintLinear>(vars, A, b, true);

    cs->add(lincon);
    //cs->add(lincon2); // TODO: continue here tomorrow!
    //cs->add(lincon3);

    ConstraintSetPtr cs2 = std::make_shared<ConstraintSet>(*cs); // Shallow copy

    ConstraintSetPtr cs4 = std::make_shared<ConstraintSet>();
    cs4->add(cs); // Shallow copy

    ConstraintSetPtr cs5 = std::make_shared<ConstraintSet>(*cs, true); // Deep copy

    DenseVector x = DenseVector::Zero(3);
    x(0) = 1;
    x(1) = 1;
    x(2) = 1;

    DenseVector y;
    DenseVector dy;
    DenseVector ddy;

    {
        cout << "Evaluating cs" << endl;

        y = cs->eval(x);
        dy = cs->evalJacobian(x);
        ddy = cs->evalHessian(x);
        cout << "y:" << endl;
        cout << y << endl;
        cout << "dy" << endl;
        cout << dy << endl;
        cout << "ddy" << endl;
        cout << ddy << endl;
    }

    ConstraintPtr cs3 = cs->clone(); // Deep copy

    {
        cout << "Evaluating cs without x" << endl;

        std::vector<VariablePtr> vars = cs->getVariables();
        vars.at(0)->setValue(2);

        y = cs->eval();
        dy = cs->evalJacobian();
        ddy = cs->evalHessian();
        cout << "y:" << endl;
        cout << y << endl;
        cout << "dy" << endl;
        cout << dy << endl;
        cout << "ddy" << endl;
        cout << ddy << endl;
    }

    {
        cout << "Evaluating cs2 (shallow copy)." << endl;

        y = cs2->eval();
        dy = cs2->evalJacobian();
        ddy = cs2->evalHessian();
        cout << "y:" << endl;
        cout << y << endl;
        cout << "dy" << endl;
        cout << dy << endl;
        cout << "ddy" << endl;
        cout << ddy << endl;
    }

    {
        cout << "Evaluating cs3 (deep copy)." << endl;

        y = cs3->eval();
        dy = cs3->evalJacobian();
        ddy = cs3->evalHessian();
        cout << "y:" << endl;
        cout << y << endl;
        cout << "dy" << endl;
        cout << dy << endl;
        cout << "ddy" << endl;
        cout << ddy << endl;
    }

    {
        cout << "Evaluating cs4 (shallow composite)." << endl;

        y = cs4->eval();
        dy = cs4->evalJacobian();
        ddy = cs4->evalHessian();
        cout << "y:" << endl;
        cout << y << endl;
        cout << "dy" << endl;
        cout << dy << endl;
        cout << "ddy" << endl;
        cout << ddy << endl;

    }

    {
        cout << "Evaluating cs5 (deep copy)." << endl;

        y = cs5->eval();
        dy = cs5->evalJacobian();
        ddy = cs5->evalHessian();
        cout << "y:" << endl;
        cout << y << endl;
        cout << "dy" << endl;
        cout << dy << endl;
        cout << "ddy" << endl;
        cout << ddy << endl;
    }
}

void testConvexRelaxations()
{
    DenseMatrix A(1,3);
    A << 1, 2, 3;

    DenseVector b(1);
    b(0) = 1;

    std::vector<VariablePtr> vars =
    {
        std::make_shared<Variable>(1,0,6),
        std::make_shared<Variable>(2,0,4),
        std::make_shared<Variable>(3,0,2)
    };

    ConstraintPtr lincon1 = std::make_shared<ConstraintLinear>(vars, A, b, true);

    A << 3, 2, 1;
    ConstraintPtr lincon2 = std::make_shared<ConstraintLinear>(vars, A, b, true);

    Splinter::DataTable data;
    data.addSample(0,0);
    data.addSample(0.5,0.5);
    data.addSample(1,1);

    Splinter::BSpline bs(data, Splinter::BSplineType::LINEAR);
    auto bsvars = {vars.at(0), vars.at(1)};
    ConstraintPtr bscon = std::make_shared<ConstraintBSpline>(bsvars, bs, true);

    ConstraintSetPtr cs = std::make_shared<ConstraintSet>();
    cs->add(lincon1);
    cs->add(lincon2);
    cs->add(bscon);

    cout << cs->getVariables() << endl;

    cout << "Relaxed" << endl;

    ConstraintPtr csr = cs->getConvexRelaxation();

    cout << csr->getVariables() << endl;
}
