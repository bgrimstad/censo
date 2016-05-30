project(CENSO CXX)
cmake_minimum_required(VERSION 2.8)
string(TOLOWER ${PROJECT_NAME} PROJECT_NAME_LOWER)

# Default configuration values
set(DEFAULT_BUILD_TYPE "Release")
set(DEFAULT_ARCH "x86-64")
set(DEFAULT_HEADER_INSTALL_DIRECTORY "include")
set(DEFAULT_LIBRARY_INSTALL_DIRECTORY "lib")
set(BITNESS "64")

set(DEFAULT_EIGEN_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/Thirdparty/Eigen)

# Set CMake settings
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m${BITNESS} -std=c++11 -Werror=return-type")
#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wno-long-long")

# Set build type
set(CMAKE_BUILD_TYPE ${DEFAULT_BUILD_TYPE})

# Include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/Thirdparty/Eigen)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/Thirdparty/splinter/include)
include_directories(/home/bjarne/Solvers/Gurobi/gurobi605/linux64/include)
include_directories(/home/bjarne/Solvers/Bonmin/Bonmin-1.7.5/build/include/coin)

# Bonmin and Ipopt link directories
link_directories(/home/bjarne/Solvers/Bonmin/Bonmin-1.7.5/build/lib)
link_directories(/usr/lib/gcc/x86_64-linux-gnu/4.8)
link_directories(/usr/x86_64-linux-gnu)
link_directories(/usr/lib)
link_directories(/lib)

# Gurobi and SPLINTER link directories
link_directories(/home/bjarne/Solvers/Gurobi/gurobi605/linux64/lib)
link_directories(${CMAKE_CURRENT_SOURCE_DIR}/Thirdparty/splinter/)

# These are the sources we need for compilation of the library
set(SRC_LIST
    Utils/definitions.h
    Utils/metrics.h
    Utils/timer.h
    Utils/linearsolvers.h
    Utils/linearinterpolation.h
    Utils/bsplinepoly.h
    OptimizationProblem/variable.h
    OptimizationProblem/constraint.h
    OptimizationProblem/constraintlinear.h
    OptimizationProblem/constraintsine.h
    OptimizationProblem/constraintset.h
    OptimizationProblem/constraintquadratic.h
    OptimizationProblem/constraintbilinear.h
    OptimizationProblem/constraintpolynomial.h
    OptimizationProblem/constraintexponential.h
    OptimizationProblem/constraintbspline.h
    SolverInterface/solver.h
    SolverInterface/solveripopt.h
    SolverInterface/interfaceipopt.h
    SolverInterface/solverbonmin.h
    SolverInterface/interfacebonmin.h
    SolverInterface/solvergurobi.h
    BranchAndBound/node.h
    BranchAndBound/nodelist.h
    BranchAndBound/branchandbound.h
    BranchAndBound/bbutils.h
    BranchAndBound/Branching/brancher.h
    BranchAndBound/Branching/brancherinteger.h
    BranchAndBound/Branching/branchercontinuous.h
    BranchAndBound/Branching/branchercontinuousrandomvariable.h
    BranchAndBound/Branching/branchercontinuousmostpromising.h
    BranchAndBound/Branching/branchercontinuouslongestinterval.h
    BranchAndBound/Branching/brancherintegermostfractional.h
    BranchAndBound/Branching/brancheradmin.h
    BranchAndBound/BoundsTightening/boundstightener.h
    BranchAndBound/BoundsTightening/obbt.h
    BranchAndBound/BoundsTightening/fbbt.h
    Projects/JOGO/jogoexamples.h
    Projects/FlowNetwork/graph.h
    Projects/FlowNetwork/flownetwork.h
    Projects/FlowNetwork/flownetworkdata.h
    TestProblems/testproblem.h
    TestProblems/alltestproblems.h
    TestProblems/pop01.h
    TestProblems/pop02.h
    TestProblems/pop03.h
    TestProblems/pop04.h
    TestProblems/pop05.h
    TestProblems/pop06.h
    TestProblems/pop07.h
    TestProblems/pop08.h
    TestProblems/pop09.h
    TestProblems/pop10.h
    TestProblems/pop11.h
    TestProblems/pop12.h
    TestProblems/pop13.h
    TestProblems/michalewicz.h
    TestProblems/lp01.h
    TestProblems/lp02.h
    TestProblems/lp03.h
    TestProblems/qp01.h
    TestProblems/sixhumpcamel.h
    TestProblems/bilinearrelaxationtest.h
    TestProblems/asaadi1.h
    TestProblems/linearintervalanalysistest.h
    Test/unittests.h
    Utils/definitions.cpp
    Utils/metrics.cpp
    Utils/linearinterpolation.cpp
    Utils/bsplinepoly.cpp
    OptimizationProblem/constraint.cpp
    OptimizationProblem/constraintlinear.cpp
    OptimizationProblem/constraintsine.cpp
    OptimizationProblem/constraintset.cpp
    OptimizationProblem/constraintquadratic.cpp
    OptimizationProblem/constraintbilinear.cpp
    OptimizationProblem/constraintpolynomial.cpp
    OptimizationProblem/constraintexponential.cpp
    OptimizationProblem/constraintbspline.cpp
    OptimizationProblem/writetogams.cpp
    SolverInterface/solveripopt.cpp
    SolverInterface/interfaceipopt.cpp
    SolverInterface/solverbonmin.cpp
    SolverInterface/interfacebonmin.cpp
    SolverInterface/solvergurobi.cpp
    BranchAndBound/node.cpp
    BranchAndBound/nodelist.cpp
    BranchAndBound/branchandbound.cpp
    BranchAndBound/bbutils.cpp
    BranchAndBound/Branching/brancher.cpp
    BranchAndBound/Branching/brancherinteger.cpp
    BranchAndBound/Branching/branchercontinuous.cpp
    BranchAndBound/Branching/branchercontinuousrandomvariable.cpp
    BranchAndBound/Branching/branchercontinuousmostpromising.cpp
    BranchAndBound/Branching/branchercontinuouslongestinterval.cpp
    BranchAndBound/Branching/brancherintegermostfractional.cpp
    BranchAndBound/Branching/brancheradmin.cpp
    BranchAndBound/BoundsTightening/obbt.cpp
    BranchAndBound/BoundsTightening/fbbt.cpp
    Projects/JOGO/jogoexamples.cpp
    Projects/FlowNetwork/flownetwork.cpp
    TestProblems/testproblem.cpp
    TestProblems/pop01.cpp
    TestProblems/pop02.cpp
    TestProblems/pop03.cpp
    TestProblems/pop04.cpp
    TestProblems/pop05.cpp
    TestProblems/pop06.cpp
    TestProblems/pop07.cpp
    TestProblems/pop08.cpp
    TestProblems/pop09.cpp
    TestProblems/pop10.cpp
    TestProblems/pop11.cpp
    TestProblems/pop12.cpp
    TestProblems/pop13.cpp
    TestProblems/michalewicz.cpp
    TestProblems/lp01.cpp
    TestProblems/lp02.cpp
    TestProblems/lp03.cpp
    TestProblems/qp01.cpp
    TestProblems/sixhumpcamel.cpp
    TestProblems/bilinearrelaxationtest.cpp
    TestProblems/asaadi1.cpp
    TestProblems/linearintervalanalysistest.cpp
    Test/testproblems.cpp
    Test/unittests.cpp
    main.cpp
)

# Testing executable
set(TEST "${PROJECT_NAME_LOWER}-test")
add_executable(${TEST} ${SRC_LIST})
target_link_libraries(${TEST} bonmin ipopt Cbc Cgl OsiClp Osi Clp CoinUtils dl coinhsl coinlapack coinblas gfortran m quadmath z gurobi_c++ gurobi60 splinter-1-4)
