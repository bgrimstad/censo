# Bonmin 1.7.4 libs
#LIBS += -L/home/bjarne/Solvers/Bonmin/Bonmin-1.7.5/build/lib -lbonmin -lipopt -lCbc -lCgl -lOsiClp -lOsi -lClp -lCoinUtils -lz
LIBS += -L/home/bjarne/Solvers/Bonmin/Bonmin-1.7.5/build/lib -L/usr/lib/gcc/x86_64-linux-gnu/4.8 -L/usr/lib/gcc/x86_64-linux-gnu/4.8/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/4.8/../../../../lib -L/lib/../lib -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/4.8/../../.. -lbonmin -lipopt -lCbc -lCgl -lOsiClp -lOsi -lClp -lCoinUtils -ldl -lcoinhsl -lcoinlapack -lcoinblas -lgfortran -lm -lquadmath -lz
INCLUDEPATH += /home/bjarne/Solvers/Bonmin/Bonmin-1.7.5/build/include/coin

# Ipopt 3.11.8 libs
#LIBS += -L/home/bjarne/Solvers/Ipopt/Ipopt-3.11.8/build/lib -L/usr/lib/gcc/x86_64-linux-gnu/4.8 -L/usr/lib/gcc/x86_64-linux-gnu/4.8/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/4.8/../../../../lib -L/lib/../lib -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/4.8/../../.. -lipopt -ldl -lcoinhsl -lcoinlapack -lcoinblas -lgfortran -lm -lquadmath
#INCLUDEPATH += /home/bjarne/Solvers/Ipopt/Ipopt-3.11.8/build/include/coin

# Gurobi 5.6.0 libs
LIBS += -L/home/bjarne/Solvers/Gurobi/gurobi600/linux64/lib -lgurobi_c++ -lgurobi60
INCLUDEPATH += /home/bjarne/Solvers/Gurobi/gurobi600/linux64/include

# Eigen library
INCLUDEPATH += /home/bjarne/C++/libs/eigen-3.2.4

# Splinter library
LIBS += -L/home/bjarne/C++/splinter/splinter/build/lib -lsplinter-1-2
INCLUDEPATH += /home/bjarne/C++/splinter/splinter/build/include

QMAKE_CXXFLAGS_RELEASE += -O2 # level of optimization
QMAKE_CXXFLAGS += -std=c++11 # use c++11
QMAKE_CXXFLAGS_WARN_ON = ""
QMAKE_CXXFLAGS += -W
QMAKE_CXXFLAGS += -Wall
QMAKE_CXXFLAGS += -Wno-unused-variable
QMAKE_CXXFLAGS += -Wno-unused-parameter
