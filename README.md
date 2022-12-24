## Logic System
A fuzzy logic system implementation in C++.

# Requirements
To compile the project you will need Makefile, clang++ and bear. However, you can copy the FuzzyLogic does not have any dependencies and you can copy and paste them in your project.

# Compilation
- make test: Compile and run the tests
- make run: Compile and run the main.cpp file, which contains an example of use of the logic system for an agent in a videogame.

# Data Types
- FuzzySetId: an unsigned integer representing a fuzzy set in the system.
- FuzzySetTriangular: a triangular membership function with left bound, peak, and right bound values.
- FuzzySetTrapezoidal: a trapezoidal membership function with left bound, left peak, right peak, and right bound values.
- FuzzySetGaussian: a gaussian membership function with mean and standard deviation values.
- MembershipFunction: a variant type containing either a FuzzySetTriangular, FuzzySetTrapezoidal, or FuzzySetGaussian membership function.
- FuzzySetType: an enum representing the type of membership function, with values Triangular, Trapezoidal, or Gaussian.
- FuzzyOperator: an enum representing the type of fuzzy operator, with values Not, And, Or, or If.
- FuzzyNotOp: a struct representing a fuzzy Not operator with a single operand.
- FuzzyAndOp: a struct representing a fuzzy And operator with two operands.
- FuzzyOrOp: a struct representing a fuzzy Or operator with two operands.
- FuzzyIfOp: a struct representing a fuzzy If operator with two operands.
- FuzzyOp: a variant type containing either a FuzzySetId, a FuzzyNotOp, a FuzzyAndOp, a FuzzyOrOp, or a FuzzyIfOp.
- FuzzySystem: a struct containing vectors of MembershipFunction, FuzzyIfOp, FuzzySetId, and float values representing the fuzzy sets, rules, output sets, and membership values in the system.

# Function Prototypes
- FuzzyOp clone(const FuzzyOp& fuzzyOp): creates a deep copy of a FuzzyOp object.
- void printFuzzyOp(const FuzzyOp& op, bool root): prints a string representation of a FuzzyOp object.
- float computeMembership(float x, FuzzySetTriangular mf): computes the membership value of a given value x in a triangular membership function mf.
- float computeMembership(float x, FuzzySetTrapezoidal mf): computes the membership value of a given value x in a trapezoidal membership function mf.
- float computeMembership(float x, FuzzySetGaussian mf): computes the membership value of a given value x in a gaussian membership function mf.
- float computeMembership(float x, MembershipFunction mf): computes the membership value of a given value x in a membership function mf of any type.
- std::vector<MembershipFunction> createTriangularFuzzySets(int fuzzySetCount, float overlappingAmount): creates a vector of triangular membership functions divided into fuzzySetCount sets with a given overlappingAmount of overlap.
