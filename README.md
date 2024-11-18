# Fuzzy Logic System
A fuzzy logic system implementation in C++17. The objective of this project is to create a flexible, easy to use interface so it can be used without much boilerplate.

## Requirements
- CMake +3.10
- A C++17 compiler

## Features
- Arbitrary logical expressions. You are able to use the operators NOT, AND, OR, to build any arbitrary logical expression in the body of the rules.
- Overloading of the operators "!", "&&", "||", ">>", so you can write the fuzzy rules in a natural way. The rules are of the form "(BODY) >> HEAD", where the HEAD must the identifier of a fuzzy set and BODY is an arbitrary logical expression of fuzzy set identifiers. It is important to note that the BODY must be between parenthesis because C++ does not allow me to change the priority of the operator ">>" (the IF operator) to be the lowest.
- Transitive fuzzy rules. This means that you can have for example the following 3 fuzzy rules:
	- a -> b
	- b -> c
	- c -> d <br>
The fuzzy logic system is able to perform a topological sort of the rules so they are computed in the correct order. This means that you must be aware that there can't be any cycles in the rules. The dependency between the fuzzy sets must be able to be represented as a Directed Acyclyc Graph (DAG).

## Compilation
```bash
# Create and navigate to the build directory
mkdir build
cd build

# Run CMake to generate build files
cmake ../src

# Compile the project
make

# Run tests
ctest

# Run simple example
./example/FuzzyLogicExample
```

## Files
- src/FuzzyLogic.hpp: Header for the fuzzy logic system with all data types and function prototypes.
- src/FuzzyLogic.cpp: Implementation of the fuzzy logic system.
- src/example/fuzzyLogicExample.cpp: An example of use of the fuzzy logic system for an agent in a videogame. 
- src/test/test.cpp: The tests that checks if the fuzzy logic system works.

## Data Types
- FuzzySetId: an unsigned integer representing a fuzzy set in the system.
- FuzzySetTriangular: a triangular membership function with left bound, peak, and right bound values.
- FuzzySetTrapezoidal: a trapezoidal membership function with left bound, left peak, right peak, and right bound values.
- FuzzySetGaussian: a gaussian membership function with mean and standard deviation values.
- MembershipFunction: a variant type containing either a FuzzySetTriangular, FuzzySetTrapezoidal, or FuzzySetGaussian membership function.
- FuzzyNotOp: a struct representing a fuzzy Not operator with a single operand.
- FuzzyAndOp: a struct representing a fuzzy And operator with two operands.
- FuzzyOrOp: a struct representing a fuzzy Or operator with two operands.
- FuzzyIfOp: a struct representing a fuzzy If operator with two operands.
- FuzzyOp: a variant type containing either a FuzzySetId, a FuzzyNotOp, a FuzzyAndOp, a FuzzyOrOp, or a FuzzyIfOp.
- FuzzySystem: a struct containing vectors of MembershipFunction, FuzzyIfOp, FuzzySetId, and float values representing the fuzzy sets, rules, output sets, and membership values in the system.

## Function Prototypes
- FuzzyOp clone(const FuzzyOp& fuzzyOp): creates a deep copy of a FuzzyOp object.
- void printFuzzyOp(const FuzzyOp& op, bool root): prints a string representation of a FuzzyOp object to the console.
- float computeMembership(float x, FuzzySetTriangular mf): computes the membership value of a given value x in a triangular membership function mf.
- float computeMembership(float x, FuzzySetTrapezoidal mf): computes the membership value of a given value x in a trapezoidal membership function mf.
- float computeMembership(float x, FuzzySetGaussian mf): computes the membership value of a given value x in a gaussian membership function mf.
- float computeMembership(float x, MembershipFunction mf): computes the membership value of a given value x in a membership function mf of any type.
- std::vector<MembershipFunction> createTriangularFuzzySets(int fuzzySetCount, float overlappingAmount): creates a vector of triangular membership functions divided into fuzzySetCount sets with a given overlappingAmount of overlap, which must be in the range [0-1].
- std::vector<MembershipFunction> createTrapezoidalFuzzySets(int fuzzySetCount, float overlapingAmmount); creates a vector of trapezoidal membership functions divided into fuzzySetCount sets with a given overlappingAmount of overlap, which must be in the range [0-1].
- std::vector<MembershipFunction> createGaussianFuzzySets(int fuzzySetCount, float overlapingAmmount); creates a vector of gaussian membership functions divided into fuzzySetCount sets with a given overlappingAmount of overlap, which must be in the range [0-1].
- FuzzyOp addFuzzySet(FuzzySystem& fs, const MembershipFunction& mf): Adds a fuzzy set to the fuzzy system and return its id.
- std::vector<FuzzyOp> addFuzzySet(FuzzySystem& fs, const std::vector<MembershipFunction>& mfVector): Adds a vector of fuzzy sets to the fuzzy system and returns their id.
- void fuzzySystemAddRule(FuzzySystem& fs, const FuzzyOp& node): Adds a fuzzy rule to the fuzzy set. Must be of the form "a -> b", an implication.
- void fuzzySystemAddRule(FuzzySystem& fs, const FuzzyOp&& node): Adds a fuzzy rule to the fuzzy set. Must be of the form "a -> b", an implication.
- void fuzzyRulesTopologicalSort(FuzzySystem& fs): Performs a topological sort of all the fuzzy rules. Must be called after all rules has been inserted and before to insert input values to the fuzzy sets.
- void setFuzzyInput(FuzzySystem& fs, const FuzzyOp& fuzzySetId, float value): Inserts the given crips value to the given fuzzy set.
- void setFuzzyInput(FuzzySystem& fs, const std::vector<FuzzyOp>& fuzzySetVector, float value): Inserts the given crips value to the given vector of fuzzy sets.
- void computeFuzzyRules(FuzzySystem& fs): Computes the membership value of all the output fuzzy sets.
- float defuzzify(const FuzzySystem& fs, const FuzzyOp& fuzzySetVector): Computes the crips value from the given fuzzy set.
- float defuzzify(const FuzzySystem& fs, const std::vector<FuzzyOp>& fuzzySetVector): Computes the crips value from the given vector of fuzzy sets.


## TODO
- Change the code so it is possible to compute in parallel multiple fuzzy inferences with the same rules set.
- The functions that create membership functions over the range [0-1] have the limitation that all the functions have the same area. An improvment is to be able to tell the functions which membership functions should have more are than the others.
- Add the capability of understand predicates and rules with a syntax similar to clingo or prolog.
- The logical expressions are represented as a tree where the memory allocation is in random places in the heap. An optimization would be to compact all of that on a single vector, so all logical expressions and logic rules are on the same cache line, making the inference faster without cache misses.
- Add the capability to work as a fuzzy production system with a working memory.
