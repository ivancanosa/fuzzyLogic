#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <iostream>
#include <memory>
#include <numeric>
#include <queue>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <variant>
#include <vector>

// ******************** Data type definitions *******************8

using FuzzySetId = unsigned int;

struct FuzzySetTriangular {
    float left_bound, peak, right_bound;
};

struct FuzzySetTrapezoidal {
    float left_bound, left_peak, right_peak, right_bound;
};

struct FuzzySetGaussian {
    float mean, standard_deviation;
};

// The variant type for a membership function
using MembershipFunction =
    std::variant<FuzzySetTriangular, FuzzySetTrapezoidal, FuzzySetGaussian>;

struct FuzzyNotOp;
struct FuzzyAndOp;
struct FuzzyOrOp;
struct FuzzyIfOp;
using FuzzyOp =
    std::variant<FuzzySetId, std::unique_ptr<FuzzyNotOp>,
                 std::unique_ptr<FuzzyAndOp>, std::unique_ptr<FuzzyOrOp>,
                 std::unique_ptr<FuzzyIfOp>>;

struct FuzzyNotOp {
    FuzzyOp lhs;
};

struct FuzzyAndOp {
    FuzzyOp lhs;
    FuzzyOp rhs;
};

struct FuzzyOrOp {
    FuzzyOp lhs;
    FuzzyOp rhs;
};

struct FuzzyIfOp {
    FuzzyOp lhs;
    FuzzyOp rhs;
};

struct FuzzySystem {
    std::vector<MembershipFunction> fuzzySets{};
    std::vector<FuzzyIfOp> fuzzyRules{};

    // Internal data
    std::vector<FuzzySetId> outputFuzzySets{};
    std::vector<float> fuzzyMembershipValues{};
};

// ***************** Function prototypes definitions *******************

FuzzyOp clone(const FuzzyOp& fuzzyOp);
void printFuzzyOp(const FuzzyOp& op, bool root);

// Fuzzizy a value to a membership set
float computeMembership(float x, FuzzySetTriangular mf);
float computeMembership(float x, FuzzySetTrapezoidal mf);
float computeMembership(float x, FuzzySetGaussian mf);
float computeMembership(float x, MembershipFunction mf);

// Divide a range [0-1] in multiple fuzzy sets of maximum height 1
std::vector<MembershipFunction>
createTriangularFuzzySets(int fuzzySetCount, float overlapingAmmount);
std::vector<MembershipFunction>
createTrapezoidalFuzzySets(int fuzzySetCount, float overlapingAmmount);
std::vector<MembershipFunction>
createGaussianFuzzySets(int fuzzySetCount, float overlapingAmmount);

// Add fuzzy sets to the fuzzy system
FuzzyOp addFuzzySet(FuzzySystem& fs, const MembershipFunction& mf);
std::vector<FuzzyOp>
addFuzzySet(FuzzySystem& fs, const std::vector<MembershipFunction>& mfVector);

// Set up rules
void fuzzySystemAddRule(FuzzySystem& fs, const FuzzyOp& node);
void fuzzySystemAddRule(FuzzySystem& fs, const FuzzyOp&& node);

// Topological sort
std::vector<std::unordered_set<FuzzySetId>>
buildDependencyGraph(const FuzzySystem& fs);
std::vector<FuzzySetId> fuzzySetsSort(const FuzzySystem& fs);
void fuzzyRulesTopologicalSort(FuzzySystem& fs);

// Insert the input to the fuzzy system
void setFuzzyInput(FuzzySystem& fs, const FuzzyOp& fuzzySetId, float value);
void setFuzzyInput(FuzzySystem& fs, const std::vector<FuzzyOp>& fuzzySetVector,
                   float value);

// Compute fuzzy rules
float computeFuzzyRule(const FuzzySystem& fs, const FuzzyOp& fuzzyOp);
void computeFuzzyRules(FuzzySystem& fs);

// Defuzzify
float defuzzify(const FuzzySystem& fs, const FuzzyOp& fuzzySetVector);
float defuzzify(const FuzzySystem& fs,
                const std::vector<FuzzyOp>& fuzzySetVector);

// Overload of all logic operations to compute an expression tree
bool operator==(const FuzzyOp& lhs, const FuzzyOp& rhs);
bool operator==(const FuzzyOp&& lhs, const FuzzyOp& rhs);
bool operator==(const FuzzyOp& lhs, const FuzzyOp&& rhs);
bool operator==(const FuzzyOp&& lhs, const FuzzyOp&& rhs);
FuzzyOp operator!(FuzzyOp& lhs);
FuzzyOp operator!(FuzzyOp&& lhs);
FuzzyOp operator&&(FuzzyOp& lhs, FuzzyOp& rhs);
FuzzyOp operator&&(FuzzyOp&& lhs, FuzzyOp& rhs);
FuzzyOp operator&&(FuzzyOp& lhs, FuzzyOp&& rhs);
FuzzyOp operator&&(FuzzyOp&& lhs, FuzzyOp&& rhs);
FuzzyOp operator||(FuzzyOp& lhs, FuzzyOp& rhs);
FuzzyOp operator||(FuzzyOp&& lhs, FuzzyOp& rhs);
FuzzyOp operator||(FuzzyOp& lhs, FuzzyOp&& rhs);
FuzzyOp operator||(FuzzyOp&& lhs, FuzzyOp&& rhs);
FuzzyOp operator>>(FuzzyOp& lhs, FuzzyOp& rhs);
FuzzyOp operator>>(FuzzyOp&& lhs, FuzzyOp& rhs);
FuzzyOp operator>>(FuzzyOp& lhs, FuzzyOp&& rhs);
FuzzyOp operator>>(FuzzyOp&& lhs, FuzzyOp&& rhs);
