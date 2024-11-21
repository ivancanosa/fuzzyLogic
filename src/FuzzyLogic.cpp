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

#include "FuzzyLogic.hpp"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

//******************* Topological sort ************************

// Helper function that constructs a dependency graph of the fuzzy rules in the
// FuzzySystem
std::vector<std::unordered_set<FuzzySetId>>
buildDependencyGraph(const FuzzySystem& fs) {
    std::vector<std::unordered_set<FuzzySetId>> graph(fs.fuzzySets.size());

    // Recursive lambda that traverses a rule to get all the sets
    std::function<void(const FuzzyOp&, FuzzySetId)> computeDependencies;
    computeDependencies = [&computeDependencies, &graph](auto&& fuzzyOp,
                                                         FuzzySetId nodeId) {
        std::visit(
            [&computeDependencies, &nodeId, &graph](auto&& arg) {
                using T = std::decay_t<decltype(arg)>;
                if constexpr (std::is_same_v<T, FuzzySetId>) {
                    graph[nodeId].insert(arg);
                } else if constexpr (std::is_same_v<
                                         T, std::unique_ptr<FuzzyNotOp>>) {
                    graph[nodeId].insert(std::get<FuzzySetId>(arg->lhs));
                } else if constexpr (std::is_same_v<
                                         T, std::unique_ptr<FuzzyAndOp>>) {
                    computeDependencies(arg->lhs, nodeId);
                    computeDependencies(arg->rhs, nodeId);
                } else if constexpr (std::is_same_v<
                                         T, std::unique_ptr<FuzzyOrOp>>) {
                    computeDependencies(arg->lhs, nodeId);
                    computeDependencies(arg->rhs, nodeId);
                } else {
                    assert(false && "Fuzzy operation not allowed");
                }
            },
            fuzzyOp);
    };

    for (auto&& rule : fs.fuzzyRules) {
        // Compute the output fuzzy set ID for the current rule

        FuzzySetId outputFuzzySet = std::get<FuzzySetId>(rule.rhs);
        computeDependencies(rule.lhs, outputFuzzySet);
    }
    return graph;
}

std::vector<FuzzySetId> fuzzySetsSort(FuzzySystem& fs) {
    if (fs.fuzzyRules.empty()) {
        return {};
    }

    // Build the dependency graph
    auto graph = buildDependencyGraph(fs);

    // Create the result vector
    std::vector<FuzzySetId> result;
    result.reserve(fs.fuzzyRules.size());

    // Create a vector to keep track of the in-degree for each node
    std::vector<int> inDegrees(graph.size());
    std::fill(inDegrees.begin(), inDegrees.end(), 0);
    for (int i = 0; i < graph.size(); i++) {
        for (auto neighbor : graph[i]) {
            inDegrees[neighbor]++;
        }
    }

    // Create a queue for the nodes with in-degree 0
    // Get all the sets that are output
    std::queue<FuzzySetId> queue;
    fs.outputFuzzySets.clear();
    for (int i = 0; i < graph.size(); i++) {
        if (inDegrees[i] == 0) {
            queue.push(i);
            fs.outputFuzzySets.push_back(i);
        } else {
        }
    }

    // Perform the topological sort
    while (!queue.empty()) {
        // Get the next node with in-degree 0
        auto node = queue.front();
        queue.pop();

        // Add the node to the result
        result.push_back(node);

        // Decrement the in-degree for all the neighbors of the node
        for (auto neighbor : graph[node]) {
            inDegrees[neighbor]--;

            // If the in-degree of the neighbor is now 0, add it to the queue
            if (inDegrees[neighbor] == 0) {
                queue.push(neighbor);
            }
        }
    }

    // Check if there are any nodes with non-zero in-degree, indicating a cycle
    // in the graph
    for (auto degree : inDegrees) {
        if (degree != 0) {
            // There is a cycle in the graph, return an empty vector
            return {};
        }
    }

    return result;
}

void fuzzyRulesTopologicalSort(FuzzySystem& fs) {
	//Initialize all the membership values of the sets to 0
	for(int i=0; i<fs.fuzzyMembershipValues.size(); i++){
		fs.fuzzyMembershipValues[i] = 0.;
	}

    auto sortedSets = fuzzySetsSort(fs);
    // Create a map from fuzzy set id to position in the sorted sets vector
    std::vector<int> setPositions(sortedSets.size());
    for (int i = 0; i < sortedSets.size(); i++) {
        setPositions[sortedSets[i]] = i;
    }

    // Sort the rules using the set positions as a comparison key
    auto cmp = [&setPositions](const FuzzyIfOp& lhs, const FuzzyIfOp& rhs) {
        auto first  = setPositions[std::get<FuzzySetId>(lhs.rhs)];
        auto second = setPositions[std::get<FuzzySetId>(rhs.rhs)];
        return first > second;
    };
    std::sort(fs.fuzzyRules.begin(), fs.fuzzyRules.end(), cmp);
}

//******************* Fuzzy logic inference ************************

float computeFuzzyRule(const FuzzySystem& fs, const FuzzyOp& fuzzyOp) {
    float value = 0.;
    std::visit(
        [&fuzzyOp, &value, &fs](auto&& arg) {
            using T = std::decay_t<decltype(arg)>;
            if constexpr (std::is_same_v<T, FuzzySetId>) {
                value = fs.fuzzyMembershipValues[arg];
            } else if constexpr (std::is_same_v<T,
                                                std::unique_ptr<FuzzyNotOp>>) {
                float lhs = computeFuzzyRule(fs, arg->lhs);
                value     = 1. - lhs;
            } else if constexpr (std::is_same_v<T,
                                                std::unique_ptr<FuzzyAndOp>>) {
                float lhs = computeFuzzyRule(fs, arg->lhs);
                float rhs = computeFuzzyRule(fs, arg->rhs);
                value     = std::min(lhs, rhs);
            } else if constexpr (std::is_same_v<T,
                                                std::unique_ptr<FuzzyOrOp>>) {
                float lhs = computeFuzzyRule(fs, arg->lhs);
                float rhs = computeFuzzyRule(fs, arg->rhs);
                value     = std::max(lhs, rhs);
            } else {
                assert(false && "Forbidden fuzzy operation");
            }
        },
        fuzzyOp);
    return value;
}

void computeFuzzyRules(FuzzySystem& fs) {
    // Set all output fuzzy sets to 0
    int aux = 0;
    for (FuzzySetId id : fs.outputFuzzySets) {
        fs.fuzzyMembershipValues[id] = 0.;
        aux += 1;
    }

    // Iterate over the fuzzy rules and apply each rule to the membership values
    // in the fuzzyMembershipValues member
    for (auto&& rule : fs.fuzzyRules) {
        // Compute the output of the fuzzy rule
        float output = computeFuzzyRule(fs, rule.lhs);

        // Update the membership value of the output fuzzy set
        // with the computed output of the fuzzy rule
        FuzzySetId outputFuzzySet = std::get<FuzzySetId>(rule.rhs);
        fs.fuzzyMembershipValues[outputFuzzySet] =
            std::max(fs.fuzzyMembershipValues[outputFuzzySet], output);
    }
}

void fuzzySystemAddRule(FuzzySystem& fs, const FuzzyOp& node) {
    auto op = clone(node);
    std::unique_ptr<FuzzyIfOp> a =
        std::move(std::get<std::unique_ptr<FuzzyIfOp>>(op));
    FuzzyIfOp result;
    result.lhs = std::move(a->lhs);
    result.rhs = std::move(a->rhs);
    fs.fuzzyRules.push_back(std::move(result));
}

void fuzzySystemAddRule(FuzzySystem& fs, const FuzzyOp&& node) {
    auto op = clone(node);
    std::unique_ptr<FuzzyIfOp> a =
        std::move(std::get<std::unique_ptr<FuzzyIfOp>>(op));
    FuzzyIfOp result;
    result.lhs = std::move(a->lhs);
    result.rhs = std::move(a->rhs);
    fs.fuzzyRules.push_back(std::move(result));
}

//**************************************************************************
//*********************** Implementation ***********************************
//**************************************************************************

void setFuzzyInput(FuzzySystem& fs, const FuzzyOp& fuzzySetId, float value) {
    int setId = std::get<FuzzySetId>(fuzzySetId);
    fs.fuzzyMembershipValues[setId] =
        computeMembership(value, fs.fuzzySets[setId]);
}

void setFuzzyInput(FuzzySystem& fs, const std::vector<FuzzyOp>& fuzzySetVector,
                   float value) {
    for (const auto& fuzzySetId : fuzzySetVector) {
        setFuzzyInput(fs, fuzzySetId, value);
    }
}

FuzzyOp addFuzzySet(FuzzySystem& fs, const MembershipFunction& mf) {
    // Add the membership function to the fuzzy system's list of fuzzy sets
    fs.fuzzySets.push_back(mf);

    // Add a value to the fuzzyMembershipValues field with the default value 0.0
    fs.fuzzyMembershipValues.push_back(0.0);

    // Return a FuzzyOp representing the added fuzzy set
    return FuzzySetId(fs.fuzzySets.size() - 1);
}

std::vector<FuzzyOp>
addFuzzySet(FuzzySystem& fs, const std::vector<MembershipFunction>& mfVector) {
    std::vector<FuzzyOp> fuzzyOps;

    for (auto& mf : mfVector) {
        fuzzyOps.push_back(addFuzzySet(fs, mf));
    }

    return fuzzyOps;
}

std::vector<MembershipFunction>
createTriangularFuzzySets(int fuzzySetCount, float overlappingAmount) {
    std::vector<MembershipFunction> sets;
    float pointCount = fuzzySetCount * 2 - 2;
    float peak;
    float leftBound;
    float rightBound;
    for (int i = 0; i < fuzzySetCount; i++) {
        peak       = i * 2 / pointCount;
        leftBound  = (i * 2 - 1) / pointCount - overlappingAmount / 2.;
        rightBound = (i * 2 + 1) / pointCount + overlappingAmount / 2.;
        sets.emplace_back(FuzzySetTriangular{leftBound, peak, rightBound});
    }
    return sets;
}

std::vector<MembershipFunction>
createTrapezoidalFuzzySets(int fuzzySetCount, float overlappingAmount) {
    std::vector<MembershipFunction> sets;
    float pointCount = fuzzySetCount * 2 - 2;
    float leftPeak;
    float rightPeak;
    float leftBound;
    float rightBound;
    float overlap = overlappingAmount * 1. / (fuzzySetCount * 2.);
    for (int i = 0; i < fuzzySetCount; i++) {
        leftBound  = (i * 2 - 1) / pointCount - overlap;
        leftPeak   = (i * 2 - 1) / pointCount + overlap;
        rightBound = (i * 2 + 1) / pointCount + overlap;
        rightPeak  = (i * 2 + 1) / pointCount - overlap;

        sets.emplace_back(
            FuzzySetTrapezoidal{leftBound, leftPeak, rightPeak, rightBound});
    }
    return sets;
}

std::vector<MembershipFunction>
createGaussianFuzzySets(int fuzzySetCount, float overlapingAmmount) {
    std::vector<MembershipFunction> fuzzySets;

    // Determine the width of each fuzzy set
    float setWidth = (1.0 - overlapingAmmount) / fuzzySetCount;

    // Create each fuzzy set
    for (int i = 0; i < fuzzySetCount; ++i) {
        FuzzySetGaussian mf;
        mf.mean               = i * setWidth + setWidth / 2;
        mf.standard_deviation = setWidth / 3;

        fuzzySets.push_back(mf);
    }

    return fuzzySets;
}

float computeMembership(float x, FuzzySetTriangular mf) {
    if (x <= mf.left_bound || x >= mf.right_bound) {
        return 0.0;
    }
    if (x < mf.peak) {
        return (x - mf.left_bound) / (mf.peak - mf.left_bound);
    }
    return (mf.right_bound - x) / (mf.right_bound - mf.peak);
}

float computeMembership(float x, FuzzySetTrapezoidal mf) {
    if (x <= mf.left_bound || x >= mf.right_bound) {
        return 0.0;
    }
    if (x < mf.left_peak) {
        return (x - mf.left_bound) / (mf.left_peak - mf.left_bound);
    }
    if (x < mf.right_peak) {
        return 1.0;
    }
    return (mf.right_bound - x) / (mf.right_bound - mf.right_peak);
}

float computeMembership(float x, FuzzySetGaussian mf) {
    return std::exp(-0.5 * std::pow((x - mf.mean) / mf.standard_deviation, 2)) /
           (mf.standard_deviation * std::sqrt(2 * M_PI));
}

float computeMembership(float x, MembershipFunction mf) {
    return std::visit(
        [&x](auto&& arg) -> float { return computeMembership(x, arg); }, mf);
}

// *************** Rules operator overload *******************

FuzzyOp clone(const FuzzyOp& fuzzyOp) {
    FuzzyOp result;
    std::visit(
        [&result](const auto& op) {
            using T = std::decay_t<decltype(op)>;
            if constexpr (std::is_same_v<T, FuzzySetId>) {
                result = op;
            } else if constexpr (std::is_same_v<T,
                                                std::unique_ptr<FuzzyNotOp>>) {
                auto aux = std::make_unique<FuzzyNotOp>();
                aux->lhs = std::move(clone(op->lhs));
                result   = std::move(aux);
            } else if constexpr (std::is_same_v<T,
                                                std::unique_ptr<FuzzyAndOp>>) {
                auto aux = std::make_unique<FuzzyAndOp>();
                aux->lhs = std::move(clone(op->lhs));
                aux->rhs = std::move(clone(op->rhs));
                result   = std::move(aux);
            } else if constexpr (std::is_same_v<T,
                                                std::unique_ptr<FuzzyOrOp>>) {
                auto aux = std::make_unique<FuzzyOrOp>();
                aux->lhs = std::move(clone(op->lhs));
                aux->rhs = std::move(clone(op->rhs));
                result   = std::move(aux);
            } else if constexpr (std::is_same_v<T,
                                                std::unique_ptr<FuzzyIfOp>>) {
                auto aux = std::make_unique<FuzzyIfOp>();
                aux->lhs = std::move(clone(op->lhs));
                aux->rhs = std::move(clone(op->rhs));
                result   = std::move(aux);
            }
        },
        fuzzyOp);
    return std::move(result);
}

bool operator==(const FuzzyOp& lhs, const FuzzyOp& rhs) {
    if (lhs.index() != rhs.index()) {
        return false; // the variants hold different types
    }

    // Check if both variants hold a FuzzySetId
    if (std::holds_alternative<FuzzySetId>(lhs) &&
        std::holds_alternative<FuzzySetId>(rhs)) {
        return std::get<FuzzySetId>(lhs) == std::get<FuzzySetId>(rhs);
    }

    // Check if both variants hold a FuzzyNotOp
    if (std::holds_alternative<std::unique_ptr<FuzzyNotOp>>(lhs) &&
        std::holds_alternative<std::unique_ptr<FuzzyNotOp>>(rhs)) {
        const auto& lhsNotOp = std::get<std::unique_ptr<FuzzyNotOp>>(lhs);
        const auto& rhsNotOp = std::get<std::unique_ptr<FuzzyNotOp>>(rhs);
        return lhsNotOp->lhs == rhsNotOp->lhs;
    }

    // Check if both variants hold a FuzzyAndOp
    if (std::holds_alternative<std::unique_ptr<FuzzyAndOp>>(lhs) &&
        std::holds_alternative<std::unique_ptr<FuzzyAndOp>>(rhs)) {
        const auto& lhsAndOp = std::get<std::unique_ptr<FuzzyAndOp>>(lhs);
        const auto& rhsAndOp = std::get<std::unique_ptr<FuzzyAndOp>>(rhs);
        return lhsAndOp->lhs == rhsAndOp->lhs && lhsAndOp->rhs == rhsAndOp->rhs;
    }

    // Check if both variants hold a FuzzyOrOp
    if (std::holds_alternative<std::unique_ptr<FuzzyOrOp>>(lhs) &&
        std::holds_alternative<std::unique_ptr<FuzzyOrOp>>(rhs)) {
        const auto& lhsOrOp = std::get<std::unique_ptr<FuzzyOrOp>>(lhs);
        const auto& rhsOrOp = std::get<std::unique_ptr<FuzzyOrOp>>(rhs);
        return lhsOrOp->lhs == rhsOrOp->lhs && lhsOrOp->rhs == rhsOrOp->rhs;
    }

    if (std::holds_alternative<std::unique_ptr<FuzzyIfOp>>(lhs) &&
        std::holds_alternative<std::unique_ptr<FuzzyIfOp>>(rhs)) {
        const auto& lhsIfOp = std::get<std::unique_ptr<FuzzyIfOp>>(lhs);
        const auto& rhsIfOp = std::get<std::unique_ptr<FuzzyIfOp>>(rhs);
        return lhsIfOp->lhs == rhsIfOp->lhs && lhsIfOp->rhs == rhsIfOp->rhs;
    }

    return false;
}

bool operator==(const FuzzyOp&& lhs, const FuzzyOp& rhs) {
    if (lhs.index() != rhs.index()) {
        return false; // the variants hold different types
    }

    // Check if both variants hold a FuzzySetId
    if (std::holds_alternative<FuzzySetId>(lhs) &&
        std::holds_alternative<FuzzySetId>(rhs)) {
        return std::get<FuzzySetId>(lhs) == std::get<FuzzySetId>(rhs);
    }

    // Check if both variants hold a FuzzyNotOp
    if (std::holds_alternative<std::unique_ptr<FuzzyNotOp>>(lhs) &&
        std::holds_alternative<std::unique_ptr<FuzzyNotOp>>(rhs)) {
        const auto& lhsNotOp = std::get<std::unique_ptr<FuzzyNotOp>>(lhs);
        const auto& rhsNotOp = std::get<std::unique_ptr<FuzzyNotOp>>(rhs);
        return lhsNotOp->lhs == rhsNotOp->lhs;
    }

    // Check if both variants hold a FuzzyAndOp
    if (std::holds_alternative<std::unique_ptr<FuzzyAndOp>>(lhs) &&
        std::holds_alternative<std::unique_ptr<FuzzyAndOp>>(rhs)) {
        const auto& lhsAndOp = std::get<std::unique_ptr<FuzzyAndOp>>(lhs);
        const auto& rhsAndOp = std::get<std::unique_ptr<FuzzyAndOp>>(rhs);
        return lhsAndOp->lhs == rhsAndOp->lhs && lhsAndOp->rhs == rhsAndOp->rhs;
    }

    // Check if both variants hold a FuzzyOrOp
    if (std::holds_alternative<std::unique_ptr<FuzzyOrOp>>(lhs) &&
        std::holds_alternative<std::unique_ptr<FuzzyOrOp>>(rhs)) {
        const auto& lhsOrOp = std::get<std::unique_ptr<FuzzyOrOp>>(lhs);
        const auto& rhsOrOp = std::get<std::unique_ptr<FuzzyOrOp>>(rhs);
        return lhsOrOp->lhs == rhsOrOp->lhs && lhsOrOp->rhs == rhsOrOp->rhs;
    }

    if (std::holds_alternative<std::unique_ptr<FuzzyIfOp>>(lhs) &&
        std::holds_alternative<std::unique_ptr<FuzzyIfOp>>(rhs)) {
        const auto& lhsIfOp = std::get<std::unique_ptr<FuzzyIfOp>>(lhs);
        const auto& rhsIfOp = std::get<std::unique_ptr<FuzzyIfOp>>(rhs);
        return lhsIfOp->lhs == rhsIfOp->lhs && lhsIfOp->rhs == rhsIfOp->rhs;
    }

    return false;
}

bool operator==(const FuzzyOp& lhs, const FuzzyOp&& rhs) {
    if (lhs.index() != rhs.index()) {
        return false; // the variants hold different types
    }

    // Check if both variants hold a FuzzySetId
    if (std::holds_alternative<FuzzySetId>(lhs) &&
        std::holds_alternative<FuzzySetId>(rhs)) {
        return std::get<FuzzySetId>(lhs) == std::get<FuzzySetId>(rhs);
    }

    // Check if both variants hold a FuzzyNotOp
    if (std::holds_alternative<std::unique_ptr<FuzzyNotOp>>(lhs) &&
        std::holds_alternative<std::unique_ptr<FuzzyNotOp>>(rhs)) {
        const auto& lhsNotOp = std::get<std::unique_ptr<FuzzyNotOp>>(lhs);
        const auto& rhsNotOp = std::get<std::unique_ptr<FuzzyNotOp>>(rhs);
        return lhsNotOp->lhs == rhsNotOp->lhs;
    }

    // Check if both variants hold a FuzzyAndOp
    if (std::holds_alternative<std::unique_ptr<FuzzyAndOp>>(lhs) &&
        std::holds_alternative<std::unique_ptr<FuzzyAndOp>>(rhs)) {
        const auto& lhsAndOp = std::get<std::unique_ptr<FuzzyAndOp>>(lhs);
        const auto& rhsAndOp = std::get<std::unique_ptr<FuzzyAndOp>>(rhs);
        return lhsAndOp->lhs == rhsAndOp->lhs && lhsAndOp->rhs == rhsAndOp->rhs;
    }

    // Check if both variants hold a FuzzyOrOp
    if (std::holds_alternative<std::unique_ptr<FuzzyOrOp>>(lhs) &&
        std::holds_alternative<std::unique_ptr<FuzzyOrOp>>(rhs)) {
        const auto& lhsOrOp = std::get<std::unique_ptr<FuzzyOrOp>>(lhs);
        const auto& rhsOrOp = std::get<std::unique_ptr<FuzzyOrOp>>(rhs);
        return lhsOrOp->lhs == rhsOrOp->lhs && lhsOrOp->rhs == rhsOrOp->rhs;
    }

    if (std::holds_alternative<std::unique_ptr<FuzzyIfOp>>(lhs) &&
        std::holds_alternative<std::unique_ptr<FuzzyIfOp>>(rhs)) {
        const auto& lhsIfOp = std::get<std::unique_ptr<FuzzyIfOp>>(lhs);
        const auto& rhsIfOp = std::get<std::unique_ptr<FuzzyIfOp>>(rhs);
        return lhsIfOp->lhs == rhsIfOp->lhs && lhsIfOp->rhs == rhsIfOp->rhs;
    }

    return false;
}

bool operator==(const FuzzyOp&& lhs, const FuzzyOp&& rhs) {
    if (lhs.index() != rhs.index()) {
        return false; // the variants hold different types
    }

    // Check if both variants hold a FuzzySetId
    if (std::holds_alternative<FuzzySetId>(lhs) &&
        std::holds_alternative<FuzzySetId>(rhs)) {
        return std::get<FuzzySetId>(lhs) == std::get<FuzzySetId>(rhs);
    }

    // Check if both variants hold a FuzzyNotOp
    if (std::holds_alternative<std::unique_ptr<FuzzyNotOp>>(lhs) &&
        std::holds_alternative<std::unique_ptr<FuzzyNotOp>>(rhs)) {
        const auto& lhsNotOp = std::get<std::unique_ptr<FuzzyNotOp>>(lhs);
        const auto& rhsNotOp = std::get<std::unique_ptr<FuzzyNotOp>>(rhs);
        return lhsNotOp->lhs == rhsNotOp->lhs;
    }

    // Check if both variants hold a FuzzyAndOp
    if (std::holds_alternative<std::unique_ptr<FuzzyAndOp>>(lhs) &&
        std::holds_alternative<std::unique_ptr<FuzzyAndOp>>(rhs)) {
        const auto& lhsAndOp = std::get<std::unique_ptr<FuzzyAndOp>>(lhs);
        const auto& rhsAndOp = std::get<std::unique_ptr<FuzzyAndOp>>(rhs);
        return lhsAndOp->lhs == rhsAndOp->lhs && lhsAndOp->rhs == rhsAndOp->rhs;
    }

    // Check if both variants hold a FuzzyOrOp
    if (std::holds_alternative<std::unique_ptr<FuzzyOrOp>>(lhs) &&
        std::holds_alternative<std::unique_ptr<FuzzyOrOp>>(rhs)) {
        const auto& lhsOrOp = std::get<std::unique_ptr<FuzzyOrOp>>(lhs);
        const auto& rhsOrOp = std::get<std::unique_ptr<FuzzyOrOp>>(rhs);
        return lhsOrOp->lhs == rhsOrOp->lhs && lhsOrOp->rhs == rhsOrOp->rhs;
    }

    if (std::holds_alternative<std::unique_ptr<FuzzyIfOp>>(lhs) &&
        std::holds_alternative<std::unique_ptr<FuzzyIfOp>>(rhs)) {
        const auto& lhsIfOp = std::get<std::unique_ptr<FuzzyIfOp>>(lhs);
        const auto& rhsIfOp = std::get<std::unique_ptr<FuzzyIfOp>>(rhs);
        return lhsIfOp->lhs == rhsIfOp->lhs && lhsIfOp->rhs == rhsIfOp->rhs;
    }

    return false;
}

float defuzzify(const FuzzySystem& fs, const FuzzyOp& fuzzySet) {
    if (auto fuzzySetId = std::get_if<FuzzySetId>(&fuzzySet)) {
        float membershipValue = fs.fuzzyMembershipValues[*fuzzySetId];
        MembershipFunction membershipFunction = fs.fuzzySets[*fuzzySetId];
        float center                          = 0;
        if (auto triangular =
                std::get_if<FuzzySetTriangular>(&membershipFunction)) {
            center = triangular->peak;
        } else if (auto trapezoidal =
                       std::get_if<FuzzySetTrapezoidal>(&membershipFunction)) {
            center = (trapezoidal->left_peak + trapezoidal->right_peak) / 2;
        } else {
            auto gaussian = std::get_if<FuzzySetGaussian>(&membershipFunction);
            center        = gaussian->mean;
        }

        return center;
    }
    return 0.;
}

// Defuzzify using the weighted average method
float defuzzify(const FuzzySystem& fs,
                const std::vector<FuzzyOp>& fuzzySetVector) {
    float numerator   = 0.;
    float denominator = 0.;

    for (const auto& fuzzySet : fuzzySetVector) {
        if (auto fuzzySetId = std::get_if<FuzzySetId>(&fuzzySet)) {
            float membershipValue   = fs.fuzzyMembershipValues[*fuzzySetId];
            auto membershipFunction = fs.fuzzySets[*fuzzySetId];
            float center            = 0;
            if (auto triangular =
                    std::get_if<FuzzySetTriangular>(&membershipFunction)) {
                center = triangular->peak;
            } else if (auto trapezoidal = std::get_if<FuzzySetTrapezoidal>(
                           &membershipFunction)) {
                center =
                    (trapezoidal->left_peak + trapezoidal->right_peak) / 2.;
            } else {
                auto gaussian =
                    std::get_if<FuzzySetGaussian>(&membershipFunction);
                center = gaussian->mean;
            }

            numerator += membershipValue * center;
            denominator += membershipValue;
        }
    }

    // Check if the denominator is zero
    if (denominator == 0.) {
        return 0.;
    } else {
        return numerator / denominator;
    }
}

void printFuzzyOp(const FuzzyOp& op, bool root = true) {
    std::visit(
        [](const auto& op) {
            using T = std::decay_t<decltype(op)>;
            if constexpr (std::is_same_v<T, FuzzySetId>) {
                // For a FuzzySetId, simply print the id
                std::cout << " " << op << " ";
            } else if constexpr (std::is_same_v<T,
                                                std::unique_ptr<FuzzyNotOp>>) {
                // For a FuzzyNotOp, print "not" followed by the operand
                // recursively
                std::cout << "~(";
                printFuzzyOp(op->lhs, false);
                std::cout << ")";
            } else if constexpr (std::is_same_v<T,
                                                std::unique_ptr<FuzzyAndOp>>) {
                // For a FuzzyAndOp, print the operands recursively surrounded
                // by "and"
                std::cout << "(";
                printFuzzyOp(op->lhs, false);
                std::cout << ") && (";
                printFuzzyOp(op->rhs, false);
                std::cout << ")";
            } else if constexpr (std::is_same_v<T,
                                                std::unique_ptr<FuzzyOrOp>>) {
                // For a FuzzyOrOp, print the operands recursively surrounded by
                // "or"
                std::cout << "(";
                printFuzzyOp(op->lhs, false);
                std::cout << ") || (";
                printFuzzyOp(op->rhs, false);
                std::cout << ")";
            } else if constexpr (std::is_same_v<T,
                                                std::unique_ptr<FuzzyIfOp>>) {
                // For a FuzzyIfOp, print the operands recursively surrounded by
                // "if" and "then"
                std::cout << "(";
                printFuzzyOp(op->lhs, false);
                std::cout << ") -> (";
                printFuzzyOp(op->rhs, false);
                std::cout << ")";
            }
        },
        op);
    if (root) {
        std::cout << std::endl;
    }
}

// Define the 'not' operator as a unary operator on FuzzyOpNode
FuzzyOp operator!(FuzzyOp& lhs) {
    // Create a FuzzyNotOp with the given FuzzyOp as the left-hand side
    auto notOp = std::make_unique<FuzzyNotOp>();
    notOp->lhs = std::move(lhs);

    // Return the FuzzyNotOp as a FuzzyOp
    return notOp;
}

// Define the 'not' operator as a unary operator on FuzzyOpNode
FuzzyOp operator!(FuzzyOp&& lhs) {
    // Create a FuzzyNotOp with the given FuzzyOp as the left-hand side
    auto notOp = std::make_unique<FuzzyNotOp>();
    notOp->lhs = std::move(lhs);

    // Return the FuzzyNotOp as a FuzzyOp
    return notOp;
}

// Define the 'and' operator as a binary operator on FuzzyOp
FuzzyOp operator&&(FuzzyOp& lhs, FuzzyOp& rhs) {
    // Create a FuzzyAndOp with the given FuzzyOp as the left-hand and
    // right-hand sides
    auto andOp = std::make_unique<FuzzyAndOp>();
    andOp->lhs = std::move(lhs);
    andOp->rhs = std::move(rhs);

    // Return the FuzzyAndOp as a FuzzyOp
    return andOp;
}

// Define the 'and' operator as a binary operator on FuzzyOp
FuzzyOp operator&&(FuzzyOp&& lhs, FuzzyOp& rhs) {
    // Create a FuzzyAndOp with the given FuzzyOp as the left-hand and
    // right-hand sides
    auto andOp = std::make_unique<FuzzyAndOp>();
    andOp->lhs = std::move(lhs);
    andOp->rhs = std::move(rhs);

    // Return the FuzzyAndOp as a FuzzyOp
    return andOp;
}

// Define the 'and' operator as a binary operator on FuzzyOp
FuzzyOp operator&&(FuzzyOp& lhs, FuzzyOp&& rhs) {
    // Create a FuzzyAndOp with the given FuzzyOp as the left-hand and
    // right-hand sides
    auto andOp = std::make_unique<FuzzyAndOp>();
    andOp->lhs = std::move(lhs);
    andOp->rhs = std::move(rhs);

    // Return the FuzzyAndOp as a FuzzyOp
    return andOp;
}

// Define the 'and' operator as a binary operator on FuzzyOp
FuzzyOp operator&&(FuzzyOp&& lhs, FuzzyOp&& rhs) {
    // Create a FuzzyAndOp with the given FuzzyOp as the left-hand and
    // right-hand sides
    auto andOp = std::make_unique<FuzzyAndOp>();
    andOp->lhs = std::move(lhs);
    andOp->rhs = std::move(rhs);

    // Return the FuzzyAndOp as a FuzzyOp
    return andOp;
}

// Define the 'or' operator as a binary operator on FuzzyOp
FuzzyOp operator||(FuzzyOp& lhs, FuzzyOp& rhs) {
    // Create a FuzzyOrOp with the given FuzzyOp as the left-hand and right-hand
    // sides
    auto orOp = std::make_unique<FuzzyOrOp>();

    orOp->lhs = std::move(lhs);
    orOp->rhs = std::move(rhs);
    // Return the FuzzyOrOp as a FuzzyOp
    return orOp;
}

// Define the 'or' operator as a binary operator on FuzzyOp
FuzzyOp operator||(FuzzyOp&& lhs, FuzzyOp& rhs) {
    // Create a FuzzyOrOp with the given FuzzyOp as the left-hand and right-hand
    // sides
    auto orOp = std::make_unique<FuzzyOrOp>();

    orOp->lhs = std::move(lhs);
    orOp->rhs = std::move(rhs);
    // Return the FuzzyOrOp as a FuzzyOp
    return orOp;
}

// Define the 'or' operator as a binary operator on FuzzyOp
FuzzyOp operator||(FuzzyOp& lhs, FuzzyOp&& rhs) {
    // Create a FuzzyOrOp with the given FuzzyOp as the left-hand and right-hand
    // sides
    auto orOp = std::make_unique<FuzzyOrOp>();

    orOp->lhs = std::move(lhs);
    orOp->rhs = std::move(rhs);
    // Return the FuzzyOrOp as a FuzzyOp
    return orOp;
}

// Define the 'or' operator as a binary operator on FuzzyOp
FuzzyOp operator||(FuzzyOp&& lhs, FuzzyOp&& rhs) {
    // Create a FuzzyOrOp with the given FuzzyOp as the left-hand and right-hand
    // sides
    auto orOp = std::make_unique<FuzzyOrOp>();

    orOp->lhs = std::move(lhs);
    orOp->rhs = std::move(rhs);
    // Return the FuzzyOrOp as a FuzzyOp
    return orOp;
}

// Define the 'if' operator as a binary operator on FuzzyOp
FuzzyOp operator>>(FuzzyOp& lhs, FuzzyOp& rhs) {
    // Create a FuzzyIfOp with the given FuzzyOp as the left-hand and right-hand
    // sides
    auto ifOp = std::make_unique<FuzzyIfOp>();
    ifOp->lhs = std::move(lhs);
    ifOp->rhs = std::move(rhs);

    // Return the FuzzyIfOp as a FuzzyOp
    return ifOp;
}

// Define the 'if' operator as a binary operator on FuzzyOp
FuzzyOp operator>>(FuzzyOp&& lhs, FuzzyOp& rhs) {
    // Create a FuzzyIfOp with the given FuzzyOp as the left-hand and right-hand
    // sides
    auto ifOp = std::make_unique<FuzzyIfOp>();
    ifOp->lhs = std::move(lhs);
    ifOp->rhs = std::move(rhs);

    // Return the FuzzyIfOp as a FuzzyOp
    return ifOp;
}

// Define the 'if' operator as a binary operator on FuzzyOp
FuzzyOp operator>>(FuzzyOp& lhs, FuzzyOp&& rhs) {
    // Create a FuzzyIfOp with the given FuzzyOp as the left-hand and right-hand
    // sides
    auto ifOp = std::make_unique<FuzzyIfOp>();
    ifOp->lhs = std::move(lhs);
    ifOp->rhs = std::move(rhs);

    // Return the FuzzyIfOp as a FuzzyOp
    return ifOp;
}

// Define the 'if' operator as a binary operator on FuzzyOp
FuzzyOp operator>>(FuzzyOp&& lhs, FuzzyOp&& rhs) {
    // Create a FuzzyIfOp with the given FuzzyOp as the left-hand and right-hand
    // sides
    auto ifOp = std::make_unique<FuzzyIfOp>();
    ifOp->lhs = std::move(lhs);
    ifOp->rhs = std::move(rhs);

    // Return the FuzzyIfOp as a FuzzyOp
    return ifOp;
}
