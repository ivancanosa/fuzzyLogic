#include "FuzzyLogic.hpp"
#include <iostream>

bool floatEquals(float a, float b, float epsilon = 0.001) {
    return std::abs(a - b) < epsilon;
}

void testComputeMembership() {
    // Test triangular membership function
    FuzzySetTriangular triangularMF = {0, 5, 10};
    assert(std::abs(computeMembership(0, triangularMF) - 0) < 0.001);
    assert(std::abs(computeMembership(5, triangularMF) - 1) < 0.001);
    assert(std::abs(computeMembership(10, triangularMF) - 0) < 0.001);
    assert(std::abs(computeMembership(2.5, triangularMF) - 0.5) < 0.001);
    assert(std::abs(computeMembership(7.5, triangularMF) - 0.5) < 0.001);

    // Test trapezoidal membership function
    FuzzySetTrapezoidal trapezoidalMF = {0, 1, 9, 10};
    assert(std::abs(computeMembership(0, trapezoidalMF) - 0) < 0.001);
    assert(std::abs(computeMembership(1, trapezoidalMF) - 1) < 0.001);
    assert(std::abs(computeMembership(5, trapezoidalMF) - 1) < 0.001);
    assert(std::abs(computeMembership(9, trapezoidalMF) - 1) < 0.001);
    assert(std::abs(computeMembership(10, trapezoidalMF) - 0) < 0.001);
    assert(std::abs(computeMembership(0.5, trapezoidalMF) - 0.5) < 0.001);
    assert(std::abs(computeMembership(9.5, trapezoidalMF) - 0.5) < 0.001);

    // Test Gaussian membership function
    FuzzySetGaussian gaussianMF = {5, 2};
    assert(std::abs(computeMembership(0, gaussianMF) - 0.008764) < 0.001);
}

void testCreateTriangularFuzzySets() {
    auto sets = createTriangularFuzzySets(3, 0.1);

    assert(sets.size() == 3);
    assert(floatEquals(std::get<FuzzySetTriangular>(sets[0]).left_bound, -0.3));
    assert(floatEquals(std::get<FuzzySetTriangular>(sets[0]).peak, 0));
    assert(floatEquals(std::get<FuzzySetTriangular>(sets[0]).right_bound, 0.3));
    assert(floatEquals(std::get<FuzzySetTriangular>(sets[1]).left_bound, 0.2));
    assert(floatEquals(std::get<FuzzySetTriangular>(sets[1]).peak, 0.5));
    assert(floatEquals(std::get<FuzzySetTriangular>(sets[1]).right_bound, 0.8));
    assert(floatEquals(std::get<FuzzySetTriangular>(sets[2]).left_bound, 0.7));
    assert(floatEquals(std::get<FuzzySetTriangular>(sets[2]).peak, 1.));
    assert(floatEquals(std::get<FuzzySetTriangular>(sets[2]).right_bound, 1.3));
}

void testCreateTrapezoidalFuzzySets() {
    std::vector<MembershipFunction> sets = createTrapezoidalFuzzySets(3, 0.1);

    assert(sets.size() == 3);
    assert(
        floatEquals(std::get<FuzzySetTrapezoidal>(sets[0]).left_bound, -0.266));
    assert(
        floatEquals(std::get<FuzzySetTrapezoidal>(sets[0]).left_peak, -0.233));
    assert(
        floatEquals(std::get<FuzzySetTrapezoidal>(sets[0]).right_peak, 0.233));
    assert(floatEquals(std::get<FuzzySetTrapezoidal>(sets[0]).right_bound,
                       0.26666));
    assert(
        floatEquals(std::get<FuzzySetTrapezoidal>(sets[1]).left_bound, 0.2333));
    assert(
        floatEquals(std::get<FuzzySetTrapezoidal>(sets[1]).left_peak, 0.2666));
    assert(floatEquals(std::get<FuzzySetTrapezoidal>(sets[1]).right_peak,
                       0.73333));
    assert(floatEquals(std::get<FuzzySetTrapezoidal>(sets[1]).right_bound,
                       0.76666));
    assert(
        floatEquals(std::get<FuzzySetTrapezoidal>(sets[2]).left_bound, 0.7333));
    assert(
        floatEquals(std::get<FuzzySetTrapezoidal>(sets[2]).left_peak, 0.7666));
    assert(
        floatEquals(std::get<FuzzySetTrapezoidal>(sets[2]).right_peak, 1.2333));
    assert(floatEquals(std::get<FuzzySetTrapezoidal>(sets[2]).right_bound,
                       1.2666));
}

void testFuzzySystemAddRule() {
    // create a fuzzy system
    FuzzySystem fs;

    // create a fuzzy set
    MembershipFunction mf = FuzzySetTriangular{0.0f, 0.5f, 1.0f};
    FuzzyOp op            = addFuzzySet(fs, mf);

    // add the rule to the fuzzy system
    fuzzySystemAddRule(fs, op >> op);

    // check if the rule was added to the system
    assert(fs.fuzzyRules.size() == 1);
    assert(std::holds_alternative<FuzzySetId>(fs.fuzzyRules[0].lhs));
    assert(std::holds_alternative<FuzzySetId>(fs.fuzzyRules[0].rhs));
}

void testComputeFuzzyRules() {
    MembershipFunction mf0 = FuzzySetTriangular{0., 0.5, 1.};
    MembershipFunction mf1 = FuzzySetTriangular{0.25, 0.5, 1.};
    MembershipFunction mf2 = FuzzySetTriangular{0.25, 0.5, 1.};
    // Not operator test
    {
        FuzzySystem fs;
        auto input0 = addFuzzySet(fs, mf0);
        auto output = addFuzzySet(fs, mf2);

        fuzzySystemAddRule(fs, (!input0 >> output));
        setFuzzyInput(fs, input0, 1.0);
        computeFuzzyRules(fs);
        float result = fs.fuzzyMembershipValues[std::get<FuzzySetId>(output)];
        assert(floatEquals(result, 1. - computeMembership(1., mf0)));
    }
    // And operator test
    {
        FuzzySystem fs;
        auto input0 = addFuzzySet(fs, mf0);
        auto input1 = addFuzzySet(fs, mf1);
        auto output = addFuzzySet(fs, mf2);

        fuzzySystemAddRule(fs, ((input0 && input1) >> output));
        setFuzzyInput(fs, input0, 1.0);
        setFuzzyInput(fs, input1, 0.5);
        computeFuzzyRules(fs);
        float result = fs.fuzzyMembershipValues[std::get<FuzzySetId>(output)];
        assert(floatEquals(result, computeMembership(1., mf0)));
    }
    // Or operator test
    {
        FuzzySystem fs;
        auto input0 = addFuzzySet(fs, mf0);
        auto input1 = addFuzzySet(fs, mf1);
        auto output = addFuzzySet(fs, mf2);

        fuzzySystemAddRule(fs, ((input0 || input1) >> output));
        setFuzzyInput(fs, input0, 1.0);
        setFuzzyInput(fs, input1, 0.5);
        computeFuzzyRules(fs);
        float result = fs.fuzzyMembershipValues[std::get<FuzzySetId>(output)];
        assert(floatEquals(result, computeMembership(0.5, mf1)));
    }
}

void testTopologicalSort() {
    // Create a FuzzySystem with 2 fuzzy sets
    MembershipFunction mf0 = FuzzySetTriangular{1.0f, 2.0f, 3.0f};
    MembershipFunction mf1 = FuzzySetTriangular{1.0f, 2.0f, 3.0f};
    MembershipFunction mf2 = FuzzySetTriangular{1.0f, 2.0f, 3.0f};
    MembershipFunction mf3 = FuzzySetTriangular{1.0f, 2.0f, 3.0f};
    MembershipFunction mf4 = FuzzySetTriangular{1.0f, 2.0f, 3.0f};

    {
        // Test with already sorted rules
        FuzzySystem fs;
        auto set0 = addFuzzySet(fs, mf0);
        auto set1 = addFuzzySet(fs, mf1);
        auto set2 = addFuzzySet(fs, mf2);

        fuzzySystemAddRule(fs, set0 >> set1);
        fuzzySystemAddRule(fs, set1 >> set2);

        fuzzyRulesTopologicalSort(fs);
        const auto& rule0Head = fs.fuzzyRules[0].rhs;
        const auto& rule1Head = fs.fuzzyRules[1].rhs;
        assert(set1 == rule0Head);
        assert(set2 == rule1Head);
    }
    {
        // Test with simple transitive rules
        FuzzySystem fs;
        auto set0 = addFuzzySet(fs, mf0);
        auto set1 = addFuzzySet(fs, mf1);
        auto set2 = addFuzzySet(fs, mf2);
        auto set3 = addFuzzySet(fs, mf3);
        auto set4 = addFuzzySet(fs, mf4);

        fuzzySystemAddRule(fs, set2 >> set3);
        fuzzySystemAddRule(fs, set1 >> set2);
        fuzzySystemAddRule(fs, set3 >> set4);
        fuzzySystemAddRule(fs, set0 >> set1);

        fuzzyRulesTopologicalSort(fs);

        assert(set1 == fs.fuzzyRules[0].rhs);
        assert(set2 == fs.fuzzyRules[1].rhs);
        assert(set3 == fs.fuzzyRules[2].rhs);
        assert(set4 == fs.fuzzyRules[3].rhs);
    }
    {
        // Test with complex expressions in the left side of the rules
        FuzzySystem fs;
        auto set0 = addFuzzySet(fs, mf0);
        auto set1 = addFuzzySet(fs, mf1);
        auto set2 = addFuzzySet(fs, mf2);
        auto set3 = addFuzzySet(fs, mf3);
        auto set4 = addFuzzySet(fs, mf4);

        fuzzySystemAddRule(fs, (set1 || set2) >> set3);
        fuzzySystemAddRule(fs, set0 >> set1);
        fuzzySystemAddRule(fs, (set1 || set2 || set3) >> set4);
        fuzzySystemAddRule(fs, set0 >> set2);
        fuzzySystemAddRule(fs, (set1 && set2) >> set3);

        fuzzyRulesTopologicalSort(fs);

        assert((set1 == fs.fuzzyRules[0].rhs) ||
               (set2 == fs.fuzzyRules[0].rhs));
        assert((set1 == fs.fuzzyRules[1].rhs) ||
               (set2 == fs.fuzzyRules[1].rhs));
        assert(set3 == fs.fuzzyRules[2].rhs);
        assert(set3 == fs.fuzzyRules[3].rhs);
        assert(set4 == fs.fuzzyRules[4].rhs);
    }
}

void testDefuzzify() {
    // Set up fuzzy system
    FuzzySystem fs;
    auto mfVector       = createTriangularFuzzySets(3, 0.5);
    auto fuzzySetVector = addFuzzySet(fs, mfVector);
    float defuzzifiedValue;

    fs.fuzzyMembershipValues[0] = 0.;
    fs.fuzzyMembershipValues[1] = 0.3;
    fs.fuzzyMembershipValues[2] = 1.;
    defuzzifiedValue            = defuzzify(fs, fuzzySetVector);
    assert(defuzzifiedValue >= 0.8);

    fs.fuzzyMembershipValues[0] = 0.8;
    fs.fuzzyMembershipValues[1] = 0.3;
    fs.fuzzyMembershipValues[2] = 0.;
    defuzzifiedValue            = defuzzify(fs, fuzzySetVector);
    assert(defuzzifiedValue <= 0.3);

    fs.fuzzyMembershipValues[0] = 0.2;
    fs.fuzzyMembershipValues[1] = 0.8;
    fs.fuzzyMembershipValues[2] = 0.2;
    defuzzifiedValue            = defuzzify(fs, fuzzySetVector);
    assert(defuzzifiedValue >= 0.4);
    assert(defuzzifiedValue <= 0.6);
}

//******************* MAIN *********************

int main() {
    testComputeMembership();
    testCreateTriangularFuzzySets();
    testCreateTrapezoidalFuzzySets();
    testFuzzySystemAddRule();
    testComputeFuzzyRules();
    testTopologicalSort();
    testDefuzzify();

    std::cout << "All tests passed!" << std::endl;

    return 0;
}
