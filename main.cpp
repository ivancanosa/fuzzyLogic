#include "FuzzyLogic.hpp"


/* This example simulates the AI of a knight in a videogame. For each input/output,
 * it creates 3 sets (low, medium, high) to represent it.
 * Input Sets:
 * 		-Health: Actual health of the knight in the range [0-1]
 * 		-Melee damage: Melee damage of the knight
 * 		-Range damage: Long distance damage of the knight
 * Output Sets:
 * 		- Retreat: Action for the knight to retreat from battle
 * 		- Heal: Action for the knight to heal himself.
 * 		- Attack melee: Action for the knight to attack at melee range
 * 		- Attack range: Action for the knight to attack at long distance 
 */

void exampleAgentKnightInVideogame() {
    // Create a fuzzy system
    FuzzySystem fs;
    // Define 3 membership trapezoidal functions in the range [0-1].
	// They will be low, medium, and high respectively by their position in the return vector
    auto mfVector3 = createTrapezoidalFuzzySets(3, 1.);

    // Create the input fuzzy sets
    auto healthSet      = addFuzzySet(fs, mfVector3);
    auto meleDamageSet  = addFuzzySet(fs, mfVector3);
    auto rangeDamageSet = addFuzzySet(fs, mfVector3);

    // Create the output fuzzy sets for each action
    auto retreatSet     = addFuzzySet(fs, mfVector3);
    auto healSet        = addFuzzySet(fs, mfVector3);
    auto attackRangeSet = addFuzzySet(fs, mfVector3);
    auto attackMeleeSet = addFuzzySet(fs, mfVector3);

    // Definition of the rules. You must be aware that every fuzzy set must be
	// present at least once in some fuzzy rule.
	fuzzySystemAddRule(fs, healthSet[0] >> healSet[2]);
	fuzzySystemAddRule(fs, healthSet[1] >> healSet[1]);
	fuzzySystemAddRule(fs, healthSet[2] >> healSet[0]);

	fuzzySystemAddRule(fs, (healthSet[2] || meleDamageSet[2]) >> attackMeleeSet[2]);
	fuzzySystemAddRule(fs, (healthSet[2] || meleDamageSet[1]) >> attackMeleeSet[1]);
	fuzzySystemAddRule(fs, (healthSet[2] && meleDamageSet[0]) >> attackMeleeSet[0]);

	fuzzySystemAddRule(fs, (rangeDamageSet[0] || meleDamageSet[0]) >> retreatSet[2]);
	fuzzySystemAddRule(fs, (rangeDamageSet[0] || meleDamageSet[1]) >> retreatSet[1]);
	fuzzySystemAddRule(fs, (rangeDamageSet[1] && meleDamageSet[2]) >> retreatSet[0]);

	fuzzySystemAddRule(fs, (rangeDamageSet[2] || meleDamageSet[0]) >> attackRangeSet[0]);
	fuzzySystemAddRule(fs, (rangeDamageSet[2] || meleDamageSet[1]) >> attackRangeSet[1]);
	fuzzySystemAddRule(fs, (rangeDamageSet[2] && meleDamageSet[2]) >> attackRangeSet[2]);

    // Sort the rules once they are all inserted
    fuzzyRulesTopologicalSort(fs);

    // Input the sensor values
    float health      = 0.4; // Knight's health
    float meleDamage  = 0.7; // Knight's melee damage
    float rangeDamage = 0.3; // Knight's range damage

	// Insert the values in the input fuzzy sets. 
    setFuzzyInput(fs, healthSet, health);
    setFuzzyInput(fs, meleDamageSet, meleDamage);
    setFuzzyInput(fs, rangeDamageSet, rangeDamage);

	// Compute the membership values for the output sets
    computeFuzzyRules(fs);

    // Defuzzify the output fuzzy sets to get the crips action values
    float retreatActionValue     = defuzzify(fs, retreatSet);
    float healActionValue        = defuzzify(fs, healSet);
    float attackRangeActionValue = defuzzify(fs, attackRangeSet);
    float attackMeleeActionValue = defuzzify(fs, attackMeleeSet);
    std::cout << "retreat value: " << retreatActionValue  << std::endl;
    std::cout << "heal value: " << healActionValue << std::endl;
    std::cout << "melee value: " << attackMeleeActionValue << std::endl;
    std::cout << "range value: " << attackRangeActionValue << std::endl;

    // Select the best action based on the action values
    if (retreatActionValue > healActionValue &&
        retreatActionValue > attackRangeActionValue &&
        retreatActionValue > attackMeleeActionValue) {
        std::cout << "The best action is to retreat." << std::endl;
    } else if (healActionValue > retreatActionValue &&
               healActionValue > attackRangeActionValue &&
               healActionValue > attackMeleeActionValue) {
        std::cout << "The best action is to heal." << std::endl;
    } else if (attackRangeActionValue > retreatActionValue &&
               attackRangeActionValue > healActionValue &&
               attackRangeActionValue > attackMeleeActionValue) {
        std::cout << "The best action is to attack range." << std::endl;
    } else {
        std::cout << "The best action is to attack melee." << std::endl;
    }
}


int main(){
    exampleAgentKnightInVideogame();
	return 0;
}
