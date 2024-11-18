add_executable(FuzzyLogicTest test.cpp)
target_link_libraries(FuzzyLogicTest fuzzyLogic)
add_test(NAME FuzzyLogicTest COMMAND FuzzyLogicTest)
