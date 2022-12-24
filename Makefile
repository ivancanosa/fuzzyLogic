main: main.cpp FuzzyLogic.cpp FuzzyLogic.hpp
	bear -- clang++ -g -stdlib=libc++ -std=c++17 main.cpp FuzzyLogic.cpp -o main 

testOutput: test.cpp FuzzyLogic.cpp FuzzyLogic.hpp
	bear -- clang++ -g -stdlib=libc++ -std=c++17 test.cpp FuzzyLogic.cpp -o testOutput 

test: testOutput
	./testOutput

run: main
	./main

clean:
	rm -rf main testOutput log
