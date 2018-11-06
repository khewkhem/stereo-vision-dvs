process_csv: process_csv.cpp
	g++ process_csv.cpp -o process_csv -DARMA_DONT_USE_WRAPPER -lopenblas -llapack -std=c++14

debug: process_csv.cpp
	g++ -g process_csv.cpp -o debug -DARMA_DONT_USE_WRAPPER -lopenblas -llapack -std=c++14
