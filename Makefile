all : main.cpp mesh.cpp
	g++ main.cpp mesh.cpp -I/usr/include/eigen3 -std=c++14 -o CalcPtArea
clean : 
	rm CalcPtArea
