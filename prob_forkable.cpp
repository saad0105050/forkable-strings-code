#include "ReachAndMargin.h"
#include <iostream>

void prob_forkable(const int N, const int R, double advFrac) {
	auto advFracArr = { advFrac };
	ReachAndMargin* rm = new ReachAndMargin(N, R, advFracArr);
	for (int t = 1; t <= N; t++) {
		rm->evolve();
	}
	std::cout << "\tPr[forkable] <= " << rm->forkableProbability() << std::endl;
	DELETE(rm);

}


void main(){
	int N, R;	
	double eps;

	std::cout <<
		"N is the length of execution (integer, at least 1)" << std::endl <<
		"R is the maximum reach of the prefix (between 0 and " << ReachAndMargin::RMAX << " inclusive)" << std::endl <<
		"eps is the bias (at least 0 but less than 1)" << std::endl << std::endl;

	std::cout << "Enter N, R, and eps (separate by space): "; 
	std::cin >> N;
	std::cin >> R;
	std::cin >> eps;
	if (N < 1) {
		std::cerr << "N must be at least 1";
		exit(-1);
	}
	if (R < 0 || R > ReachAndMargin::RMAX) {
		std::cerr << "R must be between 0 and " << ReachAndMargin::RMAX << "  (inclusive)";
		exit(-1);
	}
	if (eps < 0 || eps > 1) {
		std::cerr << "eps must be between 0 and 1";
		exit(-1);
	}
	std::cout
		<< "Received N: " << N
		<< " R: " << R
		<< " eps: " << eps
		<< std::endl;
	
	double advFrac = (1 - eps) / 2.0;
	prob_forkable(N, R, advFrac);

	system("pause");
}