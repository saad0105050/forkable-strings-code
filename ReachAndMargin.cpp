#include "ReachAndMargin.h"
#include <cassert>
#include <iomanip>




void ReachAndMargin::mult(double* dest, const double* one, const double* two) {
	for (auto t = 0; t < NFRACS; t++)
		dest[t] = one[t] * two[t];
}

void ReachAndMargin::add(double* dest, const double* one, const double* two) {
	for (auto t = 0; t < NFRACS; t++)
		dest[t] = one[t] + two[t];
}

// stationary probability for the event "reach == r"
double stationaryRho(const unsigned int r, const double adversarialFraction) {
	double p = 1 - adversarialFraction;
	double beta = (1 - p) / p;
	return (1 - beta) * pow(beta, r);
}
// stationary probability for the event "reach >= r"
double stationaryRhoTail(const unsigned int R, const double adversarialFraction) {
	double p = 1 - adversarialFraction;
	double beta = (1 - p) / p;
	return pow(beta, R);
}


// create instances of static const members here (i.e., in the cpp file); 
// their values are specified in the class header
const int ReachAndMargin::MAX_COLS;
const int ReachAndMargin::MAX_ROWS;
const int ReachAndMargin::MuZero;  // index of the row corresponding to mu = 0
const int ReachAndMargin::TotalCount;



ReachAndMargin::ReachAndMargin(const int N, const int R, const DoubleVector& advFractions)
	: N(N), R(R < 0 || R > RMAX ? RMAX : R),
	NFRACS(advFractions.size()),
	advFractions(advFractions)
{
	// set up special matrices
	zero = new double[NFRACS];
	prDown = new double[NFRACS];
	prUp = new double[NFRACS];
	for (int i = 0; i < NFRACS; i++) {
		zero[i] = 0;
		prUp[i] = advFractions[i];
		prDown[i] = 1 - prUp[i];
	}
	// set up temporary matrices
	tempA = new double[NFRACS];
	tempB = new double[NFRACS];
	tempC = new double[NFRACS];

	// set up the two main matrices 
	// one for previous state and the other for current state
	memset(matrixOneRaw, 0, sizeof(double) * TotalCount);	// current state
	memset(matrixTwoRaw, 0, sizeof(double) * TotalCount);	// previous state
	curMat = matrixOneRaw;
	prevMat = matrixTwoRaw;


	// time starts at zero
	time = 0;


	// set initial probabilities
	if (R == 0) {
		// consider only initial reach = 0
		// all mass at (0, 0)
		auto cell = get(curMat, 0, 0, 0);
		for (int alpha = 0; alpha < NFRACS; alpha++)
			cell[alpha] = 1;
	}
	else {
		// consider all initial reach <= R
		for (int rho = 0; rho <= R; rho++) {
			// considering reach = rho

			// relative margin of an empty string equals the reach of the prefix
			auto mu = rho;

			// grab the cell corresponding to reach = rho
			// its value should be the stationary probability of having reach == rho
			auto cell = get(curMat, 0, rho, mu);
			for (int alpha = 0; alpha < NFRACS; alpha++)
				cell[alpha] = stationaryRho(rho, advFractions[alpha]);
		}
	}

}



ReachAndMargin::~ReachAndMargin()
{
	DELETE_ARR(zero);
	DELETE_ARR(prDown);
	DELETE_ARR(prUp);
	DELETE_ARR(tempA);
	DELETE_ARR(tempB);
	DELETE_ARR(tempC);
}

void ReachAndMargin::swapMat() {
	MatrixType temp = prevMat;
	prevMat = curMat;
	curMat = temp;
}


void ReachAndMargin::calcNewValue(const int r, const int m) {
	// initialize to zero
	//double* newValue = curMat[r][MuZero + m];
	double* newValue = get(curMat, time, r, m);
	for (int a = 0; a < NFRACS; a++)
		newValue[a] = 0;

	//auto diagonallyDown = get(time - 1, r - 1, m - 1);
	//auto up = get(time - 1, r, m + 1);
	//auto down = get(time - 1, r, m - 1);
	//auto left = get(time - 1, r - 1, m);


	if (r > 0) {
		if (m > 0 || m <= -2) {
			// unrestricted
			// reach and margin move together
			auto diagonallyDown = get(prevMat, time - 1, r - 1, m - 1);
			auto diagonallyUp = get(prevMat, time - 1, r + 1, m + 1);

			mult(tempA, prUp, diagonallyDown);
			mult(tempB, prDown, diagonallyUp);
			add(newValue, tempA, tempB);
			// *newValue = prUp * diagonallyDown + prDown * diagonallyUp;
		}
		else if (m == 0) {
			// reach positive, margin zero
			// up => diagonally up, down => margin stays
			auto diagonallyDown = get(prevMat, time - 1, r - 1, m - 1);
			auto diagonallyUp = get(prevMat, time - 1, r + 1, m + 1);
			auto right = get(prevMat, time - 1, r + 1, m);

			add(tempC, diagonallyUp, right);
			mult(tempA, prDown, tempC);
			mult(tempB, prUp, diagonallyDown);
			add(newValue, tempA, tempB);
			// *newValue = prDown * (diagonallyUp + right) + prUp * diagonallyDown;
		}
		else if (m == -1) {
			// no incoming down
			auto diagonallyDown = get(prevMat, time - 1, r - 1, m - 1);
			mult(newValue, prUp, diagonallyDown);
			// *newValue = prUp * diagonallyDown;
		}
	}
	else {
		// r == 0
		if (m == 0) {
			// both zero
			auto right = get(prevMat, time - 1, r + 1, m);
			auto diagonallyUp = get(prevMat, time - 1, r + 1, m + 1);

			add(tempA, diagonallyUp, right);
			mult(newValue, prDown, tempA);
			// *newValue = prDown * (diagonallyUp + right);
		}
		else if (m == -1) {
			// no incoming down
			auto up = get(prevMat, time - 1, r, m + 1);
 			
			mult(newValue, prDown, up);
			// *newValue = prDown * up;
		}
		else {
			// reach zero, margin negative
			// up => diagonally up, down => reach stays
			auto up = get(prevMat, time - 1, r, m + 1);
			auto diagonallyUp = get(prevMat, time - 1, r + 1, m + 1);

			add(tempA, diagonallyUp, up);
			mult(newValue, prDown, tempA);
			// *newValue = prDown * (up + diagonallyUp);
		}
	}
}

void ReachAndMargin::evolve() {
	// prepare for the evolution
	time++;
	swapMat();

	for (int rho = 0; rho <= R + time; rho++) {			
		int below = -time;
		for (int mu = below; mu <= rho; mu++) {
			//if (!isValidCoord(time, rho, mu)) 
			//	continue;
			calcNewValue(rho, mu);			
		}
	}

}

DoubleVector ReachAndMargin::probability(const int minMu) const {
	assert(minMu >= -time && minMu <= R + time);

	DoubleVector probs(NFRACS, 0);
	// chain at the final state
	for (int rho = 0; rho <= R + time; rho++) {
		for (int mu = minMu; mu <= rho; mu++) {
			// these are the cells with nonnegative margin
			auto cell = get(curMat, time, rho, mu);
			for (int a = 0; a < NFRACS; a++)
				probs[a] += cell[a];
		}
	}

	// Now add tail probabilities
	if (R != 0)
		for (auto a = 0; a < NFRACS; a++)
			probs[a] += stationaryRhoTail(R + 1, advFractions[a]);
	// Done
	return probs;
}

DoubleVector ReachAndMargin::forkableProbability() const {
	return probability(0);
}
DoubleVector ReachAndMargin::totalProbability() const {
	return probability(-time);
}

bool ReachAndMargin::isValidCoord(const int time, const int rho, const int mu) const {
	int below = -time;

	bool goodTime = 0 <= time && time <= N;
	bool goodRho = 0 <= rho && rho <= (R + time);
	bool goodMu = (below <= mu) && (mu <= rho);

	return goodTime && goodRho && goodMu;

}


double* ReachAndMargin::get(MatrixType mat, const int time, const int rho, const int mu) const {
	if (isValidCoord(time, rho, mu) )  {
		// return mat[rho][MuZero + mu];
		return &mat[((MuZero + mu) * MAX_COLS + rho) * MAX_FRACS];
	}
	else
		return zero;
}

std::ostream& operator<<(std::ostream& out, const DoubleVector& vec) {
	out << "\t";
	int n = vec.size();
	for (int t = 0; t < n; t++){
		if (t > 0)
			out << ", ";
		out << std::setprecision(12) << vec[t];
	}
	//out << "\t";
	return out;
}


std::ostream& operator<<(std::ostream& out, const ReachAndMargin& self) {
	auto M = self.curMat;
	out << "=== Matrix at t: " << time << " ===\n";

	int below = - self.time;
	int above = self.R + self.time;
	for (int mu = above; mu >= below; mu--) {
		for (int rho = 0; rho <= above; rho++) {
			auto cell = self.get(self.curMat, self.time, rho, mu);
			out << "[";
			for (int a = 0; a < self.NFRACS; a++) {
				if (a > 0)
					out << ", ";
				out << std::setw(8) << cell[a];
			}
			out << "]";
		}
		out << "\t" << std::endl;
	}
	return out;
}



