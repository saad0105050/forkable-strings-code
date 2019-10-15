#pragma once

// Memory leak detection for windows
#ifdef _WIN32
#include "win_mem_leak.h"
#endif

#include <chrono>
#include<vector>
#include<iostream>

using UintVector = std::vector < unsigned int >;
using DoubleVector = std::vector < double >;
using ValueType = double*;
std::ostream& operator<<(std::ostream& out, const DoubleVector& vec);


class ReachAndMargin
{
public:
	enum{
		NMAX = 500,		// maximum length of evolution (i.e., number of timesteps)
		RMAX,			// maximum reach
		MAX_FRACS = 50	// maximum number of adversarial fractions to consider
	};
	const int	N,		// length of evolution
		R,		// starting reach
		NFRACS	// number of adversarial fractions to consider
		;
private:
	static const int
		MAX_COLS = NMAX + RMAX,
		MAX_ROWS = 2 * NMAX + RMAX + 1,
		MuZero = NMAX,
		TotalCount = MAX_COLS * MAX_ROWS * MAX_FRACS;

	using MatrixType = double*;
	double*** matrixOne;
	double*** matrixTwo;
	double matrixOneRaw[TotalCount];
	double matrixTwoRaw[TotalCount];
	MatrixType prevMat;
	MatrixType curMat;

	int time;


	const UintVector Ns;

	const DoubleVector advFractions;
	double *zero, *prDown, *prUp, *tempA, *tempB, * tempC;

	void calcNewValue(const int rho, const int mu);


	void mult(double* dest, const double* one, const double* two);
	void add(double* dest, const double* one, const double* two);
	DoubleVector probability(const int minMu) const;
	bool isValidCoord(const int time, const int rho, const int mu) const;
	void swapMat();
public:
	ReachAndMargin(
		const int N, 
		const int R, 
		const DoubleVector& advFractions);
	virtual ~ReachAndMargin();

	double* get(MatrixType mat, const int time, const int rho, const int mu) const;
	// you can't modify a value obtained by get

	void evolve();
	DoubleVector forkableProbability() const;
	DoubleVector totalProbability() const;
	long totalMemory() const {
		return 2 * TotalCount * sizeof(double);
	}

	friend std::ostream& operator<< (std::ostream& out, const ReachAndMargin& self);

};



