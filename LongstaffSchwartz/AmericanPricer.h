
#pragma once
#include <iostream>
#include "BlackScholesModel.h";
#include "BasketOption.h";
#include <pnl/pnl_basis.h>
class AmericanPricer  {
public:
	PnlRng *rng;
	int nbSamples;	
	int nbSSJ;
	double nbTimeSteps;
	double r;
	PnlVect* spots;
	int d;

	AmericanPricer(int d, int nbSamples, BasketOption* option, BlackScholesModel* BS, PnlRng *rng);
	void Print();
	// t désigne la date à laquelle on calcule l'epérance condtionnelle
	double compute_Conditional_Expectation( int t, double tauNext,int indexSumu);
	int compute_tau_opt(int indexSimu);
	double price();
    ~AmericanPricer();

private:
	// Vecteur de M pointeurs vers des paths
	PnlMat** simulated_paths;
	BasketOption *option;

	PnlBasis* g ;
	PnlVect* alpha ;
	PnlMat* S ;
	PnlVect* Payoff;
	double* a;
	double* b;
};