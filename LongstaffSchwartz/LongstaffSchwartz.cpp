// LongstaffSchwartz.cpp : définit le point d'entrée pour l'application console.
//
#include <iostream>
#include <ctime>
#include "BlackScholesModel.h";
#include "AmericanPricer.h";
#include <pnl/pnl_basis.h>

int main()
{
	int nbSSJ = 3;
	PnlMat *corr = pnl_mat_create_from_zero(nbSSJ, nbSSJ);
	PnlVect *spots = pnl_vect_create(nbSSJ);
	PnlVect *sigma = pnl_vect_create(nbSSJ);
	PnlVect *dividends = pnl_vect_create(nbSSJ);

	
	PnlVect* rates = pnl_vect_create(nbSSJ);
	
	


	PnlVect *trends = pnl_vect_create(nbSSJ);
	double nbTimeSteps = 10;
	double T = 3.0;


	PnlRng *rng = pnl_rng_create(PNL_RNG_MERSENNE);
	pnl_rng_sseed(rng, std::time(NULL));
	
	for (int i = 0; i < nbSSJ; i++) {
		LET(dividends, i) = 0;
		LET(rates, i) = 0.05;
		LET(sigma, i) = 0.2;
		LET(spots, i) = 100.0;
		MLET(corr, i, i) = 1;
		LET(trends, i) = 0.06;
	}
	RatesMarkets *ratesMarkets = new ConstantRatesMarkets(EUR, rates);



	BlackScholesModel* BS = new BlackScholesModel(nbSSJ, ratesMarkets, trends, dividends, sigma,
		spots, corr);

	// Nombre de samples 
	int M = 1000;

	PnlVect* lambda = pnl_vect_create(nbSSJ);
    for (int i = 0; i < nbSSJ; i++) {
		LET(lambda, i) = -1.0 / (double)nbSSJ;;
    }
	double strike = -100.0;
	BasketOption* option = new BasketOption(T, nbSSJ, nbTimeSteps, lambda, strike, EUR);
	int d = 3;
	AmericanPricer AP = AmericanPricer(d, M,option, BS, rng);

	//AP.Print();
	//PnlVect* v = pnl_vect_new();
	//AP.compute_g(v, 3, 3);
	

	double tauNext = nbTimeSteps;
	//std::cout << "conditional expectation = " <<AP.compute_Conditional_Expectation(p, d, nbTimeSteps - 1, tauNext, 0) << std::endl;
	std::cout << "price = " << AP.price() << "\n";

	system("PAUSE");

    return 0;
}

