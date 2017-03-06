#include "AmericanPricer.h";

// Factorielle de x
inline int Factorial(int x) {
	return (x == 1 ? x : x * Factorial(x - 1));
}

// Constructeur 
AmericanPricer::AmericanPricer(int d,int nbSamples, BasketOption* option, BlackScholesModel* BS, PnlRng *rng) {
	this->nbSamples = nbSamples;
	this->nbSSJ = BS->getNbSSJ();
	this->nbTimeSteps = option->nbTimeSteps;
	this->rng = rng;
	this->option = option;
	this->r = GET(BS->getRatesMarket()->getSpots(),0);
	this->spots = BS->getSpots();
	this->d = d;
	// Vecteur de M pointeurs vers des paths
    simulated_paths = new PnlMat*[nbSamples];
	for (int i = 0; i < nbSamples; i++) {
		simulated_paths[i] = pnl_mat_create(nbTimeSteps + 1, nbSSJ);
		BS->simulateUnderRiskNeutralProba(simulated_paths[i], option->maturity, 0.0, nbTimeSteps, rng);
	}
	g = pnl_basis_create_from_degree(CANONICAL, d, nbSSJ);
	alpha = pnl_vect_create(g->len_T);
	S = pnl_mat_create(nbSamples, nbSSJ);
	Payoff = pnl_vect_create(nbSamples);

	a = new double[nbSSJ];
	b = new double[nbSSJ];
	
	double T = option->maturity;
	for (int i = 0; i < nbSSJ; i++) {
		double sigma = GET(BS->volatilities, i);
		a[i] = GET(spots,i) * exp(  MIN((r - 0.5*sigma*sigma)*T,0.0) - 3*sigma*sqrt(T) );
		b[i] = GET(spots, i) * exp(MAX((r - 0.5*sigma*sigma)*T, 0.0) + 3*sigma*sqrt(T));
	}
}


// Print
void AmericanPricer::Print() {
	for (int i = 0; i < nbSamples; i++) {
		pnl_mat_print(simulated_paths[i]);
		std::cout << std::endl;
	}
}


double AmericanPricer::price() {
	double timeStep = option->maturity / nbTimeSteps;
	/*
	  Calcul du prix à partir des tau_opt des différentes simulations 
	*/
	double sum = 0.0;
	for (int i = 0; i < nbSamples; i++) {
		std::cout << i << " / " << nbSamples << std::endl;
		double tau_opt = compute_tau_opt( i);
		//std::cout << "tau_opt = " << tau_opt << "\n";
		sum += exp(-r*tau_opt*timeStep)*option->payoff(simulated_paths[i], tau_opt);
	}
	sum /= (double) nbSamples;
	return sum;
}

int AmericanPricer::compute_tau_opt( int indexSimu) {
	double timeStep = option->maturity / nbTimeSteps;
	double tauNext = nbTimeSteps;
	/*
	   Pour n = (nbTimeSteps - 1) .. 1 :
	   tau_n = tn * 1{exp(-r*tn) * Phi(Stn) > esperance cond. }  + tau_n+1 * 1{exp(-r*tn) * Phi(Stn) <= esperance cond. } 
	*/
	for (int n = nbTimeSteps - 1; n > 0; n--) {
		double ce = compute_Conditional_Expectation(n, tauNext, indexSimu);
		double t_n = n*timeStep;
		if (ce < exp(-r*t_n)*option->payoff(simulated_paths[indexSimu], n)) {
			//std::cout << " Esperance conditionnelle  = " <<  ce << std::endl;
			//std::cout << " Payoff  = " << exp(-r*t_n)*option->payoff(simulated_paths[indexSimu], n) << std::endl;
			tauNext = n;
		}
	}

	/*
	  Dernière étape : on regarde si il vaut mieux exercer en t=0 où en t = timeStep
	*/
	double esperance = 0.0;
	for (int i = 0; i < nbSamples; i++) {
		esperance += option->payoff(simulated_paths[i], 1);
	}
	esperance /= (double)nbSamples;
	//std::cout << " Esperance conditionnelle  = " << exp(-r*timeStep)* esperance << std::endl;
	//std::cout << " Payoff  = " << option->payoff(simulated_paths[indexSimu], 0) << std::endl;
	if (option->payoff(simulated_paths[indexSimu],0) > exp(-r*timeStep)* esperance) {		
		tauNext = 0.0;
	}

	return tauNext;
}

double AmericanPricer::compute_Conditional_Expectation( int t, double tauNext,int indexSimu) {

	/* 
	   Creation de la matice S de taille nbSamples * nbSSJ : chaque ligne i de S 
	   contient les valeurs des sous-jacents de la simulation i
	*/
	PnlVect* row = pnl_vect_create(nbSSJ);
	for (int i = 0; i < nbSamples; i++) {
		pnl_mat_get_row(row,simulated_paths[i],t);
		pnl_mat_set_row(S, row, i);
	}
	
	/* 
	   Normalisation des sous-jacents (pour éviter l'explosion des polynomes)
	*/
	for (int i = 0; i < nbSamples; i++) {
		for (int j = 0; j < nbSSJ; j++) {
			MLET(S, i, j) = (MGET(S, i, j) - 0.5*(a[j] + b[j]) )/ (b[j] - a[j]);
		}
	}
	pnl_vect_free(&row);

   	//pnl_mat_print(S);
	
	/*
	   Construction du vecteur Payoff des payoffs actualisés 
	*/
	double timeStep = option->maturity / nbTimeSteps;
	for (int i = 0; i < nbSamples; i++) {
		double t_tauNext = tauNext*timeStep ;
		LET(Payoff, i) = exp(-r*t_tauNext) * option->payoff(simulated_paths[i],tauNext);
	}
	
	/*
	   Construction du vecteur alpha (coefficient de l'espérance conditionnelle dans la base choisie) 
	*/
	pnl_basis_fit_ls(g,alpha,S,Payoff);
    //pnl_vect_print(alpha);

	/*
	   Construction du vecteur X (point auquel on calcule le polynome)
	*/
	double* X = new double [nbSSJ];
	for (int j = 0; j < nbSSJ; j++) {
		X[j] = (MGET(simulated_paths[indexSimu],t,j) - 0.5*(a[j] + b[j]) ) / (b[j] - a[j]);
		//std::cout << X[j] << "\n";
	}	

	/*
	   Calcul de l'espérance conditionnelle à partir d'une base de polynômes (g), des coefficients du polynôme optimal (alpha) 
	   et du point où on calcule le polynôme (X)
	*/
	double ce = pnl_basis_eval(g, alpha, X);
	delete X;

	return ce;

}

// Destructeur
AmericanPricer::~AmericanPricer() {	

	for (int i = 0; i < nbSamples; i++) {
		pnl_mat_free(&simulated_paths[i]);
	}
	delete(simulated_paths);
	pnl_mat_free(&S);
	pnl_basis_free(&g);	
	pnl_vect_free(&alpha);
	pnl_vect_free(&Payoff);
	delete a;
	delete b;

}