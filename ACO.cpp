#include "ACO.h"

#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <cmath>
#include <limits>
#include <climits>

using namespace std;

ACO::ACO (int nAnts, int nNodes, 
		double alpha, double beta, double q, double ro, double taumax,
		int initnode) {
	
	NUMBEROFANTS 	= nAnts;
	NUMBEROFNODES 	= nNodes;
	ALPHA 			= alpha;
	BETA 			= beta;
	Q 				= q;
	RO 				= ro;
	TAUMAX 			= taumax;
	INITIALnode		= initnode;

	randoms = new Randoms (21);	
}
ACO::~ACO () {
	for(int i=0; i<NUMBEROFNODES; i++) {
		delete [] GRAPH[i];
		delete [] NODES[i];
		delete [] PHEROMONES[i];
		delete [] DELTAPHEROMONES[i];
		if(i < NUMBEROFNODES - 1) {
			delete [] PROBS[i];
		}
	}
	delete [] GRAPH;
	delete [] NODES;
	delete [] PHEROMONES;
	delete [] DELTAPHEROMONES;
	delete [] PROBS;
}

void ACO::init () {
	GRAPH 			= new int*[NUMBEROFNODES];
	NODES 			= new double*[NUMBEROFNODES];
	PHEROMONES 		= new double*[NUMBEROFNODES];
	DELTAPHEROMONES = new double*[NUMBEROFNODES];
	PROBS 			= new double*[NUMBEROFNODES-1];
	for(int i=0; i<NUMBEROFNODES; i++) {
		GRAPH[i] 			= new int[NUMBEROFNODES];
		NODES[i] 			= new double[2];
		PHEROMONES[i] 		= new double[NUMBEROFNODES];
		DELTAPHEROMONES[i] 	= new double[NUMBEROFNODES];
		PROBS[i] 			= new double[2];
		for (int j=0; j<2; j++) {
			NODES[i][j] = -1.0;
			PROBS[i][j]  = -1.0;
		}
		for (int j=0; j<NUMBEROFNODES; j++) {
			GRAPH[i][j] 			= 0;
			PHEROMONES[i][j] 		= 0.0;
			DELTAPHEROMONES[i][j] 	= 0.0;
		}
	}	

	ROUTES = new int*[NUMBEROFANTS];
	for (int i=0; i<NUMBEROFANTS; i++) {
		ROUTES[i] = new int[NUMBEROFNODES];
		for (int j=0; j<NUMBEROFNODES; j++) {
			ROUTES[i][j] = -1;
		}
	}
	
	BESTcongestion = (double) INT_MAX;
	BESTROUTE  = new int[NUMBEROFNODES];
	for (int i=0; i<NUMBEROFNODES; i++) {
		BESTROUTE[i] = -1;	
	}
}


void ACO::connectNODES (int nodei, int nodej) {
	GRAPH[nodei][nodej] = 1;
	PHEROMONES[nodei][nodej] = randoms -> Uniforme() * TAUMAX;
	GRAPH[nodej][nodei] = 1;
	PHEROMONES[nodej][nodei] = PHEROMONES[nodei][nodej];
}
void ACO::setnodePOSITION (int node, double x, double y) {
	NODES[node][0] = x;
	NODES[node][1] = y;
}
void ACO::printPHEROMONES () {	
	cout << " PHEROMONES: " << endl;
	cout << "  | ";
	for (int i=0; i<NUMBEROFNODES; i++) {
		printf("%5d   ", i);
	}
	cout << endl << "- | ";
	for (int i=0; i<NUMBEROFNODES; i++) {
		cout << "--------";
	}
	cout << endl;
	for (int i=0; i<NUMBEROFNODES; i++) {
		cout << i << " | ";
		for (int j=0; j<NUMBEROFNODES; j++) {
			if (i == j) {
				printf ("%5s   ", "x");
				continue;
			}
			if (exists(i, j)) {
				printf ("%7.3f ", PHEROMONES[i][j]);
			}
			else {
				if(PHEROMONES[i][j] == 0.0) {
					printf ("%5.0f   ", PHEROMONES[i][j]);
				}
				else {
					printf ("%7.3f ", PHEROMONES[i][j]);
				}
			}
		}
		cout << endl;
	}
	cout << endl;
}


double ACO::distance (int nodei, int nodej) {
	return (double) 
		sqrt (pow (NODES[nodei][0] - NODES[nodej][0], 2) + 
 			  pow (NODES[nodei][1] - NODES[nodej][1], 2));
}
bool ACO::exists (int nodei, int nodec) {
	return (GRAPH[nodei][nodec] == 1);
}
bool ACO::vizited (int antk, int c) {
	for (int l=0; l<NUMBEROFNODES; l++) {
		if (ROUTES[antk][l] == -1) {
			break;
		}
		if (ROUTES[antk][l] == c) {
			return true;
		}
	}
	return false;
}
double ACO::PHI (int nodei, int nodej, int antk) {
	double ETAij = (double) pow (1 / distance (nodei, nodej), BETA);
	double TAUij = (double) pow (PHEROMONES[nodei][nodej],   ALPHA);

	double sum = 0.0;
	for (int c=0; c<NUMBEROFNODES; c++) {
		if (exists(nodei, c)) {
			if (!vizited(antk, c)) {
				double ETA = (double) pow (1 / distance (nodei, c), BETA);
				double TAU = (double) pow (PHEROMONES[nodei][c],   ALPHA);
				sum += ETA * TAU;
			}	
		}	
	}
	return (ETAij * TAUij) / sum;
}

double ACO::congestion (int antk) {
	double sum = 0.0;
	for (int j=0; j<NUMBEROFNODES-1; j++) {
		sum += distance (ROUTES[antk][j], ROUTES[antk][j+1]);	
	}
	return sum;
}

int ACO::node () {
	double xi = randoms -> Uniforme();
	int i = 0;
	double sum = PROBS[i][0];
	while (sum < xi) {
		i++;
		sum += PROBS[i][0];
	}
	return (int) PROBS[i][1];
}

void ACO::route (int antk) {
	ROUTES[antk][0] = INITIALnode;
	for (int i=0; i<NUMBEROFNODES-1; i++) {		
		int nodei = ROUTES[antk][i];
		int count = 0;
		for (int c=0; c<NUMBEROFNODES; c++) {
			if (nodei == c) {
				continue;	
			}
			if (exists (nodei, c)) {
				if (!vizited (antk, c)) {
					PROBS[count][0] = PHI (nodei, c, antk);
					PROBS[count][1] = (double) c;
					count++;
				}

			}
		}
		
		// deadlock
		if (0 == count) {
			return;
		}
		
		ROUTES[antk][i+1] = node();
	}
}
int ACO::valid (int antk, int iteration) {
	for(int i=0; i<NUMBEROFNODES-1; i++) {
		int nodei = ROUTES[antk][i];
		int nodej = ROUTES[antk][i+1];
		if (nodei < 0 || nodej < 0) {
			return -1;	
		}
		if (!exists(nodei, nodej)) {
			return -2;	
		}
		for (int j=0; j<i-1; j++) {
			if (ROUTES[antk][i] == ROUTES[antk][j]) {
				return -3;
			}	
		}
	}
	
	if (!exists (INITIALnode, ROUTES[antk][NUMBEROFNODES-1])) {
		return -4;
	}
	
	return 0;
}

void ACO::printGRAPH () {
	cout << " GRAPH: " << endl;
	cout << "  | ";
	for( int i=0; i<NUMBEROFNODES; i++) {
		cout << i << " ";
	}
	cout << endl << "- | ";
	for (int i=0; i<NUMBEROFNODES; i++) {
		cout << "- ";
	}
	cout << endl;
	int count = 0;
	for (int i=0; i<NUMBEROFNODES; i++) {
		cout << i << " | ";
		for (int j=0; j<NUMBEROFNODES; j++) {
			if(i == j) {
				cout << "x ";	
			}
			else {
				cout << GRAPH[i][j] << " ";	
			}
			if (GRAPH[i][j] == 1) {
				count++;	
			}
		}
		cout << endl;
	}
	cout << endl;
	cout << "Number of connections: " << count << endl << endl;
}
void ACO::printRESULTS () {
	BESTcongestion += distance (BESTROUTE[NUMBEROFNODES-1], INITIALnode);
	cout << " BEST ROUTE:" << endl;
	for (int i=0; i<NUMBEROFNODES; i++) {
		cout << BESTROUTE[i] << " ";
	}
	cout << endl << "congestion: " << BESTcongestion << endl;
	
}

void ACO::updatePHEROMONES () {
	for (int k=0; k<NUMBEROFANTS; k++) {
		double rcongestion = congestion(k);
		for (int r=0; r<NUMBEROFNODES-1; r++) {
			int nodei = ROUTES[k][r];
			int nodej = ROUTES[k][r+1];
			DELTAPHEROMONES[nodei][nodej] += Q / rcongestion;
			DELTAPHEROMONES[nodej][nodei] += Q / rcongestion;
		}
	}
	for (int i=0; i<NUMBEROFNODES; i++) {
		for (int j=0; j<NUMBEROFNODES; j++) {
			PHEROMONES[i][j] = (1 - RO) * PHEROMONES[i][j] + DELTAPHEROMONES[i][j];
			DELTAPHEROMONES[i][j] = 0.0;
		}	
	}
}


void ACO::optimize (int ITERATIONS) {
	for (int iterations=1; iterations<=ITERATIONS; iterations++) {
		cout << flush;
		cout << "ITERATION " << iterations << " HAS STARTED!" << endl << endl;

		for (int k=0; k<NUMBEROFANTS; k++) {
			cout << " : ant " << k << " has been released!" << endl;
			while (0 != valid(k, iterations)) {
				cout << "  :: releasing ant " << k << " again!" << endl;
				for (int i=0; i<NUMBEROFNODES; i++) {
					ROUTES[k][i] = -1;	
				}
				route(k);
			}
			
			for (int i=0; i<NUMBEROFNODES; i++) {
				cout << ROUTES[k][i] << " ";	
			}
			cout << endl;
			
			cout << "  :: route done" << endl;
			double rcongestion = congestion(k);

			if (rcongestion < BESTcongestion) {
				BESTcongestion = rcongestion;
				for (int i=0; i<NUMBEROFNODES; i++) {
					BESTROUTE[i] = ROUTES[k][i];
				}
			}
			cout << " : ant " << k << " has ended!" << endl;				
		}		

		cout << endl << "updating PHEROMONES . . .";
		updatePHEROMONES ();
		cout << " done!" << endl << endl;
		printPHEROMONES ();
		
		for (int i=0; i<NUMBEROFANTS; i++) {
			for (int j=0; j<NUMBEROFNODES; j++) {
				ROUTES[i][j] = -1;
			}
		}

		cout << endl << "ITERATION " << iterations << " HAS ENDED!" << endl << endl;
	}
}
