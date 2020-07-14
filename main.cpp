/*
@author:	Diogo A. B. Fernandes
@contact:	diogoabfernandes@gmail.com
@license:	see LICENSE
*/

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <climits>
#include <stdio.h>
#include "ACO.h"

using namespace std;
#define ITERATIONS		(int) 5

#define NUMBEROFANTS	(int) 6
#define NUMBEROFNODES	(int) 16

// if (ALPHA == 0) { stochastic search & sub-optimal route }
#define ALPHA			(double) 0.5
// if (BETA  == 0) { sub-optimal route }
#define BETA			(double) 0.8
// Estimation of the suspected best route.
#define Q				(double) 80
// Pheromones evaporation. 
#define RO				(double) 0.2
// Maximum pheromone random number.
#define TAUMAX			(int) 2

//#define INITIALnode		(int) 15

int main() {
	cout<<"Enter Starting node"<<endl;
	int strt;
	cin>>strt;

	ACO *ANTS = new ACO (NUMBEROFANTS, NUMBEROFNODES, 
			 			ALPHA, BETA, Q, RO, TAUMAX,
			 			strt);

	ANTS -> init();
    ANTS -> connectNODES (0, 1);
	ANTS -> connectNODES (0, 4);
	ANTS -> connectNODES (1, 2);
	ANTS -> connectNODES (1, 5);
	ANTS -> connectNODES (2, 3);
	ANTS -> connectNODES (2, 6);
	ANTS -> connectNODES (3, 7);
	ANTS -> connectNODES (4, 0);
	ANTS -> connectNODES (4, 8);
	ANTS -> connectNODES (4, 5);
	ANTS -> connectNODES(5, 6);
	ANTS -> connectNODES(5, 9); 
	ANTS -> connectNODES(6, 7);
        ANTS -> connectNODES(6, 10);
        ANTS -> connectNODES(7, 11);
	ANTS -> connectNODES(8, 9); 
        ANTS -> connectNODES(8, 4);
	ANTS -> connectNODES(8, 12);
	ANTS -> connectNODES(9, 10);
	ANTS -> connectNODES(9, 13);
        ANTS -> connectNODES(10, 11);
	ANTS -> connectNODES(10, 14);
	ANTS -> connectNODES(11, 15);
        ANTS -> connectNODES(12, 13);
        ANTS -> connectNODES(13, 14);
        ANTS -> connectNODES(14, 15);


	ANTS -> setnodePOSITION (0,  1,  1);
	ANTS -> setnodePOSITION (1, 10, 10);
	ANTS -> setnodePOSITION (2, 20, 10);
	ANTS -> setnodePOSITION (3, 10, 30);
	ANTS -> setnodePOSITION (4, 15,  5);
	ANTS -> setnodePOSITION (5, 10,  1);
	ANTS -> setnodePOSITION (6, 20, 20);
	ANTS -> setnodePOSITION (7, 20, 30);
    ANTS -> setnodePOSITION (8, 26, 10);
	ANTS -> setnodePOSITION (9, 10, 12);
	ANTS -> setnodePOSITION (10, 20, 15);
	ANTS -> setnodePOSITION (11, 16, 30);
	ANTS -> setnodePOSITION (12, 15, 25);
	ANTS -> setnodePOSITION (13, 10, 11);
	ANTS -> setnodePOSITION (14, 20, 22);
	ANTS -> setnodePOSITION (15, 22, 32);

	ANTS -> printGRAPH ();

	ANTS -> printPHEROMONES ();

	ANTS -> optimize (ITERATIONS);

	ANTS -> printRESULTS ();

	return 0;
}
