#include "Randoms.cpp"

class ACO {
public:
	ACO (int nAnts, int nNodes, 
		double alpha, double beta, double q, double ro, double taumax,
		int initnode);
	virtual ~ACO ();
	
	void init ();
	
	void connectNODES (int nodei, int nodej);
	void setnodePOSITION (int node, double x, double y);
	
	void printPHEROMONES ();
	void printGRAPH ();
	void printRESULTS ();
	
	void optimize (int ITERATIONS);

private:
	double distance (int nodei, int nodej);
	bool exists (int nodei, int nodec);
	bool vizited (int antk, int c);
	double PHI (int nodei, int nodej, int antk);
	
	double congestion (int antk);
	
	int node ();
	void route (int antk);
	int valid (int antk, int iteration);
	
	void updatePHEROMONES ();

	
	int NUMBEROFANTS, NUMBEROFNODES, INITIALnode;
	double ALPHA, BETA, Q, RO, TAUMAX;
	
	double BESTcongestion;
	int *BESTROUTE;

	int **GRAPH, **ROUTES;
	double **NODES, **PHEROMONES, **DELTAPHEROMONES, **PROBS;

	Randoms *randoms;
};

