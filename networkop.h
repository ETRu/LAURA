
#pragma once

#if defined (WIN32)
#include "igraph/igraph.h"
#else
#include "igraph.h"
#endif

#include "igraphutil.h"

#include "randomgen.h"
#include "logdebug.h"
#include <math.h>
#include <time.h>

//-- shift for drawing correctly the layout
#define SHIFT 0.8


//-------------------------------------------- GENERAL USE FUNCTIONS

//layout
void adjustlayout();

//  compute stationary state
void StationaryState();


//-------------------------------------------- NETWORK GENERATORS

//lattices generators
void generatelattice(int mylatticedim, int mylatticeside, int myistoro, int myrand, int diss, double drate, int dnode);
void generatepotlat2d(int mylatticeside, int myistoro, int diss, double drate, int dnode);

//random generators
int generaterandom1(int n, float probability, int diss, double drate, int dnode);


//clustered generators
void generateclustsym(int clusnumber, int clusdim, double interconn, double intraconn, double wall);

void generateclustsympots(int clusnumber, int clusdim, int cnodes, double interconn, double intraconn, double wall, int dissipate, double dissipation, int dnode);

//clustered gerarchic 2 layer generator

void generateclusger2(int clusnumber1, int clusnumber2, int clusdim1, int cnodes1,
                      int cnodes2, double interconn, double intraconn1, double intraconn2,
                      double wall1, double wall2, int dissipate, double dissipation, int source);



//-------------------------------------------- LOAD / SAVE NETWORKS

void loadnetwork(char *mypath);
void savenetwork(char *mypath);




















