
#pragma once

#if defined (WIN32)
#include "igraph/igraph.h"
#else
#include "igraph.h"
#endif

#include "igraphutil.h"
#include "randomgen.h"
#include "dataanalysis.h"
#include "form.h"
#include <time.h>

#define NTRYSTOP 100

void InitialStateTAS(int beginner, int useran, int maxrun);
void InitialStateRAN(int seed, int maxrun);
void InitialStateLOAD();
void InitialStateSTAT(int maxrun);
void Evolution(double dt);
void TurboRun(double dt);