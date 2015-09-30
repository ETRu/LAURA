
#pragma once

#if defined (WIN32)
#include "igraph/igraph.h"
#else
#include "igraph.h"
#endif

#include "logdebug.h"
#include "igraphutil.h"

class LAURA_Histogram_1D {
    
public:
    //histogram's extremes
    double min;
    double max;
    //number of bins
    int nbins;
    //vector storing the histogram bins value:
    igraph_vector_t bins;
    
    //mean and variance.
    double mean;
    double variance;
    
    void Clear();
    void CreateFromArray(igraph_vector_t inarray, int innbins, int normalize);
    void CreateFromArrayMinMax(igraph_vector_t inarray, int innbins, float min, float max, int normalize);
    
};


void testda();

