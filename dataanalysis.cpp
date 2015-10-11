
#include "dataanalysis.h"







/*This function destroys all the igraph variables of a LAURA_Histogram_1D object.
 */

void LAURA_Histogram_1D::Clear(){
    igraph_vector_destroy(&bins);
}




/*This function generates a 1-D histogram from an IGRAPH data array.
 The histogram is represented as an igraph array.
 * inarray: igraph array containing data.
 * innbins: number of bins.
 * normalize: boolean, 1 if the graph must be normalized, 0 otherwhise.
 */

void LAURA_Histogram_1D::CreateFromArray(igraph_vector_t inarray, int innbins, int normalize){
    
    double binwidth;
    int entries;
    
    islog=0;
    
    if (igraph_vector_empty(&inarray)!=0) {
        logDebug("\nError: empty vector. No data. Impossible to generate a Histogram\n");
        return;
    }
    
    nbins=innbins;
    min=igraph_vector_min(&inarray);
    max=igraph_vector_max(&inarray);
    entries=igraph_vector_size(&inarray);
    //adjust extremes
    max=max+max/100;
    
    //logDebug("\n NBINS=%i    min=%f   max=%f\n\n", nbins,min,max);
    binwidth=(max-min)/nbins;
    igraph_vector_init(&bins,nbins);
    igraph_vector_null(&bins);
    
    //BINNING
    for(int i=0;i<nbins;++i){
        for (int j=0; j<entries; ++j) {
            
            if(VECTOR(inarray)[j]>=(min+i*binwidth) && VECTOR(inarray)[j]<(min+(i+1)*binwidth) ){
                
                ++VECTOR(bins)[i];
                
            }
            
        }
    }
    
    //NORMALIZATION
    if (normalize==1) {
        double norm;
        norm=igraph_vector_sum(&bins);
        igraph_vector_scale(&bins,1./norm);
        
    }
    
    
    //calculate mean
    mean=0;
    mean=igraph_vector_sum(&bins);
    mean=mean/entries;
    
    //calculate variance
    /*  variance=0;
     for(int i=0;i<entries;++i){
     variance=variance+(VECTOR(bins)[i]-mean)*(VECTOR(bins)[i]-mean);
     }
     
     variance=variance/(entries-1);
     */
    
    
    
    return;
}



void LAURA_Histogram_1D::CreateFromArrayLog(igraph_vector_t inarray, int innbins, int normalize){
    
    double binwidth;
    int entries;
    
    islog=1;
    
    if (igraph_vector_empty(&inarray)!=0) {
        logDebug("\nError: empty vector. No data. Impossible to generate a Histogram\n");
        return;
    }
    
    nbins=innbins;
    min=igraph_vector_min(&inarray);
    max=igraph_vector_max(&inarray);
    entries=igraph_vector_size(&inarray);
    //adjust extremes
    max=max+max/100;
    
    //logDebug("\n NBINS=%i    min=%f   max=%f\n\n", nbins,min,max);
    binwidth=(max-min)/nbins;
    igraph_vector_init(&bins,nbins);
    igraph_vector_null(&bins);
    
    
    
    //BINNING
    for(int i=0;i<nbins;++i){
        for (int j=0; j<entries; ++j) {
            
            if(VECTOR(inarray)[j]>=(min+i*binwidth) && VECTOR(inarray)[j]<(min+(i+1)*binwidth) ){
                
                ++VECTOR(bins)[i];
                
            }
            
        }
    }
    
    
    
    //BINNING
   /* for(int i=0;i<nbins;++i){
        for (int j=0; j<entries; ++j) {
            //printf("\n BIN %i entry %i value %f exp %f ", i, j,VECTOR(inarray)[j],pow(10,VECTOR(inarray)[j]) );
            if(VECTOR(inarray)[j]>0)  {
                double tb;
                tb = VECTOR(inarray)[j];
                
                //if( log10(tb)) >= (min+(i*binwidth)) &&
                //log10(tb) <  (min+((i+1)*binwidth))
                // )
    
                
                if( pow(10,tb) >= (min+(i*binwidth)) &&
                   pow(10,tb) <  (min+((i+1)*binwidth))
                   )
                {
                    
                    ++VECTOR(bins)[i];
                    // printf(" - added");
                    
                    
                }
                
            }
            //printf("\n");
            
        }
    }
    */
    
    //NORMALIZATION
    if (normalize==1) {
        double norm;
        norm=igraph_vector_sum(&bins);
        igraph_vector_scale(&bins,1./norm);
        
    }
    
    
    //calculate mean
    mean=0;
    mean=igraph_vector_sum(&bins);
    mean=mean/entries;
    
    //calculate variance
    /*  variance=0;
     for(int i=0;i<entries;++i){
     variance=variance+(VECTOR(bins)[i]-mean)*(VECTOR(bins)[i]-mean);
     }
     
     variance=variance/(entries-1);
     */
    
    
    
    return;
}




void LAURA_Histogram_1D::CreateFromArrayMinMax(igraph_vector_t inarray, int innbins, float min, float max, int normalize){
    
    double binwidth;
    int entries;
    
    islog=0;
    
    if (igraph_vector_empty(&inarray)!=0) {
        logDebug("\nError: empty vector. No data. Impossible to generate a Histogram\n");
        return;
    }
    
    nbins=innbins;
    entries=igraph_vector_size(&inarray);
    //adjust extremes
    max=max+max/100;
    
    //logDebug("\n NBINS=%i    min=%f   max=%f\n\n", nbins,min,max);
    binwidth=(max-min)/nbins;
    igraph_vector_init(&bins,nbins);
    igraph_vector_null(&bins);
    
    //BINNING
    for(int i=0;i<nbins;++i){
        for (int j=0; j<entries; ++j) {
            
            if(VECTOR(inarray)[j]>=(min+i*binwidth) && VECTOR(inarray)[j]<(min+(i+1)*binwidth) ){
                
                ++VECTOR(bins)[i];
                
            }
            
        }
    }
    
    //NORMALIZATION
    if (normalize==1) {
        double norm;
        norm=igraph_vector_sum(&bins);
        igraph_vector_scale(&bins,1./norm);
        
    }
    
    
    
    
    //calculate mean
    mean=0;
    for(int i=0; i<nbins;++i){
        float naltrotemp;
        naltrotemp=(min+((i+1./2)*binwidth))*VECTOR(bins)[i];
        mean=mean+naltrotemp;
        
    }
    
    mean=mean/igraph_vector_sum(&bins);
    
    
    
    //calculare variance
    
    variance=0;
    for(int i=0; i<nbins;++i){
        float naltrotemp;
        naltrotemp=((min+((i+1./2)*binwidth))-mean)*((min+((i+1./2)*binwidth))-mean)*VECTOR(bins)[i];
        variance=variance+naltrotemp;
        
    }
    
    variance=sqrt(variance);
    //variance=variance/(igraph_vector_sum(&bins)-1);
    
    
    return;
    
}


void LAURA_Histogram_1D::CreateFromArrayMinMaxLog(igraph_vector_t inarray, int innbins, float min, float max, int normalize){
    
    double binwidth;
    int entries;
    
    islog=1;
    
    if (igraph_vector_empty(&inarray)!=0) {
        logDebug("\nError: empty vector. No data. Impossible to generate a Histogram\n");
        return;
    }
    
    nbins=innbins;
    entries=igraph_vector_size(&inarray);
    //adjust extremes
    max=max+max/100;
    
    //logDebug("\n NBINS=%i    min=%f   max=%f\n\n", nbins,min,max);
    binwidth=(max-min)/nbins;
    igraph_vector_init(&bins,nbins);
    igraph_vector_null(&bins);
    
    //BINNING
    
    for(int i=0;i<nbins;++i){
        for (int j=0; j<entries; ++j) {
            //printf("\n BIN %i entry %i value %f exp %f ", i, j,VECTOR(inarray)[j],pow(10,VECTOR(inarray)[j]) );
            if(VECTOR(inarray)[j]>0)  {
                double tb;
                tb = VECTOR(inarray)[j];
                
                /*if( log10(tb)) >= (min+(i*binwidth)) &&
                   log10(tb) <  (min+((i+1)*binwidth))
                   )
                */
                
                if( pow(10,tb) >= (min+(i*binwidth)) &&
                    pow(10,tb) <  (min+((i+1)*binwidth))
                    )
                {
                    
                    ++VECTOR(bins)[i];
                   // printf(" - added");
                    
                    
                }
                
            }
            //printf("\n");
            
        }
    }
    
    //NORMALIZATION
    if (normalize==1) {
        double norm;
        norm=igraph_vector_sum(&bins);
        igraph_vector_scale(&bins,1./norm);
        
    }
    
    
    
    
    //calculate mean
    mean=0;
    for(int i=0; i<nbins;++i){
        float naltrotemp;
        naltrotemp=(min+(i+1./2)*binwidth)*VECTOR(bins)[i];
        mean=mean+naltrotemp;
        
    }
    
    mean=mean/igraph_vector_sum(&bins);
    
    
    
    //calculare variance
    
    variance=0;
    for(int i=0; i<nbins;++i){
        float naltrotemp;
        naltrotemp=((min+((i+1./2)*binwidth))-mean)*((min+((i+1./2)*binwidth))-mean)*VECTOR(bins)[i];
        variance=variance+naltrotemp;
        
    }
    
    variance=sqrt(variance);
    //variance=variance/(igraph_vector_sum(&bins)-1);
    
    
    return;
    
}





/////////***********

void testda(){
    
    
    LAURA_Histogram_1D thisto;
    
    igraph_vector_t myarray;
    
    printf("\n");
    logDebug(" TEST \n");
    
    
    igraph_vector_init(&myarray,10);
    igraph_vector_fill(&myarray,1.2);
    VECTOR(myarray)[0]=0.1;
    VECTOR(myarray)[1]=4.3;
    
    print_vector_line(&myarray, stdout);
    
    thisto.CreateFromArray(myarray,4,0);
    
    print_vector_line(&(thisto.bins),stdout);
    
    thisto.Clear();
    
    thisto.CreateFromArray(myarray,4,1);
    
    print_vector_line(&(thisto.bins),stdout);
    
    printf("\n");
    logDebug("END TEST \n");
    
    
    thisto.Clear();
    igraph_vector_destroy(&myarray);
    
    return;
    
}










