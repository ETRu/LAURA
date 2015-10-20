#include <math.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>

#if defined (WIN32)
#include <io.h>
#else
#include <unistd.h>
#endif


using namespace std;

                                                #define N 1000000      // <<-------  MAXIMUM TREATABLE DATA





void analysys(int window){
    
    
    char osrm[100]; //output stream for running mean
    char osfl[100]; //output stream for fluctuuations
    char osflm[100]; //output stream for fluctuations on running mean
    char oshfl[100]; //output stream for fluctuuations' histogram
    char oshflm[100]; //output stream for fluctuations' on running mean histogram

    sprintf(osrm, "rm%i.txt",window);
    sprintf(osfl, "fluct.txt");
    sprintf(osflm, "rmfluct%i.txt",window);
    sprintf(oshfl, "hfluct.txt");
    sprintf(oshflm, "hrmfluct%i.txt",window);
    
    ifstream indata;
                                                /* INPUT FILE  */ indata.open("test/th9007.txt");
    
    
    FILE *outactimean;
    
                                                /* OUTPUT FILE - MEAN ACTIVATION*/
                                                outactimean=fopen(osrm,"w");
    
    FILE *outfluct;
    
                                                /* OUTPUT FILE - FLUCTUATIONS*/
                                                outfluct=fopen(osfl,"w");
    
    FILE *outwfluct;
    
                                                /* OUTPUT FILE - MEAN FLUCTUATIONS*/
                                                outwfluct=fopen(osflm,"w");
    
    FILE *outflucthist;
    
                                                /* OUTPUT FILE - FLUCTUATIONS HISTOGRAM*/
                                                outflucthist=fopen(oshfl,"w");
    
    FILE *outwflucthist;
    
                                                /* OUTPUT FILE - MEAN FLUCTUATIONS HISTOGRAM*/
                                                outwflucthist=fopen(oshflm,"w");
    
    
    
    double data[N];
    double histo[10000];
    for(int i=0;i<N;++i){data[i]=0;}
    
    double temp;
    int ndata;
    int cont;
    
    int nbins=50;
    double min, max;
    double binw;
    
    double mean;
    
    double tmean;
    
    
    printf("\n:::::::::::::::::::::::::::::: RUNNING ANALYSYS - WINDOW USED: %i \n::::\n::::\n::::\n::::\n::::",window);
    
    //GET DATA
    indata >> temp;
    for (ndata=0; ndata<N; ++ndata) {
        indata >> data[ndata];
        if(indata.eof()){break;}
        //printf(" %f ", data[ndata]);
    }
    
    printf("\n::::::: (%i-1)=%i DATA READEN",ndata, ndata-1);
    
    if(ndata==0) abort();
    
    
    //evaluate total mean
    tmean=0;
    for(int i=0; i<ndata;++i){
        tmean=tmean+data[i];
    }
    tmean=tmean/ndata;
    
    printf("\n::::::: TOTAL MEAN: %f",tmean);
    
    
    
    
    //EVALUATE FLUCTUATIONS --------------------------------------------------------------------------------------
    
    printf("\n::::::: Evaluating FLUCTUATIONS");
    for(int i=0; i<ndata; ++i){
        fprintf(outfluct,"%f\n",fabs( tmean-data[i]));
    }
    

    
    // EVALUATE WINDOW MEAN AND FLUCTUATION ----------------------------------------------------------------------
    printf("\n::::::: RUNNING MEAN \n ------------ Window dimension = %i, ndata = %i", window,ndata);
    
    if(ndata<window) {printf(" ---> NOT COMPUTABLE \n \n ------ STOP. \n\n"); abort();}
    printf(" ---> COMPUTABLE ");
    
    cont=0;
    do{
     
        mean=0;
        for(int i=0; i<window;++i){
            mean=mean+data[cont+i];
        }
        mean=mean/window;
     
        fprintf(outactimean,"%i %f\n",(window/2)+cont,mean);
        fprintf(outwfluct, "%f\n", fabs(mean-tmean));
        
       ++cont;
    
    } while(cont<ndata-window);
    
    
    
    
    
    
    indata.close();
    fclose(outactimean);
    fclose(outfluct);
    fclose(outwfluct);
    
    
    
    printf("\n::::::: COMPUTING HISTOGRAMS :::::::::::::::::::::: ");
    
    
    // FLUCTUATION HISTOGRAM -----------------

    printf("\n::::::: Normal Fluctuations");
    
    //read from the file
    indata.open(osfl);
    for (ndata=0; ndata<N; ++ndata) {
        indata >> data[ndata];
        if(indata.eof()){break;}
        //printf(" %f ", data[ndata]);
    }
    
    indata.close();
    
    
    printf(" ndata= %i", ndata);
    
    printf("\n::::::: Binning");
    
    //set min
    min=1000000;
    for (int i=0; i<ndata; ++i) {
        if(data[i]<min){min=data[i];}
    }
    
    
    //set max
    max=-1;
    for (int i=0; i<ndata; ++i) {
        if(data[i]>max){max=data[i];}
    }
    max=max+max/10000;
    
    //set binwidth
    binw=(max-min)/nbins;
    
    //clean histogram
    for (int i=0; i<nbins; ++i) {
        histo[i]=0;
    }
    
    //binning
    for (int j=0; j<ndata; ++j) {
        for (int i=0; i<nbins; ++i) {
            if( min+(i*binw) <= data[j] && data[j] < min+((i+1)*binw)) {++histo[i];
                //printf("\n dato %i (%f) in bin %i (%f - %f)", j, data[j], i, min+(i*binw),min+((i+1)*binw));
            }
        }
    }
    //normalize
    for (int i=0; i<nbins; ++i) {
        histo[i]=histo[i]/ndata;
    }
    
    
    
    //print data
    for (int i=0; i<ndata; ++i) {
        if(histo[i]!=0) {fprintf(outflucthist,"%f %f\n", min+ (i+0.5)*binw, histo[i]);}
    }
    

    
    
    
    
    
    
    
    
    
    // RUNNING MEAN FLUCTUATION HISTOGRAM -----------------
    
    
    
    printf("\n::::::: Running Mean Fluctuations");
    
    //read from the file
    indata.open(osflm);
    for (ndata=0; ndata<N; ++ndata) {
        indata >> data[ndata];
        if(indata.eof()){break;}
        //printf(" %f ", data[ndata]);
    }
    
    indata.close();
    
    
    printf(" ndata= %i", ndata);
 
    printf("\n::::::: Binning");

    //set min
    min=1000000;
    for (int i=0; i<ndata; ++i) {
        if(data[i]<min){min=data[i];}
    }
    
    //set max
    max=-1;
    for (int i=0; i<ndata; ++i) {
        if(data[i]>max){max=data[i];}
    }
    max=max+max/10000;
    
    //set binwidth
    binw=(max-min)/nbins;
    
    //clean histogram
    for (int i=0; i<nbins; ++i) {
        histo[i]=0;
    }
    
    //binning
    for (int j=0; j<ndata; ++j) {
        for (int i=0; i<nbins; ++i) {
            if( min+(i*binw) <= data[j] && data[j] < min+((i+1)*binw)) {++histo[i];
                //printf("\n dato %i (%f) in bin %i (%f - %f)", j, data[j], i, min+(i*binw),min+((i+1)*binw));
            }
        }
    }
    //normalize
    for (int i=0; i<nbins; ++i) {
        histo[i]=histo[i]/ndata;
    }
    
    
    
    //print data
    for (int i=0; i<ndata; ++i) {
        if(histo[i]!=0) {fprintf(outwflucthist,"%f %f\n", min+ (i+0.5)*binw, histo[i]);}
    }
    

    printf("\n\n");
    
    fclose(outflucthist);
    fclose(outwflucthist);

}




int main(){

    analysys(5);    //2
  
    analysys(10);   //3
    /*
    analysys(20);   //4
    analysys(30);   //5
    analysys(40);   //6
    analysys(50);   //7
    analysys(75);   //8
    analysys(100);  //9
    analysys(125);  //10
    */
    
    return 0;}





