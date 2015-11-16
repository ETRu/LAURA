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

                                                #define N 10000000      // <<-------  MAXIMUM TREATABLE DATA


double data[N];
double var[N];

int  peak[N];
int  plen[N];


char folder[100];

void analysys(int window, int d, int t, int s, double *eavg, double *efluct){
    
    
    char osrm[100]; //output stream for running mean
    char osfl[100]; //output stream for fluctuuations
    char osflm[100]; //output stream for fluctuations on running mean
    char oshp[100];  //output stream for PEAKS HISTOGRAM
    char oshfl[100]; //output stream for fluctuuations' histogram
    char oshflm[100]; //output stream for fluctuations' on running mean histogram
    char infile[100];

    sprintf(osrm, "%s/rm%i-d%i-s%i.txt",folder,window,d,s);
    sprintf(osfl, "%s/fluct.txt",folder);
    sprintf(osflm, "%s/rmfluct%i-d%i-s%i.txt",folder, window,d,s);
    sprintf(oshp, "%s/hp%i-d%i-s%i.txt",folder,window,d,s);
    sprintf(oshfl, "%s/hfluct.txt",folder);
    sprintf(oshflm, "%s/hrmfluct%i.txt",folder,window);
    
    //sprintf(infile, "data/clus-small/s10/-p%it%is%i-7.txt",p,t,s);
    
    sprintf(infile, "%s/-d%it%is%i-7.txt",folder,d,t,s);
    
    printf("---->  INFILE:  %s   \n",infile);
    
   
    ifstream indata;
                                                /* INPUT FILE  */ indata.open(infile);
    
    
    FILE *outactimean;
    
                                                /* OUTPUT FILE - MEAN ACTIVATION*/
                                                outactimean=fopen(osrm,"w");
    
    FILE *outfluct;
    
                                                /* OUTPUT FILE - FLUCTUATIONS*/
                                                outfluct=fopen(osfl,"w");
    
    FILE *outwfluct;
    
                                                /* OUTPUT FILE - MEAN FLUCTUATIONS*/
                                                outwfluct=fopen(osflm,"w");
    
    
    FILE *outpeakhist;
    
    /* OUTPUT FILE - FLUCTUATIONS HISTOGRAM*/
    outpeakhist=fopen(oshp,"w");
    
    
    FILE *outflucthist;
    
                                                /* OUTPUT FILE - FLUCTUATIONS HISTOGRAM*/
                                                outflucthist=fopen(oshfl,"w");
    
    FILE *outwflucthist;
    
                                                /* OUTPUT FILE - MEAN FLUCTUATIONS HISTOGRAM*/
                                                outwflucthist=fopen(oshflm,"w");
    
    
    
    double histo[10000];
    for(int i=0;i<N;++i){data[i]=0; var[i]=0; peak[i]=0; plen[i]=0;}
    
    double temp;
    int ndata;
    int cont;
    
    int nbins=50;
    double min, max;
    double binw;
    
    double mean;
    double f;
    
    double tmean;
    
    int ifrom;
    
    
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
    
    ifrom=1000;
    
    //evaluate total mean
    tmean=0;
    for(int i=ifrom; i<ndata;++i){
        tmean=tmean+data[i];
    }
    tmean=tmean/(ndata-ifrom);
    
    printf("\n::::::: TOTAL MEAN: %f",tmean);
    
    
    
    
    //EVALUATE FLUCTUATIONS --------------------------------------------------------------------------------------
    
    printf("\n::::::: Evaluating FLUCTUATIONS ---- L2 NORM");
    for(int i=0; i<ndata; ++i){
        var[i]=(tmean-data[i])*(tmean-data[i]);
        fprintf(outfluct,"%f\n",var[i]);
    }
    
    double meanfluct=0;
    for (int i=ifrom; i<ndata; ++i) {
        
        meanfluct=meanfluct+var[i];
        
    }
    
    meanfluct=sqrt(meanfluct/(ndata-1-ifrom));
    
    
    //evaluate peaks
    printf("\n::::::: Evaluating PEAKS  ");
    for (int i=ifrom; i<ndata; ++i) {
        if(sqrt(var[i])>meanfluct) peak[i]=1;
        
         //printf(" \n %i MF= %f  F(i)=%f %i ",i, meanfluct, sqrt(var[i]), peak[i]);
        
    }
    
 
    
    printf("\n::::::: Evaluating PEAKS DISTRIB  \n");
    
    int peaksnumber=0;
    int inpeak=0;
    int peakdim=0;
    
    if (peak[ifrom]==1) {
        ++peaksnumber;
        inpeak=1;
        ++peakdim;
    }
    
    for (int i=ifrom+1; i<ndata; ++i) {
        
        if (inpeak==1 && peak[i]==1) {
            ++peakdim;
        }
        
        else if (inpeak==1 && peak[i]==0){
            inpeak=0;
            plen[peaksnumber-1]=peakdim;
            peakdim=0;
            
        }
        
        else if (inpeak==0 && peak[i]==1){
            inpeak=1;
            ++peakdim;
            ++peaksnumber;
        }
        
        
    }
    
   for (int i=0; i<peaksnumber; ++i) {printf(" %i ", plen[i]);  }
    
    
    
    //BIN PEAKS HISTROGRAM
    
    
    
    // FLUCTUATION HISTOGRAM -----------------
    
    printf("\n::::::: Peak Histrogram");
    
    nbins=50;

    printf("\n::::::: Binning");
    
    //set min
    min=1000000;
    for (int i=0; i<peaksnumber; ++i) {
        if(plen[i]<min){min=plen[i];}
    }
    
    
    
    //set max
    max=-1;
    for (int i=0; i<peaksnumber; ++i) {
        if(plen[i]>max){max=plen[i];}
    }
    max=max+1;
    
    //set binwidth
    binw=(max-min)/nbins;
    
    printf("     MIN=%lf MAX=%lf  nbins=%i  binw=%lf", min,max,nbins,binw);
    
    //clean histogram
    for (int i=0; i<nbins; ++i) {
        histo[i]=0;
    }
    
    //binning
    for (int j=0; j<peaksnumber; ++j) {
        for (int i=0; i<nbins; ++i) {
            if( min+(i*binw) <= plen[j] && plen[j] < min+((i+1)*binw)) {++histo[i];
                //printf("\n dato %i (%f) in bin %i (%f - %f)", j, data[j], i, min+(i*binw),min+((i+1)*binw));
            }
        }
    }
    //normalize
    for (int i=0; i<nbins; ++i) {
        histo[i]=histo[i]/peaksnumber;
    }
    
    
    
    //print data
    for (int i=0; i<nbins; ++i) {
        //if(histo[i]!=0)
        {fprintf(outpeakhist,"%f %f\n", min+ (i+0.5)*binw, histo[i]);}
    }
    
    
    
    

    
    
    fclose(outpeakhist);
    
    

    
    // EVALUATE WINDOW MEAN AND FLUCTUATION ----------------------------------------------------------------------
    printf("\n::::::: RUNNING MEAN \n ------------ Window dimension = %i, ndata = %i", window,ndata);
    
    if(ndata<window) {printf(" ---> NOT COMPUTABLE \n \n ------ STOP. \n\n"); abort();}
    printf(" ---> COMPUTABLE ");
    
    cont= (ndata-window)-100;
    cont=0;
    
    do{
     
        mean=0;
        f=0;
        for(int i=0; i<window;++i){
            mean=mean+data[cont+i];
            f=f+var[cont+i];
        }
        mean=mean/window;
        f=f/window;
        f=sqrt(f);
     
        fprintf(outactimean,"%i %f\n",(window/2)+cont,mean);
        //fprintf(outwfluct, "%f\n", fabs(mean-tmean));
       
        fprintf(outwfluct,"%f \n",f);
        
       ++cont;
    
    } while(cont<ndata-window);
    
    printf(" ---- LAST MEAN ----  LAST FLUCTUATION:   %f \n \n",mean, f);
    printf("%i %f %f\n\n", window*2,mean,f);
    
    *eavg=mean;
    *efluct=f;
    
    
    
    indata.close();
    fclose(outactimean);
    fclose(outfluct);
    fclose(outwfluct);
    
    
    /*
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
     
     */
    
    fclose(outflucthist);
    fclose(outwflucthist);

}




int main(){
    
    double avg, f;
    
    int win=1000;
    
    FILE * totout1;
    FILE * totout2;
    
    char st[100];

    int maxs=1   ;
    int maxd=1   ;
    int maxt=1   ;
    
    double sr0=0.25000; double ds=0.0450;
    double d0=0.001; double dd=0.00045;
    
    
    
    
    
    int t=0;
    
    sprintf(folder,"data/D1");
    
    for (int s=0; s<maxs; ++s) {
        
        sprintf(st,"%s/lastavg-s-%i.txt",folder,s);
        totout1=fopen(st,"w");
        
        for (int d=0; d<maxd; ++d) {
            
            analysys(win, d, t, s, &avg, &f);
            
            fprintf(totout1, "%lf %i %lf %lf %lf \n", d0+dd*d, t, sr0+ds*s, avg, f);
            
            printf(" [ d = %i   |  t = %i  | s = %i ]\n ",d,t,s);
            
            
        }
        
        
        fclose(totout1);
   
    }
    
    
    
    

    for (int d=0; d<maxd; ++d)
    {
        sprintf(st,"%s/lastavg-d-%i.txt",folder,d);
        totout1=fopen(st,"w");
            for (int s=0; s<maxs; ++s) {
            
            analysys(win, d, t, s, &avg, &f);
            
            fprintf(totout1, "%lf %i %lf %lf %lf \n", d0+dd*d, t, sr0+ds*s, avg, f);
            
                printf(" [ d = %i   |  t = %i  | s = %i ]\n ",d,t,s);
            
            
        }
        
        
        fclose(totout1);
        
    }
    
    
    return 0;}





