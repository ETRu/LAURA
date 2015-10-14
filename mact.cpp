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

#define TNODES 26


void analysys(int n, int eps, int th){
    
    char inpath[100];
    char outpath[100];
    
    sprintf(inpath, "lattice/N%ieps0%ith%i1.txt",n,eps,th);
    sprintf(outpath, "lattice/ssavg.txt");
    
    ifstream indata;                    /* INPUT FILE  */ indata.open(inpath);
    
    FILE * outdata;
    outdata=fopen(outpath,"a");
    
    
    //calculate average
    
    int cont=0;
    double avg=0;
    double temp;
    
    
    do{
        double prevsum=0;
        for (int i=0; i<TNODES; ++i) {
            indata >> temp;
            //printf(" %i %f \n", i, temp);
            if(i!=0){prevsum=prevsum+temp;}
            if(indata.eof()){break;}
            
        }
        
        indata >> temp;
        
        avg=avg+temp;
        
        for (int i=100; i<102; ++i) {
            indata >> temp;
            //printf(" %i %f \n", i, temp);
            if(indata.eof()){break;}
            
        }
        
        
        // printf("\n prevsum = %f , avg = %f\n\n", prevsum, avg);
        
        
        ++cont;}
    while(!indata.eof());
    
    
    indata.close();
    
    avg=avg/cont;
    
    fprintf(outdata,"%i %i %f %i \n", n, eps, avg, th);
    
    
    fclose(outdata);
    
    
}




int main(){
    
    
    
    analysys(100,1,0);
    analysys(100,2,0);
    analysys(100,3,0);
    analysys(100,4,0);
    analysys(100,5,0);
    analysys(100,6,0);
    analysys(100,7,0);
    analysys(100,8,0);
    analysys(100,9,0);
    analysys(100,10,0);
    
    
    
    analysys(100,1,1);
    analysys(100,2,1);
    analysys(100,3,1);
    analysys(100,4,1);
    analysys(100,5,1);
    analysys(100,6,1);
    analysys(100,7,1);
    analysys(100,8,1);
    analysys(100,9,1);
    analysys(100,10,1);
    
    
    return 0;}





