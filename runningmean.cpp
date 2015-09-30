#include <math.h>
#include <iostream>
#include <fstream>

#if defined (WIN32)
#include <io.h>
#else
#include <unistd.h>
#endif


using namespace std;

                                                #define N 1000      // <<-------  WINDOW'S DIMENSION

int main(){

    FILE *indata;
    FILE *outdata;
    
                                                /* OUTPUT FILE */ outdata=fopen("pmedia.txt","w");
    
                                                /* INPUT FILE  */ indata=fopen("data/p7.txt","r");
    
    double data[N]; double temp;
    int cont;
    double mean;
    
    

    
    
    for(int i=0;i<N;++i){
       fscanf(indata,"%lf",&temp);
        data[i]=temp;
        
        //printf("\n %i %f ", i,temp);
    }
    
    mean=0;
    for(int i=0; i<N;++i){
        mean=mean+data[i];
    }
    mean=mean/N;

    fprintf(outdata,"%i %f\n",(N/2),mean);
    
 //   for (int i=0; i<N; ++i) {printf("%f ",data[i]);} printf("\n");
    
    cont=0;
    while (fscanf(indata, "%lf", &temp) == 1){
    
        for(int i=0; i<(N-1);++i){
            data[i]=data[i+1];
        }
        data[N-1]=temp;
        
       // for (int i=0; i<N; ++i) {printf("%f ",data[i]);}   printf("\n");
        
        mean=0;
        for(int i=0; i<N;++i){
            mean=mean+data[i];
        }
        mean=mean/N;
        ++cont;
        fprintf(outdata,"%i %f\n",(N/2)+cont,mean);
    
    
    }
    
    
    
    
    fclose(indata);
    fclose(outdata);


}









