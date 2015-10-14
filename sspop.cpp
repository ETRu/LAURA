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
    
    char inpath[100]; char inpath1[100];
    char outpath[100]; char outpath1[100];
    
    sprintf(inpath, "lattice/N%ieps0%ith%i1.txt",n,eps,th);
    sprintf(inpath1, "lattice/N%ieps0%ith%i7.txt",n,eps,th);
    sprintf(outpath, "lattice/ssavg.txt");
    sprintf(outpath1, "lattice/actm.txt");
    
    ifstream indata; indata.open(inpath);
    ifstream indata1; indata1.open(inpath1);
    
    FILE * outdata;
    outdata=fopen(outpath,"a");
    
    FILE * outdata1;
    outdata1=fopen(outpath1,"a");
    

    
    int cont=0;
    double avg=0; double avg1=0;
    double temp;
    
    //---------------------------------------------------------------------------------------- sinksourcenode
    
    
    avg=0; cont=0;
    
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
    

    
    
    //---------------------------------------------------------------------------------------- activation
    
    
    avg1=0; cont =0;
    
    do{

       
        indata1 >> temp;
        
        avg1=avg1+temp;
        
        ++cont;}
    while(!indata1.eof());

    avg1=avg1/cont;
    
        indata1.close();
    
    
    
    
        fprintf(outdata,"%i %i %f %f %i \n", n, eps, avg, avg1, th);
    
    
    fclose(outdata);
    fclose(outdata1);
    
}




int main(){
    
    /*
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
*/
    
    analysys(25,1,0);
    analysys(50,1,0);
    analysys(75,1,0);
        analysys(100,1,0);
        analysys(125,1,0);
        analysys(150,1,0);
        analysys(175,1,0);
        analysys(200,1,0);
        analysys(250,1,0);
        analysys(300,1,0);
    
    analysys(25,1,1);
    analysys(50,1,1);
    analysys(75,1,1);
    analysys(100,1,1);
    analysys(125,1,1);
    analysys(150,1,1);
    analysys(175,1,1);
    analysys(200,1,1);
    analysys(250,1,1);
    analysys(300,1,1);
    
    return 0;}





