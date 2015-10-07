/*
 
 @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 @                                  @
 @  LAURA                           @
 @                  v 0.1.2         @
 @                                  @
 @     Last Update  23 JULY 2015    @
 @                                  @
 @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 
 * By Elena Tea RUSSO
 * Physics MSc student @ UNIBO, Italy – Physics Department.
 
 * For bug reporting, help and complaint: "elentatea.russo at gmail.com"
 
 Creative Common Licensed.
 LAURA is free and opensource.
 
 
 ******************* SOME COMPILING / DEPENDENCES INFO ******************************************

 LAURA usess "FLTK" and "IGRAPH" C/C++ libraries.
 LAURA has been developed using MACOSX Mavericks and also tested on WINDOWS 8 with Visual Studio. Thus, it should compile and work properly under those OSs / IDE. Most common issues are due to libraries referncing problems of igraph or FLTK.
 LAURA is compiled by its MAKEFILE or Visual Studio project file. Please check libraries' locations before compiling.

 
 ************************************************************************************************
 
 
 *
 */


#include <FL/Fl.H>

#include "form.h"
#include "Frame.h"
#include "logdebug.h"
#include "igraphutil.h"
#include "randomgen.h"
#include "dinamica.h"
#include "dataanalysis.h"
#include "networkop.h"

#if defined (WIN32)
#include <io.h>
#include "igraph/igraph.h"
#else
#include <unistd.h>
#include "igraph.h"
#endif

#include <FL/Fl_Text_Display.H>
#include <FL/Fl_Light_Button.H>
#include <FL/Fl_Check_Button.H>
#include <FL/Fl_Box.H>

// FLTK FRAMES / WINDOWS / other
extern Frame   *scene;
extern DataFrame   *datascene;
extern FluctFrame   *fluctscene;
extern int framexmin, framexmax, frameymin, frameymax;
extern Fl_Light_Button *runbutton;



    /*
    /*  MAIN VARIABLES, MOSTLY USED AS "EXTERNAL" FROM THE REST OF THE MODULES
    /*
     */

//GRAPHS:
igraph_t graph;     //this is the graph where I work
igraph_t sgraph;    //this is a SIMPLIFIED graph; used by the dynamical routines

//LAYOUT
igraph_matrix_t layout;

//DYNAMICAL MATRICES/VECTORS
                            //(THOSE ARE MATRICES BECAUSE OF POSSIBLE MULTIPLE RUN)
igraph_matrix_t state;      //the state of the system (this is an array of actual integers, representing particles population)
igraph_matrix_t density, densityold;    //density-state of the system
igraph_matrix_t statenew;   //new state (computed during the evolution function
igraph_matrix_t loadedstate;//state loaded when LAURAFILE is loaded
igraph_matrix_t flux;       //fluxes on the edges
igraph_matrix_t bridgeslinks;//links
igraph_vector_t gain; igraph_matrix_t loss; //"gain"&"loss" vector (how many particles the i-th node gains/loses).
igraph_vector_t statstate;  //stationary state of the network (eigenvec of eigenval 1)

igraph_matrix_t dissipation;


igraph_matrix_t admatrix;   //ADJACENCE MATRIX

igraph_matrix_t estates;    //EIGENSTATES OF ADMATRIX

int nodesnumber;

int ticks;      //time elapsed
int tickstep;   //time espales in a single "step" a steprun
double deltat;  //delta t (makes it faster or slower if "ispii" (i.e. loops) is activated

//OUTPUTS
FILE * output1;
FILE * output2;
FILE * output3;
FILE * output4;
FILE * output5;
FILE * output6;
FILE * output7;


     /*
     /*  BOOLEANS VARIABLES, MOSTLY USED AS "EXTERNAL" FROM THE REST OF THE MODULES
     /*
      */

int drawingcontrol; //draw or not the scenes
int drawhisto;      //draw or not the histogram
int drnodes;        //draw nodes
int drlinks;        //draw links
int drfluxes;       //draw fluxes on links

int runningcontrol; // 1-> simulation is running, 0-> simulation is not running
int amstepping;     //stepping (smal steprun)
int step;           //how long is a steprun
int cleared;        //everything is cleared ( see void clear() in form.h )
int graphisloaded;  //do I have a graph?
int newgraph;       //i am working on a new graph?
int islattice;      //is lattice graph?
int istoro;         //is (lattice) toroidal?
int islatrand;      //is lattice with random weights?
int islatunif;      //is lattice with uniform weights?
int islatpot;       //is lattice made from a 2D potential?
int israndomER1;    //is random graph? (erdosh reny)
int isclustered;    //is a clustered graph?
int isclustger2l;   //is a clustered gerarchic graph (2 layers)
int usesteady;      //use (i.e. i have computed) steady state
int haveloadedstate;//have i got a loaded state?
int ispii;          //is pii (i.e., i have a p_ii=1/2*dt (prob of self loops)
int havecstop;      //do i have a stop condition?
int isrelaxed;      //has my system relaxed according to the stop condition?
int isdissipating;  //do I have a dissipation node?

int havepath;       //do I have a path for my job?
int haveout;        //do I have a working output file?

int rewrite;        //mainidlecb (fltk event routine) has to rewrite something?
int error;          //some error occured

int viewingdata;    //is the fluctuation window open?



    /*
    /* STRINGS, MOSTLY USED AS "EXTERNAL" FROM THE REST OF THE MODULES
    /*
     */

char path[300];         // path of my job (to load/save)
char errorstring[100];  //string containing the ERROR MESSAGE dispayed



    /*
    /* DATA ANALYSIS (from dataanalysys)
    /*
    */

LAURA_Histogram_1D histdist;
double dm, dv;




    /*
    /*  BUFFERS AND OTHER, extern from form.h
    /*
     */

extern Fl_Text_Buffer      *databuff;
extern Fl_Text_Buffer      *pathbuff;
extern Fl_Text_Buffer      *tickbuff;

extern Fl_Light_Button *drawbutton;
extern Fl_Check_Button *printdatabutton;
extern Fl_Check_Button *printcorrbutton;


/*
 /*  OTHER VARIABLES, MOSTLY USED AS "EXTERNAL" FROM THE REST OF THE MODULES
 /*
 */

//lattice characteristics
int latticedim;
int latticeside;


//presets
int vincolo=10000000;
int threshold=1;
int particles=10000;
int beginner;
int totrun=1;
int maxtime=1000;





/*__________________________________________________________________________ MAIN CYCLE (FLTK CYCLE) ____*/
void mainidle_cb(void*){    //this routine updates the program.
                            //thus, it computes the EVOLUTION
    
    double shooted;
    double dens, err;
    double totdens, toterr, totimerr;
    char s[100];
    
    
    
    
    
    
    // ---- running controls AND PRINTING
    if(
       (amstepping==0 && runningcontrol==1 && graphisloaded==1 && ticks<=maxtime ) ||
       (amstepping==1 && runningcontrol==1 && graphisloaded==1 && ticks<=maxtime &&  tickstep<=step-1)
       )
    
    {
  
        //EVOLUTION STEP
        Evolution(deltat);
        
        
        //PRINTS
        if((int)printdatabutton->value()==1){
          
            
            
            
            
            igraph_matrix_t activation;
            igraph_matrix_init(&activation,nodesnumber,totrun);
            igraph_matrix_null(&activation);
            
          
            //---------------------------------------------------------------------- if have steady state
            if(usesteady==1){
                
                fprintf(output1,"%i ", ticks);
                fprintf(output2,"%i ", ticks);
                fprintf(output5,"%i ", ticks);

                
                totdens=0; toterr=0;
                for(int i=0;i<nodesnumber;++i){
                    shooted=0;
                    dens=0;
                    err=0;
                    for(int j=0; j<totrun; ++j){
                        dens=dens+MATRIX(density,i,j);
                        err=err+((VECTOR(statstate)[i]-MATRIX(density,i,j))*(VECTOR(statstate)[i]-MATRIX(density,i,j)));
                        shooted=shooted+MATRIX(loss,i,j);

                        if(MATRIX(loss,i,j)!=0 && MATRIX(loss,i,j)!=MATRIX(dissipation,i,j)){++MATRIX(activation,i,j);}
                    }
                    
                    dens=dens/totrun;
                    err=sqrt(err)/totrun;
                    shooted=shooted/totrun;
                    totdens=totdens+dens;
                    toterr=toterr+err;
                    
                    fprintf(output1,"%f ",dens);
                    fprintf(output2,"%f ",err);
                    fprintf(output5,"%f ",shooted);
                }
                
                totdens=totdens/nodesnumber;
                toterr=toterr/nodesnumber;
                fprintf(output1,"%f ",totdens);
                fprintf(output2,"%f ",toterr);
                
              
                
                
                
                
                //calculate error on run
                
                
                igraph_matrix_t distl1; igraph_matrix_t distimel1;
                igraph_matrix_init(&distl1,nodesnumber,totrun); igraph_matrix_init(&distimel1,nodesnumber,totrun);
                igraph_matrix_null(&distl1); igraph_matrix_null(&distimel1);
                
                igraph_vector_t rundistl1; igraph_vector_t rundistimel1;
                igraph_vector_init(&rundistl1,totrun); igraph_vector_init(&rundistimel1,totrun);
                igraph_vector_null(&rundistl1); igraph_vector_null(&rundistimel1);
                
                toterr=0; totimerr=0;
                //for every run
                for(int j=0;j<totrun;++j){
                    
                    //i evaluate the distance between the state and the stationary state (toterr) and the distance between old and new density (totimerr)
                    for(int i=0; i<nodesnumber; ++i)
                    {
                        //L1 DISTANCE WRT STATSTATE & DENSITY
                        MATRIX(distl1,i,j)=fabs(VECTOR(statstate)[i]-MATRIX(density,i,j));
                        
                        //L1 DISTANCE WRT OLD DENSITY & DENSITY
                        MATRIX(distimel1,i,j)=fabs(MATRIX(densityold,i,j)-MATRIX(density,i,j));
                    }
                    
                }
                
                igraph_matrix_rowsum(&distl1,&rundistl1); igraph_matrix_rowsum(&distimel1,&rundistimel1);
                igraph_vector_scale(&rundistl1,(1./nodesnumber)); igraph_vector_scale(&rundistimel1,(1./nodesnumber));
                
                toterr=  (double)( igraph_vector_sum(&rundistl1) )  /   (double)totrun ;
                totimerr= (double)( igraph_vector_sum(&rundistimel1)) / (double)totrun;
                
                igraph_vector_destroy(&rundistl1); igraph_vector_destroy(&rundistimel1);
                igraph_matrix_destroy(&distl1); igraph_matrix_destroy(&distimel1);
                
                fprintf(output2,"%f %f",toterr, totimerr);
                
                fprintf(output1,"\n");
                fprintf(output2,"\n");
                fprintf(output5,"\n");

                
                
                
                //if i have BRIDGES ("clustered" graph), I print the traffic on the BRIDGES, using "output3" file
                if(isclustered==1){
                    
                    //for each bridge
                    fprintf(output3,"%i ",ticks);
                    for(int nbri=0; nbri<igraph_matrix_ncol(&bridgeslinks); ++nbri){
                        for(int i=0; i<igraph_matrix_nrow(&bridgeslinks);++i){
                            double tfl=0;
                            for(int j=0; j<totrun; ++j){
                                int beid;
                                beid=(int)MATRIX(bridgeslinks,i,nbri);
                                tfl=tfl+MATRIX(flux,beid,j);
                            }
                            fprintf(output3,"%f ",tfl);
                        }
                        
                        
                    }
                    fprintf(output3,"\n");
                }
                
                
                
                
                
            }
            
            //if i HAVENT STEADY STATE
            else {
                
                fprintf(output1,"%i ", ticks);
                fprintf(output5,"%i ", ticks);
                
                for(int i=0;i<nodesnumber;++i){
                    shooted=0;
                    dens=0;
                    for(int j=0; j<totrun; ++j){
                        dens=dens+MATRIX(density,i,j);
                        shooted=shooted+MATRIX(loss,i,j);
                        if(MATRIX(loss,i,j)!=0 && MATRIX(loss,i,j)!=MATRIX(dissipation,i,j)){++MATRIX(activation,i,j);}
                    }
                    dens=dens/totrun;
                    shooted=shooted/totrun;
                    
                    fprintf(output1,"%f " ,dens);
                    fprintf(output5,"%f " ,shooted);
                }
                
                fprintf(output1,"\n");
                fprintf(output5,"\n");
                
            }
            
            
            
            // ---------------------------------- CORRELATION ---
            
  
            if((int)printcorrbutton->value()==1) {
                
                
               // print_matrix_ur(&activation,stdout);
                
                
            fprintf(output6,"%i ", ticks);
            
                
            igraph_vector_t correlation;
            igraph_vector_init(&correlation,(nodesnumber*nodesnumber));
            igraph_vector_null(&correlation);
            
            igraph_vector_t meanactivation;
            igraph_vector_init(&meanactivation,nodesnumber);
            igraph_vector_null(&meanactivation);
                
                
            //calculate mean activation
            for (int j=0; j<totrun; ++j) {
                for (int i=0; i<nodesnumber; ++i) {
                    VECTOR(meanactivation)[i]=VECTOR(meanactivation)[i]+MATRIX(activation,i,j);
                }
            }
            igraph_vector_scale(&meanactivation,1./totrun);
            
            //calculate actual correlation
            for (int x=0; x<nodesnumber ; ++x) {
                for(int y=0; y<nodesnumber; ++y){
                    
                    double prod=0;
                    for (int j=0; j<totrun; ++j) {
                        prod=prod+ ( MATRIX(activation,x,j)*MATRIX(activation,y,j) );
                    }
                    prod=prod/totrun;
                    
                    VECTOR(correlation)[(x*nodesnumber+y)] = prod - (VECTOR(meanactivation)[x]*VECTOR(meanactivation)[y]);
                    
                }
            }
            
            igraph_vector_destroy(&meanactivation);
            
            for (int i=0; i<(nodesnumber*nodesnumber); ++i) {
                fprintf(output6,"%f ",VECTOR(correlation)[i]);
            }
            
            
            fprintf(output6,"\n");
            
            
            
            
            igraph_vector_destroy(&correlation);
            }
            
            
            
            
            
            igraph_matrix_destroy(&activation);
            
            
            
            
        }
    }
    
    
    if((ticks==maxtime || tickstep==step) && runningcontrol==1){
        run();
    }
    
    
    // ---- no graph loaded
    if(graphisloaded==0 && rewrite==1) {
        runbutton->deactivate();
        sprintf(s,"No\nnetwork\nloaded");
        databuff->text(s);
    }
    
    
    // ---- graph loaded
    else {
        runbutton->activate();
        
        if(islattice==1 && rewrite==1){
            if(istoro==1){
                sprintf(s,"Nodes=%i\nToroidal\nLattice\n%iD Side=%i",nodesnumber, latticedim, latticeside);
            }
            else{
                sprintf(s,"Nodes=%i\nLattice\n%iD Side=%i",nodesnumber, latticedim, latticeside);
            }
            
           databuff->text(s);
            
        }
        else if(rewrite==1){
            sprintf(s,"Nodes=%i",nodesnumber);
           databuff->text(s);
        }
        
    }
    
    
    
    //have path
    if(havepath==1 && rewrite==1){pathbuff->text(path); }
    else if(havepath==0 && rewrite==1){pathbuff->text("No Path");}
    
    
    if(error==1 && rewrite==1){
        sprintf(s,errorstring);
       databuff->text(s);
    }
    
    
    if (ticks<=maxtime){
    scene->redraw();
    datascene->redraw();
    if(viewingdata==1){fluctscene->redraw();}
    }
    
 
    rewrite=0;
    
    
    
}


//-------------------------------------------------------------------------------------------------
int main(int argc, char **argv) {
    
    igraph_set_error_handler(igraph_error_handler_printignore);
        igraph_i_set_attribute_table(&igraph_cattribute_table);
    
    
    
    drawingcontrol=1;
	runningcontrol=0;
    cleared=1;
    graphisloaded=0;
	newgraph=1;
    usesteady=0;
    ispii=0;
    haveout=0;
    isrelaxed=0;
    havecstop=0;
    isdissipating=0;
    viewingdata=0;
    
    
    //RANDOM NUMBER GENERATOR INITIALISE
    init_genrand(0); //srand(3);
    
    
    //preloaded graph:
    {
        latticedim=2;
        latticeside=5;
        istoro=1;
        islattice=1;
        isclustered=0;
        israndomER1=0;
        generatelattice(2,5,1,0,0,0,0);
        //inizializations
        igraph_matrix_init(&state, 0, 0);
        igraph_matrix_init(&density, 0, 0);
        igraph_matrix_init(&densityold,0,0);
        igraph_matrix_init(&statenew, 0, 0);
        igraph_matrix_init(&flux, 0, 0);
        igraph_matrix_init(&loss, 0, 0);
        
        igraph_matrix_init(&dissipation,0,0);
        
        beginner=(nodesnumber/2);
        InitialStateTAS(beginner,0,totrun);
        igraph_vector_init(&gain,nodesnumber);
        graphisloaded=1;
        newgraph=1;
        haveloadedstate=0;
        
        
        drnodes=1;
        drlinks=1;
        drfluxes=0;
    }
    
    
    
    
    havepath=0;
    sprintf(path,"No path");
  
	
    logDebug("\n DEFAULT ADJ MATRIX:\n");
    print_matrix_ur(&admatrix,stdout);
    printf("\n");
	printf("\n");
    rewrite=1;
    
    CreateMyWindow();
    DataWindow();
    
    
  Fl::add_idle(mainidle_cb, 0);
  Fl::run();
    
    
    histdist.Clear();
    
    igraph_destroy(&graph);
    igraph_destroy(&sgraph);
    igraph_matrix_destroy(&admatrix);
    igraph_matrix_destroy(&layout);
    igraph_matrix_destroy(&density);
    igraph_matrix_destroy(&densityold);
    igraph_matrix_destroy(&state);
    igraph_matrix_destroy(&loadedstate);
    igraph_matrix_destroy(&statenew);
    igraph_matrix_destroy(&flux);
    igraph_vector_destroy(&gain);
    igraph_matrix_destroy(&loss);
    igraph_matrix_destroy(&dissipation);
    igraph_vector_destroy(&statstate);
    igraph_matrix_destroy(&estates);
    
    closeout();
    
  return 0;
}
//-------------------------------------------------------------------------------------------------
