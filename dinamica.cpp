

#include "dinamica.h"


extern igraph_t graph;
extern igraph_t sgraph;
extern igraph_matrix_t admatrix;
extern igraph_matrix_t density; extern igraph_matrix_t densityold;
extern igraph_matrix_t state;
extern igraph_matrix_t statenew, loadedstate;
extern igraph_matrix_t flux;
extern igraph_matrix_t bridgeslinks;
extern igraph_vector_t gain;
extern igraph_matrix_t loss;
extern igraph_vector_t statstate;
extern int ticks;
extern int tickstep;


extern int nodesnumber;

extern int vincolo;
extern int threshold;
extern int particles;
extern int maxtime;
extern int ispii;
extern int isrelaxed;
extern int havecstop;
extern int haveloadedstate;
extern int usesteady;
extern int isdissipating;

extern int isclustered;

extern int totrun;

extern FILE * output1;
extern FILE * output2;
extern FILE * output3;
extern FILE * output4;
extern FILE * output5;
extern FILE * output6;
extern FILE * output7;

extern LAURA_Histogram_1D histdist;
extern double dm, dv;
extern double hism, hisv;
 extern Fl_Round_Button *roundstop;
 extern Fl_Value_Input *inmeanstop;
 extern Fl_Value_Input *invarstop;
 extern Fl_Value_Input *inmeanerstop;
 extern Fl_Value_Input *invarerstop;
 
 
//------------------------------------------------------------------------------
void InitialStateTAS(int beginner, int useran, int maxrun){
    //inf useran==1 use random beginners
    int run;
    int ednum;
    
    
    
    nodesnumber=igraph_matrix_nrow(&admatrix);
    
    
    igraph_matrix_resize(&state, nodesnumber, maxrun);
    igraph_matrix_null(&state);
    igraph_matrix_resize(&loss, nodesnumber, maxrun);
    igraph_matrix_null(&loss);
    
    ednum=igraph_ecount(&sgraph);
    igraph_matrix_resize(&flux, ednum, maxrun);
    igraph_matrix_null(&flux);
    if(useran==1){init_genrand(time(NULL));}
    
    for(run=0;run<maxrun;++run ){
        if(useran==1){beginner=genrand_int31()%nodesnumber;}
        MATRIX(state,beginner,run)=particles;
        
    }
    
    
    //print_matrix(&state, stdout);
    
    igraph_matrix_update(&density,&state);
    igraph_matrix_scale(&density,1./particles);
    igraph_matrix_update(&densityold,&state);
    igraph_matrix_scale(&densityold,1./particles);
    
    igraph_matrix_update(&statenew,&state);
    
    ticks=0;
    
    //print_matrix(&state, stdout);printf("\nSTATENEW:\n"); print_matrix(&state, stdout); printf(" \n\n");
    
    
    init_genrand(time(NULL));
    
    haveloadedstate=0;
    
    
}




//------------------------------------------------------------------------------
void InitialStateRAN(int seed, int maxrun){
    
    int res;
    int selected;
    int run;
    int ednum;
    
    
    igraph_vector_t tempvector;
    
    init_genrand(seed);
    
    nodesnumber=igraph_matrix_nrow(&admatrix);
    
    
    igraph_matrix_resize(&state, nodesnumber, maxrun);
    igraph_matrix_null(&state);
    igraph_matrix_resize(&loss, nodesnumber, maxrun);
    igraph_matrix_null(&loss);
    
    ednum=igraph_ecount(&sgraph);
    igraph_matrix_resize(&flux, ednum, maxrun);
    igraph_matrix_null(&flux);
    
    
    
    res=particles;
    run=0;
    
    while (res>0) {
        selected=genrand_int31()%nodesnumber;
        ++MATRIX(state,selected,run);
        --res;
    }
    
    
    if(maxrun>1){
        
        igraph_vector_init(&tempvector,0);
        igraph_matrix_get_col(&state,&tempvector,0);
        
        for (run=1; run<maxrun; ++run) {
            igraph_matrix_set_col(&state,&tempvector,run);
        }
        
        igraph_vector_destroy(&tempvector);
        
    }
    
    
    
    //print_matrix(&state, stdout);
    
    igraph_matrix_update(&density,&state);
    igraph_matrix_scale(&density,1./particles);
    igraph_matrix_update(&densityold,&state);
    igraph_matrix_scale(&densityold,1./particles);
    
    igraph_matrix_update(&statenew,&state);
    
    ticks=0;
    
    
    // print_matrix(&state, stdout);printf("\nSTATENEW:\n"); print_matrix(&state, stdout); printf(" \n\n");
    
    init_genrand(time(NULL));
    
    haveloadedstate=0;
    
}


//------------------------------------------------------------------------------
void InitialStateLOAD(){
    
    int ednum;
    
    igraph_matrix_update(&state,&loadedstate);
    igraph_matrix_resize(&loss, nodesnumber, totrun);
    igraph_matrix_null(&loss);
    
    ednum=igraph_ecount(&sgraph);
    igraph_matrix_resize(&flux, ednum, totrun);
    igraph_matrix_null(&flux);
    
    particles=igraph_matrix_sum(&state)/totrun;
    
    igraph_matrix_update(&density,&state);
    igraph_matrix_scale(&density,1./particles);
    igraph_matrix_update(&densityold,&state);
    igraph_matrix_scale(&densityold,1./particles);
    
    igraph_matrix_update(&statenew,&state);
    
    ticks=0;
    
    if(usesteady==1){dm=0; dv=0;}
    
    init_genrand(time(NULL));
    
}


//------------------------------------------------------------------------------
void InitialStateSTAT(int maxrun){
    
    int run;
    int ednum;
    
    init_genrand(0);
    
    nodesnumber=igraph_matrix_nrow(&admatrix);
    
    
    igraph_matrix_resize(&state, nodesnumber, maxrun);
    igraph_matrix_null(&state);
    igraph_matrix_resize(&loss, nodesnumber, maxrun);
    igraph_matrix_null(&loss);
    
    ednum=igraph_ecount(&sgraph);
    igraph_matrix_resize(&flux, ednum, maxrun);
    igraph_matrix_null(&flux);
    
    
    for(run=0;run<maxrun;++run ){
        
        for (int p=0; p<particles; ++p) {
            double coin; int selnode;
            do{
                coin=genrand_real1();
                selnode=genrand_int31()%nodesnumber;
            }
            while(VECTOR(statstate)[selnode]<coin);
            
            ++MATRIX(state,selnode,run);
            
        }
        
        
    }
    
    
    //print_matrix(&state, stdout);
    
    igraph_matrix_update(&density,&state);
    igraph_matrix_scale(&density,1./particles);
    igraph_matrix_update(&densityold,&state);
    igraph_matrix_scale(&densityold,1./particles);
    
    igraph_matrix_update(&statenew,&state);
    
    ticks=0;
    
    //print_matrix(&state, stdout);printf("\nSTATENEW:\n"); print_matrix(&state, stdout); printf(" \n\n");
    
    
    init_genrand(time(NULL));
    
    haveloadedstate=0;
    
    
}





//------------------------------------------------------------------------------


void Evolution(double dt){
    
    igraph_vector_t neis;
    igraph_vector_init(&neis, 0);
    
    igraph_matrix_t normstate;
    igraph_matrix_init(&normstate, 0,0);
    igraph_vector_t avgnormstate;
    igraph_vector_init(&avgnormstate, 0);
    igraph_vector_t anothertemp;
    igraph_vector_init(&anothertemp, 0);
    
    int maxwheel=10000000;
    
    int run, i, j, k, repete, pinstate, checkparticles;
    int eid;
    int ederror;
    int tempint;
    
    double tempreal, dontmove;
    
    
    if (ticks!=0){
        nodesnumber=igraph_matrix_nrow(&admatrix);
        
        igraph_matrix_null(&flux);
        igraph_matrix_null(&loss);
        
        for(run=0;run<totrun;++run){
            igraph_vector_null(&gain);
            
            for(j=0;j<nodesnumber;++j){
                pinstate=MATRIX(state,j,run);
                
                if(pinstate!=0){
                    
                    igraph_neighbors(&sgraph, &neis, j, IGRAPH_ALL);
                    
                    
                    for(i=0;i<pinstate && i<vincolo;++i){
                        
                        if(ispii==1)dontmove=genrand_real1();
                        else dontmove=0;
                        
                        if(dontmove<1-(0.5*dt)){
                            
                            
                            
                            //°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°   WHEEEEEEEEEEEEL °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
                            
                            for(repete=0;repete<maxwheel;++repete)
                                
                            {
                                int pippo;
                                //_____seleziona potenziale nodo di arrivo
                                tempint = genrand_int31();
                                k=tempint%(igraph_vector_size(&neis));
                                
                                //printf("\n tempint=%i nodesnumber=%i, k=%i", tempint, nodesnumber,k);
                                
                                //_____lancia una moneta per vedere se arriva
                                tempreal=genrand_real1();
                                //printf("\n tempreal=%f, MATRIX(admatrix,j,k)=%f \n", tempreal, MATRIX(admatrix,j,k));
                                //printf("no \n");
                                pippo=VECTOR(neis)[k];
                                if ( tempreal<(MATRIX(admatrix,j,pippo)*dt)) {
                                    
                                    if(isdissipating==1){
                                        if(pippo==(nodesnumber-1) || j==(nodesnumber-1) || pinstate>threshold){//add to state vector
                                    ++MATRIX(loss,j,run); ++VECTOR(gain)[pippo];
                                    
                                    //get eids for flux
                                    igraph_get_eid(&sgraph, &eid,j,pippo,IGRAPH_UNDIRECTED,1);
                                    //printf("\n ederror:%i, %i->%i (%i) || ", ederror,j,pippo,eid);
                                    //printf("%f = %f + ((%i-%i)/abs(%i-%i))   ",MATRIX(flux,eid,run),MATRIX(flux,eid,run),j,pippo,j,pippo);
                                    MATRIX(flux,eid,run) = MATRIX(flux,eid,run) + (double)( (j-pippo) / abs(j-pippo) );
                                    
                                    //printf(" [%i]        %f    ", (j-pippo) / abs(j-pippo),MATRIX(flux,eid,run));
                                    
                                        break;}
                                        else {break;}
                                        
                                    }
                                    
                                    else if(pinstate>threshold){//add to state vector
                                        ++MATRIX(loss,j,run); ++VECTOR(gain)[pippo];
                                        
                                        //get eids for flux
                                        igraph_get_eid(&sgraph, &eid,j,pippo,IGRAPH_UNDIRECTED,1);
                                        //printf("\n ederror:%i, %i->%i (%i) || ", ederror,j,pippo,eid);
                                        //printf("%f = %f + ((%i-%i)/abs(%i-%i))   ",MATRIX(flux,eid,run),MATRIX(flux,eid,run),j,pippo,j,pippo);
                                        MATRIX(flux,eid,run) = MATRIX(flux,eid,run) + (double)( (j-pippo) / abs(j-pippo) );
                                        
                                        //printf(" [%i]        %f    ", (j-pippo) / abs(j-pippo),MATRIX(flux,eid,run));
                                        
                                        break;}
                                    
                                    else {break;}
                                    
                                }
                                
                            }
                            
                            
                            //°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
                        }
                        
                    }
                    
                }
                
                
            }
            
            for (j=0;j<nodesnumber;++j){
                MATRIX(statenew,j,run)=MATRIX(state,j,run)+VECTOR(gain)[j]-MATRIX(loss,j,run);
            }
            
            
            //CONTROLLO DI NON AVER PERSO PARTICELLE PER STRADA
            checkparticles=0;
            for (j=0;j<nodesnumber;++j){checkparticles=checkparticles+MATRIX(statenew,j,run);}
            if (checkparticles!=particles) { printf("\n \n   HO PERSO PARTICELLE PER STRADA! \n \n \n ");
                printf("RUN %i ---- # PARTICLES %i (* %i *) \n ", run, checkparticles, particles);
                printf(" \n\n state: \n\n "); print_matrix_ur(&state,stdout);
                printf(" \n\n loss: \n\n "); print_matrix_ur(&loss,stdout);
                printf(" \n\n statenew: \n\n "); print_matrix_ur(&statenew,stdout);
            
                abort();}
            
            
        }
        
        
        //salvo vecchia densità
        igraph_matrix_update(&densityold,&density);
        //stato nuovo ------> stato; azzero stato nuovo
        igraph_matrix_update(&state,&statenew);
        igraph_matrix_fill(&statenew, 0);
        
        
        //stato densità
        igraph_matrix_update(&density,&state);
        igraph_matrix_scale(&density,1./particles);
        
        
        //---------------------------------------------- ISTOGRAMMA!
        
        if(usesteady==1){
            igraph_matrix_update(&normstate,&state);
            
            //normalization
            igraph_matrix_scale(&normstate,1./particles);
            
            igraph_matrix_get_col(&normstate,&avgnormstate,0);
            igraph_vector_null(&avgnormstate);
            for(int i=0;i<totrun;++i){
                igraph_matrix_get_col(&normstate,&anothertemp,i);
                igraph_vector_div(&anothertemp,&statstate);
                igraph_vector_add(&avgnormstate,&anothertemp);
            }
            igraph_vector_scale(&avgnormstate,1./totrun);
            
            
            
            histdist.CreateFromArrayMinMax(anothertemp, 40,0.,4., 1);
            
           
        }
   
        
        
        
        
        
        
        //print_matrix(&state, stdout); printf(" \n\n");
    }
    
    
    ++ticks;
    ++tickstep;
    
    
     if(roundstop->value()==1
     && hisv < (invarstop->value()+invarerstop->value()) && hisv > (invarstop->value()- invarerstop->value() )
     && hism < (inmeanstop->value()+inmeanerstop->value()) && hism > (inmeanstop->value() -inmeanerstop->value()))
     {
         ++isrelaxed;
    
     
     }
    
    if (isrelaxed>NTRYSTOP) {
         maxtime=ticks;
    }
     
    igraph_vector_destroy(&neis);
    igraph_matrix_destroy(&normstate);
    igraph_vector_destroy(&anothertemp);
}
















void TurboRun(double dt){
    
    igraph_vector_t neis;
    igraph_vector_init(&neis, 0);
    
    double shooted;
    double dens, err;
    double totdens, toterr, totimerr;
    double totreldens, totrelerr, densrel, errrel;
    
    int maxwheel=10000000;
    
    int run, i, j, k, repete, pinstate, checkparticles;
    int eid;
    int ederror;
    int tempint;
    
    int time;
    
    double totactivity;
    
    double tempreal, dontmove;
    
    
    nodesnumber=igraph_matrix_nrow(&admatrix);
    
    
    
    igraph_matrix_t dissipation;
    igraph_matrix_init(&dissipation,nodesnumber,totrun);
    
    
    if(havecstop==1){isrelaxed=0;}
    
    
    for(time=0;time<maxtime;++time){
        
        
        igraph_matrix_null(&loss);
        igraph_matrix_null(&flux);
        
        igraph_matrix_null(&dissipation);
        
        
        for(run=0;run<totrun;++run){
            
            igraph_vector_null(&gain);
            
            for(j=0;j<nodesnumber;++j){
                
                pinstate=MATRIX(state,j,run);
                
                if(pinstate!=0){
                    
                    igraph_neighbors(&sgraph, &neis, j, IGRAPH_ALL);
                    
                    for(i=0;i<pinstate && i<vincolo;++i){
                        
                        if(ispii==1)dontmove=genrand_real1();
                        else dontmove=0;
                        
                        if(dontmove<1-(0.5*dt)){
                        
                            
                            
                            //°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°   WHEEEEEEEEEEEEL °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
                            
                            for(repete=0;repete<maxwheel;++repete)
                                
                            {
                                int pippo;
                                //_____seleziona potenziale nodo di arrivo
                                tempint = genrand_int31();
                                k=tempint%(igraph_vector_size(&neis));
                                
                                //printf("\n tempint=%i nodesnumber=%i, k=%i", tempint, nodesnumber,k);
                                
                                //_____lancia una moneta per vedere se arriva
                                tempreal=genrand_real1();
                                //printf("\n tempreal=%f, MATRIX(admatrix,j,k)=%f \n", tempreal, MATRIX(admatrix,j,k));
                                //printf("no \n");
                                pippo=VECTOR(neis)[k];
                                if ( tempreal<(MATRIX(admatrix,j,pippo)*dt)) {
                                    
                                    //dissipative case
                                    if( isdissipating==1 ){
                                        
                                        if(pippo==(nodesnumber-1) || j==(nodesnumber-1)){
                                        //add to state vector
                                        ++MATRIX(loss,j,run); ++MATRIX(dissipation,j,run);
                                        ++VECTOR(gain)[pippo];
                                    //get eids for flux
                                    igraph_get_eid(&sgraph, &eid,j,pippo,IGRAPH_UNDIRECTED,1);
                                    //printf("\n ederror:%i, %i->%i (%i) || ", ederror,j,pippo,eid);
                                    //printf("%f = %f + ((%i-%i)/abs(%i-%i))   ",MATRIX(flux,eid,run),MATRIX(flux,eid,run),j,pippo,j,pippo);
                                    MATRIX(flux,eid,run) = MATRIX(flux,eid,run) + (double)( (j-pippo) / abs(j-pippo) );
                                    //printf(" [%i]        %f    ", (j-pippo) / abs(j-pippo),MATRIX(flux,eid,run));
                                    break;
                                    }
                                        
                                        else if(pinstate>threshold){
                                            //add to state vector
                                            ++MATRIX(loss,j,run); ++VECTOR(gain)[pippo];
                                            //get eids for flux
                                            igraph_get_eid(&sgraph, &eid,j,pippo,IGRAPH_UNDIRECTED,1);
                                            //printf("\n ederror:%i, %i->%i (%i) || ", ederror,j,pippo,eid);
                                            //printf("%f = %f + ((%i-%i)/abs(%i-%i))   ",MATRIX(flux,eid,run),MATRIX(flux,eid,run),j,pippo,j,pippo);
                                            MATRIX(flux,eid,run) = MATRIX(flux,eid,run) + (double)( (j-pippo) / abs(j-pippo) );
                                            //printf(" [%i]        %f    ", (j-pippo) / abs(j-pippo),MATRIX(flux,eid,run));
                                            break;
                                        }
                                        
                                        
                                        else {break;}
                                    }
                                    
                                    else if( pinstate>threshold ){
                                        //add to state vector
                                            ++MATRIX(loss,j,run); ++VECTOR(gain)[pippo];
                                        //get eids for flux
                                        igraph_get_eid(&sgraph, &eid,j,pippo,IGRAPH_UNDIRECTED,1);
                                        //printf("\n ederror:%i, %i->%i (%i) || ", ederror,j,pippo,eid);
                                        //printf("%f = %f + ((%i-%i)/abs(%i-%i))   ",MATRIX(flux,eid,run),MATRIX(flux,eid,run),j,pippo,j,pippo);
                                        MATRIX(flux,eid,run) = MATRIX(flux,eid,run) + (double)( (j-pippo) / abs(j-pippo) );
                                        //printf(" [%i]        %f    ", (j-pippo) / abs(j-pippo),MATRIX(flux,eid,run));
                                        break;
                                        
                                    }
                                    
                                    else {break;}

                                    
                                }
                                
                            }
                            
                            
                            //°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
                        }
                        
                    }
                    
                }
                
                
            }
            
            for (j=0;j<nodesnumber;++j){
                MATRIX(statenew,j,run)=MATRIX(state,j,run)+VECTOR(gain)[j]-MATRIX(loss,j,run);
            }
            
            
            //CONTROLLO DI NON AVER PERSO PARTICELLE PER STRADA
            checkparticles=0;
            for (j=0;j<nodesnumber;++j){checkparticles=checkparticles+MATRIX(statenew,j,run);}
            if (checkparticles!=particles) { printf("\n \n   HO PERSO PARTICELLE PER STRADA! \n \n \n "); abort();}
            //printf("RUN %i ---- # PARTICLES %i \n ", run, checkparticles);
            
            
        }
        
        
        
        //salvo vecchia densità
        igraph_matrix_update(&densityold,&density);
        
        //stato nuovo ------> stato; azzero stato nuovo
        igraph_matrix_update(&state,&statenew);
        igraph_matrix_fill(&statenew, 0);
        
        //stato densità
        igraph_matrix_update(&density,&state);
        igraph_matrix_scale(&density,1./particles);
        
        
        
        //print_matrix(&state, stdout); printf(" \n\n");
        
        
        //if have steady state
        if(usesteady==1){
            
            igraph_matrix_t activation;
            igraph_matrix_init(&activation,nodesnumber,totrun);
            igraph_matrix_null(&activation);
            
            igraph_vector_t correlation;
            igraph_vector_init(&correlation,(nodesnumber*nodesnumber)); igraph_vector_null(&correlation);

            
            fprintf(output1,"%i ", time);
            fprintf(output2,"%i ", time);
            fprintf(output5,"%i ", time);
            fprintf(output6,"%i ", time);
            
            totdens=0; toterr=0;   totactivity=0;
            for(int i=0;i<nodesnumber;++i){
    
                shooted=0;
                dens=0; densrel=0;
                err=0; errrel=0;
                
                for(int j=0; j<totrun; ++j){
                    dens=dens+MATRIX(density,i,j);
                    densrel=densrel+MATRIX(density,i,j)/VECTOR(statstate)[i];
                    
                    err=err+((VECTOR(statstate)[i]-MATRIX(density,i,j))*(VECTOR(statstate)[i]-MATRIX(density,i,j)));
                    errrel=errrel+( 1 - MATRIX(density,i,j)/VECTOR(statstate)[i] )*( 1 - MATRIX(density,i,j)/VECTOR(statstate)[i] );
                    
                    shooted=shooted+MATRIX(loss,i,j);
                    
                    if(MATRIX(loss,i,j)!=0 && MATRIX(loss,i,j)!=MATRIX(dissipation,i,j)){++MATRIX(activation,i,j);}
                    
                }
                
                dens=dens/totrun; densrel=densrel/totrun;
                err=sqrt(err)/totrun; errrel=errrel/totrun;
                
                totdens=totdens+dens; totreldens=totreldens+densrel;
                toterr=toterr+err; totrelerr=totrelerr+errrel;
                
                shooted=shooted/totrun;
                
                
                fprintf(output1,"%f ",dens);
                fprintf(output2,"%f ",err);
                fprintf(output5,"%f ",shooted);
                
            }
            
            
            
            totdens=totdens/nodesnumber;
            toterr=toterr/nodesnumber;
            totreldens=totreldens/nodesnumber;
            totrelerr=totrelerr/nodesnumber;
            
            totactivity=igraph_matrix_sum(&activation);
            totactivity=totactivity/(nodesnumber*totrun);
            
            fprintf(output1,"%f %f ",totdens, totreldens);
            fprintf(output2,"%f %f ",toterr, totrelerr);
            fprintf(output7,"%f",totactivity);
            
            
            // ---- CORRELATION ---
            
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
                //fprintf(output6,"%f ",VECTOR(correlation)[i]);
            }
            
            
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
                
                
            
                //i evaluate the distance (err)
                for(int i=0; i<nodesnumber; ++i)
                {
                    //L1 DISTANCE
                    MATRIX(distl1,i,j)=fabs((VECTOR(statstate)[i]-MATRIX(density,i,j)));
                    //L1 DISTANCE WRT OLD DENSITY & DENSITY
                    MATRIX(distimel1,i,j)=fabs(MATRIX(densityold,i,j)-MATRIX(density,i,j));

                }
                
            }
            
            igraph_matrix_rowsum(&distl1,&rundistl1); igraph_matrix_rowsum(&distimel1,&rundistimel1);
            igraph_vector_scale(&rundistl1,(1./nodesnumber)); igraph_vector_scale(&rundistimel1,(1./nodesnumber));
            
            toterr=  (double)( igraph_vector_sum(&rundistl1) )  /   (double)totrun ;
            totimerr= (double)(igraph_vector_sum(&rundistimel1))/  (double)totrun;
            
            igraph_vector_destroy(&rundistl1); igraph_vector_destroy(&rundistimel1);
            igraph_matrix_destroy(&distl1); igraph_matrix_destroy(&distimel1);
            
            fprintf(output2,"%f %f",toterr, totimerr);
            
            fprintf(output1,"\n");
            fprintf(output2,"\n");
            fprintf(output5,"\n");
            fprintf(output6,"\n");
            fprintf(output7,"\n");
            
            
            
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
            
            
            
            
            //--------- IF I HAVE THE STOP CONDITION
            if(havecstop==1){
            
                if(roundstop->value()==1
                   && totrelerr < (invarstop->value()+invarerstop->value()) && totrelerr > (invarstop->value()- invarerstop->value() )
                   && totreldens < (inmeanstop->value()+inmeanerstop->value()) && totreldens > (inmeanstop->value() -inmeanerstop->value()))
                {
                    ++isrelaxed;
                    
                    
                }
                
                if (isrelaxed>NTRYSTOP) {
                    
                    printf("\n\n\n\n ---- has relaxed ----");
                    break;
                }

            
            }
            
            igraph_matrix_destroy(&activation);
            igraph_vector_destroy(&correlation);
            
        }
        
        //if i HAVENT STEADY STATE
        else {
            
            fprintf(output1,"%i ", time);
            fprintf(output5,"%i ", time);
            
            for(int i=0;i<nodesnumber;++i){
                dens=0;
                shooted=0;
                for(int j=0; j<totrun; ++j){
                    dens=dens+MATRIX(density,i,j);
                    shooted=shooted+MATRIX(loss,i,j);
                }
                dens=dens/totrun;
                shooted=shooted/totrun;
                
                fprintf(output1,"%f " ,dens);
                fprintf(output5,"%i ", shooted);
            }
            
            fprintf(output1,"\n");
            fprintf(output5,"\n");
            
        }
        
        
        
    }
    
    igraph_vector_destroy(&neis);
    
    igraph_matrix_destroy(&dissipation);
    
    printf("       ISRELAXED #=%i \n \n \n", isrelaxed);
    
}





