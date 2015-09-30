#include "form.h"


extern int nodesnumber;
extern igraph_matrix_t admatrix;
extern igraph_t graph;
extern igraph_matrix_t layout;
extern igraph_matrix_t density, densityold, state, statenew, loadedstate;
extern igraph_vector_t gain;
extern igraph_matrix_t loss;

extern int latticedim;
extern int latticeside;

extern int vincolo;
extern int threshold;
extern int particles;
extern int beginner;
extern int totrun;
extern int maxtime;
extern double deltat;

extern int ticks;
extern int tickstep;

extern FILE * output1;
extern FILE * output2;
extern FILE * output3;
extern FILE * output4;
extern FILE * output5;
extern FILE * output6;
extern FILE * output7;

extern int framexmin, framexmax, frameymin, frameymax;

//----- CONTROL VARIABLES -  FROM MAIN.CPP

extern int drawingcontrol;
extern int drawhisto;
extern int drnodes;
extern int drlinks;
extern int drfluxes;

extern int runningcontrol;
extern int amstepping; extern int step;
extern int graphisloaded;
extern int cleared;
extern int islattice;
extern int istoro;
extern int islatrand;
extern int islatunif;
extern int islatpot;
extern int israndomER1;
extern int isclustered;
extern int isclustger2l;
extern int newgraph;
extern int usesteady;
extern int ispii;
extern int havecstop;
extern int isrelaxed;
extern int haveloadedstate;

extern int rewrite;

extern int havepath;
extern int haveout;
extern char path[100];


extern int error;
extern char errorstring[100];


//main window defines-----------------------
#define BORDER1 20
#define SCREEN_WIDTH  500
#define SCREEN_HEIGHT 500
#define COLWIDTH 120
#define COLWIDTHLARGE 200
#define LEFT_SPACE COLWIDTH*2+2*BORDER1
#define RIGHT_SPACE COLWIDTH+COLWIDTHLARGE+3*BORDER1

#define BUTTON_WL 100
#define BUTTON_WR 100
#define BUTTON_L 200
#define BUTTON_H 50
#define BUTTON_H1 30


Fl_Window    *form;
Fl_Group    *mainwindow;
Frame        *scene;

igraph_matrix_t oldlayout;


// button in left column 1
Fl_Group *settingsgroup;
Fl_Input *setoutname;
Fl_Value_Input *in1;
Fl_Value_Input *in2;
Fl_Value_Input *inthreshold;
Fl_Value_Input *in3;
Fl_Value_Input *indt;
Fl_Value_Input *in4;
Fl_Light_Button *activatesteady;
Fl_Round_Button *pii;
Fl_Group *initstategroup;
Fl_Round_Button *roundTAS;
Fl_Check_Button *roundTASRAN;
Fl_Round_Button *roundRAN;
Fl_Round_Button *roundSTAT;
Fl_Value_Input *in5;
Fl_Value_Input *in6;

// button in left column 2
Fl_Text_Buffer  *tickbuff;
Fl_Text_Display *tickdisp;
Fl_Light_Button *runbutton;
Fl_Button *runstepbutton;
Fl_Value_Input *instep;
Fl_Group *drawgroup;
Fl_Light_Button *drawbutton;
Fl_Check_Button *drawnodes;
Fl_Check_Button *drawlinks;
Fl_Check_Button *drawfluxes;
Fl_Button *clearbutton;
Fl_Check_Button *printdatabutton;
Fl_Button *turbobutton;

// button in right column 1
/*Fl_Button *bnewlattice;
 Fl_Button *bnewrandom;
 Fl_Button *bnewclustered;*/
Fl_Button *bnewnet;
Fl_Button *bload; Fl_File_Chooser *loadchooser ;
Fl_Button *bsave; Fl_File_Chooser *savechooser ;
Fl_Button *bsaveas;
Fl_Button *bsetlayout;
Fl_Button *buttonexit;

// button in right column 2
DataFrame        *datascene;
Fl_Text_Buffer      *databuff;
Fl_Text_Display *datadisp;

Fl_Text_Buffer      *varbuff;
Fl_Text_Buffer      *meanbuff;
Fl_Text_Buffer      *dvarbuff;
Fl_Text_Buffer      *dmeanbuff;

Fl_Group *histogroup;
Fl_Check_Button *drawhistobutton;
Fl_Text_Display *meandisp;
Fl_Text_Display *vardisp;
Fl_Round_Button *roundstop;
Fl_Value_Input *inmeanstop;
Fl_Value_Input *invarstop;
Fl_Value_Input *inmeanerstop;
Fl_Value_Input *invarerstop;
double hism, hisv;
Fl_Text_Display *dmeandisp;
Fl_Text_Display *dvardisp;
Fl_Button *printbutton1;
Fl_Button *cstopbutton;

Fl_Button *test;


//bottom buttons
Fl_Group *bottomgroup;
Fl_Round_Button *round1;
Fl_Round_Button *round2;
Fl_Round_Button *round3;
Fl_Group *palettegroup;
Fl_Round_Button *palette0;
Fl_Round_Button *palette1;
Fl_Round_Button *palette2;
Fl_Text_Buffer *pathbuff;
Fl_Text_Display *pathdisp;
extern igraph_vector_t statstate;

//--------------------- dialogues  ---------------------
#define DIAL_W 200
#define DIAL_WT 400
#define DIAL_H 500

#define BUTTON_WD 100



Fl_Window    *dial0;

//new network dial
Fl_Window *dialnewnet;
Fl_Button *bnewlattice;
Fl_Button *bnewrandom;
Fl_Button *bnewclustered;
Fl_Button *bnewclusger2;

// new lattice dial
Fl_Window    *dial1;
Fl_Counter *dimcounter;
Fl_Counter *sidecounter;
Fl_Group *nlatgroup;
Fl_Check_Button *toroidal;
Fl_Round_Button *unifconn;
Fl_Round_Button *randomconn;
Fl_Round_Button *potentialconn;
Fl_Button *setpotentialbutton;
Fl_Button *done1;
Fl_Check_Button *dissipatelatcheck;
Fl_Value_Input *disslat;
Fl_Value_Input *dnlat;

//new random dial
Fl_Window    *dial2;
Fl_Value_Input *innodesnumber;
Fl_Value_Input *inprob;
Fl_Value_Input *inedges;
Fl_Button *done2;
Fl_Check_Button *dissipaterancheck;
Fl_Value_Input *dissran;
Fl_Value_Input *dnran;

//new clustered dial
Fl_Window    *dial3;
Fl_Value_Input *ininterconn;
Fl_Value_Input *inintraconn;
Fl_Value_Input *inbridge;
Fl_Value_Input *inclus;
Fl_Value_Input *inclusdim;
Fl_Value_Input *inwall;
Fl_Value_Input *instim;
Fl_Button *done3;
Fl_Check_Button *dissipatecluscheck;
Fl_Value_Input *dissclus;
Fl_Value_Input *dnclus;

//new clustered gerarchic 2 layers dial
Fl_Window    *dial4;
Fl_Value_Input *ininterconn2l;
Fl_Value_Input *inintraconn2l1;
Fl_Value_Input *inintraconn2l2;
Fl_Value_Input *inbridge2l1;
Fl_Value_Input *inbridge2l2;
Fl_Value_Input *inclus2l1;
Fl_Value_Input *inclus2l2;
Fl_Value_Input *inclusdim2l;
Fl_Value_Input *inwall2l1;
Fl_Value_Input *inwall2l2;
Fl_Check_Button *dissipateclus2lcheck;
Fl_Value_Input *dissclus2l;
Fl_Value_Input *dnclus2l;


Fl_Button *done4;


#define DEFINTER 1
#define DEFINTRA 1
#define DEFBRIDGE 1
#define DEFWALL 5
#define DEFDISS 1
#define DEFCN 2
#define DEFCS 10



#define DEFCS2L 10
#define DEFCN2L1 2
#define DEFINTER2L 1
#define DEFINTRA2L1 1
#define DEFBRIDGE2L1 1
#define DEFWALL2L1 5

#define DEFCN2L2 2
#define DEFINTRA2L2 1
#define DEFBRIDGE2L2 1
#define DEFWALL2L2 5

#define DEFDISS2L 1



//turbo dial
Fl_Window    *dialturbo;
Fl_Button *turbogo;



int initdrawbut=0;








//---------------------------------------- TURBO --------------------------------------------------------
void TurboDial(void) {
    
    int w_est,h_est;
    int ypos, xpos, widgw, widgh;
    
    char text[1000];

    
    w_est= DIAL_WT;   h_est= 260;
    
    dialturbo = new Fl_Window(w_est,h_est,"Turbo Run");
    dialturbo->set_modal();
    
    ypos=BORDER1;
    xpos=BORDER1;
    
    
    widgh=200;
    widgw=(DIAL_WT-(2*BORDER1)-(2*BORDER1));
    Fl_Box *turbotext = new Fl_Box(xpos,ypos,widgw,widgh, "");
    //turbotext->labelfont(FL_BOLD);
    turbotext->align(FL_ALIGN_INSIDE + FL_ALIGN_TOP_LEFT );
    

    turbotext->label("WARNING  \n You won\'t be able to access LAURA during TURBO RUN.\n\nPlease check --all-- your settings before TURBO RUN.\n\nIf you have to kill a TURBO RUN, kill the main program.");
    
    ypos=ypos+widgh;
    
    
    
    widgh=BUTTON_H;
    widgw=(DIAL_WT-(2*BORDER1)-(2*BORDER1))/2;
    Fl_Button *turbogo = new Fl_Button(xpos, ypos, widgw, widgh, "Turbo Run");
    turbogo->labelfont(FL_BOLD);
    xpos=xpos+widgw+2*BORDER1;
    
    widgh=BUTTON_H;
    Fl_Button *cancelturbo = new Fl_Button(xpos, ypos, widgw, widgh, "Cancel");
    ypos=ypos+widgh;
    
    ypos=ypos+10;
    h_est=ypos;
    dialturbo->size(w_est,h_est);
    
    cancelturbo->callback(exitdialturbocb,0);
    turbogo->callback(goturbocb,0);
    
    dialturbo->end();
    dialturbo->show();
    
    
}



//--------------------------------------------
void exitdialturbocb(Fl_Widget *, void *) {
    dialturbo->hide();
}



void goturborun(){
    
    clear();
    
    //set parameters
    particles=(int)in1->value();
    vincolo=(int)in2->value();
    threshold=(int)inthreshold->value();
    totrun=(int)in3->value();
    maxtime=(int)in4->value();
    ispii=(int)pii->value();
   
    if ( (int)roundTAS->value()==1 ){
        logDebug("\nTAS initial state.\n");
        InitialStateTAS((int)in5->value(),(int)roundTASRAN->value(),(int)in3->value());
        
    }
    else if((int)roundRAN->value()==1){
        logDebug("\nRAN initial state.\n");
        InitialStateRAN((int)in6->value(), (int)in3->value());
        
    }
    else if((int)roundSTAT->value()==1){
        logDebug("\nSTAT initial state.\n");
        InitialStateSTAT((int)in3->value());
        
    }
    
    else if(havepath==0){
        InitialStateTAS(0,0,(int)in3->value()); roundTAS->value(1);
    }
    else if(haveloadedstate==1){
        InitialStateLOAD();
        in1->value(particles);
        in3->value(totrun);
        maxtime=(int)in4->value();
    }
    
    
    //set OUTPUT
    openout();
    
    TurboRun((double)indt->value());
    
    closeout();
}


void goturbocb(Fl_Widget *, void *) {
    
    dialturbo->hide();
    logDebug("Turbo run Started\n");
    goturborun();
    logDebug("Turbo run Completed. (maxtime=%i) \n",maxtime);
    
}



//----------------------------------------CHANGE LAYOUT--------------------------------------------------------
void ChangeLayout(void) {
    
    int w_est,h_est;
    int ypos, xpos, widgw, widgh;
    
    
    igraph_matrix_init(&oldlayout,0,0);
    
    igraph_matrix_update(&oldlayout,&layout);
    
    w_est= DIAL_W;   h_est= 260;
    
    dial0 = new Fl_Window(w_est,h_est,"Change Layout");
    dial0->set_modal();
    
    ypos=BORDER1;
    xpos=BORDER1;
    
    widgh=BUTTON_H1; widgw=DIAL_W-(2*xpos);
    Fl_Button *gridlaybut = new Fl_Button(xpos, ypos, widgw, widgh, "Grid");
    ypos=ypos+widgh;
    
    widgh=BUTTON_H1; widgw=DIAL_W-(2*xpos);
    Fl_Button *starlaybut = new Fl_Button(xpos, ypos, widgw, widgh, "Star");
    ypos=ypos+widgh;
    
    widgh=BUTTON_H1; widgw=DIAL_W-(2*xpos);
    Fl_Button *circlelaybut = new Fl_Button(xpos, ypos, widgw, widgh, "Circle");
    ypos=ypos+widgh;
    
    widgh=BUTTON_H1; widgw=DIAL_W-(2*xpos);
    Fl_Button *mdslaybut = new Fl_Button(xpos, ypos, widgw, widgh, "MDS by Torgerson");
    ypos=ypos+widgh;
    
    widgh=BUTTON_H1; widgw=DIAL_W-(2*xpos);
    Fl_Button *forceunwlaybut = new Fl_Button(xpos, ypos, widgw, widgh, "Force UnWeighted");
    ypos=ypos+widgh;
    
    widgh=BUTTON_H1; widgw=DIAL_W-(2*xpos);
    Fl_Button *forcewlaybut = new Fl_Button(xpos, ypos, widgw, widgh, "Force Weighted");
    ypos=ypos+widgh;
    
    
    
    widgh=BUTTON_H1;
    ypos=ypos+30;
    Fl_Button *reset = new Fl_Button(xpos, ypos, widgw, widgh, "Reset Layout");
    ypos=ypos+widgh;
    
    widgh=BUTTON_H1;
    ypos=ypos+10;
    Fl_Button *done0 = new Fl_Button(20, ypos, DIAL_W-40, widgh, "Close and Keep Layout");
    ypos=ypos+widgh;
    
    
    gridlaybut->callback(gridcb,0);
    starlaybut->callback(starcb,0);
    circlelaybut->callback(circlecb,0);
    mdslaybut->callback(mdscb,0);
    forceunwlaybut->callback(forceunwlaycb,0);
    forcewlaybut->callback(forcewlaycb,0);
    
    reset->callback(resetlayoutcb,0);
    done0->callback(exitdial0cb,0);
    
    ypos=ypos+10;
    h_est=ypos;
    dial0->size(w_est,h_est);
    
    
    dial0->end();
    dial0->show();
    
    
}




void gridcb(Fl_Widget *, void *) {
    
    igraph_layout_grid(&graph, &layout, sqrt((double)nodesnumber));
    
    adjustlayout();
}

void starcb(Fl_Widget *, void *) {
    
    int central;
    
    igraph_vector_t res;
    igraph_vector_init(&res,0);
    
    igraph_degree(&graph, &res, igraph_vss_all(), IGRAPH_ALL, IGRAPH_NO_LOOPS);
    
    central=igraph_vector_which_max(&res);
    
    igraph_vector_destroy(&res);
    
    igraph_layout_star(&graph, &layout,central,0);
    
    adjustlayout();
    
    
    
}

void circlecb(Fl_Widget *, void *) {
    
    
    igraph_layout_circle(&graph, &layout);
    
    adjustlayout();
    
    
    
}


void mdscb(Fl_Widget *, void *) {
    
    igraph_matrix_t dist;
    
    igraph_matrix_init(&dist,0,0);
    igraph_matrix_update(&dist,&admatrix);
    igraph_matrix_scale(&dist,-1);
    igraph_matrix_add_constant(&dist,1);
    
    igraph_layout_mds(&graph, &layout, &dist, 2, 0);
    
    igraph_matrix_destroy(&dist);
    
    adjustlayout();
    
    
}

void forcewlaycb(Fl_Widget *, void *) {
    
    
    igraph_t tempgraph;
    igraph_vector_t weight;
    
    igraph_weighted_adjacency(&tempgraph, &admatrix,IGRAPH_ADJ_DIRECTED, "w",0);
    
    igraph_vector_init(&weight,0);
    
    EANV(&tempgraph,"w",&weight);
    
    
    
    igraph_destroy(&tempgraph);
    
    igraph_vector_scale(&weight,10.);
    
    igraph_layout_fruchterman_reingold(&graph, &layout,500, (float)nodesnumber,
                                       sqrt((double)nodesnumber), 1.5,
                                       (float)(nodesnumber*nodesnumber), 0,
                                       &weight,
                                       NULL,NULL,NULL,NULL);
    
    igraph_vector_destroy(&weight);
    
    adjustlayout();
    
    
}

void forceunwlaycb(Fl_Widget *, void *) {
    
    
    
    
    igraph_layout_fruchterman_reingold(&graph, &layout,500, (float)nodesnumber,
                                       sqrt((double)nodesnumber), 1.5,
                                       (float)(nodesnumber*nodesnumber), 0,
                                       NULL,
                                       NULL,NULL,NULL,NULL);
    
    
    
    adjustlayout();
    
    
}


void resetlayoutcb(Fl_Widget *, void *) {
    
    igraph_matrix_update(&layout,&oldlayout);
}




//--------------------------------------------
void exitdial0cb(Fl_Widget *, void *) {
    
    igraph_matrix_destroy(&oldlayout);
    dial0->hide();
}






//----------------------------------------NEW NETWORK (GENERAL)---------------------------------------------------------
void DialogueNewNet(void) {
    
    int w_est,h_est;
    int ypos, xpos, widgh, widgw;
    
    w_est= DIAL_W;   h_est= 260;
    
    dialnewnet =  new Fl_Window(w_est,h_est,"New Network");
    dialnewnet->set_modal();
    
    ypos=BORDER1; xpos=BORDER1/2;
    widgw=DIAL_W-20;
    
    
    widgh=BUTTON_H1;
    bnewlattice = new Fl_Button(xpos, ypos, widgw, widgh, "Lattice");
    ypos=ypos+widgh;
    
    ypos=ypos+10;
    widgh=BUTTON_H1;
    bnewrandom = new Fl_Button(xpos, ypos, widgw, widgh, "Random");
    ypos=ypos+widgh;
    
    ypos=ypos+10;
    widgh=BUTTON_H1;
    bnewclustered = new Fl_Button(xpos, ypos, widgw, widgh, "Clustered");
    ypos=ypos+widgh;
    
    ypos=ypos+10;
    widgh=BUTTON_H1;
    bnewclusger2 = new Fl_Button(xpos, ypos, widgw, widgh, "Clust - Gerar (2L)");
    ypos=ypos+widgh;
    
    bnewlattice->callback(newlatcb,0);
    bnewrandom->callback(newrandcb,0);
    bnewclustered->callback(newcluscb,0);
    bnewclusger2->callback(newclusger2cb,0);
    
    
    widgh=BUTTON_H1;
    ypos=ypos+30;
    Fl_Button *cancel = new Fl_Button(20, ypos, DIAL_W-40, widgh, "CANCEL");
    ypos=ypos+widgh;
    
    
    cancel->callback(exitdialnncb,0);
    
    ypos=ypos+10;
    h_est=ypos;
    dialnewnet->size(w_est,h_est);
    
    
    dialnewnet->end();
    dialnewnet->show();
    
    
    
}


//--------------------------------------------
void exitdialnncb(Fl_Widget *, void *) {
    dialnewnet->hide();
    graphisloaded=1; error=0;
}





//--------------------------------------------
void newlatcb(Fl_Widget *, void *) {
    DialogueNewLat();
}





//--------------------------------------------
void newrandcb(Fl_Widget *, void *) {
    DialogueNewRand();
}



//--------------------------------------------
void newcluscb(Fl_Widget *, void *) {
    DialogueNewClus();

}


//--------------------------------------------
void newclusger2cb(Fl_Widget *, void *) {
    DialogueNewClusGer2();
}




//----------------------------------------NEW LATTICE---------------------------------------------------------
void DialogueNewLat(void) {
    
    int w_est,h_est;
    int ypos=0, xpos=0, widgh, widgw;
    
    w_est= DIAL_W;   h_est= 260;
    
    dial1 =  new Fl_Window(w_est,h_est,"New Lattice");
    dial1->set_modal();
    
    ypos=BORDER1; xpos=BORDER1/2;
    widgw=DIAL_W-20;
    
    widgh=20;
    dimcounter = new Fl_Counter(xpos, ypos, widgw, widgh, "Dimensions (1D, 2D, 3D...)");
    dimcounter->type(FL_SIMPLE_COUNTER);
    dimcounter->step(1);
    dimcounter->minimum(1);
    dimcounter->value(2);
    ypos=ypos+widgh;
    
    widgh=20;
    ypos=ypos+30;
    sidecounter = new Fl_Counter(xpos, ypos, widgw, widgh, "Lattice Side");
    sidecounter->step(1, 5);
    sidecounter->minimum(2);
    sidecounter->value(5);
    ypos=ypos+widgh;
    
    widgh=20;
    ypos=ypos+30;
    toroidal = new Fl_Check_Button(xpos, ypos, widgw, widgh, "Toroidal");
    toroidal->set();
    ypos=ypos+widgh;
    
    ypos=ypos+30;
    
    nlatgroup = new Fl_Group(xpos-5,ypos-5,widgw+10,20+widgh*7, "Set Connectivities");
    nlatgroup->align(FL_ALIGN_TOP_LEFT);
    nlatgroup->labelfont(FL_BOLD);
    nlatgroup->box(FL_BORDER_BOX);
    
    
    ypos=ypos+10;
    widgh=30;
    unifconn = new Fl_Round_Button(10, ypos, DIAL_W-20, widgh, "Uniform Connectivities");
    unifconn->type(FL_RADIO_BUTTON);
    
    ypos=ypos+widgh;
    
    ypos=ypos+10;
    widgh=30;
    randomconn = new Fl_Round_Button(10, ypos, DIAL_W-20, widgh, "Random Connectivities");
    randomconn->type(FL_RADIO_BUTTON);
    
    ypos=ypos+widgh;
    
    ypos=ypos+10;
    widgh=30;
    potentialconn = new Fl_Round_Button(10, ypos, DIAL_W-20, widgh, "Use Potential (2D Only)");
    potentialconn->type(FL_RADIO_BUTTON);
    potentialconn->setonly();
    ypos=ypos+widgh;
    
    widgh=BUTTON_H1;
    setpotentialbutton = new Fl_Button(20, ypos, DIAL_W-40, widgh, "Set Potential");
    setpotentialbutton->deactivate();
    ypos=ypos+widgh;
    
    
    nlatgroup->end();
    
    xpos=BORDER1/2;
    ypos=ypos+20;
    widgh=BUTTON_H1;
    dissipatelatcheck = new Fl_Check_Button(xpos, ypos, DIAL_W-40, widgh, "Dissipate");
    ypos=ypos+widgh;
    ypos=ypos+20;
    
    widgh=20;
    disslat = new Fl_Value_Input(xpos,ypos,DIAL_W-40,widgh, "Dissipation Rate:");
    disslat->align(FL_ALIGN_TOP_LEFT);
    disslat->step(0.001);
    disslat->minimum(0);
    disslat->maximum(1);
    disslat->value(0.1);
    ypos=ypos+widgh;

    ypos=ypos+20;
    dnlat = new Fl_Value_Input(xpos,ypos,DIAL_W-40,widgh, "Source Node:");
    dnlat->align(FL_ALIGN_TOP_LEFT);
    dnlat->step(1);
    dnlat->minimum(0);
    dnlat->maximum(10000);
    dnlat->value(0);
    ypos=ypos+widgh;
    
    ypos=ypos+widgh;
    
    
    widgh=BUTTON_H1;
    ypos=ypos+20;
    done1 = new Fl_Button(20, ypos, DIAL_W-40, widgh, "Create");
    ypos=ypos+widgh;
    
    widgh=BUTTON_H1;
    ypos=ypos+10;
    Fl_Button *cancel = new Fl_Button(20, ypos, DIAL_W-40, widgh, "CANCEL");
    ypos=ypos+widgh;
    
    done1->callback(generatelatticecb,0);
    cancel->callback(exitdial1cb,0);
    
    ypos=ypos+10;
    h_est=ypos;
    dial1->size(w_est,h_est);
    
    
    dial1->end();
    dial1->show();
    
    
    
}




//--------------------------------------------
void generatelatticecb(Fl_Widget *, void *) {
    
    
    graphisloaded=0;
    
    
    if((int)randomconn->value()==1 || (int)unifconn->value()==1){
        
        latticedim=(int)dimcounter->value();
        latticeside=(int) sidecounter->value();
        islatrand=(int)randomconn->value();
        islatunif=(int)unifconn->value();
        islatpot=(int)potentialconn->value();
        istoro=(int)toroidal->value();
        usesteady=0;
        
        generatelattice( (int)dimcounter->value(), (int) sidecounter->value(), (int)toroidal->value(), (int)randomconn->value(), (int)dissipatelatcheck->value(), (double)disslat->value(), (int)dnlat->value());
        
    }
    
    else if ((int)potentialconn->value()==1){
        
        if(dimcounter->value()==2){
            
            latticedim=(int)dimcounter->value();
            latticeside=(int) sidecounter->value();
            islatrand=(int)randomconn->value();
            islatunif=(int)unifconn->value();
            islatpot=(int)potentialconn->value();
            istoro=(int)toroidal->value();
            
            generatepotlat2d(latticeside, istoro,(int)dissipatelatcheck->value(), (double)disslat->value(), (int)dnlat->value());
            
            
        }
        
        else {
            
            printf("\n ERRROR, NO 2D LATTICE. \n");
            sprintf(errorstring,"ERROR\nGENERATING\nPOTENTIAL\n\nNO 2D\nLATTICE");
            error=1;
            rewrite=1;
            return;
            
            
        }
        
    }
    
    
    
    //inizializations
    igraph_matrix_init(&state, 0, 0);
    igraph_matrix_init(&density, 0, 0); igraph_matrix_init(&densityold,0,0);
    igraph_matrix_init(&statenew, 0, 0);
    igraph_matrix_init(&loss, 0, 0);
    
    beginner=(nodesnumber/2); in5->value(beginner);
    InitialStateTAS(beginner,0,totrun);
    
    igraph_vector_init(&gain,nodesnumber);
    
    
    // set global control variables
    
    if(islatpot==0){
        usesteady=0;
        activatesteady->value(0);
        round1->setonly(); bottomgroup->deactivate();
        rewrite=1;
        roundSTAT->deactivate();
    }
    else if(islatpot==1){
        activatesteady->value(0);
        bottomgroup->activate();
        histogroup->activate();
        usesteady=0;
        round1->setonly();
        rewrite=1;
        roundSTAT->deactivate();
        
        
    }
    
    
    pii->value(0); ispii=pii->value();
    
    islattice=1;  rewrite=1;
    israndomER1=0; isclustered=0;
    graphisloaded=1;
    newgraph=1;
    
    dial1->hide();    dialnewnet->hide();
    
    clear();
    
    roundTAS->setonly();
    in5->maximum(igraph_vcount(&graph)-1);
    
    havepath=0;
    
}


//--------------------------------------------
void exitdial1cb(Fl_Widget *, void *) {
    dial1->hide();
    graphisloaded=1; error=0;
}










//----------------------------------------NEW RANDOM---------------------------------------------------------



//-------------------------------------------------------------------------------------------------
void DialogueNewRand(void) {
    int w_est,h_est;
    
    int xpos, ypos, widgh, widw;
    
    w_est= DIAL_W;   h_est= 200;
    
    dial2=   new Fl_Window(w_est,h_est,"New E-R Random Graph");
    dial2->set_modal();
    
    ypos=20;
    
    widgh=20;
    ypos=ypos;
    innodesnumber = new Fl_Value_Input(10,ypos,w_est-20,widgh, "Nodes (n):");
    innodesnumber->align(FL_ALIGN_TOP_LEFT);
    innodesnumber->step(1);
    innodesnumber->minimum(2);
    innodesnumber->maximum(100000000000);
    innodesnumber->value(25);
    ypos=ypos+widgh;
    
    widgh=20;
    ypos=ypos+30;
    inprob = new Fl_Value_Input(10,ypos,w_est-20,widgh, "Probability (p):");
    inprob->align(FL_ALIGN_TOP_LEFT);
    inprob->step(0.001);
    inprob->minimum(0.01);
    inprob->maximum(1);
    inprob->value(0.500);
    ypos=ypos+widgh;
    
    xpos=BORDER1/2;
    ypos=ypos+20;
    widgh=BUTTON_H1;
    dissipaterancheck = new Fl_Check_Button(xpos, ypos, DIAL_W-40, widgh, "Dissipate");
    ypos=ypos+widgh;
    ypos=ypos+20;
    
    widgh=20;
    dissran = new Fl_Value_Input(xpos,ypos,DIAL_W-40,widgh, "Dissipation Rate:");
    dissran->align(FL_ALIGN_TOP_LEFT);
    dissran->step(0.001);
    dissran->minimum(0);
    dissran->maximum(1);
    dissran->value(0.1);
    ypos=ypos+widgh;
    
    ypos=ypos+20;
    dnran = new Fl_Value_Input(xpos,ypos,DIAL_W-40,widgh, "Source Node:");
    dnran->align(FL_ALIGN_TOP_LEFT);
    dnran->step(1);
    dnran->minimum(0);
    dnran->maximum(10000);
    dnran->value(0);
    ypos=ypos+widgh;
    
    ypos=ypos+widgh;
    
    widgh=BUTTON_H1;
    ypos=ypos+30;
    done2 = new Fl_Button(20, ypos, w_est-40, widgh, "Create E-R G(n,p)");
    ypos=ypos+widgh;
    
    ypos=ypos+10;
    Fl_Button *cancel = new Fl_Button(20, ypos, DIAL_W-40, widgh, "CANCEL");
    ypos=ypos+widgh;
    
    done2->callback(generaterandom1cb,0);
    cancel->callback(exitdial2cb,0);
    
    
    ypos=ypos+10;
    h_est=ypos;
    dial2->size(w_est,h_est);
    dial2->end();
    dial2->show();
    
    
    
}




//--------------------------------------------
void generaterandom1cb(Fl_Widget *, void *) {
    
    graphisloaded=0;
    
    beginner=generaterandom1((int)innodesnumber->value(), (float)inprob->value(), (int)dissipaterancheck->value(), (double)dissran->value(), (int)dnran->value());
    
    if(graphisloaded==1){
        
        
        
        //inizializations
        igraph_matrix_init(&state, 0, 0);
        
        igraph_matrix_init(&density, 0, 0);
        igraph_matrix_init(&densityold,0,0);
        igraph_matrix_init(&statenew, 0, 0);
        igraph_matrix_init(&loss, 0, 0);
        
        //beginner=0;
        in5->value(beginner);
        InitialStateTAS(beginner,0,totrun);
        
        igraph_vector_init(&gain,nodesnumber);
        
        //  StationaryState();
        
        usesteady=0; roundSTAT->deactivate(); activatesteady->value(0); round1->setonly();bottomgroup->deactivate();
        pii->value(0); ispii=pii->value();
        islattice=0; israndomER1=1;  rewrite=1; isclustered=0;
        newgraph=1;
        
        dial2->hide();    dialnewnet->hide();
        
        clear();
        
        roundTAS->setonly();
        in5->maximum(igraph_vcount(&graph)-1);
        
        havepath=0;
    }
    graphisloaded=1;
}


//--------------------------------------------
void exitdial2cb(Fl_Widget *, void *) {
    dial2->hide();
    graphisloaded=1; error=0;
}













//---------------------------------------- NEW clustered ---------------------------------------------------------
void DialogueNewClus(void) {
    
    int w_est; int h_est=10;
    int ypos, xpos, widgh, widgw;
    
    w_est= DIAL_W;
    
    dial3 =  new Fl_Window(w_est,h_est,"New clustered");
    dial3->set_modal();
    
    ypos=BORDER1; xpos=BORDER1/2;
    
    widgw=DIAL_W-BORDER1;
    
    widgh=20;
    inclus = new Fl_Value_Input(10,ypos,w_est-20,widgh, "# of clusters:");
    inclus->align(FL_ALIGN_TOP_LEFT);
    inclus->step(1);
    inclus->minimum(1);
    inclus->maximum(10);
    inclus->value(DEFCN);
    ypos=ypos+widgh;
    
    
    widgh=20;
    ypos=ypos+30;
    inclusdim = new Fl_Value_Input(10,ypos,w_est-20,widgh, "Clusters' dimension:");
    inclusdim->align(FL_ALIGN_TOP_LEFT);
    inclusdim->step(1);
    inclusdim->minimum(1);
    inclusdim->maximum(1000);
    inclusdim->value(DEFCS);
    ypos=ypos+widgh;
    
    widgh=20;
    ypos=ypos+30;
    ininterconn = new Fl_Value_Input(10,ypos,w_est-20,widgh, "Inter-connectivity:");
    ininterconn->align(FL_ALIGN_TOP_LEFT);
    ininterconn->step(0.01);
    ininterconn->minimum(0.01);
    ininterconn->maximum(1);
    ininterconn->value(DEFINTER);
    ypos=ypos+widgh;
    
    widgh=20;
    ypos=ypos+30;
    inintraconn = new Fl_Value_Input(10,ypos,w_est-20,widgh, "Intra-connectivity:");
    inintraconn->align(FL_ALIGN_TOP_LEFT);
    inintraconn->step(0.0001);
    inintraconn->minimum(0.01);
    inintraconn->maximum(1);
    inintraconn->value(DEFINTRA);
    ypos=ypos+widgh;
    
    widgh=20;
    ypos=ypos+30;
    inbridge = new Fl_Value_Input(10,ypos,w_est-20,widgh, "Nodes on Bridge:");
    inbridge->align(FL_ALIGN_TOP_LEFT);
    inbridge->step(1);
    inbridge->minimum(1);
    inbridge->maximum(100000);
    inbridge->value(DEFBRIDGE);
    ypos=ypos+widgh;
    
    widgh=20;
    ypos=ypos+30;
    inwall = new Fl_Value_Input(10,ypos,w_est-20,widgh, "Wall:");
    inwall->align(FL_ALIGN_TOP_LEFT);
    inwall->step(0.01);
    inwall->minimum(0.);
    inwall->maximum(1000000);
    inwall->value(DEFWALL);
    ypos=ypos+widgh;
    
    
    xpos=BORDER1/2;
    ypos=ypos+20;
    widgh=BUTTON_H1;
    dissipatecluscheck = new Fl_Check_Button(xpos, ypos, DIAL_W-40, widgh, "Dissipate");
    ypos=ypos+widgh;
    ypos=ypos+20;
    
    widgh=20;
    dissclus= new Fl_Value_Input(xpos,ypos,DIAL_W-40,widgh, "Dissipation Rate:");
    dissclus->align(FL_ALIGN_TOP_LEFT);
    dissclus->step(0.001);
    dissclus->minimum(0);
    dissclus->maximum(1);
    dissclus->value(0.1);
    ypos=ypos+widgh;
    
    ypos=ypos+20;
    dnclus= new Fl_Value_Input(xpos,ypos,DIAL_W-40,widgh, "Source Node:");
    dnclus->align(FL_ALIGN_TOP_LEFT);
    dnclus->step(1);
    dnclus->minimum(0);
    dnclus->maximum(10000);
    dnclus->value(0);
    ypos=ypos+widgh;
    
    ypos=ypos+widgh;
    
    
    
    
    widgh=BUTTON_H1;
    ypos=ypos+10;
    Fl_Button *done3 = new Fl_Button(xpos, ypos, widgw, widgh, "Done");
    ypos=ypos+widgh;
    
    
    widgh=BUTTON_H1;
    ypos=ypos+10;
    Fl_Button *cancel = new Fl_Button(xpos, ypos, widgw, widgh, "CANCEL");
    ypos=ypos+widgh;
    
    
    
    ypos=ypos+10;
    h_est=ypos;
    dial3->size(w_est,h_est);
    
    cancel->callback(exitdial3cb,0);
    done3->callback(generateclustsymcb,0);
    
    dial3->end();
    dial3->show();
    
    
    
}





//--------------------------------------------
void generateclustsymcb(Fl_Widget *, void *) {
    
    graphisloaded=0;
    
    
    generateclustsympots((int)inclus->value(),
                         (int)inclusdim->value(),
                         (int)inbridge->value(),
                         (double)ininterconn->value(),
                         (double)inintraconn->value(),
                         (double)inwall->value(),
                         (int)dissipatecluscheck->value(),
                         (double)dissclus->value(),
                         (int)dnclus->value());
    
    // generateclustsym(2, 10, 1, 0.5, 0);
    
    
    
    
    if(graphisloaded==1){
        
        //inizializations
        igraph_matrix_init(&state, 0, 0);
        
        igraph_matrix_init(&density, 0, 0);
        igraph_matrix_init(&densityold,0,0);
        igraph_matrix_init(&statenew, 0, 0);
        igraph_matrix_init(&loss, 0, 0);
        
        igraph_vector_init(&gain,nodesnumber);
        
        //beginner=0; set initial state as TAS state
        InitialStateRAN(1, totrun);
        
        usesteady=0; roundSTAT->deactivate(); activatesteady->value(0); round1->setonly();bottomgroup->deactivate();
        pii->value(0); ispii=pii->value();
        islattice=0; israndomER1=0;  rewrite=1; isclustered=0;
        newgraph=1;
        
        dial3->hide();    dialnewnet->hide();
        
        clear();
        
        roundRAN->setonly();
        in5->maximum(igraph_vcount(&graph)-1);
        
        havepath=0;
    }
    graphisloaded=1;
    isclustered=1;
    
    
    print_matrix_ur(&admatrix,stdout);
    
    
    
}



//--------------------------------------------
void exitdial3cb(Fl_Widget *, void *) {
    dial3->hide();
    graphisloaded=1; error=0;
}










//---------------------------------------- NEW clustered gerarchic 2 layers --------------------------------------------------
void DialogueNewClusGer2(void) {
    
    int w_est; int h_est=10;
    int ypos, xpos, widgh, widgw;
    
    w_est= DIAL_W;
    
    dial4 =  new Fl_Window(w_est,h_est,"New Clus-GERARCHIC - 2 Layers");
    dial4->set_modal();
    
    ypos=BORDER1/2;
    xpos=BORDER1/2;
    
    widgw=DIAL_W-BORDER1;
    
    
    widgh=20;
    ypos=ypos+30;
    inclusdim2l = new Fl_Value_Input(xpos,ypos,widgw,widgh, "L1 Clusters' dimension:");
    inclusdim2l->align(FL_ALIGN_TOP_LEFT);
    inclusdim2l->step(1);
    inclusdim2l->minimum(1);
    inclusdim2l->maximum(1000);
    inclusdim2l->value(DEFCS2L);
    ypos=ypos+widgh;
    
    
    
    
    widgh=20;
    ypos=ypos+30;
    inclus2l1 = new Fl_Value_Input(xpos,ypos,widgw,widgh, "# of clusters (L1):");
    inclus2l1->align(FL_ALIGN_TOP_LEFT);
    inclus2l1->step(1);
    inclus2l1->minimum(1);
    inclus2l1->maximum(10);
    inclus2l1->value(DEFCN2L1);
    ypos=ypos+widgh;
    
    
    
    widgh=20;
    ypos=ypos+30;
    ininterconn2l = new Fl_Value_Input(xpos,ypos,widgw,widgh, "Inter-connectivity: (L1)");
    ininterconn2l->align(FL_ALIGN_TOP_LEFT);
    ininterconn2l->step(0.01);
    ininterconn2l->minimum(0.01);
    ininterconn2l->maximum(1);
    ininterconn2l->value(DEFINTER2L);
    ypos=ypos+widgh;
    
    widgh=20;
    ypos=ypos+30;
    inintraconn2l1 = new Fl_Value_Input(xpos,ypos,widgw,widgh, "Intra-connectivity: (L1)");
    inintraconn2l1->align(FL_ALIGN_TOP_LEFT);
    inintraconn2l1->step(0.0001);
    inintraconn2l1->minimum(0.01);
    inintraconn2l1->maximum(1);
    inintraconn2l1->value(DEFINTRA2L1);
    ypos=ypos+widgh;
    
    
    widgh=20;
    ypos=ypos+30;
    inbridge2l1 = new Fl_Value_Input(xpos,ypos,widgw,widgh, "Nodes on Bridge: (L1)");
    inbridge2l1->align(FL_ALIGN_TOP_LEFT);
    inbridge2l1->step(1);
    inbridge2l1->minimum(1);
    inbridge2l1->maximum(100000);
    inbridge2l1->value(DEFBRIDGE2L1);
    ypos=ypos+widgh;
    
    
    widgh=20;
    ypos=ypos+30;
    inwall2l1 = new Fl_Value_Input(xpos,ypos,widgw,widgh, "Wall (L1):");
    inwall2l1->align(FL_ALIGN_TOP_LEFT);
    inwall2l1->step(0.01);
    inwall2l1->minimum(0.);
    inwall2l1->maximum(1000000);
    inwall2l1->value(DEFWALL2L1);
    ypos=ypos+widgh;
    
    xpos=xpos+widgw+BORDER1/2;
    ypos=BORDER1/2;
    
    ypos=ypos+widgh+30;
    
    widgh=20;
    ypos=ypos+30;
    inclus2l2 = new Fl_Value_Input(xpos,ypos,widgw,widgh, "# of clusters (L2):");
    inclus2l2->align(FL_ALIGN_TOP_LEFT);
    inclus2l2->step(1);
    inclus2l2->minimum(1);
    inclus2l2->maximum(10);
    inclus2l2->value(DEFCN2L2);
    ypos=ypos+widgh;
    
    ypos=ypos+widgh+30;
    
    widgh=20;
    ypos=ypos+30;
    inintraconn2l2 = new Fl_Value_Input(xpos,ypos,widgw,widgh, "Intra-connectivity: (L2)");
    inintraconn2l2->align(FL_ALIGN_TOP_LEFT);
    inintraconn2l2->step(0.01);
    inintraconn2l2->minimum(0.01);
    inintraconn2l2->maximum(1);
    inintraconn2l2->value(DEFINTRA2L2);
    ypos=ypos+widgh;
    
    
    widgh=20;
    ypos=ypos+30;
    inbridge2l2 = new Fl_Value_Input(xpos,ypos,widgw,widgh, "Nodes on Bridge: (L2)");
    inbridge2l2->align(FL_ALIGN_TOP_LEFT);
    inbridge2l2->step(1);
    inbridge2l2->minimum(1);
    inbridge2l2->maximum(100000);
    inbridge2l2->value(DEFBRIDGE2L2);
    ypos=ypos+widgh;
    
    widgh=20;
    ypos=ypos+30;
    inwall2l2 = new Fl_Value_Input(xpos,ypos,widgw,widgh, "Wall (L2):");
    inwall2l2->align(FL_ALIGN_TOP_LEFT);
    inwall2l2->step(0.01);
    inwall2l2->minimum(0.);
    inwall2l2->maximum(1000000);
    inwall2l2->value(DEFWALL2L2);
    ypos=ypos+widgh;
    
    
    ypos=ypos+100;
    widgh=30;
    dissipateclus2lcheck = new Fl_Check_Button(xpos,ypos,widgw,widgh, "Dissipate");
    ypos=ypos+widgh;
    
    widgh=20;
    ypos=ypos+10;
    dissclus2l= new Fl_Value_Input(xpos,ypos,widgw,widgh, "Dissipation rate:");
    dissclus2l->align(FL_ALIGN_TOP_LEFT);
    dissclus2l->step(0.001);
    dissclus2l->minimum(0);
    dissclus2l->maximum(1);
    dissclus2l->value(DEFDISS2L);
    ypos=ypos+widgh;
    
    widgh=20;
    ypos=ypos+30;
    dnclus2l = new Fl_Value_Input(xpos,ypos,widgw,widgh, "Source node:");
    dnclus2l->align(FL_ALIGN_TOP_LEFT);
    dnclus2l->step(1);
    dnclus2l->minimum(0);
    dnclus2l->maximum(10000);
    dnclus2l->value(0);
    ypos=ypos+widgh;
    
    
    
    
    widgh=BUTTON_H1;
    ypos=ypos+30;
    done4 = new Fl_Button(xpos, ypos, widgw, widgh, "Done");
    ypos=ypos+widgh;
    
    
    widgh=BUTTON_H1;
    ypos=ypos+10;
    Fl_Button *cancel = new Fl_Button(xpos, ypos, widgw, widgh, "CANCEL");
    ypos=ypos+widgh;
    
    
    
    ypos=ypos+10;
    h_est=ypos;
    w_est=xpos+widgw+BORDER1/2;
    dial4->size(w_est,h_est);
    
    cancel->callback(exitdial4cb,0);
    done4->callback(generateclusger2cb,0);
    
    dial4->end();
    dial4->show();
    
    
    
}





//--------------------------------------------
void generateclusger2cb(Fl_Widget *, void *) {
    
    graphisloaded=0;
    
   /*  generateclusger2(2//int clusnumber1
                          ,2//int clusnumber2
                          ,4//int clusdim1
                          ,1//int cnodes1
                          ,1//int cnodes2
                          ,1.0//double interconn
                          ,1.0//double intraconn1
                          ,1.0// double intraconn2
                          ,5.// double wall1
                          ,5.// double wall2
                          ,0// int dissipate
                          ,0// double dissipation
                          );
*/
    
    
    generateclusger2((int)inclus2l1->value() //int clusnumber1
                     ,(int)inclus2l2->value()//int clusnumber2
                     ,(int)inclusdim2l->value()//int clusdim1
                     ,(int)inbridge2l1->value()//int cnodes1
                     ,(int)inbridge2l2->value()//int cnodes2
                     ,(double)ininterconn2l->value()//double interconn
                     ,(double)inintraconn2l1->value()//double intraconn1
                     ,(double)inintraconn2l2->value()// double intraconn2
                     ,(double)inwall2l1->value()// double wall1
                     ,(double)inwall2l2->value()// double wall2
                     ,(int)dissipateclus2lcheck->value()// int dissipate
                     ,(double)dissclus2l->value()// double dissipation
                     ,(int)dnclus2l->value());
    
    
    
    if(graphisloaded==1){
        
        printf("\n\n graph done. Inizializating dyinamics...\n"); fflush(stdout);
        
        //inizializations
        igraph_matrix_init(&state, 0, 0);
        printf("init state, "); fflush(stdout);
        igraph_matrix_init(&density, 0, 0); igraph_matrix_init(&densityold,0,0);
        printf("init density, "); fflush(stdout);
        igraph_matrix_init(&statenew, 0, 0);
        printf("init statenew, "); fflush(stdout);;
        igraph_matrix_init(&loss, 0, 0);
        printf("init loss, "); fflush(stdout);
        
        igraph_vector_init(&gain,nodesnumber);
        printf("init gain. "); fflush(stdout);
        
        //beginner=0; set initial state as TAS state
        InitialStateTAS(1,0,totrun);
        printf("TAS "); fflush(stdout);
        
        usesteady=0;
        roundSTAT->deactivate();
        activatesteady->value(0);
        round1->setonly();
        bottomgroup->deactivate();
        
        pii->value(0); ispii=pii->value();
        islattice=0; israndomER1=0;  rewrite=1; isclustered=0; isclustger2l=1;
        newgraph=1;
        
        dial4->hide();
        dialnewnet->hide();
        
        clear();
        
        roundRAN->setonly();
        in5->maximum(igraph_vcount(&graph)-1);
        
        havepath=0;
    }
    graphisloaded=1;
    isclustered=1;
    
    
    printf("\n\n\n\n");
    print_matrix_ur(&admatrix,stdout);
    
    printf("\n\n\n\n");
    
    
}



//--------------------------------------------
void exitdial4cb(Fl_Widget *, void *) {
    dial4->hide();
    graphisloaded=1; error=0;
}













//---------------- MAIN WINDOWS BUTTONS CALLBACKS -------------------------------------



//--------------------------------------------
void exitcb(Fl_Widget *, void *) {
    exit(0);
}




//--------------------------------------------
void newnetcb(Fl_Widget *, void *) {
    DialogueNewNet();
}


//--------------------------------------------
void activatesteadycb(Fl_Widget *, void *) {
    if(usesteady==0){
        bottomgroup->activate();
        histogroup->activate();
        roundSTAT->activate();
        StationaryState();
        usesteady=1;}
    else{
        histogroup->deactivate();
        bottomgroup->deactivate();
        roundSTAT->deactivate();
        usesteady=0;
    }
    
}






//--------------------------------------------
void piicb(Fl_Widget *, void *) {
    if(ispii==0) ispii=1;
    else    ispii=0;
    
}


//--------------------------------------------
void prevcb(Fl_Widget *, void *) {
   
    
    if((int)roundTAS->value()==1)InitialStateTAS((int)in5->value(),(int)roundTASRAN->value(),(int)in3->value());
    else if((int)roundRAN->value()==1)InitialStateRAN((int)in6->value(), (int)in3->value());
    else if((int)roundSTAT->value()==1) InitialStateSTAT((int)in3->value());
   
    
}


void openout(){
    
    char filename[100];
    
    if(haveout==0) {
        
        sprintf(filename,"%s%i.txt",setoutname->value(),1);
        
        if( (output1=fopen(filename,"w")) ==NULL) {
            logDebug("\nERROR: CANNOT OPEN OUTPUT FILE '%s'\n",filename);
            sprintf(errorstring,"ERROR\nOPENING\nOUTPUT\nFILE");
            error=1;
            rewrite=1;
            runbutton->value(0);
            return;
        }
        
        sprintf(filename,"%s%i.txt",setoutname->value(),2);
        
        if( (output2=fopen(filename,"w")) ==NULL) {
            logDebug("\nERROR: CANNOT OPEN OUTPUT FILE '%s'\n",filename);
            sprintf(errorstring,"ERROR\nOPENING\nOUTPUT\nFILE");
            error=1;
            rewrite=1;
            runbutton->value(0);
            return;
        }
        
        sprintf(filename,"%s%i.txt",setoutname->value(),3);
        
        if( (output3=fopen(filename,"w")) ==NULL) {
            logDebug("\nERROR: CANNOT OPEN OUTPUT FILE '%s'\n",filename);
            sprintf(errorstring,"ERROR\nOPENING\nOUTPUT\nFILE");
            error=1;
            rewrite=1;
            runbutton->value(0);
            return;
        }
        
        
        sprintf(filename,"%s%i.txt",setoutname->value(),5);
        
        if( (output5=fopen(filename,"w")) ==NULL) {
            logDebug("\nERROR: CANNOT OPEN OUTPUT FILE '%s'\n",filename);
            sprintf(errorstring,"ERROR\nOPENING\nOUTPUT\nFILE");
            error=1;
            rewrite=1;
            runbutton->value(0);
            return;
        }
        
        sprintf(filename,"%s%i.txt",setoutname->value(),6);
        
        if( (output6=fopen(filename,"w")) ==NULL) {
            logDebug("\nERROR: CANNOT OPEN OUTPUT FILE '%s'\n",filename);
            sprintf(errorstring,"ERROR\nOPENING\nOUTPUT\nFILE");
            error=1;
            rewrite=1;
            runbutton->value(0);
            return;
        }
        
        sprintf(filename,"%s%i.txt",setoutname->value(),7);
        
        if( (output7=fopen(filename,"w")) ==NULL) {
            logDebug("\nERROR: CANNOT OPEN OUTPUT FILE '%s'\n",filename);
            sprintf(errorstring,"ERROR\nOPENING\nOUTPUT\nFILE");
            error=1;
            rewrite=1;
            runbutton->value(0);
            return;
        }
        
        haveout=1;
    }
    
}


void closeout(){
    
    if (haveout==1) {
        fclose(output1);
        fclose(output2);
        fclose(output3);
        fclose(output5);
        fclose(output6);
        fclose(output7);
    }
    
    haveout=0;
    
}


void run(){
    
    if (runningcontrol==0) {
        
        if(cleared==1){
            
            //set parameters
            particles=(int)in1->value();
            vincolo=(int)in2->value();
            threshold=(int)inthreshold->value();
            totrun=(int)in3->value();
            maxtime=(int)in4->value();
            
            deltat=(double)indt->value();
           
            if ( (int)roundTAS->value()==1 ){
                logDebug("\nTAS initial state.\n");
                InitialStateTAS((int)in5->value(),(int)roundTASRAN->value(),(int)in3->value());
                haveloadedstate=0;
                in3->activate();
                in1->activate();
            }
            else if((int)roundRAN->value()==1){
                logDebug("\nRAN initial state.\n");
                InitialStateRAN((int)in6->value(), (int)in3->value());
                haveloadedstate=0;
                in3->activate();
                in1->activate();
            }
            
            else if((int)roundSTAT->value()==1){
                logDebug("\nSTAT initial state.\n");
                InitialStateSTAT((int)in3->value());
                haveloadedstate=0;
                in3->activate();
                in1->activate();
            }
            
            
            else if(havepath==0){
                InitialStateTAS(0,0,(int)in3->value()); roundTAS->value(1);
                haveloadedstate=0;
                in3->activate();
                in1->activate();
            }
            else if(haveloadedstate==1){
                InitialStateLOAD();
                in1->value(particles);
                in3->value(totrun);
                maxtime=(int)in4->value();
            }
            
            
            
            //set OUTPUT
            openout();
            
            cleared=0;}
        
        runningcontrol=1;
        bload->deactivate();
        bsave->deactivate();
        bsaveas->deactivate();
        bsetlayout->deactivate();
        
        /*in1->deactivate();
         in2->deactivate();
         in3->deactivate();
         in4->deactivate();
         activatesteady->deactivate();
         pii->deactivate();
         initstategroup->deactivate();*/
        
        settingsgroup->deactivate();
        clearbutton->deactivate();
        turbobutton->deactivate();
    }
    else {
        runningcontrol=0;
        
        bload->activate();
        if(havepath==1) bsave->activate();
        bsaveas->activate();
        bsetlayout->activate();
        
        clearbutton->activate();
    }
    
}



//--------------------------------------------
void runcontrolcb(Fl_Widget *, void *) {
    
    amstepping=0;
    run();
    
}


//--------------------------------------------
void runstepecb(Fl_Widget *, void *) {
    
    amstepping=1;
    step=instep->value();
    
    tickstep=0;
    run();
    
    
    
    
}



//--------------------------------------------
void clearcb(Fl_Widget *, void *) {
    
    clear();
}




void clear(){
    
    settingsgroup->activate();
    turbobutton->activate();
    
    runbutton->value(0);
    runningcontrol=runbutton->value();
    
    pii->value(0);
    ispii=pii->value();
    
    in5->maximum(igraph_vcount(&graph)-1);
    
    deltat=(double)indt->value();
    
 
    if ( (int)roundTAS->value()==1 ){
        InitialStateTAS((int)in5->value(),(int)roundTASRAN->value(),(int)in3->value());
        haveloadedstate=0;
        in3->activate();
        in1->activate();
    }
    else if((int)roundRAN->value()==1){
        InitialStateRAN((int)in6->value(), (int)in3->value());
        haveloadedstate=0;
        in3->activate();
        in1->activate();
    }
    
    else if((int)roundSTAT->value()==1){
        InitialStateSTAT((int)in3->value());
        haveloadedstate=0;
        in3->activate();
        in1->activate();
    }
    
    
    else if(havepath==0){
        InitialStateTAS(0,0,(int)in3->value()); roundTAS->value(1);
        haveloadedstate=0;
        in3->activate();
        in1->activate();
    }
    else if(haveloadedstate==1){
        InitialStateLOAD();
        in1->value(particles);
        in3->value(totrun);
        maxtime=(int)in4->value();
    }
    
    
    if(haveout==1){closeout(); haveout=0;}
    
    ticks=0;
    isrelaxed=0;
    havecstop=0;
    cleared=1;
    
    
}

//--------------------------------------------
void drawcb(Fl_Widget *, void *) {
    if (drawingcontrol==0) {drawingcontrol=1; drawgroup->activate();}
    else {drawingcontrol=0;drawgroup->deactivate();}
    
}


//--------------------------------------------
void drawnodescb(Fl_Widget *, void *) {
    if (drnodes==0) drnodes=1;
    else drnodes=0;
    
}

//--------------------------------------------
void drawlinkscb(Fl_Widget *, void *) {
    if (drlinks==0) drlinks=1;
    else drlinks=0;
    
}

//--------------------------------------------
void drawfluxescb(Fl_Widget *, void *) {
    if (drfluxes==0) drfluxes=1;
    else drfluxes=0;
    
}





//--------------------------------------------
void loadcb(Fl_Widget *, void *) {
    
    char inpath[5000];
    
    graphisloaded=0;
    
    if ( ! loadchooser ) {
        loadchooser = new Fl_File_Chooser("", "", Fl_File_Chooser::SINGLE, "");
    }
    
    loadchooser->directory("./");
    loadchooser->filter(NULL);
    loadchooser->label("Load Network");
    loadchooser->preview(0);
    loadchooser->show();
    
    // Block until user picks something.
    while(loadchooser->shown()) Fl::wait();
    
    // Print the results
    if (  loadchooser->value() == NULL ) {
        logDebug("\nNo file selected: network not loaded.\n");
        graphisloaded=1;
        return;
    }
    
    
    sprintf(inpath,"%s",loadchooser->value());
    
    //open graph;
    loadnetwork(inpath);
    
    
    if(graphisloaded==1){
        
        logDebug("\nFile correctly Loaded.\n");
        
        sprintf(path,"%s",loadchooser->value());
        havepath=1; rewrite=1;

        igraph_matrix_init(&loss, 0, 0);
        igraph_vector_init(&gain,nodesnumber);
        
        
        activatesteady->value(usesteady);
        round1->setonly();
        in3->value(totrun);
        if(usesteady==0){bottomgroup->deactivate();}
        else {bottomgroup->activate();}
        
        newgraph=1;
        
        roundTAS->value(0);
        roundRAN->value(0);
        roundSTAT->value(0);
        
        bsave->activate();
        
        if (haveloadedstate==1) {
            InitialStateLOAD();
            in3->deactivate();
            in1->deactivate();
        }
        
        else {roundTAS->value(1);beginner=0; particles=100;}
        
        
        in1->value(particles);
        
        clear();
        
        
    }
    
    
    
    
}


//--------------------------------------------

void saveascb(Fl_Widget *, void *) {
    
    char filepath[100];
    
    if ( ! savechooser ) {
        savechooser = new Fl_File_Chooser("", "*", Fl_File_Chooser::CREATE, "");
    }
    
    savechooser->directory("./");
    savechooser->filter(NULL);
    savechooser->preview(0);
    savechooser->label("Save As");
    
    savechooser->show();
    
    // Block until user picks something.
    while(savechooser->shown()) Fl::wait();
    
    // Print the results
    if (  savechooser->value() == NULL ) {
        logDebug("\nNo path selected.\n");
        return;
    }
    
    sprintf(filepath,"%s",savechooser->value());
    
    //save graph;
    savenetwork(filepath);
    
    if(havepath==1){
        logDebug("\nDone. Saved as ' %s '.\n", filepath);
        sprintf(path,filepath);
        rewrite=1;
        bsave->activate();
    }
    
    
}



//--------------------------------------------
void savecb(Fl_Widget *, void *) {
    
    savenetwork(path);
    logDebug("\nSaved.\n");
    
    
}



void changelayoutcb(Fl_Widget *, void *){
    
    ChangeLayout();
    
}


void turbocb(Fl_Widget *, void *){
    
    printdatabutton->value(1);
    TurboDial();
    
}



void drawhistocb(Fl_Widget *, void *) {
    if (drawhisto==0) {drawhisto=1; datascene->activate();}
    else {drawhisto=0;datascene->deactivate();}
    
}



void printhistodatacb(Fl_Widget *,void *){
    
    char filename[100];
    
    sprintf(filename,"%s%i.txt",setoutname->value(),4);
    
    if( (output4=fopen(filename,"a")) ==NULL) {
        logDebug("\nERROR: CANNOT OPEN OUTPUT FILE '%s'\n",filename);
        sprintf(errorstring,"ERROR\nOPENING\nOUTPUT\nFILE");
        error=1;
        rewrite=1;
        return;
    }
    
    
    fprintf(output4,"%i %i %f %f %f\n", particles, nodesnumber, atof(meanbuff->text()), atof(varbuff->text()), atof(meanbuff->text())/sqrt((float)nodesnumber) );
    
    fclose(output4);
    
}



void cstopcb(Fl_Widget *, void *){
    
    particles=in1->value();
    totrun=in3->value();
    cstop();
    
}

void cstop(){

    int ntry=NTRYSTOP;
    
    double resmean, resvar, resvarvar;
    
    igraph_vector_t trymean, tryvar;
    igraph_vector_init(&trymean,ntry); igraph_vector_init(&tryvar,ntry);
    igraph_vector_null(&trymean); igraph_vector_null(&tryvar);
    
    
    //build the ntrys steady states:
    for (int i=0; i<ntry; ++i) {
        
    igraph_matrix_t trystate;
    igraph_matrix_init(&trystate,nodesnumber,totrun);
    igraph_matrix_null(&trystate);
        
        //simulate the state
        for(int r=0;r<totrun;++r ){
            for (int p=0; p<particles; ++p) {
                double coin;
                int selnode;
                do{
                    coin=genrand_real1();
                    selnode=genrand_int31()%nodesnumber;
                }
                while(VECTOR(statstate)[selnode]<coin);
                
                ++MATRIX(trystate,selnode,r);
                
            }
        }
        
        
        //fix and generate histogram
        
        igraph_matrix_scale(&trystate,1./particles);
        
        igraph_vector_t avgtrystate, anothertemp;
        igraph_vector_init(&avgtrystate,0); igraph_vector_init(&anothertemp,0);
        igraph_matrix_get_col(&trystate,&avgtrystate,0);
        igraph_vector_null(&avgtrystate);
        
        for(int i=0;i<totrun;++i){
            igraph_matrix_get_col(&trystate,&anothertemp,i);
            igraph_vector_div(&anothertemp,&statstate);
            igraph_vector_add(&avgtrystate,&anothertemp);
        }
        igraph_vector_scale(&avgtrystate,1./totrun);
        
        //histogram
        LAURA_Histogram_1D tryhist;
        tryhist.CreateFromArrayMinMax(anothertemp, 40,0.,4., 1);
        
        //get from histogram
        VECTOR(trymean)[i]=tryhist.mean;
        VECTOR(tryvar)[i]=tryhist.variance;
               
        tryhist.Clear();
        //delete &tryhist;
        
        igraph_vector_destroy(&avgtrystate);
        igraph_vector_destroy(&anothertemp);
        igraph_matrix_destroy(&trystate);
    }
    
    
    resmean=igraph_vector_sum(&trymean)/ntry;
    resvar=igraph_vector_sum(&tryvar)/ntry;
    
    resvarvar=0;
    for (int i=0; i<ntry; ++i) {
        resvarvar= resvarvar + (VECTOR(tryvar)[i]-resvar)*(VECTOR(tryvar)[i]-resvar);
    }
    resvarvar=sqrt(resvarvar)/ntry;

    igraph_vector_destroy(&trymean);
    igraph_vector_destroy(&tryvar);
    
    inmeanstop->value(resmean);
    invarstop->value(resvar);
    inmeanerstop->value((3*resvar)/2);
    invarerstop->value((3*resvarvar)/2);
    
    

}





void testcb(Fl_Widget *, void *){
    
    testda();
    
}



/************************************** MAIN WINDOWD DRAW **********************************/


//-------------------------------------------------------------------------------------------------
void CreateMyWindow(void) {
    
    char *s;
    s= new char[100];
    
    int w_est,h_est;
    
    int ypos=0;
    int xpos=0;
    int widgh=0;
    int widgw=0;
    
    w_est= LEFT_SPACE+SCREEN_WIDTH+RIGHT_SPACE;
    h_est= 23+SCREEN_HEIGHT+150;
    
    form = new Fl_Window(w_est,h_est,"LAURA");
    
    mainwindow = new Fl_Group(0,0,w_est,h_est,"");
    
    
    
    new Fl_Box(FL_DOWN_FRAME,LEFT_SPACE,60,SCREEN_WIDTH+6,SCREEN_HEIGHT+6,"");
    
    scene =  new Frame(LEFT_SPACE+3,63,SCREEN_WIDTH,SCREEN_HEIGHT, 0);
    
    
    
    
    
    
    
    
    
    
    
    //---------- LEFT COLUMN I ----------------------------------------------------------------
    ypos=10;
    xpos=20;
    
    
    widgh=BUTTON_H1;
    Fl_Box *dynbox = new Fl_Box(xpos,ypos, BUTTON_WL,widgh, " Dynamics ");
    dynbox->labelsize(16);
    dynbox->labelfont(FL_BOLD);
    dynbox->box(FL_ROUNDED_BOX);
    ypos=ypos+widgh;
    
    ypos=ypos+30;
    settingsgroup = new Fl_Group(xpos-5,ypos-5,BUTTON_WL+10,600, "Settings");
    settingsgroup->align(FL_ALIGN_TOP_LEFT);
    settingsgroup->labelfont(FL_BOLD);
    settingsgroup->box(FL_BORDER_BOX);
    settingsgroup->color(FL_LIGHT2);

    widgh=20;
    ypos=ypos+15;
    in1 = new Fl_Value_Input(xpos,ypos,BUTTON_WL,widgh, "Particles:");
    in1->align(FL_ALIGN_TOP_LEFT);
    in1->step(1);
    in1->minimum(1);
    in1->value(particles);
    ypos=ypos+widgh;
    
    widgh=20;
    ypos=ypos+25;
    in2 = new Fl_Value_Input(xpos,ypos,BUTTON_WL,widgh, "Constraint:");
    in2->align(FL_ALIGN_TOP_LEFT);
    in2->step(1);
    in2->minimum(1);
    in2->value(vincolo);
    ypos=ypos+widgh;
    
    widgh=20;
    ypos=ypos+25;
    inthreshold = new Fl_Value_Input(xpos,ypos,BUTTON_WL,widgh, "Threshold:");
    inthreshold->align(FL_ALIGN_TOP_LEFT);
    inthreshold->step(1);
    inthreshold->minimum(0);
    inthreshold->value(threshold);
    ypos=ypos+widgh;
    
    
    widgh=20;
    ypos=ypos+25;
    in3 = new Fl_Value_Input(xpos,ypos,BUTTON_WL,widgh, "Contemp. runs:");
    in3->align(FL_ALIGN_TOP_LEFT);
    in3->step(1);
    in3->minimum(1);
    in3->value(totrun);
    ypos=ypos+widgh;
    
    widgh=20;
    ypos=ypos+25;
    indt = new Fl_Value_Input(xpos,ypos,BUTTON_WL,widgh, "Delta t:");
    indt->align(FL_ALIGN_TOP_LEFT);
    indt->step(0.001);
    indt->minimum(0.0001);
    indt->maximum(2);
    indt->value(1);
    ypos=ypos+widgh;
    
    widgh=20;
    ypos=ypos+25;
    in4 = new Fl_Value_Input(xpos,ypos,BUTTON_WL,widgh, "Max Time:");
    in4->align(FL_ALIGN_TOP_LEFT);
    in4->step(1);
    in4->minimum(1);
    in4->value(maxtime);
    ypos=ypos+widgh;
    
    ypos=ypos+5;
    widgh=BUTTON_H1;
    activatesteady = new Fl_Light_Button(xpos, ypos, BUTTON_WL, widgh, "St. State");
    activatesteady->color(FL_LIGHT2);
    ypos=ypos+widgh;
    
    
    widgh=BUTTON_H1;
    pii = new Fl_Round_Button(xpos, ypos, BUTTON_WL, widgh, "P_ii = 1/2");
    pii->value(0);
    ypos=ypos+widgh;
    
    ypos=ypos+15;
    initstategroup = new Fl_Group(xpos,ypos,BUTTON_WL,200, "Initial State:");
    initstategroup->align(FL_ALIGN_TOP_LEFT);
    initstategroup->labelfont(FL_BOLD);
    //initstategroup->box(FL_FLAT_BOX);
    
    widgh=20;
    roundTAS = new Fl_Round_Button(xpos+10,ypos,BUTTON_WL-10,widgh,"TAS");
    roundTAS->type(FL_RADIO_BUTTON);
    roundTAS->setonly();
    ypos=ypos+widgh;
    
    roundTAS->setonly();
    
    widgh=20;
    ypos=ypos+10;
    in5 = new Fl_Value_Input(xpos+20,ypos,BUTTON_WL-20,widgh, "Initial node");
    in5->align(FL_ALIGN_TOP_LEFT);
    in5->step(1);
    in5->minimum(0);
    in5->maximum(igraph_vcount(&graph)-1);
    in5->value(beginner);
    ypos=ypos+widgh;
    
    widgh=20;
    roundTASRAN = new Fl_Check_Button(xpos+20,ypos,BUTTON_WL-10,widgh,"Random");
    ypos=ypos+widgh;
    
    
    widgh=20;
    ypos=ypos+10;
    roundRAN = new Fl_Round_Button(xpos+10,ypos,BUTTON_WL-10,widgh,"Random");
    roundRAN->type(FL_RADIO_BUTTON);
    roundRAN->setonly();
    ypos=ypos+widgh;
    
    
    widgh=20;
    ypos=ypos+10;
    in6 = new Fl_Value_Input(xpos+20,ypos,BUTTON_WL-20,widgh, "Rand. Seed");
    in6->align(FL_ALIGN_TOP_LEFT);
    in6->step(1);
    in6->minimum(1);
    in6->value(0);
    ypos=ypos+widgh;
    
    
    widgh=20;
    ypos=ypos+10;
    roundSTAT = new Fl_Round_Button(xpos+10,ypos,BUTTON_WL-10,widgh,"Stat State");
    roundSTAT->type(FL_RADIO_BUTTON);
    roundSTAT->deactivate();
    ypos=ypos+widgh;
    
    
    widgh=20;
    ypos=ypos+10;
    Fl_Button *prev  = new Fl_Button(xpos+10, ypos, BUTTON_WL-10, widgh, "Preview");
    prev->color(FL_LIGHT2);
    ypos=ypos+widgh;
    
    
    initstategroup->end();
    
    
    
    
    widgh=20;
    ypos=ypos+40;
    setoutname = new Fl_Input(xpos,ypos,BUTTON_WL,widgh, "Output file:");
    setoutname->align(FL_ALIGN_TOP_LEFT);
    setoutname->labelfont(FL_BOLD);
    setoutname->value("myout");
    ypos=ypos+widgh;
    
    settingsgroup->end();
    
    
    
    
    
    //---------- LEFT COLUMN II ----------------------------------------------------------------
    ypos=10;
    xpos=20+BUTTON_WL+40;
    
    
    widgh=BUTTON_H1;
    Fl_Box *controlbox = new Fl_Box(xpos,ypos, BUTTON_WL,widgh, " Controls ");
    controlbox->labelsize(16);
    controlbox->labelfont(FL_BOLD);
    controlbox->box(FL_ROUNDED_BOX);
    ypos=ypos+widgh;
    
    
    ypos=ypos+20;
    widgh=25;
    tickbuff = new Fl_Text_Buffer();
    tickdisp = new Fl_Text_Display(xpos,ypos, BUTTON_WR,widgh, "Ticks");
    tickdisp->align(FL_ALIGN_TOP_LEFT);
    tickdisp->buffer(tickbuff);
    tickdisp->textsize(12);
    
    ypos=ypos+widgh;
    
    ypos=ypos+10;
    
    widgh=BUTTON_H;
    runbutton  = new Fl_Light_Button(xpos, ypos, BUTTON_WL, widgh, "RUN");
    ypos=ypos+widgh;
    
    ypos=ypos+5;
    
    widgh=BUTTON_H1;
    runstepbutton  = new Fl_Button(xpos, ypos, BUTTON_WL, widgh, "Run Step");
    ypos=ypos+widgh;
    
    widgh=20;
    ypos=ypos+20;
    instep = new Fl_Value_Input(xpos+20,ypos,BUTTON_WL-20,widgh, "Step:");
    instep->align(FL_ALIGN_TOP_LEFT);
    instep->step(1);
    instep->minimum(0);
    instep->maximum(maxtime-ticks);
    instep->value(1);
    ypos=ypos+widgh;
    
    
    ypos=ypos+10;
    
    widgh=BUTTON_H1;
    printdatabutton = new Fl_Check_Button(xpos,ypos,BUTTON_WL,widgh,"Print Data");
    printdatabutton->value(1);
    ypos=ypos+widgh;
    
    widgh=BUTTON_H;
    ypos=ypos+5;
    clearbutton  = new Fl_Button(xpos, ypos, BUTTON_WL, widgh, "CLEAR");
    ypos=ypos+widgh;
    
    
    ypos=ypos+10;
    
    widgh=BUTTON_H;
    drawbutton  = new Fl_Light_Button(xpos, ypos, BUTTON_WL, widgh, "DRAW");
    if (initdrawbut==0) {initdrawbut=1; if(drawingcontrol==1) drawbutton->set();}
    ypos=ypos+widgh;
    
    
    ypos=ypos+20;
    drawgroup = new Fl_Group(xpos-5,ypos,BUTTON_WL+10,205, "Draw Options");
    drawgroup->align(FL_ALIGN_TOP_LEFT);
    drawgroup->labelfont(FL_BOLD);
    drawgroup->box(FL_BORDER_BOX);
    drawgroup->color(FL_DARK1);
    
    ypos=ypos+5;
    widgh=BUTTON_H1;
    drawnodes = new Fl_Check_Button(xpos,ypos,BUTTON_WL,widgh,"Draw Nodes");
    drawnodes->value(1);
    if(drnodes==1) drawnodes->set();
    
    ypos=ypos+widgh;
    widgh=BUTTON_H1;
    drawlinks = new Fl_Check_Button(xpos,ypos,BUTTON_WL,widgh,"Draw Links");
    drawlinks->value(1);
    if(drlinks==1) drawlinks->set();
    ypos=ypos+widgh;
    
    widgh=BUTTON_H1;
    drawfluxes = new Fl_Check_Button(xpos,ypos,BUTTON_WL,widgh,"Draw Fluxes");
    drawfluxes->value(0);
    ypos=ypos+widgh;
    
    ypos=ypos+20;
    
    palettegroup = new Fl_Group(xpos,ypos,BUTTON_WL,3*widgh, "Palette");
    palettegroup->align(FL_ALIGN_TOP_LEFT);
    palettegroup->labelfont(FL_BOLD);
    
    widgh=BUTTON_H1;
    palette0 = new Fl_Round_Button(xpos,ypos,BUTTON_WL,widgh,"Standard");
    palette0->type(FL_RADIO_BUTTON);
    palette0->setonly();
    ypos=ypos+widgh;
    
    widgh=BUTTON_H1;
    palette1 = new Fl_Round_Button(xpos,ypos,BUTTON_WL,widgh,"Sqrt");
    palette1->type(FL_RADIO_BUTTON);
    ypos=ypos+widgh;
    
    widgh=BUTTON_H1;
    palette2 = new Fl_Round_Button(xpos,ypos,BUTTON_WL,widgh,"DoubleSqrt");
    palette2->type(FL_RADIO_BUTTON);
    ypos=ypos+widgh;
    
    palettegroup->end();
    
    drawgroup->end();
    
    
    
    
    

    widgh=BUTTON_H;
    widgw=BUTTON_WL;
    ypos=h_est-widgh-10;
    turbobutton  = new Fl_Button(xpos, ypos, widgw, widgh, "TURBO");
    turbobutton->labelfont(FL_BOLD);
    turbobutton->labelcolor(FL_DARK_RED);
    ypos=ypos+widgh;
    
    
    
    
    
    
    //---------- RIGHT COLUMN I ----   R.C.1 ------------------------------------------------------------
    
    ypos=10;
    xpos=LEFT_SPACE+SCREEN_WIDTH+20;
    
    
    widgh=BUTTON_H1;
    Fl_Box *netbox = new Fl_Box(xpos,ypos,BUTTON_WR,widgh, " Network ");
    netbox->labelsize(16);
    netbox->labelfont(FL_BOLD);
    netbox->box(FL_ROUNDED_BOX);
    ypos=ypos+widgh;
    
    widgh=3*BUTTON_H;
    
    ypos=ypos+10;
    
    databuff = new Fl_Text_Buffer();
    datadisp = new Fl_Text_Display(xpos,ypos, BUTTON_WR,widgh, "");
    datadisp->buffer(databuff);
    ypos=ypos+widgh;

    
    ypos=ypos+40;
    
    widgh=BUTTON_H;
    bsetlayout = new Fl_Button(xpos, ypos, BUTTON_WR, widgh, "Layout");
    ypos=ypos+widgh;
    
    /* widgh=BUTTON_H1;
     ypos=ypos+10;
     Fl_Box *box = new Fl_Box(xpos, ypos, BUTTON_WR, widgh, "New Network:");
     box->labelsize(16);
     ypos=ypos+widgh;
     
     widgh=BUTTON_H1;
     bnewlattice = new Fl_Button(xpos, ypos, BUTTON_WR, widgh, "New Lattice");
     ypos=ypos+widgh;
     
     widgh=BUTTON_H1;
     bnewrandom = new Fl_Button(xpos, ypos, BUTTON_WR, widgh, "New Random");
     ypos=ypos+widgh;
     
     ypos=ypos+5;
     widgh=BUTTON_H1;
     bnewclustered = new Fl_Button(xpos, ypos, BUTTON_WR, widgh, "clustered");
     ypos=ypos+widgh;   */
    
    
    ypos=ypos+50;
    widgh=BUTTON_H;
    bnewnet = new Fl_Button(xpos, ypos, BUTTON_WR, widgh, "New Network");
    ypos=ypos+widgh;
    
    /*widgh=BUTTON_H1;
    ypos=ypos+10;
    Fl_Box *box1 = new Fl_Box(xpos, ypos, BUTTON_WR, widgh, "Load/Save Net.:");
    box1->labelsize(16);
    ypos=ypos+widgh;
    */
    ypos=ypos+10;
    widgh=BUTTON_H1;
    bload = new Fl_Button(xpos, ypos, BUTTON_WR, widgh, "Load Network");
    ypos=ypos+widgh;

    
    widgh=BUTTON_H1;
    bsave = new Fl_Button(xpos, ypos, BUTTON_WR, widgh, "Save Network");
    bsave->deactivate();
    ypos=ypos+widgh;
    
    
    widgh=BUTTON_H1;
    bsaveas = new Fl_Button(xpos, ypos, BUTTON_WR, widgh, "Save Net. As");
    ypos=ypos+widgh;
    


    
    
    
    widgh=BUTTON_H;
    widgw=BUTTON_WR;
    ypos=h_est-widgh-10;
    buttonexit = new Fl_Button(xpos, ypos,widgw, widgh, "EXIT");
    buttonexit->labelfont(FL_BOLD);
    ypos=ypos+widgh;
    
    
    
    
    
    //---------- RIGHT COLUMN II -------- R.C.2 --------------------------------------------------------
    ypos=10;
    xpos=LEFT_SPACE+SCREEN_WIDTH+20+COLWIDTH+10;
    
    
    widgh=BUTTON_H1;
    Fl_Box *mydatabox = new Fl_Box(xpos,ypos,BUTTON_L,widgh, " Data ");
    mydatabox->labelsize(16);
    mydatabox->labelfont(FL_BOLD);
    mydatabox->box(FL_ROUNDED_BOX);
    ypos=ypos+widgh;
    
    ypos=ypos+10;
    
    
    histogroup = new Fl_Group(xpos-5,ypos-5,BUTTON_L+10,530, "");
    histogroup->align(FL_ALIGN_TOP_LEFT);
    histogroup->labelfont(FL_BOLD);
    histogroup->box(FL_BORDER_BOX);
    histogroup->color(FL_LIGHT2);
    
    ypos=ypos+5;
    
    widgh=BUTTON_L; widgw=BUTTON_L;
    datascene =  new DataFrame(xpos,ypos,widgw,widgh, 0);
    ypos=ypos+widgh;
    
    ypos=ypos+5;
    widgh=BUTTON_H1;  widgw=BUTTON_L;
    drawhistobutton = new Fl_Check_Button(xpos,ypos,BUTTON_L,widgh,"Draw Histo");
    drawhistobutton->value(1);    drawhisto=1;
    ypos=ypos+widgh;
    
    ypos=ypos+10;
    widgw=BUTTON_L-10;
    meanbuff = new Fl_Text_Buffer();
    meandisp = new Fl_Text_Display(xpos, ypos, widgw, 25, "Mean");
    meandisp->buffer(meanbuff);
    meandisp->textsize(10);
    ypos=ypos+widgh;
    
    ypos=ypos+10;
    widgw=BUTTON_L-10;
    varbuff = new Fl_Text_Buffer();
    vardisp = new Fl_Text_Display(xpos, ypos, widgw, 25, "Variance");
    vardisp->buffer(varbuff);
    vardisp->textsize(10);
    ypos=ypos+widgh;
    
    ypos=ypos+5;
    widgh=BUTTON_H1;
    widgw=BUTTON_WR-40;
    printbutton1 = new Fl_Button(xpos+(COLWIDTHLARGE-widgw),ypos,widgw,widgh, "Print ^");
    ypos=ypos+widgh;
    
    
   ypos=ypos+10;
    widgh=30;
    xpos=xpos-5;
    roundstop = new Fl_Round_Button(xpos, ypos, widgw, widgh, "STOP WHEN STEADY");
    ypos=ypos+widgh;
    
    xpos=xpos+20;
    widgh=20;
    widgw=BUTTON_WL-20;
    ypos=ypos+20;
    inmeanstop = new Fl_Value_Input(xpos,ypos,widgw,widgh, "STOP-Mean:");
    inmeanstop->align(FL_ALIGN_TOP_LEFT);
    inmeanstop->step(0.0001);
    inmeanstop->minimum(0);
    inmeanstop->maximum(10);
    inmeanstop->value(1.);
    
    
    widgh=20;
    xpos=xpos+10+widgw;
    inmeanerstop = new Fl_Value_Input(xpos,ypos,widgw,widgh, " (error):");
    inmeanerstop->align(FL_ALIGN_TOP_LEFT);
    inmeanerstop->step(0.0001);
    inmeanerstop->minimum(0);
    inmeanerstop->maximum(10);
    inmeanerstop->value(0.005);
    ypos=ypos+widgh;
    xpos=xpos-10-widgw;
    
    widgh=20;
    ypos=ypos+20;
    invarstop = new Fl_Value_Input(xpos,ypos,widgw,widgh, "STOP-Var:");
    invarstop->align(FL_ALIGN_TOP_LEFT);
    invarstop->step(0.0001);
    invarstop->minimum(0);
    invarstop->maximum(10);
    invarstop->value(0.1);
    
    
    widgh=20;
    xpos=xpos+10+widgw;
    invarerstop = new Fl_Value_Input(xpos,ypos,widgw,widgh, " (error):");
    invarerstop->align(FL_ALIGN_TOP_LEFT);
    invarerstop->step(0.0001);
    invarerstop->minimum(0);
    invarerstop->maximum(10);
    invarerstop->value(0.005);
    ypos=ypos+widgh;
    xpos=xpos-10-widgw;
    
    ypos=ypos+10;
    widgh=BUTTON_H1;
    widgw=widgw+widgw+10;
    cstopbutton = new Fl_Button(xpos,ypos,widgw,widgh, "Calculate stop cond.");
    ypos=ypos+widgh;
    
    
     /* ypos=ypos+20;
     widgw=BUTTON_L-10;
     dmeanbuff = new Fl_Text_Buffer();
     dmeandisp = new Fl_Text_Display(xpos, ypos, widgw, 25, "Delta Mean");
     dmeandisp->buffer(dmeanbuff);
     dmeandisp->textsize(10);
     ypos=ypos+widgh;
    
     ypos=ypos+20;
     widgw=BUTTON_L-10;
     dvarbuff = new Fl_Text_Buffer();
     dvardisp = new Fl_Text_Display(xpos, ypos, widgw, 25, "Delta Var");
     dvardisp->buffer(dvarbuff);
     dvardisp->textsize(10);
     ypos=ypos+widgh;*/
    
    
    
    histogroup->end();
    if(usesteady==0){histogroup->deactivate();}
    
    
    
    ypos=ypos+20;
    widgh=BUTTON_H1;
    widgw=BUTTON_L-20;
    Fl_Button * fluctwinbutton = new Fl_Button(xpos,ypos,widgw,widgh, "Fluctuations");
    ypos=ypos+widgh;
    
    
    
    /*
    
    ypos=ypos+40;
    widgh=BUTTON_H1;
    test = new Fl_Button(xpos,ypos,BUTTON_WR,widgh, " TEST ");
    ypos=ypos+widgh;
     
     test->callback(testcb,0);
     
    */
    
    drawhistobutton->callback(drawhistocb,0);
    
    
    
    
    
    
    
    
    
    //----------------- BOTTOM ----------------
    
    
    widgh=BUTTON_H;
    ypos=63+SCREEN_HEIGHT+10;
    xpos=LEFT_SPACE;
    
    widgh=BUTTON_H1;
    ypos=ypos+25;
    
    bottomgroup = new Fl_Group(xpos,ypos,SCREEN_WIDTH,widgh, "Display:");
    bottomgroup->align(FL_ALIGN_TOP_LEFT);
    bottomgroup->labelfont(FL_BOLD);
    bottomgroup->box(FL_BORDER_BOX);
    bottomgroup->color(FL_LIGHT2);
    
    widgw=BUTTON_WL;
    xpos=xpos+60;
    round1 = new Fl_Round_Button(xpos,ypos,widgw,widgh,"Actual State");
    round1->type(FL_RADIO_BUTTON);
    round1->setonly();
    xpos=xpos+widgw;
    
    widgw=BUTTON_WL;
    xpos=xpos+5;
    round2 = new Fl_Round_Button(xpos,ypos,widgw,widgh,"Steady State");
    round2->type(FL_RADIO_BUTTON);
    xpos=xpos+widgw;
    
    widgw=BUTTON_WL;
    xpos=xpos+5;
    round3 = new Fl_Round_Button(xpos,ypos,widgw,widgh,"Difference");
    round3->type(FL_RADIO_BUTTON);
    xpos=xpos+widgw;
    
    
    bottomgroup->end();
    bottomgroup->deactivate();
    
    
    ypos=30;
    xpos=LEFT_SPACE+1;
    pathbuff = new Fl_Text_Buffer();
    pathdisp = new Fl_Text_Display(xpos, ypos, SCREEN_WIDTH, 20, "Network file path:");
    pathdisp->buffer(pathbuff);
    pathdisp->textsize(10);
    
    
    buttonexit->callback(exitcb,0);
    runbutton->callback(runcontrolcb,0);
    runstepbutton->callback(runstepecb,0);
    clearbutton->callback(clearcb,0);
    
    drawbutton->callback(drawcb,0);
    drawnodes->callback(drawnodescb,0);
    drawlinks->callback(drawlinkscb,0);
    drawfluxes->callback(drawfluxescb,0);
    
    activatesteady->callback(activatesteadycb,0);
    pii->callback(piicb,0);
    prev->callback(prevcb,0);
    
    turbobutton->callback(turbocb,0);
    
    /* bnewlattice->callback(newlatcb,0);
     bnewrandom->callback(newrandcb,0);
     bnewclustered->callback(newcluscb,0);*/
    bnewnet->callback(newnetcb,0);
    bload->callback(loadcb,0);
    bsaveas->callback(saveascb,0);
    bsave->callback(savecb,0);
    
    bsetlayout->callback(changelayoutcb,0);
    
    printbutton1->callback(printhistodatacb,0);
    
    cstopbutton->callback(cstopcb,0);
    
    mainwindow->end();
    
    form->end();
    form->show();
    scene->show();
    datascene->show();
}
//-------------------------------------------------------------------------------------------------


