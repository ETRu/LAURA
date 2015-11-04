#pragma once

#if defined (WIN32)
#include "igraph/igraph.h"
#else
#include "igraph.h"
#endif

#include "draw.h"

#include <string>
#include <iostream>

#include <FL/Fl.H>
#include <FL/Fl_Window.H>
#include <FL/Fl_Text_Display.H>
#include <FL/Fl_Group.H>
#include <FL/Fl_Box.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Input.H>
#include <FL/Fl_Light_Button.H>
#include <FL/Fl_Check_Button.H>
#include <FL/Fl_Round_Button.H>
#include <FL/Fl_File_Chooser.H>
#include "Frame.h"
#include "dinamica.h"
#include "igraphutil.h"
#include "networkop.h"
#include "logdebug.h"
#include "dataanalysis.h"

#include <FL/Fl_Counter.H>
#include <FL/Fl_Input.H>
#include <FL/Fl_Value_Input.H>



//CHANGE LAYOUT dial
void ChangeLayout(void);
void gridcb(Fl_Widget *, void *);
void starcb(Fl_Widget *, void *);
void circlecb(Fl_Widget *, void *);
void mdscb(Fl_Widget *, void *);
void forcewlaycb(Fl_Widget *, void *);
void forceunwlaycb(Fl_Widget *, void *);
void resetlayoutcb(Fl_Widget *, void *);
void exitdial0cb(Fl_Widget *, void *);

//NEW NETWORK (GENERAL) dial
void newnetcb(Fl_Widget *, void *);
void DialogueNewNet(void);
void exitdialnncb(Fl_Widget *, void *);

// NEW LATTICE dial
void newlatcb(Fl_Widget *, void *);
void DialogueNewLat(void) ;
void generatelatticecb(Fl_Widget *, void *);
void exitdial1cb(Fl_Widget *, void *);

// NEW RANDOM dial
void newrandcb(Fl_Widget *, void *);
void DialogueNewRand(void) ;
void generaterandom1cb(Fl_Widget *, void *);
void exitdial2cb(Fl_Widget *, void *);


// NEW clustered dial
void newcluscb(Fl_Widget *, void *);
void DialogueNewClus(void) ;
void generateclustsymcb(Fl_Widget *, void *);
void exitdial3cb(Fl_Widget *, void *);

// NEW clustered gerarchic 2 layers dial
void newclusger2cb(Fl_Widget *, void *);
void DialogueNewClusGer2(void) ;
void generateclusger2cb(Fl_Widget *, void *);
void exitdial4cb(Fl_Widget *, void *);


// MAIN WINDOW
void exitcb(Fl_Widget *, void *);

// load - save - save as
void loadcb(Fl_Widget *, void *);
void loadstatecb(Fl_Widget *, void *);
void saveascb(Fl_Widget *, void *);
void savecb(Fl_Widget *, void *);
void savestateascb(Fl_Widget *, void *);

//runs
void run();
void runcontrolcb(Fl_Widget *, void *);
void runstepecb(Fl_Widget *, void *);

//clears
void clearcb(Fl_Widget *, void *);
void clear();

//output
void openout(char *string);
void closeout();

//drawing
void drawcb(Fl_Widget *, void *);
void drawnodescb(Fl_Widget *, void *);
void drawlinkscb(Fl_Widget *, void *);
void drawfluxescb(Fl_Widget *, void *);

//steady state
void activatesteadycb(Fl_Widget *, void *);

//preview
void prevcb(Fl_Widget *, void *);

//change layout
void changelayoutcb(Fl_Widget *, void *);

//TURBO DIAL
void TurboDial(void);
void exitdialturbocb(Fl_Widget *, void *);
void goturbocb(Fl_Widget *, void *);
void goturborun();

//TURBO JOB DIAL
void TJobDial(void);
void exitdialtjobcb(Fl_Widget *, void *);
void gotjobcb(Fl_Widget *, void *);
void gotjob();

//turbo
void turbocb(Fl_Widget *, void *);

//turbo job
void tjobcb(Fl_Widget *, void *);

//calculate stop condition
void cstopcb(Fl_Widget *, void *);
void cstop();

//data window
void DataWindow(void);
void datacb(Fl_Widget *, void *);
void exitdatacb(Fl_Widget *, void *) ;


//----- GENERATE MAIN WINDOW
void CreateMyWindow(void);


