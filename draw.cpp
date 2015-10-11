
#include "draw.h"

#define CARTESIO  1
#define QUADRO    2
#define LINEA     3
#define HUEINDEX   4
#define MYARROW   5

#ifndef PI
#define PI 3.141593
#endif

#if defined (WIN32)
#define NOANIM
#else
#define NOANIM gl_draw("NO ANIMATION",5,8);
#endif

extern Fl_Text_Buffer  *tickbuff;

extern Fl_Text_Buffer  *meanbuff;
extern Fl_Text_Buffer  *varbuff;
extern Fl_Text_Buffer  *dmeanbuff;
extern Fl_Text_Buffer  *dvarbuff;


#define PRINTTIME tickbuff->text(textstring);


extern int drawingcontrol;
extern int drawhisto;
extern int drnodes;
extern int drlinks;
extern int drfluxes;

extern int ticks;
extern int maxtime;
extern int totrun;

extern int rewrite;

extern int particles;

extern int usesteady;

extern int runningcontrol;

extern int tactdisp;
extern int windowrm;

int ticked;         //has time passed?
int lastick;



using namespace std;


extern igraph_t graph;
extern igraph_t sgraph;
extern igraph_matrix_t layout;
extern igraph_matrix_t density; extern igraph_matrix_t densityold;
extern igraph_matrix_t flux;
extern igraph_vector_t statstate;



extern Fl_Round_Button *round1;
extern Fl_Round_Button *round2;
extern Fl_Round_Button *round3;
extern Fl_Round_Button *palette0;
extern Fl_Round_Button *palette1;
extern Fl_Round_Button *palette2;


//data analysis

extern LAURA_Histogram_1D histdist;
extern double hism, hisv;

extern igraph_vector_t tactvect;
extern igraph_vector_t fluctvect;
float tactcoord[MAXTACTREG*3];
extern LAURA_Histogram_1D histtact;
extern LAURA_Histogram_1D histfluct;


//Color (hue) Functions

double RedHue(double color){
    double hue;
    hue = ((-2.)*(color-1.0)*(color-1.0))+1.0;
    if(hue<0) hue=0;
    return hue;
}

double GreenHue(double color){
    double hue;
    hue=((-7.)*(color-0.4)*(color-0.4))+1.0;
    if(hue<0) hue=0;
    return hue;
}

double BlueHue(double color){
    double hue;
    //hue=((-100.)*(color-0.1)*(color-0.1))+1.0;
    hue=((-40.)*(color-0.15)*(color-0.15))+1.0;
    
    if(hue<0) hue=0;
    return hue;
}





/**********************************************                     ****************************************
 **********************************************                     ****************************************
 **********************************************     NETWORK DRAW    ****************************************
 **********************************************                     ****************************************
 **********************************************                     ****************************************
 */


void draw_init(void){
    float L=0, LT;
    L=(10./igraph_vcount(&graph));
    if(L<0.05) L=0.05;
    
    LT=L*0.5;
    
    glNewList(QUADRO, GL_COMPILE);
    glBegin(GL_QUADS);
    glVertex3f( L, L,0.0);
		  glVertex3f(-L, L,0.0);
		  glVertex3f(-L,-L,0.0);
		  glVertex3f( L,-L,0.0);
    glEnd();
    glEndList();
    
    L=0.5;
    glNewList(LINEA, GL_COMPILE);
    glBegin(GL_LINES);
    glColor3f(0.5,0.5,0.5);
    glVertex3f( -L, 0.0,0.0); glVertex3f( L,0.0,0.0);
    glColor3f(0.5,0.5,0.5);
    glEnd();
    glEndList();
    
    
    L=0.1;
    glNewList(HUEINDEX, GL_COMPILE);
    glBegin(GL_QUADS);
    glVertex3f( L, L,0.0);
		  glVertex3f(-L, L,0.0);
		  glVertex3f(-L,-L,0.0);
		  glVertex3f( L,-L,0.0);
    glEnd();
    glEndList();
    
    
    
    glNewList(MYARROW, GL_COMPILE);
    glBegin(GL_TRIANGLES);
    glVertex3f( LT, 0,0.0);
    glVertex3f(-LT, LT,0.0);
    glVertex3f(-LT, -LT,0.0);
    glEnd();
    glEndList();
    
}





void draw_scene(void){
    
    double val;
    char textstring[100];
    int ednum, vnum;
    int from, to;
    
    float arrx, arry, angle, dx, dy;
    int versus;
    
    double red, green, blue;
    
    /*****************        DRAW THE NETWORK      ********************/
    if(drawingcontrol==1){
        ednum=igraph_ecount(&sgraph);
        vnum=igraph_vcount(&sgraph);
        
        //glClear(GL_DEPTH_BUFFER_BIT);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        
        
        
        //------------------ DRAW HUE INDEX
        for (int i=0; i<100; ++i) {
            
            double pluto;
            
            pluto=(double)i/100.0;
            
            if((int)palette1->value()==1) pluto=sqrt(pluto);
            if((int)palette2->value()==1) pluto=pow(pluto,1./4.);
            else pluto=pluto;
            
            red= RedHue(pluto);
            green= GreenHue(pluto);
            blue= BlueHue(pluto);
            
            glPushMatrix();
            glTranslated(i*0.05,0.2,0);
            glColor3f(red, green, blue);
            glCallList(HUEINDEX);
            glPopMatrix();
        }
        
        
        //-------------- DRAW VERTICES
        
        if(drnodes==1){
            
            for (int i=0; i<vnum; ++i) {
                
                
                // -- view state
                if (round1->value()==1){
                    val=0;
                    for(int j=0; j<totrun; ++j){
                        val=val+(double)MATRIX(density,i,j);
                    }
                    val=val/totrun;
                }
                
                if(usesteady==1){ // -- view stationary state
                    if (round2->value()==1)
                        val=VECTOR(statstate)[i];
                    
                    // -- view difference
                    if (round3->value()==1){
                        val=0;
                        for(int j=0; j<totrun; ++j){
                            val=val+((VECTOR(statstate)[i]-MATRIX(density,i,j))*(VECTOR(statstate)[i]-MATRIX(density,i,j)));
                        }
                        val=sqrt(val)/totrun; }
                    
                }
                
                if((int)palette1->value()==1) val=sqrt(val);
                if((int)palette2->value()==1) val=sqrt(sqrt(val));
                else val=val;
                
                
                //set color
                red= RedHue(val);
                green= GreenHue(val);
                blue= BlueHue(val);
                
                
                
                
                //draw
                glPushMatrix();
                glTranslated(MATRIX(layout,i,0),MATRIX(layout,i,1),0);
                glColor3f(red, green, blue);
                glCallList(QUADRO);
                glPopMatrix();
                
            }
            
            
        }
        
        
        
        // glColor3f(0.5,0.5,0.5);
        
        //------ DRAW ARROWS
        
        if(drfluxes==1){
            for (int i=0; i<ednum; ++i) {
                
                igraph_edge(&sgraph, i, &from, &to);
                
                arrx=(MATRIX(layout,from,0)+MATRIX(layout,to,0))*0.5;
                arry=(MATRIX(layout,from,1)+MATRIX(layout,to,1))*0.5;
                
                dy=(MATRIX(layout,from,1)-MATRIX(layout,to,1));
                dx=(MATRIX(layout,from,0)-MATRIX(layout,to,0));
                
                //IF dx<<1 the link is more or less vertical (tan divergers)
                if((dx*dx)<0.001){angle=PI*0.5;}
                else {angle=tan(dy/dx);}
                
                //printf("\nangle=%f dx=%f/dy=%f", angle, dx, dy);
                
                
                val=0;
                for(int run=0;run<totrun;++run){
                    val=val+MATRIX(flux,i,run);
                }
                val=val/totrun;
                
                versus=0;
                if(val<0)versus=1;
                
                
                val=val/particles;
                val=fabs(val);
                
                if(val>0.00001){
                    
                    glPushMatrix();
                    glTranslated(arrx,arry,0);
                    //set direction
                    if(from>to) glRotated((versus*PI+angle)*(180.0/PI),0,0,1.);
                    if(from<to) glRotated((versus*PI+angle+PI)*(180.0/PI),0,0,1.);
                    //set color
                    val=sqrt(val);
                    red= RedHue(val);
                    green= GreenHue(val);
                    blue= BlueHue(val);
                    glColor3f(red,green,blue);
                    glCallList(MYARROW);
                    glPopMatrix();}
                
                
            }
            
        }
        
        
        //------ DRAW LINKS
        if(drlinks==1){
            glPushMatrix();
            glColor3f(0.5,0.5,0.5);
            
            for (int i=0; i<ednum; ++i) {
                
                igraph_edge(&sgraph, i, &from, &to);
                
                glBegin(GL_LINES);
                glVertex3d(MATRIX(layout,from,0),MATRIX(layout,from,1),0.0);
                glVertex3d(MATRIX(layout,to,0),MATRIX(layout,to,1),0.0);
                glEnd();
                
            }
            
            
            glPopMatrix();
            
        }
        
        
        
        // TICKS
        {sprintf(textstring, "%i", ticks);
        PRINTTIME}
        
        
    }
    
    
    
    /*****************        DON'T DRAW THE NETWORK      ********************/
    else{
        
        
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        
        glColor3f(1,1,1);
        gl_font(FL_HELVETICA,10);
       gl_draw("NO ANIMATION",5,8);
        
        
        if(rewrite==1){sprintf(textstring, "%i", ticks);
            PRINTTIME;}
        
    }
    
}






/**********************************************             ****************************************
 **********************************************             ****************************************
 **********************************************  DATA DRAW  ****************************************
 **********************************************             ****************************************
 **********************************************             ****************************************
 */




void draw_datainit(void){
    
    
}





void draw_datascene(void){
    
    int nbins=histdist.nbins;
    float binwidth=9./nbins;
    
    char textstring[100];
    char firstchar;
    
    float xtl, xtr, xbl, xbr, ytl, ytr, ybl, ybr;
    //maxmum is 8 points in height.
    
    
    if(drawhisto==1 && usesteady==1){
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
         glClearColor(1., 1., 1., 1);
    
    glPushMatrix();
    glTranslated(0,0,0);
        
        
        
        
        
        glColor3f(0,0,0);
        glBegin(GL_QUADS);                      // Draw A Quad
        glVertex3f(0, 11, 0.0f);              // Top Left
        glVertex3f(11, 11, 0.0f);              // Top Right
        glVertex3f(11,0, 0.0f);              // Bottom Right
        glVertex3f(0,0, 0.0f);              // Bottom Left
        glEnd();                            // Done Drawing The Quad
        
        
        glColor3f(1.,1.,1.);
        glBegin(GL_QUADS);                      // Draw A Quad
        glVertex3f(0.05, 10.95, 0.1f);              // Top Left
        glVertex3f(10.95, 10.95, 0.1f);              // Top Right
        glVertex3f(10.95,0.05, 0.1f);              // Bottom Right
        glVertex3f(0.05,0.05, 0.1f);              // Bottom Left
        glEnd();                            // Done Drawing The Quad
        
    
    
    
    for(int i=0;i<nbins;++i){
        
        xtl=1+i*binwidth; ytl=(VECTOR((histdist.bins))[i])*8+1;
        xtr=1+(1+i)*binwidth; ytr=ytl;
        xbl=xtl; ybl=1;
        xbr=xtr; ybr=ybl;
        
        //BORDER
       if(ytl>1){
        glColor3f(0,0,0);
        glBegin(GL_QUADS);                      // Draw A Quad
        glVertex3f(xtl, ytl, 0.2f);              // Top Left
        glVertex3f(xtr, ytr, 0.2f);              // Top Right
        glVertex3f(xbr,ybr, 0.2f);              // Bottom Right
        glVertex3f(xbl,ybl, 0.2f);              // Bottom Left
        glEnd();                            // Done Drawing The Quad
        
        //inner
        glColor3f(1,0,0);
        glBegin(GL_QUADS);
        glVertex3f(xtl+0.05, ytl-0.05, 0.3f);
        glVertex3f(xtr-0.05, ytr-0.05, 0.3f);
        glVertex3f(xbr-0.05, ybr+0.05, 0.3f);
        glVertex3f(xbl+0.05, ybl+0.05, 0.3f);
        glEnd();
        }
    }
        
        
        //draw mean
        
        glColor3f(0,1,0);
        glBegin(GL_LINES);
        glVertex3d(1+(9./4.)*(histdist.mean),10.,0.4f);
        glVertex3d(1+(9./4.)*(histdist.mean),0.6,0.4f);
        glEnd();
        gl_font(FL_TIMES,12);
        gl_draw("Mean",(float)(0+(9./4.)*(histdist.mean)),(float)10.);
    
        
    
    glPopMatrix();
    
    
    //chart
    glColor3f(0,0,0);
    glPushMatrix();
    glTranslated(0,0,0);
   
    //X AXIS
    glBegin(GL_LINES);
    glVertex3d(0.5,1.,0.5f);
    glVertex3d(10.5,1.,0.5f);
    glEnd();
        
         gl_font(FL_TIMES,12);
        
        glBegin(GL_LINES);
        glVertex3d(1+(9./4.),10.,0.6f);
        glVertex3d(1+(9./4.),0.6,0.6f);
        glEnd();
        gl_draw("1",(float)(0.8+(9./4.)),(float)0.3);
        
        glBegin(GL_LINES);
        glVertex3d(1+(9./4.)*2,1.,0.6f);
        glVertex3d(1+(9./4.)*2,0.8,0.6f);
        glEnd();
        gl_draw("2",(float)(0.8+(9./4.)*2),(float)0.3);
        
        glBegin(GL_LINES);
        glVertex3d(1+(9./4.)*3,1.,0.6f);
        glVertex3d(1+(9./4.)*3,0.8,0.6f);
        glEnd();
        gl_draw("3",(float)(0.8+(9./4.)*3),(float)0.3);
        
        glBegin(GL_LINES);
        glVertex3d(1+(9./4.)*4,1.,0.6f);
        glVertex3d(1+(9./4.)*4,0.8,0.6f);
        glEnd();
        gl_draw("4",(float)(0.8+(9./4.)*4),(float)0.3);
        

        
        
    //Y AXIS
    glBegin(GL_LINES);
    glVertex3d(1,1.,0.7f);
    glVertex3d(1,10.,0.7f);
    glEnd();

        
        
    glPopMatrix();
        
        if(histdist.mean!=atof(meanbuff->text())){
            /*if(runningcontrol==1){ dm=fabs(atof(meanbuff->text())-histdist.mean);
                sprintf(textstring,"%f", dm);
                dmeanbuff->text(textstring);}*/
            hism=histdist.mean;
        sprintf(textstring,"%f", hism);
        meanbuff->text(textstring);
            
        }
        
        if(histdist.variance!=atof(varbuff->text())){
           /* if(runningcontrol==1){
                dv=fabs(atof(varbuff->text())-histdist.variance);
        sprintf(textstring,"%f", dv);
                dvarbuff->text(textstring);}*/
            hisv=histdist.variance;
        sprintf(textstring,"%f", hisv);
        varbuff->text(textstring);}
        
        
    
    }
    
    else {
    
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        
        if(*(meanbuff->text())!='n'){
        sprintf(textstring,"not computable");
        meanbuff->text(textstring);
        }
        
        
        if(*(varbuff->text())!='n'){
        sprintf(textstring,"not computable");
        varbuff->text(textstring);
        }
        
        /*if(*(dmeanbuff->text())!='n'){
            sprintf(textstring,"not computable");
            dmeanbuff->text(textstring);
        }
        
        if(*(dvarbuff->text())!='n'){
            sprintf(textstring,"not computable");
            dvarbuff->text(textstring);
        }*/
        
        
        
        
        glClearColor(0., 0., 0., 1);
    
    }
    

    
    
}









/**********************************************                         ****************************************
 **********************************************                         ****************************************
 **********************************************  ACTIVATION CHART DRAW  ****************************************
 **********************************************                         ****************************************
 **********************************************                         ****************************************
 */




void draw_tactinit(void){
    
    
}





void draw_tactscene(void){
    
    
    char textstring[100];
    char firstchar;
    
    double xtl, xtr, xbl, xbr, ytl, ytr, ybl, ybr;
    //maxmum is 8 points in height.
    
    double span=20;
    double hone=8;
    double dist;
    
    int tickstodraw;
    
    tickstodraw=ticks;
    if(ticks > tactdisp){tickstodraw=tactdisp;}
    
    
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glClearColor(1., 1., 1., 1);
    
    
    //chart
    glColor3f(0,0,0);
    glPushMatrix();
    glTranslated(0,0,0);
    
    
    glColor3f(0,0,0);
    glBegin(GL_QUADS);                      // Draw A Quad
    glVertex3f(0, 22, 0.0f);              // Top Left
    glVertex3f(22, 22, 0.0f);              // Top Right
    glVertex3f(22,0, 0.0f);              // Bottom Right
    glVertex3f(0,0, 0.0f);              // Bottom Left
    glEnd();                            // Done Drawing The Quad
    
    
    glColor3f(1.,1.,1.);
    glBegin(GL_QUADS);                      // Draw A Quad
    glVertex3f(0.05, 10.95, 0.05f);              // Top Left
    glVertex3f(21.95, 10.95, 0.05f);              // Top Right
    glVertex3f(21.95,0.05, 0.05f);              // Bottom Right
    glVertex3f(0.05,0.05, 0.05f);              // Bottom Left
    glEnd();                            // Done Drawing The Quad
    
    glColor3f(0.,0.,0.);
    
    
    //X AXIS
    glBegin(GL_LINES);
    glVertex3d(1,1.,0.1f);
    glVertex3d(span+2*0.5,1.,0.1f);
    glEnd();
    
    gl_font(FL_TIMES,10);
    
    // print last x ticks
    glBegin(GL_LINES);
    glVertex3d(1+span,1.,0.4f);
    glVertex3d(1+span,0.8,0.4f);
    glEnd();
    sprintf(textstring,"%i",ticks);
    gl_draw(textstring,(float)(0.8+span),(float)0.3);
    
    //print half x ticks
    glBegin(GL_LINES);
    glVertex3d(1+span/2,1.,0.4f);
    glVertex3d(1+span/2,0.8,0.4f);
    glEnd();
    sprintf(textstring,"%.1f",(float)(ticks-tickstodraw+ticks)/2);
    gl_draw(textstring,(float)(0.8+span/2),(float)0.3);
    
    // print 0
    sprintf(textstring,"%i",(ticks-tickstodraw));
    gl_draw(textstring,(float)(0.8),(float)0.3);
    
    
    
    //Y AXIS
    glBegin(GL_LINES);
    glVertex3d(1,1.,0.5f);
    glVertex3d(1,10.,0.5f);
    glEnd();
    gl_draw("0",(float)0.5,(float)1.);

    
    glBegin(GL_LINES);
    glVertex3d(1,1.+hone, 0.1f);
    glVertex3d(1+span,1.+hone, 0.1f);
    glEnd();
    gl_draw("1",(float)0.5,(float)1.+hone);
    
    glBegin(GL_LINES);
    glVertex3d(1,1.+hone/2, 0.1f);
    glVertex3d(1+span,1.+hone/2, 0.1f);
    glEnd();
    gl_draw("0.5",(float)0.2,(float)1.+hone/2);
    
    
    
    //DRAW LINE
    
    
    
    if(lastick!=ticks){
        
        dist=(float)span/(float)tickstodraw;
        
        if(ticks <= tactdisp && ticks > 0){
            
            //x
            for (int i=0; i<=tickstodraw; ++i) {
                tactcoord[3*i]=1+i*dist;
            }
            tactcoord[3*ticks]=1+tickstodraw*dist;
            //y
            tactcoord[(3*ticks)+1]=1+VECTOR(tactvect)[ticks-1]*hone;
            //z
            tactcoord[(3*ticks)+2]=1;
        }
        
        else{
            // SHIFT y of 1 coord (i.e. of 3)
            for (int i=0; i<=(tickstodraw-1); ++i) {
                tactcoord[3*i+1]=tactcoord[3*(i+1)+1];
            }
            
            tactcoord[3*(tickstodraw)]=1+tickstodraw*dist;
            tactcoord[(3*(tickstodraw))+1]=1+VECTOR(tactvect)[ticks-1]*hone;
            
        }
        
        
        
        lastick=ticks;
        
        /*
         printf("\n\n dist = %f   ", dist);
         for (int i=0; i<ticks; ++i) {
         printf("%i %f %f (%f) %f \n", i, tactcoord[3*i],tactcoord[3*i+1], VECTOR(tactvect)[i] ,tactcoord[3*i+2]);
         }
         
         printf("\n\n");
         printf("\n\n");
         */
        
    }
    
    glColor3f(1,0,0);
    glEnableClientState(GL_VERTEX_ARRAY);
    glVertexPointer(3, GL_FLOAT, 0, tactcoord);
    glDrawArrays(GL_LINE_STRIP, 1, tickstodraw);
    glDisableClientState(GL_VERTEX_ARRAY);
    
    
    
    
    
    
    glPopMatrix();
    
    
    
    
}





/**********************************************                         ****************************************
 **********************************************                         ****************************************
 **********************************************  ACTIVATION  HIST DRAW  ****************************************
 **********************************************                         ****************************************
 **********************************************                         ****************************************
 */




void draw_htactinit(void){
    
    
}





void draw_htactscene(void){
    
    
    int nbins=histtact.nbins;
    float binwidth=9./nbins;
    
    char textstring[100];
    char firstchar;
    
    float xtl, xtr, xbl, xbr, ytl, ytr, ybl, ybr;
    //maxmum is 8 points in height.
    
    
    if(drawhisto==1){
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glClearColor(1., 1., 1., 1);
        
        glPushMatrix();
        glTranslated(0,0,0);
        
        
        
        
        
        glColor3f(0,0,0);
        glBegin(GL_QUADS);                      // Draw A Quad
        glVertex3f(0, 11, 0.0f);              // Top Left
        glVertex3f(11, 11, 0.0f);              // Top Right
        glVertex3f(11,0, 0.0f);              // Bottom Right
        glVertex3f(0,0, 0.0f);              // Bottom Left
        glEnd();                            // Done Drawing The Quad
        
        
        glColor3f(1.,1.,1.);
        glBegin(GL_QUADS);                      // Draw A Quad
        glVertex3f(0.05, 10.95, 0.1f);              // Top Left
        glVertex3f(10.95, 10.95, 0.1f);              // Top Right
        glVertex3f(10.95,0.05, 0.1f);              // Bottom Right
        glVertex3f(0.05,0.05, 0.1f);              // Bottom Left
        glEnd();                            // Done Drawing The Quad
      
        
        
        
        for(int i=0;i<nbins;++i){
            
            xtl=1+i*binwidth; ytl=(VECTOR((histtact.bins))[i])*8+1;
            xtr=1+(1+i)*binwidth; ytr=ytl;
            xbl=xtl; ybl=1;
            xbr=xtr; ybr=ybl;
            
            //BORDER
            if(ytl>1){
                glColor3f(0,0,0);
                glBegin(GL_QUADS);                      // Draw A Quad
                glVertex3f(xtl, ytl, 0.2f);              // Top Left
                glVertex3f(xtr, ytr, 0.2f);              // Top Right
                glVertex3f(xbr,ybr, 0.2f);              // Bottom Right
                glVertex3f(xbl,ybl, 0.2f);              // Bottom Left
                glEnd();                            // Done Drawing The Quad
                
                //inner
                glColor3f(1,0,0);
                glBegin(GL_QUADS);
                glVertex3f(xtl+0.05, ytl-0.05, 0.3f);
                glVertex3f(xtr-0.05, ytr-0.05, 0.3f);
                glVertex3f(xbr-0.05, ybr+0.05, 0.3f);
                glVertex3f(xbl+0.05, ybl+0.05, 0.3f);
                glEnd();
            }
        }
        
        
        //draw mean
        
        glColor3f(0,1,0);
        glBegin(GL_LINES);
        glVertex3d(1+9*(histtact.mean),10.,0.4f);
        glVertex3d(1+9*(histtact.mean),0.6,0.4f);
        glEnd();
        gl_font(FL_TIMES,12);
        gl_draw("Mean",(float)(0+9*(histtact.mean)),(float)10.);
        
        
        
        glPopMatrix();
        
        
        //chart
        glColor3f(0,0,0);
        glPushMatrix();
        glTranslated(0,0,0);
        
        //X AXIS
        glBegin(GL_LINES);
        glVertex3d(0.5,1.,0.5f);
        glVertex3d(10.5,1.,0.5f);
        glEnd();
        
        gl_font(FL_TIMES,12);
        
        glBegin(GL_LINES);
        glVertex3d(1+(9./4.),1.,0.6f);
        glVertex3d(1+(9./4.),0.8,0.6f);
        glEnd();
        gl_draw("0.25",(float)(0.8+(9./4.)),(float)0.3);
        
        glBegin(GL_LINES);
        glVertex3d(1+(9./4.)*2,1.,0.6f);
        glVertex3d(1+(9./4.)*2,0.8,0.6f);
        glEnd();
        gl_draw("0.5",(float)(0.8+(9./4.)*2),(float)0.3);
        
        glBegin(GL_LINES);
        glVertex3d(1+(9./4.)*3,1.,0.6f);
        glVertex3d(1+(9./4.)*3,0.8,0.6f);
        glEnd();
        gl_draw("0.75",(float)(0.8+(9./4.)*3),(float)0.3);
        
        glBegin(GL_LINES);
        glVertex3d(1+(9./4.)*4,1.,0.6f);
        glVertex3d(1+(9./4.)*4,0.6,0.6f);
        glEnd();
        gl_draw("1",(float)(0.8+(9./4.)*4),(float)0.3);
        
        
        
        
        //Y AXIS
        glBegin(GL_LINES);
        glVertex3d(1,1.,0.7f);
        glVertex3d(1,10.,0.7f);
        glEnd();
        
        
        
        glPopMatrix();
        
      
        
        
        
    }
    
    else {
        
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        
        
        
        
        glClearColor(0., 0., 0., 1);
        
    }
    
    
    
    
    
}












/**********************************************                         ****************************************
 **********************************************                         ****************************************
 **********************************************  FLUCTATIONS HIST DRAW  ****************************************
 **********************************************                         ****************************************
 **********************************************                         ****************************************
 */




void draw_hfluctinit(void){
    
    
}





void draw_hfluctscene(void){
    
    
    int nbins=histfluct.nbins;
    float binwidth=9./nbins;
    
    char textstring[100];
    char firstchar;

    
    float xtl, xtr, xbl, xbr, ytl, ytr, ybl, ybr;
    //maxmum is 8 points in height.
    
    
    if(drawhisto==1)
    
    {
        

        
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glClearColor(1., 1., 1., 1);
        
        glPushMatrix();
        glTranslated(0,0,0);
        
        
        
        
        
        glColor3f(0,0,0);
        glBegin(GL_QUADS);                      // Draw A Quad
        glVertex3f(0, 11, 0.0f);              // Top Left
        glVertex3f(11, 11, 0.0f);              // Top Right
        glVertex3f(11,0, 0.0f);              // Bottom Right
        glVertex3f(0,0, 0.0f);              // Bottom Left
        glEnd();                            // Done Drawing The Quad
        
        
        glColor3f(1.,1.,1.);
        glBegin(GL_QUADS);                      // Draw A Quad
        glVertex3f(0.05, 10.95, 0.1f);              // Top Left
        glVertex3f(10.95, 10.95, 0.1f);              // Top Right
        glVertex3f(10.95,0.05, 0.1f);              // Bottom Right
        glVertex3f(0.05,0.05, 0.1f);              // Bottom Left
        glEnd();                            // Done Drawing The Quad
        
        
        
        
        for(int i=0;i<nbins;++i){
            
            xtl=1+(i*binwidth); ytl=(VECTOR((histfluct.bins))[i])*8+1;
            xtr=1+(1+i)*binwidth; ytr=ytl;
            xbl=xtl; ybl=1;
            xbr=xtr; ybr=ybl;
            
            //BORDER
            if(ytl>1){
                glColor3f(0,0,0);
                glBegin(GL_QUADS);                      // Draw A Quad
                glVertex3f(xtl, ytl, 0.2f);              // Top Left
                glVertex3f(xtr, ytr, 0.2f);              // Top Right
                glVertex3f(xbr,ybr, 0.2f);              // Bottom Right
                glVertex3f(xbl,ybl, 0.2f);              // Bottom Left
                glEnd();                            // Done Drawing The Quad
                
                //inner
                glColor3f(1,0,0);
                glBegin(GL_QUADS);
                glVertex3f(xtl+0.05, ytl-0.05, 0.3f);
                glVertex3f(xtr-0.05, ytr-0.05, 0.3f);
                glVertex3f(xbr-0.05, ybr+0.05, 0.3f);
                glVertex3f(xbl+0.05, ybl+0.05, 0.3f);
                glEnd();
            }
        }
        
        
        //draw mean
        
        
        glColor3f(0,1,0);
        glBegin(GL_LINES);
        glVertex3d(1+9*(histfluct.mean),10.,0.4f);
        glVertex3d(1+9*(histfluct.mean),0.6,0.4f);
        glEnd();
        gl_font(FL_TIMES,12);
        gl_draw("Mean",(float)(0+9*(histfluct.mean)),(float)10.);
        
        
        
        glPopMatrix();
        
        
        //chart
        glColor3f(0,0,0);
        glPushMatrix();
        glTranslated(0,0,0);
        
        //X AXIS
        glBegin(GL_LINES);
        glVertex3d(0.5,1.,0.5f);
        glVertex3d(10.5,1.,0.5f);
        glEnd();
        
        gl_font(FL_TIMES,12);
        
        glBegin(GL_LINES);
        glVertex3d(1+(9./4.),1.,0.6f);
        glVertex3d(1+(9./4.),0.8,0.6f);
        glEnd();
        sprintf(textstring,"%.3f",histfluct.max / 4.);
        gl_draw(textstring,(float)(0.4+(9./4.)),(float)0.3);
        
        glBegin(GL_LINES);
        glVertex3d(1+(9./4.)*2,1.,0.6f);
        glVertex3d(1+(9./4.)*2,0.8,0.6f);
        glEnd();
        sprintf(textstring,"%.3f",histfluct.max /2.);
        gl_draw(textstring,(float)(0.4+(9./4.)*2),(float)0.3);
        
        glBegin(GL_LINES);
        glVertex3d(1+(9./4.)*3,1.,0.6f);
        glVertex3d(1+(9./4.)*3,0.8,0.6f);
        glEnd();
        sprintf(textstring,"%.3f",histfluct.max*(3./4.));
        gl_draw(textstring,(float)(0.4+(9./4.)*3),(float)0.3);
        
        glBegin(GL_LINES);
        glVertex3d(1+(9./4.)*4,1.,0.6f);
        glVertex3d(1+(9./4.)*4,0.6,0.6f);
        glEnd();
        sprintf(textstring,"%.3f",histfluct.max);
        gl_draw(textstring,(float)(0.4+(9./4.)*4),(float)0.3);
        
        
        
        
        //Y AXIS
        glBegin(GL_LINES);
        glVertex3d(1,1.,0.7f);
        glVertex3d(1,10.,0.7f);
        glEnd();
        
        
        
        glPopMatrix();
        
        
        
        
        
    }
    
    else {
        
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        
        
        
        
        glClearColor(0., 0., 0., 1);
        
    }
    
    
    
}
















