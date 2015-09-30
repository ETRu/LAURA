
#include "Frame.h"

extern int newgraph;

int framexmin=0, framexmax=11, frameymin=0, frameymax=11;


//------------------------------------------------------------------------------
void Frame::draw() {
  if (!valid()) {
    glClearColor(0.0, 0.0, 0.0, 1);                        // Turn the background color black
    glViewport(0,0,w(),h());                               // Make our viewport the whole window
    glMatrixMode(GL_PROJECTION);                           // Select The Projection Matrix
    glLoadIdentity();                                      // Reset The Projection Matrix
    gluOrtho2D(framexmin,framexmax,frameymin,frameymax);                             // (xmin,xmax,ymin,ymax)
    //gluPerspective(45.0f,w()/h(), 1 ,150.0);
    glMatrixMode(GL_MODELVIEW);                            // Select The Modelview Matrix
    glLoadIdentity();                                      // Reset The Modelview Matrix
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);    // Clear The Screen And The Depth Buffer
    glLoadIdentity();                                      // Reset The View
    //gluLookAt( 0, 0, 10,     0, 0, 0,     0, 1, 0);         // Position - View  - Up Vector
    glEnable(GL_DEPTH_TEST);
    
	  
    valid(1);
  }
	
	if(newgraph==1){draw_init(); newgraph=0;}

  draw_scene();
}
//------------------------------------------------------------------------------



//------------------------------------------------------------------------------
void DataFrame::draw() {
    if (!valid()) {
        glClearColor(1., 1.0, 1.0, 1);                        // Turn the background color white
        glViewport(0,0,w(),h());                               // Make our viewport the whole window
        glMatrixMode(GL_PROJECTION);                           // Select The Projection Matrix
        glLoadIdentity();                                      // Reset The Projection Matrix
        gluOrtho2D(framexmin,framexmax,frameymin,frameymax);                             // (xmin,xmax,ymin,ymax)
        //gluPerspective(45.0f,w()/h(), 1 ,150.0);
        glMatrixMode(GL_MODELVIEW);                            // Select The Modelview Matrix
        glLoadIdentity();                                      // Reset The Modelview Matrix
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);    // Clear The Screen And The Depth Buffer
        glLoadIdentity();                                      // Reset The View
        //gluLookAt( 0, 0, 10,     0, 0, 0,     0, 1, 0);         // Position - View  - Up Vector
        glEnable(GL_DEPTH_TEST);
        
        
        valid(1);
    }
    
    draw_datainit();
    draw_datascene();
}
//------------------------------------------------------------------------------
