
#pragma once

#if defined (WIN32)
#include <FL\Gl.h>
#else
#include <FL/gl.h>
#endif

#if defined (WIN32)
#include <FL\Glu.h>
#else
#include <FL/glu.h>
#endif

#include "draw.h"

#include <FL/Fl_Gl_Window.H>


class Frame : public Fl_Gl_Window {
  void draw();
public:
  Frame(int x,int y,int w,int h,const char *l=0) : Fl_Gl_Window(x,y,w,h,l) {}
};

class DataFrame : public Fl_Gl_Window {
    void draw();
public:
    DataFrame(int x,int y,int w,int h,const char *l=0) : Fl_Gl_Window(x,y,w,h,l) {}
};


class TactFrame : public Fl_Gl_Window {
    void draw();
public:
    TactFrame(int x,int y,int w,int h,const char *l=0) : Fl_Gl_Window(x,y,w,h,l) {}
};


class HTactFrame : public Fl_Gl_Window {
    void draw();
public:
    HTactFrame(int x,int y,int w,int h,const char *l=0) : Fl_Gl_Window(x,y,w,h,l) {}
};


class HFluctFrame : public Fl_Gl_Window {
    void draw();
public:
    HFluctFrame(int x,int y,int w,int h,const char *l=0) : Fl_Gl_Window(x,y,w,h,l) {}
};




