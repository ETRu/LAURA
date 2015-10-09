#pragma once

#if defined (WIN32)
#include <FL\Gl.h>
#else
#include <FL/gl.h>
#endif

#if defined (WIN32)
#include "igraph/igraph.h"
#else
#include "igraph.h"
#endif

#include <math.h>
#include "dataanalysis.h"

#include <FL/Fl_Round_Button.H>
#include <FL/Fl_Text_Display.H>


double RedHue(double color);
double GreenHue(double color);
double BlueHue(double color);

void draw_scene(void);
void draw_init(void);

void draw_datascene(void);
void draw_datainit(void);

void draw_tactscene(void);
void draw_tactinit(void);

void draw_htactscene(void);
void draw_htactinit(void);

void draw_hfluctscene(void);
void draw_hfluctinit(void);

#define MAXTACTREG 10000

