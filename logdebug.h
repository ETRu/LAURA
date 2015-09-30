
#pragma once

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <ctype.h>
#include <locale.h>
#include <time.h>


typedef unsigned char OCTECT;
typedef unsigned int FOURBYTES;


void asciiDateTime(FOURBYTES tstamp, char *string, int maxlen, char special);
void logDebug(char *format, ...);
