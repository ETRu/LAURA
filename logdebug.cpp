
#include "logdebug.h"


typedef unsigned char OCTECT;
typedef unsigned int FOURBYTES;

time_t now;
time_t basetime=0;
int counter4debug=0;

void asciiDateTime(FOURBYTES tstamp, char *string, int maxlen, char special)
/*  returns a string "yyyy-mm-dd.hh_mm_ss"
 ( '_' is specialChar )
 corresponding to tstamp
 (is 19 chars long up to year 10000)
 */
{  struct tm timeStruct;
    int todayYear, year, month, day, hh, mm, ss;
    time_t years68 = 1<<32, adjtstamp;
    
    if ( maxlen < 19 + 1 )
    {printf ("asciiDateTime - called with too short a maxlen: %d\n", maxlen); abort();}
    
    /* currente year, e.g. 2075 */
    
    timeStruct = *gmtime(&now);
    todayYear = timeStruct.tm_year + 1900;
    
    /* initialize structure to adjusted_tstamp = tstamp + basetime */
    adjtstamp = tstamp + basetime;
    timeStruct = *gmtime(&adjtstamp);
    /* until time_t is a signed long integer,
     this will work only up to 2037 !!! */
    
    year = timeStruct.tm_year + 1900;
    
    
    month    = timeStruct.tm_mon + 1;
    day      = timeStruct.tm_mday;
    hh       = timeStruct.tm_hour;
    mm       = timeStruct.tm_min;
    ss       = timeStruct.tm_sec;
    
    sprintf(string, "%d-%02d-%02d.%02d%c%02d%c%02d", year, month, day, hh, special, mm, special, ss);
    return;
}



void logDebug(char *format, ...)
{
    
    counter4debug++;
    
    char asciiDate[19 + 1];
    va_list arguments;
    va_start(arguments, format);
    time( &now );
    asciiDateTime(now-basetime, asciiDate, sizeof(asciiDate), ':');
    /* fprintf(stdout, "%s - ", asciiDate); */
    
    {
        fprintf(stdout, "%i - %s - ", counter4debug, asciiDate);
        
    }
    /* vfprintf(stdout, format, arguments); */
    
    
    {
        vfprintf(stdout, format, arguments);
        va_end(arguments);
        fflush (stdout);
    }
    /* fflush (stdout); */
    
    return;
}


