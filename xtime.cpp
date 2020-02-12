#include <string.h>
#include <fstream>
#include "xtime.h"

timeval clockstart, cpustart, cputotal;

void startTimer()
{
    gettimeofday(&clockstart,NULL);
    gettimeofday(&cpustart,NULL);
    memset(&cputotal,0,sizeof(timeval));
}

void outtime(const struct timeval t, FILE *filePointer)
{
    char buf[1000];
    sprintf(buf,"%ld.%06ld ",t.tv_sec,t.tv_usec);
    printf("%s",buf);
    fwrite(buf,strlen(buf),1,filePointer);
    fwrite("\n",strlen("\n"),1,filePointer);
}
void timeval_add(const struct timeval *add, struct timeval *sum)
{
    sum->tv_usec += add->tv_usec;
    int nsec = sum->tv_usec / 1000000;
    sum->tv_usec %= 1000000;
    sum->tv_sec  += add->tv_sec+nsec;
}
int timeval_subtract (struct timeval *result, const struct timeval *x, const struct timeval *y)
{
    int y_tv_sec  = y->tv_sec;
    int y_tv_usec = y->tv_usec;

    /* Perform the carry for the later subtraction by updating copy of y. */
    if (x->tv_usec < y_tv_usec)
    {
        int nsec = (y_tv_usec - x->tv_usec) / 1000000 + 1;
        y_tv_usec -= 1000000 * nsec;
        y_tv_sec  += nsec;
    }
    if (x->tv_usec - y_tv_usec > 1000000)
    {
        int nsec = (y_tv_usec - x->tv_usec) / 1000000;
        y_tv_usec += 1000000 * nsec;
        y_tv_sec  -= nsec;
    }
    /* Compute the time remaining to wait.
       tv_usec is certainly positive. */
    result->tv_sec  = x->tv_sec  - y_tv_sec;
    result->tv_usec = x->tv_usec - y_tv_usec;

    /* Return 1 if result is negative. */
    return x->tv_sec < y_tv_sec;
}
void elapsed0(const timeval* start, timeval* diff)
{
    timeval now;
    gettimeofday(&now,NULL);
    timeval_subtract (diff, &now, start);
}
void elapsed1(timeval* diff)
{
    timeval now;
    gettimeofday(&now,NULL);
    timeval_subtract (diff, &now, &clockstart);
}
void elapsed(FILE *filePointer)
{
    timeval diff;
    elapsed1(&diff);
    outtime(diff,filePointer);
}
void elapsedTime()
{
    timeval diff;
    elapsed0(&cpustart,&diff);
    timeval_add(&diff,&cputotal);
}
