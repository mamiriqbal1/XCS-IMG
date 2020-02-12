#include <sys/time.h>
#include <fstream>

void startTimer();
void outtime(const timeval t, FILE *filePointer);
void timeval_add(const timeval *add, timeval *sum);
int timeval_subtract (timeval *result, const timeval *x, const timeval *y);
void elapsed0(const timeval* start, timeval* diff);
void elapsed1(timeval* diff);
void elapsed(FILE *filePointer);
void elapsedTime();
