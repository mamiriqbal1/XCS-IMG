#include <fstream>

long seeds[30] = {114665,134296,176806,247157,262025,311756,336922,337104,344465,362234,379332,425485,470006,490758,538402,583115,584068,614795,678991,710953,715062,715911,784943,787483,790558,810292,824057,846845,968381,972820};
FILE *filePerformance=NULL, *testPerformance=NULL, *fileClassifierPopulation=NULL, *cfWritingFilePointer=NULL, *cfReadingFilePointer=NULL;

DataSource *input = NULL,*trainingData = NULL, *testingData = NULL;
void loadDataFile(DataSource inputArray[]);

void resetState(float state[], int inputFileint[][condLength+1]);
double executeAction(int action, int stateAction, bool &wasCorrect);

void header();
void Exit(FILE *fp);

void writePerformance(ClassifierMap &pop, int *performance, double *sysError, int exploreProbC);

void startXCS();
void doOneSingleStepExperiment(ClassifierMap &pop);

void doOneSingleStepTest(ClassifierMap &pop);
