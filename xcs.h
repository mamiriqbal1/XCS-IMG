#include <fstream>

long seeds[30] = {114665,134296,176806,247157,262025,311756,336922,337104,344465,362234,379332,425485,470006,490758,538402,583115,584068,614795,678991,710953,715062,715911,784943,787483,790558,810292,824057,846845,968381,972820};

DataSource *trainingData = NULL, *testingData = NULL;
double executeAction(int action, int stateAction, bool &wasCorrect);


void writePerformance(double performance, double sysError, int problem_count, double action_set_size,
                      std::ofstream &output_training_file);

void startXCS();
void doOneSingleStepExperiment();

void
doOneSingleStepTest(int training_problem_count, std::ofstream &output_test_file, double training_performance,
                    double training_error, std::string path_postfix);
