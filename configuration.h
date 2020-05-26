#ifndef CONFIGURATION_H
#define CONFIGURATION_H

#include <string>
#include <forward_list>
#include "xcsMacros.h"

extern bool use_kb;
extern std::string kb_file;
extern std::string output_path;
extern int numActions; // = 2; //0 or 1
// max_actions have to be set such that it is always more than the actual number of actions in an experiment
const int max_actions = 10;
extern int trainNumInstances;// = 11982;
extern int testNumInstances;// = 1984;
extern int totalNumInstances;// = trainNumInstances + testNumInstances;//4; //for review analysis

const int condLength = 784; //512; //256; //4626; //for review analysis 300, 500, 1000, 1500, 2000, 2257
const int image_width = 28;
const int image_height = 28;
extern int maxPopSize; //=  1000; //1 * totalNumInstances; //Specifies the maximal number of micro-classifiers in the population. [ (0.5, 1, 2, 5, 10/20, 50)*1000 for 6-, 11-, 20-, 37-, 70-, 135-bits MUX respectively]
extern int maxProblems;// = trainNumInstances; //50 * totalNumInstances; //1*100*1000; //training set = [ (1, 1.5, 2, 2.5, 3,)*100*1000 for 6-, 11-, 20-, 37-, 70-, 135-bits MUX respectively]
const double maxPayoff = 1000;
const int clfrCondLength = 1; // 784/8; // 64; //32; //300;//condLength/4; // condLength/2 and condLength/4 for 70mux and 135mux respectively.
const int run = 1;
//const char inputFile[] = "features1.txt";
extern std::string inputTrainingFile; //[] = "../data/mnist/3_8_train_mnist.txt";
extern std::string inputTestFile; //[] = "../data/mnist/3_8_test_mnist.txt";
//const char testFile[] = ""; extern bool Testing;
//const int currentProblemLevel = 1; //must be set coz previous level and current level file names are adjusted

const int precisionDigits = 2;
extern int testFrequency;// = trainNumInstances; // 1034;

const char outputFileName[] = "output_training.txt";
const char featureFileName[] = "feature_codefragments.txt";
const char ruleFileName[] = "rule_with_codefragements.txt";
const char resultFile[] = "result_testing.txt";

const int cfMaxDepth = 0;
const int cfMinDepth = 0;
const int cfMaxLength = 2;// pow(2,adfMaxDepth+1); //allow for endstop OPNOP
const int cfMaxArity = 2;
const int cfMaxStack = (cfMaxArity-1)*(cfMaxDepth-1)+2;
const int numLeaf = 1;
const int filter_size = 4;  // 3 for 3x3, 5 for 5x5 etc.

typedef int opType;
const int opSize = sizeof(opType);

const opType OPNOP = -100; //to be used as ending symbol
const opType OPAND = -101;
const opType OPOR = -102;
const opType OPNAND = -103;
const opType OPNOR = -104;
const opType OPNOT = -105;
//const opType OPUNITY = -106;

const int totalFunctions = 1;
const opType functionCodes[] = {OPAND};


struct Filter{
    int id=-1; // uninitialized value
    float lower_bounds[filter_size*filter_size];
    float upper_bounds[filter_size*filter_size];
};

struct FilterList{
    std::forward_list<Filter> filters;
    int gid = 0;
    int max_size_limit = N_filter_ol;
    float fitness = 1.0; // initial fitness of the filter
};

struct Leaf
{
    opType featureNumber;
    float lowerBound;
    float upperBound;
};

struct CodeFragment
{
    opType codeFragment[cfMaxLength];
    //Leaf leaf[numLeaf];
    int num_filters; // equal to number of leaves to be determined at run time
    Filter filter[numLeaf];
    int filter_id[numLeaf];
    int cfID;
};

struct Classifier
{
    int id;
    CodeFragment condition[clfrCondLength];
    int action;
    double prediction;
    double predictionError;
    double accuracy;
    double fitness;
    int experience;
    int numerosity;
    double actionSetSize;
    int timeStamp;
    int specificness; //number of specific CFs
};

struct ClassifierSet
{
    Classifier *classifier;
    ClassifierSet *next;
};

struct DataSource{
    float *state;
    int action;

    ~DataSource(){
        delete []state;
    }
};


#endif //CONFIGURATION_H
