#ifndef CONFIGURATION_H
#define CONFIGURATION_H

#include <string>
#include <forward_list>
#include <unordered_map>
#include "xcsMacros.h"
#include <vector>

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
const int run = 1;
//const char inputFile[] = "features1.txt";
extern std::string inputTrainingFile; //[] = "../data/mnist/3_8_train_mnist.txt";
extern std::string inputTestFile; //[] = "../data/mnist/3_8_test_mnist.txt";
//const char testFile[] = ""; extern bool Testing;
//const int currentProblemLevel = 1; //must be set coz previous level and current level file names are adjusted

const int precisionDigits = 2;
extern int testFrequency;// = trainNumInstances; // 1034;
extern int validation_frequency;
const int filter_list_management_frequency = 500;

const char output_training_file_name[] = "training_performance.txt";
const char output_test_file_name[] = "test_performance.txt";
const char output_classifier_file_name[] = "classifier.txt";
const char output_code_fragment_file_name[] = "code_fragment.txt";
const char output_promising_code_fragment_file_name[] = "promising_code_fragment.txt";
const char output_filter_file_name[] = "filter.txt";
const char output_promising_filter_file_name[] = "promising_filter.txt";
const char output_stats_file_name[] = "stats.txt";

extern int clfrCondMaxLength;//2 // 784/8; // 64; //32; //300;//condLength/4; // condLength/2 and condLength/4 for 70mux and 135mux respectively.
extern int cfMaxDepth;// = 2;
const int cfMinDepth = 0;
extern int cfMaxLength;// = 8;// 2^(cdfMaxDepth+1); //allow for endstop OPNOP
const int cfMaxArity = 2;
extern int cfMaxStack;// = (cfMaxArity-1)*(cfMaxDepth-1)+2;
extern int cfMaxLeaf;// = 4; // 2^cfMaxDepth

typedef int opType;
const int opSize = sizeof(opType);
/*
 * We have pairs of operators (positive and negative). They will cover all combinations without the use of NOT
 * There is only one case which require NOT and that is single filter with NOT.
 * So NOT is not used to build GP tree except when there is only only filter.
 * This can be used when growing a tree of depth zero. As soon as the tree depth increases, the NOT will be removed.
 * This will allow easy negation of the code fragment
 */
const opType OPNOP = -100; //to be used as ending symbol
const opType OPAND = -101;
const opType OPOR = -102;
const opType OPNAND = -103;
const opType OPNOR = -104;
const opType OPNOT = -105;
//const opType OPUNITY = -106;

const int totalFunctions = 4;
// OPNOT must be the last operator
const opType functionCodes[] = {OPAND, OPOR, OPNAND, OPNOR};
const int num_negative_binary_operators = 2;
const opType negative_binary_operators[] = {OPNAND, OPNOR};


const int num_filter_sizes = 4;
const int filter_sizes[] = {3, 5, 7, 9}; // filter sizes to be used
const int max_filter_size = filter_sizes[num_filter_sizes-1];  // the last one should be maximum
const bool allow_dilated_filters = true;

struct Filter{
    int id=-1; // uninitialized value
    int numerosity = 1; // initial numerosity when a filter is created
    // Fitness of a filter is the appearance of the filter in "promising classifiers"
    // a promising classifier is one whose error < 10 and experience > 10
    int fitness = 0; // Fitness of a filter is the appearance of the filter in "promising classifiers"
    int filter_size = -1;
    bool is_dilated = false; // what is the type of filter normal or dilated
    float lower_bounds[max_filter_size*max_filter_size];
    float upper_bounds[max_filter_size*max_filter_size];
};

typedef std::unordered_map<int, Filter> FilterMap;
extern FilterMap kb_filter;

struct FilterList{
    FilterMap filters;
    int gid = 0;
    int max_size_limit = N_filter_ol;
};

struct Leaf
{
    opType featureNumber;
    float lowerBound;
    float upperBound;
};

extern  int cf_gid;
struct CodeFragment
{
    std::vector<opType> reverse_polish;
    int num_filters; // equal to number of leaves/filters to be determined at run time
    std::vector<int> filter_id;  // leaves: ids of the filters included in this code fragment
    int cf_id;
    CodeFragment(){
        reverse_polish.reserve(cfMaxLength);
        reverse_polish.assign(cfMaxLength, OPNOP);
        filter_id.reserve(cfMaxLeaf);
        filter_id.assign(cfMaxLeaf, -1);
        num_filters = -1;
        cf_id = -1;
    }
};

typedef std::unordered_map<int, CodeFragment> CodeFragmentMap;
extern CodeFragmentMap kb_cf;
typedef std::vector<CodeFragment> CodeFragmentVector;

struct Classifier
{
    int id = -1;
    CodeFragmentVector cf;
    int action = -1;
    double prediction = 0;
    double predictionError = 0;
    double accuracy = 0;
    double fitness = 0;
    int experience = 0;
    int numerosity = 0;
    double actionSetSize = 0;
    int timeStamp = 0;
};

typedef std::unordered_map<int, Classifier> ClassifierMap;
typedef std::vector<int> ClassifierVector;
struct ClassifierSet{
    ClassifierVector ids;
    ClassifierMap& pop;
    ClassifierSet(int size, ClassifierMap& population) : pop(population){
        ids.reserve(size);
    }
};

struct DataSource{
    float *state;
    int action;

    ~DataSource(){
        delete []state;
    }
};


#endif //CONFIGURATION_H
