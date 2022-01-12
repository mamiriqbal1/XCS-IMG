#ifndef CONFIGURATION_H
#define CONFIGURATION_H

#include <string>
#include <forward_list>
#include <unordered_map>
#include "xcsMacros.h"
#include <vector>
#include <list>
#include <algorithm>
typedef std::vector<double> DoubleVector;
typedef std::vector<float> FloatVector;
typedef std::vector<int> IntVector;
typedef std::vector<char> CharVector;
typedef std::vector<FloatVector> FloatMatrix;
typedef std::vector<IntVector> IntMatrix;
typedef std::vector<CharVector> CharMatrix;
typedef std::unordered_map<int, int> IntMap;

extern bool normalize_data;
const float IMAGE_MAX_VALUE = 255.0;
const float IMAGE_MIN_VALUE = 0;
extern bool use_kb;
extern std::string kb_file;
extern std::string output_path;
extern int numActions; // = 2; //0 or 1
extern IntMap class_map;
// max_actions have to be set such that it is always more than the actual number of actions in an experiment
extern int trainNumInstances;// = 11982;
extern int testNumInstances;// = 1984;
extern int totalNumInstances;// = trainNumInstances + testNumInstances;//4; //for review analysis

extern int condLength; // = 784; //512; //256; //4626; //for review analysis 300, 500, 1000, 1500, 2000, 2257
extern int image_width; // = 28;
extern int image_height; // = 28;
extern int maxPopSize; //=  1000; //1 * totalNumInstances; //Specifies the maximal number of micro-classifiers in the population. [ (0.5, 1, 2, 5, 10/20, 50)*1000 for 6-, 11-, 20-, 37-, 70-, 135-bits MUX respectively]
extern int maxProblems;// = trainNumInstances; //50 * totalNumInstances; //1*100*1000; //training set = [ (1, 1.5, 2, 2.5, 3,)*100*1000 for 6-, 11-, 20-, 37-, 70-, 135-bits MUX respectively]
const double maxPayoff = 1000;
const int run = 1;
extern std::string inputTrainingFile; //[] = "../data/mnist/3_8_train_mnist.txt";
extern std::string inputTestFile; //[] = "../data/mnist/3_8_test_mnist.txt";

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
const char output_visualization_file_name[] = "visualization.txt";
const char output_done_file_name[] = "done.txt";  // file indicating completion of experiment
const char output_parameter_file_name[] = "parameter.txt";  // file indicating completion of experiment

extern bool visualization;
extern int clfrCondMaxLength;//2 // 784/8; // 64; //32; //300;//condLength/4; // condLength/2 and condLength/4 for 70mux and 135mux respectively.
extern int cfMaxDepth;// = 2;
extern int cfMinDepth; // = 0;
extern int cfMaxLength;// = 8;// 2^(cdfMaxDepth+1); //allow for endstop OPNOP
const int cfMaxArity = 2;
extern int cfMaxStack;// = (cfMaxArity-1)*(cfMaxDepth-1)+2;
extern int cfMaxLeaf;// = 4; // 2^cfMaxDepth
extern int cfMinLeaf;// = 1; // 2^cfMinDepth




// CF will now have a bounding box within which all filters will lie.
extern int cf_max_bounding_box_size;
extern int cf_min_bounding_box_size;

struct BoundingBox
{
    // bounding box related parameters; x, y coordinates and size, height
    int x = -1;
    int y = -1;
    int size_x = -1;  // size_x == pixel width, for a 1x1 bb size_x == 1, for 2x2 size_x == 2, it can not be zero
    int size_y = -1;  // end coordinate of bb end_x = x + size_x -1
};

struct Position
{
    // relative position of filter with respective to CF bounding box
    int x = -1;
    int y = -1;
};


const int NOT_INITIALIZED = -1;
const int ENABLED = 1;
const int DISABLED = 0;
const char UNKNOWN = -1;
const char MATCHED = 1;
const char NOT_MATCHED = 0;
const int NUMEROSITY_INITIALIZATION = -1;
const int FITNESS_INITIALIZATION = 0;
const float BINARY_THRESHOLD_INITIALIZATION = 0.5;

extern float min_binary_threshold;
extern float max_binary_threshold;

struct CodeFragment
{
    int cf_id = NOT_INITIALIZED;
    int numerosity = NUMEROSITY_INITIALIZATION;
    int fitness = FITNESS_INITIALIZATION; // Fitness of a code fragment is its appearance in "promising classifiers"
    BoundingBox bb;
    FloatMatrix pattern;
    float matching_threshold = NOT_INITIALIZED;
    float binary_threshold = BINARY_THRESHOLD_INITIALIZATION; // the threshold used to binarize image region before extracting pattern
};

//typedef std::unordered_map<int, CodeFragment> CodeFragmentMap;
//extern CodeFragmentMap kb_cf;
typedef std::vector<CodeFragment> CodeFragmentVector;

extern int classifier_gid;
struct Classifier
{
    int id = -1;
    std::vector<int> cf_ids;
    int action = -1;
    double prediction = 0;
    double predictionError = 0;
    double accuracy = 0;
    double fitness = 0;
    int experience = 0;
    int numerosity = 0;
    double actionSetSize = 0;
    int timeStamp = 0;
    BoundingBox bb; // a calculated value based on cf bounding boxes
    Classifier() : cf_ids(clfrCondMaxLength){
        cf_ids.assign(clfrCondMaxLength, -1);
    }
};


typedef std::vector<Classifier> ClassifierVector;
typedef std::vector<int> ClassifierIDVector;
struct ClassifierSet{
    ClassifierIDVector ids;
    ClassifierSet(int size){
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
