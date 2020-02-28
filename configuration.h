extern bool use_kb;
extern std::string kb_file;
extern std::string output_path;
extern std::string input_path;
extern int numActions; // = 2; //0 or 1
const int trainNumInstances = 11982;
const int testNumInstances = 1984;
const int totalNumInstances = trainNumInstances + testNumInstances;//4; //for review analysis

const int condLength = 784; //512; //256; //4626; //for review analysis 300, 500, 1000, 1500, 2000, 2257
const int maxPopSize=  1000; //1 * totalNumInstances; //Specifies the maximal number of micro-classifiers in the population. [ (0.5, 1, 2, 5, 10/20, 50)*1000 for 6-, 11-, 20-, 37-, 70-, 135-bits MUX respectively]
const int maxProblems = trainNumInstances; //50 * totalNumInstances; //1*100*1000; //training set = [ (1, 1.5, 2, 2.5, 3,)*100*1000 for 6-, 11-, 20-, 37-, 70-, 135-bits MUX respectively]
const double maxPayoff = 1000;
const int clfrCondLength = 784/8; // 64; //32; //300;//condLength/4; // condLength/2 and condLength/4 for 70mux and 135mux respectively.
const int run = 1;
//const char inputFile[] = "features1.txt";
extern std::string inputTrainingFile; //[] = "../data/mnist/3_8_train_mnist.txt";
extern std::string inputTestFile; //[] = "../data/mnist/3_8_test_mnist.txt";
//const char testFile[] = "";
const bool Testing = true;
//const int currentProblemLevel = 1; //must be set coz previous level and current level file names are adjusted

const int precisionDigits = 2;
const int testFrequency = trainNumInstances; // 1034;

const char outputFileName[] = "output_training.txt";
const char featureFileName[] = "feature_codefragments.txt";
const char ruleFileName[] = "rule_with_codefragements.txt";
const char resultFile[] = "result_testing.txt";

const int cfMaxDepth = 2;
const int cfMinDepth = 0;
const int cfMaxLength = 8;// pow(2,adfMaxDepth+1); //allow for endstop OPNOP
const int cfMaxArity = 2;
const int cfMaxStack = (cfMaxArity-1)*(cfMaxDepth-1)+2;

typedef int opType;
const int opSize = sizeof(opType);

const opType OPNOP = -100; //to be used as ending symbol
const opType OPAND = -101;
const opType OPOR = -102;
const opType OPNAND = -103;
const opType OPNOR = -104;
const opType OPNOT = -105;
//const opType OPUNITY = -106;

const int totalFunctions = 5;
const opType functionCodes[] = {OPAND,OPOR,OPNAND,OPNOR,OPNOT};

const int numLeaf = 4;

struct Leaf
{
    opType featureNumber;
    float lowerBound;
    float upperBound;
};

struct CodeFragment
{
    opType codeFragment[cfMaxLength];
    Leaf leaf[numLeaf];
    int cfID;
};

struct Classifier
{
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


