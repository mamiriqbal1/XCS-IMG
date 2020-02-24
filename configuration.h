const int numActions = 2; //0 or 1
/*
//  Two Classes Test
const int class1Instances = 1078; //
const int class2Instances = 871; //
*/

// Test H : Gernalization To Specializaiton- (Person(imagent)+ship(stl10)) to  (scubadiver+snorkldiving)+((aircraftcareer+battleship+destroyer)Imagent)
//const int class1Instances = 1044; //scubadiver 313
//const int class2Instances = 771; //snorkling   231
//const int class3Instances = 1300; //aircraftcareer  390
//const int class4Instances = 821; //battleship  246
//const int class5Instances = 640; //person destroyer  192


// For Two Classes
//const int totalNumInstances = class1Instances + class2Instances + 1;
////  For Five Classes
//const int totalNumInstances  = // class1Instances + class2Instances + class3Instances + class4Instances + class5Instances + 1;
////  For Ten Classes
//const int totalNumInstances  = class1Instances + class2Instances + class3Instances + class4Instances + class5Instances + class6Instances + class7Instances + class8Instances + class9Instances + class10Instances+1;
////  For Twenty Classes
//const int totalNumInstances = class1Instances + class2Instances + class3Instances + class4Instances + class5Instances + class6Instances + class7Instances + class8Instances + class9Instances + class10Instances + class11Instances + class12Instances + class13Instances + class14Instances + class15Instances + class16Instances + class17Instances + class18Instances + class19Instances + class20Instances + class21Instances + class22Instances + class23Instances + 1;

const int trainNumInstances = 11982;
const int testNumInstances = 1984;
const int totalNumInstances = trainNumInstances + testNumInstances;//4; //for review analysis

const int condLength = 784; //512; //256; //4626; //for review analysis 300, 500, 1000, 1500, 2000, 2257
//const int posBits = 2; //2, 3, 4, 5, 6, and 7 for 6-, 11-, 20-, 37, 70-, and 135-bits MUX respectively
const int maxPopSize=  1000; //1 * totalNumInstances; //Specifies the maximal number of micro-classifiers in the population. [ (0.5, 1, 2, 5, 10/20, 50)*1000 for 6-, 11-, 20-, 37-, 70-, 135-bits MUX respectively]

//const int class1TrainNumInstance = 0.7 * class1Instances;
//const int class2TrainNumInstance = 0.7 * class2Instances;
//const int class3TrainNumInstance = 0.7 * class3Instances;
//const int class4TrainNumInstance = 0.7 * class4Instances;
//const int class5TrainNumInstance = 0.7 * class5Instances;
/*const int class6TrainNumInstance = 0.7 * class6Instances;
const int class7TrainNumInstance = 0.7 * class7Instances;
const int class8TrainNumInstance = 0.7 * class8Instances;
const int class9TrainNumInstance = 0.7 * class9Instances;
const int class10TrainNumInstance = 0.7 * class10Instances;
/*const int class11TrainNumInstance = 0.7 * class11Instances;
const int class12TrainNumInstance = 0.7 * class12Instances;
const int class13TrainNumInstance = 0.7 * class13Instances;
const int class14TrainNumInstance = 0.7 * class14Instances;
const int class15TrainNumInstance = 0.7 * class15Instances;
const int class16TrainNumInstance = 0.7 * class16Instances;
const int class17TrainNumInstance = 0.7 * class17Instances;
const int class18TrainNumInstance = 0.7 * class18Instances;
const int class19TrainNumInstance = 0.7 * class19Instances;
const int class20TrainNumInstance = 0.7 * class20Instances;
const int class21TrainNumInstance = 0.7 * class21Instances;
const int class22TrainNumInstance = 0.7 * class22Instances;
const int class23TrainNumInstance = 0.7 * class23Instances;
*/

// for 2 classes
//const int trainNumInstances = class1TrainNumInstance + class2TrainNumInstance;
// For 5 Five classes
//const int trainNumInstances = class1TrainNumInstance + class2TrainNumInstance + class3TrainNumInstance + class4TrainNumInstance + class5TrainNumInstance; //0.5 * totalNumInstances-1; //class1Instances + class2Instances - 380;  // 25000;
// For 10 ten classes
//const int trainNumInstances = class1TrainNumInstance + class2TrainNumInstance + class3TrainNumInstance + class4TrainNumInstance + class5TrainNumInstance + class6TrainNumInstance + class7TrainNumInstance + class8TrainNumInstance + class9TrainNumInstance + class10TrainNumInstance; //0.5 * totalNumInstances-1; //class1Instances + class2Instances - 380;  // 25000;
// For 20 twenty classes
//const int trainNumInstances = class1TrainNumInstance + class2TrainNumInstance + class3TrainNumInstance + class4TrainNumInstance + class5TrainNumInstance + class6TrainNumInstance + class7TrainNumInstance + class8TrainNumInstance + class9TrainNumInstance + class10TrainNumInstance + class11TrainNumInstance + class12TrainNumInstance + class13TrainNumInstance + class14TrainNumInstance + class15TrainNumInstance + class16TrainNumInstance + class17TrainNumInstance + class18TrainNumInstance + class19TrainNumInstance + class20TrainNumInstance + class21TrainNumInstance + class22TrainNumInstance + class23TrainNumInstance;

//const int testNumInstances = totalNumInstances - trainNumInstances; //0.5 * totalNumInstances; //380; // 20% of data
//const int negTrainNumInstance = 0.5 * class1Instances; //664; //339; //664; //652; //755; //651; //736; //200;
//const int posTrainNumInstance = 0.5 * class2Instances; //731; //371; //731; //535; //610; //535; //682; //780; //0;

const int maxProblems = trainNumInstances; //50 * totalNumInstances; //1*100*1000; //training set = [ (1, 1.5, 2, 2.5, 3,)*100*1000 for 6-, 11-, 20-, 37-, 70-, 135-bits MUX respectively]
const double maxPayoff = 1000;
const int clfrCondLength = 784/8; // 64; //32; //300;//condLength/4; // condLength/2 and condLength/4 for 70mux and 135mux respectively.
const int run = 1;
//const char inputFile[] = "features1.txt";
const char inputTrainingFile[] = "../data/mnist/3_8_train_mnist.txt";
const char inputTestFile[] = "../data/mnist/3_8_test_mnist.txt";
//const char testFile[] = "";
const bool Testing = true;
const int currentProblemLevel = 1; //must be set coz previous level and current level file names are adjusted
const int precisionDigits = 2;
const int testFrequency = trainNumInstances; // 1034;

const char outputFileName[] = "/output_training.txt";
const char featureFileName[] = "/feature_codefragments.txt";
const char ruleFileName[] = "/rule_with_codefragements.txt";
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


