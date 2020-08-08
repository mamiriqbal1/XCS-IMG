#include <stdlib.h>
#include <sys/types.h>
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <assert.h>
#include <math.h>
#include <cstring>
#include <vector>
#include "xcsMacros.h"
#include "xtime.h"
#include "configuration.h"
#include "classifier.h"
#include "env.h"
#include "xcs.h"
#include <algorithm>
#include <execinfo.h>
#include <signal.h>

//using namespace std;


double pX;// = 0.5; // 0.8; // 0.04; //0.04; //The probability of applying crossover in an offspring classifier.
double pM;// = 0.5; //0.04; //0.8; //The probability of mutating one allele and the action in an offspring classifier.
double pM_allel;// = 0.1; // number of allels modified during mutation of a filter
double p_promising_filter;// = 0.5;  // probability of using a filter from observed list


bool use_kb = false;
bool Testing = true;
int numActions = 2;
int clfrCondMaxLength = 0; // 784/8; // 64; //32; //300;//condLength/4; // condLength/2 and condLength/4 for 70mux and 135mux respectively.
int cfMaxDepth = 0;
int cfMaxLength = 2;// 2^(cdfMaxDepth+1); //allow for endstop OPNOP
int cfMaxStack = 1;// = (cfMaxArity-1)*(cfMaxDepth-1)+2;
int cfMaxLeaf = 1;// = 4; // 2^cfMaxDepth
std::string inputTrainingFile;
std::string inputTestFile;
std::string kb_cf_file;
std::string kb_filter_file;
std::string output_path;
int trainNumInstances=0;// = 11982;
int testNumInstances=0;// = 1984;
int totalNumInstances=0;// = trainNumInstances + testNumInstances;//4; //for review analysis
int maxProblems=0;// = trainNumInstances; //50 * totalNumInstances; //1*100*1000; //training set = [ (1, 1.5, 2, 2.5, 3,)*100*1000 for 6-, 11-, 20-, 37-, 70-, 135-bits MUX respectively]
int testFrequency=0;// = trainNumInstances; // 1034;
int validation_frequency = 0; // to be initialized
int maxPopSize=  1000; //1 * totalNumInstances; //Specifies the maximal number of micro-classifiers in the population. [ (0.5, 1, 2, 5, 10/20, 50)*1000 for 6-, 11-, 20-, 37-, 70-, 135-bits MUX respectively]
bool analyze = false;
std::string analyze_path;

struct distanceInputClassifier{
  int posClassifier;
  float distance;
};

double executeAction(int action, int stateAction, bool &wasCorrect){  // Executes the action and determines the reward.

    int ret=0;
   // std::cout<<"--- "<<action<<" = " <<stateAction<<"--- ";

    if(action == stateAction)
    {
       // std::cout<<"-- "<<action << " = "<<stateAction<<"--";
        wasCorrect=true;
        ret = maxPayoff;
    }
    else
    {
        wasCorrect=false;
        ret = 0;
    }
    return (double)ret;
}

void header(){
    time_t t;
    char hostName[] = "act-PC"; //[MAXHOSTNAMELEN];
    time(&t);
    //gethostname(hostName,sizeof(hostName));
    printf("%s\n Seed: %ld %s",hostName,getSeed(),ctime(&t));
}//end header

void Exit(FILE *fp){
    printf("\nShutting down...\n");
    elapsedTime();
    header();
    printf("Elapsed Time: ");
    elapsed(fp);
    printf(" Secs\n");
}//end exit

/*##########---- Output ----##########*/
/**
 * Writes the performance of the XCS to the specified file.
 * The function writes time performance systemError actualPopulationSizeInMacroClassifiers.
 * Performance and system error are averaged over the last fifty trials.
 *
 * @param performance The performance in the last fifty exploration trials.
 * @param sysError The system error in the last fifty exploration trials.
 * @param problem_count The number of exploration trials executed so far.
 */
void writePerformance(ClassifierMap &pop, int *performance, double *sysError, int problem_count,
                      std::ofstream &output_training_file) {
    double perf=0.0, serr=0.0;
    int setSize;
    for(int i=0; i<testFrequency; i++)
    {
        perf+=performance[i];
        serr+=sysError[i];
    }
    perf/=testFrequency;
    serr/=testFrequency;
    setSize = pop.size();
    output_training_file << problem_count << " " << perf << " " << serr << " " << setSize << std::endl;
    std::cout << "Training: " << problem_count << "  accuracy: " << perf << "  error: " << serr << "  set size: " << setSize << std::endl;
}

/*******************************Write Test Performance*************************************/
void writeTestPerformance(ClassifierMap &pop, int *performance, double *sysError, int problem_count,
                          std::ofstream &output_test_file, int training_problem_count) {
    double perf=0.0, serr=0.0;

    for(int i=0; i<testNumInstances; i++)
    {
        perf+=performance[i];
        serr+=sysError[i];
    }
    perf/=testNumInstances;
    serr/=testNumInstances;
    int setSize = pop.size();
    output_test_file << training_problem_count << " " << problem_count << " " << perf << " " << serr << " " << setSize << std::endl;
    std::cout << "Validation: " << training_problem_count << " " << problem_count << "  accuracy: " << perf << "  error: " << serr << "  set size: " << setSize << std::endl;
}



// new unified single step using epsilon greedy strategy
void doOneSingleStepProblem(ClassifierMap &pop, DataSource *object, int counter, int img_id, int correct[],
                            double sysError[]) {

    bool wasCorrect = false;
    ClassifierSet match_set(maxPopSize, pop);
    ClassifierSet action_set(maxPopSize, pop);

    getMatchSet(pop, match_set, object->state, counter, object->action, img_id);
    getPredictionArray(match_set);

    int actionWinner;
    double reward;
    bool explore = false;

    //if(drand() < epsilon - ((float)counter/(float)maxProblems)/2.0){  // steadily decrease the exploration probability
    if(drand() < epsilon){
        explore = true;
    }else{
        explore = false;
    }

    if(explore) {
        actionWinner = randomActionWinner();
    }else{
        actionWinner = bestActionWinner();
    }
    getActionSet(actionWinner, match_set, action_set);
    reward = executeAction(actionWinner,object->action,wasCorrect);

    updateActionSet(action_set, 0.0, reward);

    // apply GA only with exploration step
    if(explore) {
        discoveryComponent(action_set, pop, counter, object->state);
    }

    // get best action for statistics
    actionWinner = bestActionWinner();
    reward = executeAction(actionWinner,object->action,wasCorrect);

    if(wasCorrect)
    {
        correct[counter%testFrequency]=1;
    }
    else
    {
        correct[counter%testFrequency]=0;
    }
    sysError[counter%testFrequency] = absoluteValue(reward - getBestValue());
}


/*
 * Code review notes:
 * Instead of using epsilon-greedy strategy for exploration and exploitation, the code is doing
 * only exploration. However it always performs a second pass for monitoring the performance.
 * It should be changed to use epsilon-greedy strategy and also the performance monitoring part should be
 * modified to calculate training and validation accuracies after every epoch.
 */
void doOneSingleStepExperiment(ClassifierMap &pop) {  //Executes one single-step experiment monitoring the performance.

    std::ofstream output_training_file;
    output_training_file.open(output_path + output_training_file_name);
    if(!output_training_file.is_open()){
        std::cout << "Could not open output training file";
        exit(1);
    }
    std::ofstream output_test_file;
    output_test_file.open(output_path + output_test_file_name);
    if(!output_test_file.is_open()){
        std::cout << "Could not open output test file";
        exit(1);
    }

    int correct[testFrequency];
    double sysError[testFrequency];

    DataSource *state = NULL;
    //int counter = 0;
    //int index = 0;
    for(int problem_count=0; problem_count <= maxProblems; problem_count++)
    {
        int pop_size = pop.size();
        int pop_numerosity = get_pop_numerosity(pop);
        //std::cout<<problem_count<<"/"<<maxProblems<<"  "<<pop_size<<"/"<<pop_numerosity<<"\r";
        // state = inputArray[irand(totalNumInstances)];
        //index = ;
        //state = &inputArray[irand(totalNumInstances)];
        int img_id = irand(trainNumInstances);
        state = &trainingData[img_id];

        doOneSingleStepProblem(pop, state, problem_count, img_id, correct, sysError);

       if(problem_count % testFrequency == 0 && problem_count > 0){
           writePerformance(pop, correct, sysError, problem_count, output_training_file);
        }
        if(problem_count % validation_frequency == 0 && problem_count > 0){
            doOneSingleStepTest(pop, problem_count, output_test_file);
        }
        if(problem_count % filter_list_management_frequency == 0 && problem_count > 0){
            manage_filter_list(pop);
        }
    }
    output_training_file.close();
    output_test_file.close();
}

void sortAll(distanceInputClassifier arrayToBeSorted[], int size){ // Bubble Sort Function for Ascending Order
  int i, j;
  bool flag = true;    // set flag to true to start first pass
  distanceInputClassifier temp;          // holding variable
  for(i = 1; i<=size && flag==true; i++){
    flag = false;
    for(j=0; j < size-1; j++){
      if(arrayToBeSorted[j+1].distance < arrayToBeSorted[j].distance){
	temp = arrayToBeSorted[j];  // swap elements
        arrayToBeSorted[j] = arrayToBeSorted[j+1];
        arrayToBeSorted[j+1] = temp;
        flag = true;  // indicates that a swap occurred.
      }
    }
  }
  return;   //arrays are passed to functions by address; nothing is returned
}

void sortK(distanceInputClassifier arrayToBeSorted[], int size){ // Bubble Sort Function for Ascending Order
  int i, j;
  bool flag = true;    // set flag to true to start first pass
  distanceInputClassifier temp;          // holding variable
  for(i = 1; i<=size && flag==true; i++){
    flag = false;
    for(j=0; j < size-1; j++){
      if(arrayToBeSorted[j+1].posClassifier < arrayToBeSorted[j].posClassifier){
	temp = arrayToBeSorted[j];  // swap elements
        arrayToBeSorted[j] = arrayToBeSorted[j+1];
        arrayToBeSorted[j+1] = temp;
        flag = true;  // indicates that a swap occurred.
      }
    }
  }
  return;   //arrays are passed to functions by address; nothing is returned
}

std::string NumberToString(int num){
    std::stringstream ss;
    ss << num;
    return ss.str();
}

void doOneSingleStepTest(ClassifierMap &pop, int training_problem_count, std::ofstream &output_test_file) {
	bool wasCorrect = false;
	int correctCounter = 0;
	int correct[testNumInstances];
	double sysError[testNumInstances];
	DataSource *testState = NULL;
	int tmpcorrectcounter = 0;
    int tmpnotmatched = 0;
	int TP = 0, TN = 0, FP = 0, FN = 0, TNeu = 0, FNeu = 0;
    int class1TP=0, class2TP=0, class3TP=0, class4TP=0, class5TP=0, class6TP=0, class7TP=0, class8TP=0, class9TP=0, class10TP=0;
    int class11TP=0, class12TP=0, class13TP=0, class14TP=0, class15TP=0, class16TP=0, class17TP=0, class18TP=0, class19TP=0, class20TP=0, class21TP=0, class22TP=0, class23TP=0;
	int class1FP=0, class2FP=0, class3FP=0, class4FP=0, class5FP=0, class6FP=0, class7FP=0, class8FP=0, class9FP=0, class10FP=0;
	int class11FP=0, class12FP=0, class13FP=0, class14FP=0, class15FP=0, class16FP=0, class17FP=0, class18FP=0, class19FP=0, class20FP=0, class21FP=0, class22FP=0, class23FP=0;


    int cc = 0;

	for(int t=0; t<testNumInstances; t++){
        ClassifierSet match_set(maxPopSize, pop);
        //std::cout<<t<<"/"<<testNumInstances<<"\r";
		bool isMatched = false;
		testState = &testingData[t];

        get_matching_classifiers(pop, testState->state, match_set, t, false);
        isMatched = match_set.ids.size() > 0;
        if(isMatched == false){
            cc++;
            tmpnotmatched++;
            int popSize = pop.size();
			distanceInputClassifier distanceArray[popSize];
			int i = 0;
			int k = popSize*tournamentSize;
			assert(k > 0); // make sure that k > 0

			for(auto& item : pop){
			  distanceArray[i].posClassifier = i;
			  distanceArray[i].distance = computeDistance(item.second, testState->state);
			  i++;
			}
			sortAll(distanceArray,popSize);
			sortK(distanceArray,k);
			int mK = 0, mT = 0;
			for(auto it = pop.begin(); it != pop.end() && mK < k; it++,mT++){
			  if(distanceArray[mK].posClassifier == mT){
			      match_set.ids.push_back(it->first);
			    mK++;
			  }
		   	}
		}
		getPredictionArray(match_set);
        int actionWinner = bestActionWinner();
        double reward = executeAction(actionWinner,testState->action,wasCorrect);
        sysError[t] = absoluteValue(reward - getBestValue());

		if(wasCorrect)
        {
                correct[t]=1;
                correctCounter++;
        }
        else
        {
               correct[t]=0;
        }
    }
    writeTestPerformance(pop, correct, sysError, testNumInstances, output_test_file, training_problem_count);

//    std::cout<<"N = "<<maxPopSize<<std::endl;
//    std::cout<<"P# = "<<P_dontcare<<std::endl;
//    std::cout<<"Total Number of Instances = "<<testNumInstances<<std::endl;
//    std::cout<<"Number of KNNs = "<<cc<<std::endl;
//    std::cout<<"Number of correct Instances = "<<correctCounter<<std::endl;
//    //std::cout<<"TP: "<<TP<<"--TN: "<<TN<<"--FP: "<<FP<<"--FN: "<<FN<<std::endl;
//    std::cout<<"\nAccuracy: "<< (correctCounter)*1.0/(testNumInstances) << "\n";
}



void startXCS(){
    startTimer();
    ClassifierMap pop;
    printf("\nLoading Input! Please wait ....\n");
    //inp = new ds[totalNumInstances];
    //loadData(inp);

    //input = new DataSource[totalNumInstances];
    //initializeInput(input,totalNumInstances);
    //loadData(input);
    //loadDataFile();
    //loadDataFile(input);
    trainingData = new DataSource[trainNumInstances];
    testingData = new DataSource[testNumInstances];
    initializeInput(trainingData,trainNumInstances);
    initializeInput(testingData,testNumInstances);
    loadDataFromFile(trainingData, inputTrainingFile.c_str(), trainNumInstances);
    loadDataFromFile(testingData, inputTestFile.c_str(), testNumInstances);
    updateRange(trainingData,trainNumInstances);
    updateRange(testingData,testNumInstances);
    if(use_kb) {
        load_kb(kb_cf_file, kb_filter_file);
    }

    printf("\nIt is in progress! Please wait ....\n");

    doOneSingleStepExperiment(pop);

    // clean filter list finally before writing to file
    manage_filter_list(pop);


    save_experiment_results(pop);
    //delete []input;
    delete []testingData;
    delete []trainingData;

}//end startXCS

int CountLines(const char* file)
{
    int count=0;
    // std::ifstream is RAII, i.e. no need to call close
    std::ifstream cFile (file);
    if (cFile.is_open()){
        std::string line;
        while(getline(cFile, line)){
            count++;
        }
    }
    else {
        std::string error("Error opening input file: ");
        error.append(file).append(", could not load data!");
        throw std::runtime_error(error);
    }
    return count;
}



void LoadConfig(char* file)
{
    int epochs = 0;
    bool create_output_file = true;
    std::ofstream output_config_file;
    // std::ifstream is RAII, i.e. no need to call close
    std::ifstream cFile (file);
    if (cFile.is_open())
    {
        std::string line;
        while(getline(cFile, line)){
            line.erase(std::remove_if(line.begin(), line.end(), isspace),
                       line.end());
            if(line[0] == '#' || line.empty())
                continue;
            auto delimiterPos = line.find("=");
            auto name = line.substr(0, delimiterPos);
            auto value = line.substr(delimiterPos + 1);
            if(name == "train_file_name"){
                inputTrainingFile = value;
            }else if(name == "test_file_name"){
                inputTestFile = value;
            }else if(name == "num_actions"){
                numActions = atoi(value.c_str());
            }else if(name == "max_depth"){
                cfMaxDepth = atoi(value.c_str());
                cfMaxLength = std::pow(2, cfMaxDepth+1);
                cfMaxStack = (cfMaxArity-1)*(cfMaxDepth-1) + 2;
                cfMaxLeaf = std::pow(2, cfMaxDepth);
            }else if(name == "max_condition_length"){
                clfrCondMaxLength = atoi(value.c_str());
            }else if(name == "pX"){
                pX = atof(value.c_str());
            }else if(name == "pM"){
                pM = atof(value.c_str());
            }else if(name == "pM_allel"){
                pM_allel = atof(value.c_str());
            }else if(name == "p_promising_filter"){
                p_promising_filter = atof(value.c_str());
            }else if(name == "use_kb"){
                if(value == "no"){
                    use_kb = false;
                }else if(value == "yes"){
                    use_kb = true;
                }
            }else if(name == "kb_cf_file"){
                kb_cf_file = value;
            }else if(name == "kb_filter_file"){
                kb_filter_file = value;
            }else if(name == "output_path"){
                output_path = value;
            }else if(name == "max_population_size"){
                maxPopSize = atoi(value.c_str());
            }else if(name == "epochs"){
                epochs = atoi(value.c_str());
            }else if(name == "test_frequency"){
                testFrequency= atoi(value.c_str());
            }else if(name == "max_problems"){
                maxProblems= atoi(value.c_str());
            }else if(name == "analyze"){
                if(value == "no"){
                    analyze = false;
                }else{
                    analyze = true;
                }
            }else if(name == "analyze_path"){
                analyze_path = value;
            }else if (name == "testing") {
                if (value == "no") {
                    Testing = false;
                } else if (value == "yes") {
                    Testing = true;
                }
            }
            if(create_output_file){
                create_output_file = false;
                if(output_path.empty()){
                   std::cout<<"Output path must be the first parameter in config file"<<std::endl;
                   exit(1);
                }
                mkdir(output_path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
                output_config_file.open(output_path + "experiment_config.txt");
                if(!output_config_file.is_open()){
                    std::cout << "Could not open output config file";
                    exit(1);
                }
            }
            std::cout << name << " " << value << '\n';
            output_config_file<< name << " " << value << std::endl;
        }
        output_config_file.close();
    }
    else {
        std::string error("Error opening input file: ");
        error.append(file).append(", could not load data!");
        throw std::runtime_error(error);
    }
    // calculate the number of train and test instances
    trainNumInstances = CountLines(inputTrainingFile.c_str());
    std::cout<<"Training instances: "<<trainNumInstances<<std::endl;
    testNumInstances = CountLines(inputTestFile.c_str());
    std::cout<<"Test instances: "<<testNumInstances<<std::endl;
    totalNumInstances = trainNumInstances + testNumInstances;
    if(testFrequency == 0) {
        testFrequency = trainNumInstances;
    }
    std::cout<<"Test frequency: "<<testFrequency<<std::endl;
    if(validation_frequency == 0){
        validation_frequency = trainNumInstances;
    }
    std::cout<<"Validation frequency: "<<validation_frequency<<std::endl;
    // calculate number of training examples to be presented based on epochs
    if(maxProblems == 0) {
        maxProblems = trainNumInstances * epochs;
    }

}

/*
void copy_filter_to_condition(Classifier *classifer, xt::xtensor<float, 1> filter)
{
     Leaf *leaf = classifer->code_fragment[0].leaf;

     for (int i = 0; i < cfMaxLeaf; i++) {
        leaf[i].lowerBound = filter(i*2);
        leaf[i].upperBound = filter(i*2+1);
     }
}

void count_matches_for_filters(xt::xtensor<float, 2> good_filters, xt::xtensor<float, 2> good_actions)
{

    int num_filters = good_actions.shape()[0];
    Classifier *classifier;
    // create a dummy classifier code_fragment and reuse for all filters later
    DataSource *obj = &trainingData[0];
    classifier = matchingCondAndSpecifiedAct(obj->state, 0, 1, 1);
    CodeFragment *cfs = classifier->code_fragment;

    for (int i = 0; i < num_filters; i++) {
        int match_0 = 0;
        int match_1 = 0;
        int matched = 0;

        copy_filter_to_condition(classifier, xt::view(good_filters, i, xt::all()));
        for (int j = 0; j < trainNumInstances; j++) {
            obj = &trainingData[j];
            if(isConditionMatched(cfs, obj->state, -1, -1)){
                matched++;
                if(obj->action == 0){
                    match_0++;
                }else{
                    match_1++;
                }
            }
        }

        std::cout << good_actions(i,0) << " , " << matched << " , " << match_0 << " , " << match_1 << std::endl;
    }
}
*/

void analyze_rules()
{

    trainingData = new DataSource[trainNumInstances];
    testingData = new DataSource[testNumInstances];
    initializeInput(trainingData,trainNumInstances);
    initializeInput(testingData,testNumInstances);
    loadDataFromFile(trainingData, inputTrainingFile.c_str(), trainNumInstances);
    loadDataFromFile(testingData, inputTestFile.c_str(), testNumInstances);
    updateRange(trainingData,trainNumInstances);
    updateRange(testingData,testNumInstances);

    // load good filters
//    std::ifstream filters_file;
//    filters_file.open(analyze_path + "all_filters.txt");
//    auto good_filters = xt::load_csv<float>(filters_file);
//    std::ifstream actions_file;
//    actions_file.open(analyze_path + "all_actions.txt");
//    auto good_actions = xt::load_csv<float>(actions_file);


    //count_matches_for_filters(good_filters, good_actions);



}


void handler(int sig) {
    void *array[100];
    size_t size;

    // get void*'s for all entries on the stack
    size = backtrace(array, 100);

    // print out all the frames to stderr
    fprintf(stderr, "Error: signal %d:\n", sig);
    backtrace_symbols_fd(array, size, STDERR_FILENO);
    exit(1);
}


int main(int argc, char **argv){
    signal(SIGSEGV, handler);   // install our handler for stack trace in case of a signal

    if(argc == 1){
        std::cout << "Please provide experiment config file" << std::endl;
        return -1;
    }
    LoadConfig(argv[1]);
    if(analyze){
        analyze_rules();
        return 0;
    }

    for(int j=0; j<run; j++)
    {
        mkdir(output_path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        //mkdir(path);
        setSeed(seeds[j]);

        if (use_kb)
        {
            // open kb file
        }
        startXCS();
    }

    exit(0);
}//end main
