#include <stdlib.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include "xcsMacros.h"
#include "configuration.h"
#include "classifier.h"
#include "env.h"
#include "xcs.h"
#include <algorithm>
#include <execinfo.h>
#include <signal.h>
#include "codeFragment.h"
#include "cf_list.h"
#include <chrono>
#include <ctime>
#include <iomanip>

auto start = std::chrono::system_clock::now();
float epsilon = 0.5; // the probability of exploration for epsilon greedy strategy
float epsilon_start = 0.5; // the probability of exploration for epsilon greedy strategy
float epsilon_end = 0.5; // the probability of exploration for epsilon greedy strategy
double beta = 0.2;
double beta_start = 0.2;
double beta_end = 0.2;
double pX = 0.8;// = 0.5; // 0.8; // 0.04; //0.04; //The probability of applying crossover in an offspring classifier.
double pX_start = 0.8;// = 0.5; // 0.8; // 0.04; //0.04; //The probability of applying crossover in an offspring classifier.
double pX_end = 0.8;// = 0.5; // 0.8; // 0.04; //0.04; //The probability of applying crossover in an offspring classifier.
double pM = 0.2;// = 0.5; //0.04; //0.8; //The probability of mutating one allele and the action in an offspring classifier.
double pM_start = 0.2;
double pM_end = 0.2; // parameter control during execution after every epoch
double p_kb = 0.5;
double p_kb_start = 0.5;
double p_kb_end = 0.5;
double pM_allel;// = 0.1; // number of allels modified during mutation of a filter
double p_promising;// = 0.5;  // probability of using a filter from observed list

bool visualization = false;
bool fixed_seed = true;
bool normalize_data = true;
bool use_kb = false;
bool Testing = false;
int numActions = 2;
IntMap class_map;
int image_width = NOT_INITIALIZED;
int image_height = NOT_INITIALIZED;
int condLength = NOT_INITIALIZED;
int clfrCondMaxLength = 0; // 784/8; // 64; //32; //300;//condLength/4; // condLength/2 and condLength/4 for 70mux and 135mux respectively.
int cfMaxDepth = 0;
int cfMinDepth = 0;
int cfMaxLength = 2;// 2^(cdfMaxDepth+1); //allow for endstop OPNOP
int cfMaxStack = 1;// = (cfMaxArity-1)*(cfMaxDepth-1)+2;
int cfMaxLeaf = 1;// = 4; // 2^cfMaxDepth
int cfMinLeaf = 1;// = 1; // 2^cfMinDepth
int cf_max_bounding_box_size = image_height;
int cf_min_bounding_box_size = max_filter_size * 2;
float filter_matching_threshold = 0.001; // thresh hold for filter matching

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
std::string resume_from;


double executeAction(int action, int stateAction, bool &wasCorrect){  // Executes the action and determines the reward.

    int ret=0;

    if(action == stateAction)
    {
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
void writePerformance(double performance, double sysError, int problem_count, double action_set_size,
                      std::ofstream &output_training_file) {

    int setSize = get_pop_size(false);
    int setTotal = get_pop_size(true);
    output_training_file << std::fixed;
    output_training_file.precision(4);
    output_training_file << std::setw(7) << problem_count << " " << performance;
    output_training_file.precision(1);
    output_training_file << " " << std::setw(5) << sysError;
    output_training_file << " " << std::setw(5) << setSize << " " << std::setw(5) << setTotal << " " << std::setw(5) << action_set_size << std::endl;
    std::cout << std::fixed;
    std::cout.precision(4);
    std::cout << "Training " << std::setw(7) << problem_count;
    std::cout << " accuracy "  << performance;
    std::cout.precision(1);
    std::cout << " error " << std::setw(5) << sysError;
    std::cout << " set " << std::setw(5) << setSize;
    std::cout << "/" << std::setw(5) << setTotal;
    std::cout << " action set " << std::setw(5) << action_set_size << std::endl;
}

/*******************************Write Test Performance*************************************/
void writeTestPerformance(double performance, double sysError, int problem_count, double action_set_size,
                          std::ofstream &output_test_file, int training_problem_count, double training_accuracy,
                          double training_error) {

    int setSize = get_pop_size(false);
    int setTotal = get_pop_size(true);
    output_test_file << std::fixed;
    output_test_file.precision(4);
    output_test_file << std::setw(7) << training_problem_count << " " << training_accuracy;
    output_test_file.precision(1);
    output_test_file << " " << std::setw(5) << training_error;
    output_test_file.precision(4);
    output_test_file << " " << std::setw(7) << problem_count << " " << performance;
    output_test_file.precision(1);
    output_test_file << " " << std::setw(5) << sysError;
    output_test_file << " " << std::setw(5) << setSize << " " << std::setw(5) << setTotal;
    output_test_file << " " << action_set_size << std::endl;

    std::cout << std::fixed;
    std::cout.precision(4);
    std::cout << "Validation " << training_problem_count << " " << problem_count << " accuracy " << performance;
    std::cout.precision(1);
    std::cout << " error " << std::setw(5) << sysError << " set " << setSize << "/" << setTotal << " action set " << action_set_size << std::endl;
}



// new unified single step using epsilon greedy strategy
void doOneSingleStepProblem(DataSource *object, int counter, int img_id, int &correct, double &sysError,
                            int &match_set_size) {

    bool wasCorrect = false;
    ClassifierSet match_set(maxPopSize);
    ClassifierSet action_set(maxPopSize);

    getMatchSet(match_set, object->state, counter, object->action, img_id);
    match_set_size = match_set.ids.size();
    getPredictionArray(match_set);

    int actionWinner;
    double reward;
    bool explore = false;

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

    updateActionSet(action_set, 0.0, reward);  // maxPrediction parameter is zero for single step

    // apply GA only with exploration step
    if(explore) {
        discoveryComponent(action_set, counter, object->state, actionWinner);
    }

    // get best action for statistics
    actionWinner = bestActionWinner();
    reward = executeAction(actionWinner,object->action,wasCorrect);

    if(wasCorrect){
        correct=1;
    } else {
        correct=0;
    }
    sysError = absoluteValue(reward - getBestValue());
}


void load_state_for_resume()
{
    load_cf_for_resume(output_path + resume_from + output_code_fragment_file_name);
    load_classifier(output_path + resume_from + output_classifier_file_name);
}

/*
 * Code review notes:
 * Instead of using epsilon-greedy strategy for exploration and exploitation, the code is doing
 * only exploration. However it always performs a second pass for monitoring the performance.
 * It should be changed to use epsilon-greedy strategy and also the performance monitoring part should be
 * modified to calculate training and validation accuracies after every epoch.
 */
void doOneSingleStepExperiment() {  //Executes one single-step experiment monitoring the performance.

    std::ofstream output_training_file;
    output_training_file.open(output_path + output_training_file_name, std::ofstream::app);
    if(!output_training_file.is_open()){
        std::cout << "Could not open output training file";
        exit(1);
    }
    std::ofstream output_test_file;
    output_test_file.open(output_path + output_test_file_name, std::ofstream::app);
    if(!output_test_file.is_open()){
        std::cout << "Could not open output test file";
        exit(1);
    }

    int correct = 0, correct_count = 0, epoch_correct_count = 0;
    double sysError = 0, error_sum = 0, epoch_error_sum = 0;
    int match_set_size = 0, match_set_size_sum = 0, epoch_match_set_size_sum = 0;
    DataSource *state = NULL;
    int problem_count=1;

    // load state from previous execution in case resume is require
    if(!resume_from.empty()){
        // insert one empty line before resuming performance recording
        output_training_file<<std::endl;
        output_test_file<<std::endl;
        load_state_for_resume();
        manage_filter_and_cf_list(); // update filter parameters etc
        problem_count = std::atoi(resume_from.c_str());
        problem_count++;  // resume from next problem
        // parameter control from initial value to final value
        epsilon = epsilon_start - (epsilon_start - epsilon_end) * problem_count/maxProblems;
        beta = beta_start - (beta_start - beta_end) * problem_count/maxProblems;
        pX = pX_start - (pX_start - pX_end) * problem_count/maxProblems;
        pM = pM_start - (pM_start - pM_end) * problem_count/maxProblems;
        p_kb = p_kb_start - (p_kb_start - p_kb_end) * problem_count/maxProblems;
    }
    for( ; problem_count <= maxProblems; problem_count++)
    {
        int img_id = irand(trainNumInstances);
        state = &trainingData[img_id];

        doOneSingleStepProblem(state, problem_count, img_id, correct, sysError, match_set_size);
        correct_count += correct;
        error_sum += sysError;
        epoch_correct_count += correct;
        epoch_error_sum += sysError;
        match_set_size_sum += match_set_size;
        epoch_match_set_size_sum += match_set_size;

       if(problem_count % testFrequency == 0 && problem_count > 0){
           double accuracy = correct_count/(double)testFrequency;
           double error = error_sum / testFrequency;
           double avg_action_set_size = (match_set_size_sum / (double)testFrequency) / numActions;
           writePerformance(accuracy, error, problem_count, avg_action_set_size, output_training_file);
           correct_count = 0;
           error_sum = 0;
           match_set_size_sum = 0;
        }
        if(problem_count % validation_frequency == 0 && problem_count > 0){

            double epoch_accuracy = epoch_correct_count/(double)validation_frequency;
            double epoch_error = epoch_error_sum/validation_frequency;
            double avg_epoch_action_set_size = (epoch_match_set_size_sum / (double)validation_frequency) / numActions;
            save_experiment_results(std::to_string(problem_count) + "/"); // save experiment results after every epoch
            doOneSingleStepTest(problem_count, output_test_file, epoch_accuracy, epoch_error,
                                std::to_string(problem_count) + "/");
            epoch_correct_count = 0;
            epoch_error_sum = 0;
            epoch_match_set_size_sum = 0;
        }
        if(problem_count % filter_list_management_frequency == 0 && problem_count > 0){
            manage_filter_and_cf_list();
        }
        // parameter control from initial value to final value
        epsilon = epsilon_start - (epsilon_start - epsilon_end) * problem_count/maxProblems;
        beta = beta_start - (beta_start - beta_end) * problem_count/maxProblems;
        pX = pX_start - (pX_start - pX_end) * problem_count/maxProblems;
        pM = pM_start - (pM_start - pM_end) * problem_count/maxProblems;
        p_kb = p_kb_start - (p_kb_start - p_kb_end) * problem_count/maxProblems;
    }
    output_training_file.close();
    output_test_file.close();
}


void
doOneSingleStepTest(int training_problem_count, std::ofstream &output_test_file, double training_performance,
                    double training_error, std::string path_postfix) {
	bool wasCorrect = false;
    int correct_count = 0;
    double error_sum = 0;
    int match_set_sum = 0;
	DataSource *testState = NULL;
    std::string output_full_path = output_path + path_postfix;

    std::ofstream output_visualization_file;
    if(visualization) {
        output_visualization_file.open(output_full_path + output_visualization_file_name);
        if (!output_visualization_file.is_open()) {
            std::cout << "Could not open output visualization file";
            exit(1);
        }
    }

	for(int t=0; t<testNumInstances; t++){
        ClassifierSet match_set(maxPopSize);
        ClassifierSet action_set(maxPopSize);
		bool isMatched = false;
		testState = &testingData[t];

        get_matching_classifiers(testState->state, match_set, t, false);
        isMatched = match_set.ids.size() > 0;
        match_set_sum += match_set.ids.size();
        int actionWinner=-1;
        if(isMatched) {
            getPredictionArray(match_set);
            actionWinner = bestActionWinner();
            // save visualization data
            if(visualization){
                getActionSet(actionWinner, match_set, action_set);
                // save visualization data - start with image id, actual action and predicted action
                output_visualization_file << t << " " << testState->action << " " << actionWinner << std::endl; // image id
                save_visualization_data(action_set, t, output_visualization_file);
            }
        }else{
            // if match set is empty then select an action randomly
            actionWinner = irand(numActions);
        }
        double reward = executeAction(actionWinner, testState->action, wasCorrect);
        error_sum += absoluteValue(reward - getBestValue());

		if(wasCorrect){
		    correct_count += 1;
        }
    }
	if(visualization){
	    output_visualization_file.close();
	}
    double accuracy = correct_count/(double)testNumInstances;
    double error = error_sum / testNumInstances;
    double avg_action_set_size = (match_set_sum / (double)testNumInstances) / numActions;
    writeTestPerformance(accuracy, error, testNumInstances, avg_action_set_size, output_test_file, training_problem_count,
                         training_performance, training_error);
    output_test_file.flush();
}


void startXCS(){
    printf("\nLoading Input! Please wait ....\n");

    trainingData = new DataSource[trainNumInstances];
    testingData = new DataSource[testNumInstances];
    initializeInput(trainingData,trainNumInstances);
    initializeInput(testingData,testNumInstances);
    loadDataFromFile(trainingData, inputTrainingFile.c_str(), trainNumInstances);
    loadDataFromFile(testingData, inputTestFile.c_str(), testNumInstances);
    if(normalize_data) {
        normalize_image(trainingData, trainNumInstances);
        normalize_image(testingData, testNumInstances);
    }
    if(use_kb) {
        load_kb(kb_cf_file, kb_filter_file);
    }
    printf("\nIt is in progress! Please wait ....\n");

    doOneSingleStepExperiment();

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
            }else if(name == "classes"){
                std::replace(value.begin(), value.end(), ',', ' ');
                std::stringstream ss(value);
                for(int i=0; i<numActions; i++)
                {
                    int c=0;
                    ss>>c;
                    class_map[c] = i;
                }
            }else if(name == "image_width"){
                image_width = atoi(value.c_str());
            }else if(name == "image_height"){
                image_height = atoi(value.c_str());
            }else if(name == "max_depth"){
                cfMaxDepth = atoi(value.c_str());
                cfMaxLength = std::pow(2, cfMaxDepth+1);
                cfMaxStack = (cfMaxArity-1)*(cfMaxDepth-1) + 2;
                cfMaxLeaf = std::pow(2, cfMaxDepth);
            }else if(name == "min_depth"){
                cfMinDepth = atoi(value.c_str());
                cfMinLeaf = std::pow(2, cfMinDepth);
            }else if(name == "max_condition_length"){
                clfrCondMaxLength = atoi(value.c_str());
            }else if(name == "cf_max_bounding_box_size"){
                cf_max_bounding_box_size = atoi(value.c_str());
            }else if(name == "cf_min_bounding_box_size"){
                cf_min_bounding_box_size = atoi(value.c_str());
            }else if(name == "filter_matching_threshold"){
                filter_matching_threshold = atof(value.c_str());
            }else if(name == "pX_start"){
                pX_start = atof(value.c_str());
                pX = pX_start;
            }else if(name == "pX_end"){
                pX_end = atof(value.c_str());
            }else if(name == "pM_start"){
                pM_start = atof(value.c_str());
                pM = pM_start;
            }else if(name == "pM_end"){
                pM_end = atof(value.c_str());
            }else if(name == "p_kb_start"){
                p_kb_start = atof(value.c_str());
                p_kb = p_kb_start;
            }else if(name == "p_kb_end"){
                p_kb_end = atof(value.c_str());
            }else if(name == "epsilon_start"){
                epsilon_start = atof(value.c_str());
                epsilon = epsilon_start;
            }else if(name == "epsilon_end"){
                epsilon_end = atof(value.c_str());
            }else if(name == "beta_start"){
                beta_start = atof(value.c_str());
                beta = beta_start;
            }else if(name == "beta_end"){
                beta_end = atof(value.c_str());
            }else if(name == "pM_allel"){
                pM_allel = atof(value.c_str());
            }else if(name == "p_promising"){
                p_promising = atof(value.c_str());
            }else if(name == "normalize_data"){
                if(value == "no"){
                    normalize_data = false;
                }else if(value == "yes"){
                    normalize_data = true;
                }
            }else if(name == "use_kb"){
                if(value == "no"){
                    use_kb = false;
                }else if(value == "yes"){
                    use_kb = true;
                }
            }else if(name == "fixed_seed"){
                if(value == "no"){
                    fixed_seed = false;
                }else if(value == "yes"){
                    fixed_seed = true;
                }
            }else if(name == "visualization"){
                if(value == "no"){
                    visualization = false;
                }else if(value == "yes"){
                    visualization = true;
                }
            }else if(name == "kb_cf_file"){
                kb_cf_file = value;
            }else if(name == "kb_filter_file"){
                kb_filter_file = value;
            }else if(name == "output_path"){
                output_path = value + "/";
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
                   std::cout<<"Output path not set"<<std::endl;
                   exit(1);
                }
                mkdir(output_path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
                output_config_file.open(output_path + "experiment_config.txt", std::ofstream::app);
                if(!output_config_file.is_open()){
                    std::cout << "Could not open output config file";
                    exit(1);
                }
                std::time_t start_time = std::chrono::system_clock::to_time_t(start);
                output_config_file << std::endl << "Started at " << std::ctime(&start_time) << std::endl;
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
    condLength = image_width*image_height;
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

    // Time accounting
    start = std::chrono::system_clock::now();

    if(argc != 3){
        std::cout << "Please provide experiment config file and output path" << std::endl;
        return -1;
    }
    output_path = argv[2];
    output_path = output_path + "/";
    // check if the experiment has already completed
    std::ifstream input_done_file;
    input_done_file.open(output_path + output_done_file_name);
    if(input_done_file.is_open()){
        std::cout << "The experiment has already completed, exiting ...";
        input_done_file.close();
        exit(0);
    }

    LoadConfig(argv[1]);

    // standardized random number generator
    initialize_random_number_generator(fixed_seed);
    initialize_parameters();
    initialize_env();
    initialize_population(maxPopSize + 10);
    // initialize size of cf_list such that it can accommodate enough code fragments before periodic cleanup
    initialize_cf_list(clfrCondMaxLength*(maxPopSize + filter_list_management_frequency*numActions));


    // check if resume is needed by looking for validation_performance file
    std::ifstream validation_file(output_path + output_test_file_name);
    if(validation_file.is_open()){
        // get the last line from the file
        std::string line;
        std::string last_line;
        while(getline(validation_file, line)){
            if(!line.empty()) {
                last_line = line;
            }
        }
        if(!last_line.empty()) {
            std::stringstream line1(last_line);
            int last_iteration = 0;
            line1 >> last_iteration;
            resume_from = std::to_string(last_iteration) + "/";
        }
        validation_file.close();
    }

    for(int j=0; j<run; j++)
    {
        mkdir(output_path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        setSeed(seeds[j]);

        startXCS();
    }

    // Time accounting
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);


    // write done.txt
    std::ofstream output_done_file;
    output_done_file.open(output_path + output_done_file_name);
    if(!output_done_file.is_open()){
        std::cout << "Could not open done file";
        exit(1);
    }
    output_done_file << "The experiment has completed"<<std::endl;
    output_done_file << "Finished at " << std::ctime(&end_time) << std::endl;
    output_done_file << "Elapsed time since last execution (e.g. if experiment was resumed): " << elapsed_seconds.count()/3600 << " hours" << std::endl;
    output_done_file.close();
    std::cout << "Done!" << std::endl;
    exit(0);
}//end main
