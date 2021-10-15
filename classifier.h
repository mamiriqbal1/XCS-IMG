#ifndef XCS_IMG_CL_H
#define XCS_IMG_CL_H

void setInitialVariables(Classifier &clfr, double setSize, int time);

void getMatchSet(ClassifierSet &match_set, float *state, int itTime, int action, int img_id);
int nrActionsInSet(ClassifierSet &match_set, bool *coveredActions);
bool isConditionMatched(Classifier &cl, float state[], int img_id, bool train);
void matchingCondAndSpecifiedAct(Classifier &cl, float *state, int act, int setSize, int time);
void createMatchingCondition(Classifier &cl, float *state);

void getPredictionArray(ClassifierSet &match_set);
double getBestValue();
int randomActionWinner();
int bestActionWinner();

void getActionSet(int action, ClassifierSet &match_set, ClassifierSet &action_set);
void updateActionSet(ClassifierSet &action_set, double maxPrediction, double reward);
void updateFitness(ClassifierSet &action_set);

void discoveryComponent(ClassifierSet &action_set, int itTime, float *situation, int action);
void getDiscoversSums(ClassifierSet &action_set, double *fitsum, int *setsum, int *gaitsum);
void setTimeStamps(ClassifierSet &action_set, int itTime);

void selectTwoClassifiers(Classifier &child1, Classifier &child2, int &parent1, int &parent2, ClassifierSet &action_set,
                          double fitsum, int setsum);

bool crossover(Classifier &cl1, Classifier &cl2, float *state);

bool mutation(Classifier &clfr, float *state);

bool mutateAction(Classifier& clfr);

void insertDiscoveredClassifier(Classifier *child, int *parent, ClassifierSet &action_set, int len);

void doActionSetSubsumption(ClassifierSet &action_set);
bool subsumeClassifier(Classifier &cl, Classifier &p1, Classifier &p2, ClassifierSet &action_set);
bool subsumeClassifierToPop(Classifier &cl);
bool subsumeClassifierToSet(Classifier &cl, ClassifierSet &cl_set);
bool subsumes(Classifier &cl1, Classifier & cl2);
bool isSubsumer(Classifier &cl);
bool isMoreGeneral(Classifier &clfr1, Classifier &clfr2);

int deleteStochClassifier(ClassifierVector &pop);
double getDelProp(Classifier &clfr, double meanFitness);

void save_experiment_results(std::string path_postfix);

void fprintClassifier(Classifier &classifier, std::ofstream &output_classifier_file);

double absoluteValue(double value);

void manage_filter_and_cf_list();
int get_pop_size(bool numerosity);
void get_matching_classifiers(float *state, ClassifierSet &match_set, int img_id, bool train);
bool is_promising_classifier(Classifier& cl);

void remove_classifier_cfs_from_list(Classifier &cl);
void add_classifier_cfs_to_list(Classifier &cl);
int count_classifier_cfs(const Classifier &cl);
void load_classifier(std::string classifier_file_name);
void initialize_population(int size);
void add_new_classifiers_to_population(float* state, int action, int itTime);
void write_classifier_header(std::ofstream &output_classifier_file);

#endif //XCS_IMG_CL_H
