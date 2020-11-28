void setInitialVariables(Classifier &clfr, double setSize, int time);

void getMatchSet(ClassifierVector &pop, ClassifierSet &match_set, float *state, int itTime, int action, int img_id);
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

void discoveryComponent(ClassifierSet &action_set, ClassifierVector &pop, int itTime, float *situation);
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
bool subsumeClassifierToPop(Classifier &cl, ClassifierVector &cl_set);
bool subsumeClassifierToSet(Classifier &cl, ClassifierSet &cl_set);
bool subsumes(Classifier &cl1, Classifier & cl2);
bool isSubsumer(Classifier &cl);
bool isMoreGeneral(Classifier &clfr1, Classifier &clfr2);

int deleteStochClassifier(ClassifierVector &pop);
double getDelProp(Classifier &clfr, double meanFitness);

void save_experiment_results(ClassifierVector &pop, std::string path_postfix);

void fprintClassifier(Classifier &classifier, std::ofstream &output_classifier_file,
                      std::ofstream &output_code_fragment_file, std::ofstream &output_promising_code_fragment_file,
                      std::ofstream &output_promising_filter_file);

double absoluteValue(double value);

void manage_filter_list(ClassifierVector &pop);
int get_pop_size(ClassifierVector &pop, bool numerosity);
void get_matching_classifiers(ClassifierVector& pop, float *state, ClassifierSet &match_set, int img_id, bool train);
bool is_promising_classifier(Classifier& cl);
void transfer_kb_filter(CodeFragment & cf);
void remove_classifier_cfs_from_list(Classifier &cl);
void add_classifier_cfs_to_list(Classifier &cl);
