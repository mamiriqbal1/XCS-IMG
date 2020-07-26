void setInitialVariables(Classifier &clfr, double setSize, int time);
double getAvgFitness(ClassifierMap &pop);
int getNumFitterCFs(ClassifierMap &pop, double avgFitness);

void *
getMatchSet(ClassifierMap &pop, ClassifierSet &match_set, float *state, int itTime, int action, int img_id);
int nrActionsInSet(ClassifierSet &match_set, bool *coveredActions);
bool isConditionMatched(Classifier &cl, float state[], int img_id, bool train);
void matchingCondAndSpecifiedAct(Classifier &cl, float *state, int act, int setSize, int time);
void createMatchingCondition(CodeFragment code_fragment[], float state[]);

void getPredictionArray(ClassifierSet &match_set);
double getBestValue();
int randomActionWinner();
int bestActionWinner();

void getActionSet(int action, ClassifierSet &match_set, ClassifierSet &action_set);
void updateActionSet(ClassifierSet &action_set, double maxPrediction, double reward);
void updateFitness(ClassifierSet &action_set);

void discoveryComponent(ClassifierSet &action_set, ClassifierMap &pop, int itTime, float *situation);
void getDiscoversSums(ClassifierSet &action_set, double *fitsum, int *setsum, int *gaitsum);
void setTimeStamps(ClassifierSet &action_set, int itTime);

void selectTwoClassifiers(Classifier &child1, Classifier &child2, int &parent1, int &parent2, ClassifierSet &action_set,
                          double fitsum, int setsum);

bool crossover(Classifier &cl1, Classifier &cl2, float situation[]);

bool mutation(Classifier &clfr, float *state);

bool mutateAction(Classifier& clfr);

void insertDiscoveredClassifier(Classifier *child, int *parent, ClassifierMap &pop, int len);

void doActionSetSubsumption(ClassifierSet &action_set);
bool subsumeClassifier(Classifier &cl, Classifier &p1, Classifier &p2, ClassifierMap &pop);
bool subsumeClassifierToPop(Classifier &cl, ClassifierMap &cl_set);
bool subsumes(Classifier &cl1, Classifier & cl2);
bool isSubsumer(Classifier &cl);
bool isMoreGeneral(Classifier &clfr1, Classifier &clfr2);

int deleteStochClassifier(ClassifierMap &pop);
double getDelProp(Classifier &clfr, double meanFitness);

void save_experiment_results(ClassifierMap &pop);

void fprintClassifier(Classifier &classifier, std::ofstream &output_classifier_file,
                      std::ofstream &output_code_fragment_file, std::ofstream &output_promising_code_fragment_file,
                      std::ofstream &output_promising_filter_file);

double absoluteValue(double value);
float computeDistance(CodeFragment code_fragment[], float cond[]);


void manage_filter_list(ClassifierMap &pop);

int get_pop_numerosity(ClassifierMap& pop);
void get_matching_classifiers(ClassifierMap& pop, float *state, ClassifierSet &match_set, int img_id, bool train);
bool is_promising_classifier(Classifier& cl);
