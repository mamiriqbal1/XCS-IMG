void setInitialVariables(Classifier &clfr, double setSize, int time);
void initializePopulation(delete_ClassifierSet **population, FILE *cfReadingFilePointer);//, FILE *cfWritingFilePointer);
double getAvgFitness(delete_ClassifierSet *set);
int getNumFitterCFs(delete_ClassifierSet *set, double avgFitness);

void *
getMatchSet(ClassifierMap &pop, ClassifierSet &match_set, float *state, int itTime, int action, int img_id);
int nrActionsInSet(ClassifierSet &match_set, bool *coveredActions);
//bool isConditionMatched(CodeFragment clfrCond[], float state[], int cl_id=-1, int img_id=-1);
bool isConditionMatched(Classifier &cl, float state[], int img_id, bool train);
void matchingCondAndSpecifiedAct(Classifier &cl, float *state, int act, int setSize, int time);
void createMatchingCondition(CodeFragment cond[], float state[]);

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

void selectTwoClassifiers(Classifier *cl, int *parents, ClassifierSet &action_set, double fitsum, int setsum);

bool crossover(Classifier &cl1, Classifier &cl2, float situation[]);

bool mutation(Classifier &clfr, float *state);

bool applyNicheMutation1(Classifier *clfr, float state[]);
bool applyNicheMutation2(Classifier *clfr, float state[]);
bool applyGeneralMutation(Classifier *clfr, float state[]);
bool mutateAction(Classifier& clfr);

void insertDiscoveredClassifier(Classifier *child, int *parent, ClassifierMap &pop, int len);

void doActionSetSubsumption(ClassifierSet &action_set);
bool subsumeClassifier(Classifier &cl, Classifier &p1, Classifier &p2);
bool subsumeClassifierToPop(Classifier &cl, ClassifierMap &cl_set);
bool subsumes(Classifier &cl1, Classifier & cl2);
bool isSubsumer(Classifier &cl);
bool isMoreGeneral(Classifier &clfr1, Classifier &clfr2);

int deleteStochClassifier(ClassifierMap &pop);
double getDelProp(Classifier &clfr, double meanFitness);

void freeSet(delete_ClassifierSet **cls);
void freeClassifierSet(delete_ClassifierSet **cls);
void freeClassifier(Classifier *cl);

void printClassifierSet(delete_ClassifierSet *head);
void fprintClassifierSet(FILE *fpClfr, FILE *fpCF, delete_ClassifierSet *head);
void printClassifier(Classifier *clfr);
void fprintClassifier(FILE *fp, Classifier *classifier);

delete_ClassifierSet* sortClassifierSet(delete_ClassifierSet **cls, int type);

double absoluteValue(double value);
float computeDistance(CodeFragment clfrCond[], float cond[]);


void manage_filter_list(ClassifierMap &pop);

int get_set_size(ClassifierMap& pop);
int get_pop_numerosity(ClassifierMap& pop);
void get_matching_classifiers(ClassifierMap& pop, float *state, ClassifierSet &match_set, int img_id, bool train);
