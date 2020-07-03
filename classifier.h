void setInitialVariables(Classifier *clfr, double setSize, int time);
void initializePopulation(ClassifierSet **population, FILE *cfReadingFilePointer);//, FILE *cfWritingFilePointer);
double getAvgFitness(ClassifierSet *set);
int getNumFitterCFs(ClassifierSet *set, double avgFitness);

void *
getMatchSet(ClassifierList &pop, ClassifierList &match_set, float *state, int itTime, int action, int img_id);
int nrActionsInSet(ClassifierList &set, bool *coveredActions);
//bool isConditionMatched(CodeFragment clfrCond[], float state[], int cl_id=-1, int img_id=-1);
bool isConditionMatched(Classifier &cl, float state[], int img_id, bool train);
Classifier* matchingCondAndSpecifiedAct(float state[], int act, int setSize, int time);
void createMatchingCondition(CodeFragment cond[], float state[]);

void getPredictionArray(ClassifierList &match_set);
double getBestValue();
int randomActionWinner();
int bestActionWinner();

void *getActionSet(int action, ClassifierList &match_set, ClassifierList &action_set);
void updateActionSet(ClassifierList &action_set, double maxPrediction, double reward, ClassifierList &pop);
void updateFitness(ClassifierList &action_set);

void discoveryComponent(ClassifierList &action_set, ClassifierList &pop, int itTime, float situation[]);
void getDiscoversSums(ClassifierList action_set, double *fitsum, int *setsum, int *gaitsum);
void setTimeStamps(ClassifierList action_set, int itTime);

void selectTwoClassifiers(Classifier cl[], Classifier parents[], ClassifierList &action_set, double fitsum, int setsum);

bool crossover(Classifier &cl1, Classifier &cl2, float situation[]);

bool mutation(Classifier &clfr, float *state);

bool applyNicheMutation1(Classifier *clfr, float state[]);
bool applyNicheMutation2(Classifier *clfr, float state[]);
bool applyGeneralMutation(Classifier *clfr, float state[]);
bool mutateAction(Classifier& clfr);

void insertDiscoveredClassifier(Classifier cl[], Classifier parents[], ClassifierList &action_set, ClassifierList &pop,
                                int len, float *state);

void doActionSetSubsumption(ClassifierList &action_set, ClassifierList &pop);
void subsumeClassifier(Classifier& cl, Classifier parents[], ClassifierList &action_set, ClassifierList &pop, float *state);
bool subsumeClassifierToSet(Classifier &cl, ClassifierList &cl_set, ClassifierSet *set);
bool subsumes(Classifier &cl1, Classifier & cl2);
bool isSubsumer(Classifier &cl);
bool isMoreGeneral(Classifier &clfr1, Classifier &clfr2);

int deleteStochClassifier(ClassifierList &pop);
double getDelProp(Classifier *clfr, double meanFitness);

void freeSet(ClassifierSet **cls);
void freeClassifierSet(ClassifierSet **cls);
void freeClassifier(Classifier *cl);

void printClassifierSet(ClassifierSet *head);
void fprintClassifierSet(FILE *fpClfr, FILE *fpCF, ClassifierSet *head);
void printClassifier(Classifier *clfr);
void fprintClassifier(FILE *fp, Classifier *classifier);

ClassifierSet* sortClassifierSet(ClassifierSet **cls, int type);

double absoluteValue(double value);
float computeDistance(CodeFragment clfrCond[], float cond[]);


void manage_filter_list(ClassifierList &pop);

int get_set_size(ClassifierList& pop);
int get_set_numerosity(ClassifierList& pop);
void get_matching_classifiers(ClassifierList& pop, float state[], ClassifierList& match_set, int img_id, bool train);
