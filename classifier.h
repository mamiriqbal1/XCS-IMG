void setInitialVariables(Classifier *clfr, double setSize, int time);
void initializePopulation(ClassifierSet **population, FILE *cfReadingFilePointer);//, FILE *cfWritingFilePointer);
int getNumerositySum(ClassifierSet *set);
int getSetSize(ClassifierSet *set);
double getAvgFitness(ClassifierSet *set);
int getNumFitterCFs(ClassifierSet *set, double avgFitness);

ClassifierSet* getMatchSet(ClassifierSet **population, ClassifierSet **killset, float state[], int itTime);
int nrActionsInSet(ClassifierSet *set, bool coveredActions[]);
bool isConditionMatched(CodeFragment clfrCond[], float state[]);
Classifier* matchingCondAndSpecifiedAct(float state[], int act, int setSize, int time);
void createMatchingCondition(CodeFragment cond[], float state[]);

void getPredictionArray(ClassifierSet *ms);
double getBestValue();
int randomActionWinner();
int bestActionWinner();
int rouletteActionWinner();

ClassifierSet* getActionSet(int action, ClassifierSet *ms);
void updateActionSet(ClassifierSet **aset, double maxPrediction, double reward, ClassifierSet **pop, ClassifierSet **killset);
void updateFitness(ClassifierSet *aset);

void discoveryComponent(ClassifierSet **set, ClassifierSet **pop, ClassifierSet **killset, int itTime, float situation[]);
void getDiscoversSums(ClassifierSet *set, double *fitsum, int *setsum, int *gaitsum);
void setTimeStamps(ClassifierSet *set, int itTime);

void selectTwoClassifiers(Classifier **cl, Classifier **parents, ClassifierSet *set, double fitsum, int setsum);
Classifier* selectClassifierUsingTournamentSelection(ClassifierSet *set, int setsum, Classifier *notMe);
Classifier* selectClassifierUsingRWS(ClassifierSet *set, double fitsum);

void subsumeCFs(Classifier *clfr, float state[]);

bool crossover(Classifier **cl, int crossoverType);
void uniformCrossover(Classifier **cl);
void onePointCrossover(Classifier **cl);
void twoPointCrossover(Classifier **cl);
bool mutation(Classifier *clfr, float state[]);
bool applyNicheMutation(Classifier *clfr, float state[]);
bool applyGeneralMutation(Classifier *clfr, float state[]);
bool mutateAction(Classifier *clfr);

void insertDiscoveredClassifier(Classifier **cl, Classifier **parents, ClassifierSet **set, ClassifierSet **pop, ClassifierSet **killset, int len, float state[]);

void doActionSetSubsumption(ClassifierSet **aset, ClassifierSet **pop, ClassifierSet **killset);
void subsumeClassifier(Classifier *cl, Classifier **parents, ClassifierSet *locset, ClassifierSet **pop,  float state[]);
bool subsumeClassifierToSet(Classifier *cl, ClassifierSet *set);
bool subsumes(Classifier *cl1, Classifier * cl2);
bool isSubsumer(Classifier *cl);
bool isMoreGeneral(Classifier *clfr1, Classifier *clfr2);
bool compareDontcares(Classifier *clfr1, Classifier *clfr2);
bool checkNonDontcares(CodeFragment cond1[], CodeFragment cond2[]);

bool addClassifierToPointerSet(Classifier *cl, ClassifierSet **pointerset);
bool addClassifierToSet(Classifier *cl, ClassifierSet **clSet);
void addNewClassifierToSet(Classifier *cl, ClassifierSet **clSet);
bool equals(Classifier *clfr1, Classifier *clfr2);

Classifier* deleteStochClassifier(ClassifierSet **pop);
double getDelProp(Classifier *clfr, double meanFitness);
Classifier* deleteTypeOfClassifier(ClassifierSet *setp, ClassifierSet *setpl, ClassifierSet **pop);
bool updateSet(ClassifierSet **uset, ClassifierSet *killset);
bool deleteClassifierPointerFromSet(ClassifierSet **set, Classifier *clp);

void freeSet(ClassifierSet **cls);
void freeClassifierSet(ClassifierSet **cls);
void freeClassifier(Classifier *cl);

void printClassifierSet(ClassifierSet *head);
void fprintClassifierSet(FILE *fpClfr, FILE *fpCF, ClassifierSet *head);
void printClassifier(Classifier *clfr);
void fprintClassifier(FILE *fp, Classifier *classifier);

ClassifierSet* sortClassifierSet(ClassifierSet **cls, int type);
void simplifyPopulation(ClassifierSet **population);

double absoluteValue(double value);
float computeDistance(CodeFragment clfrCond[], float cond[]);
