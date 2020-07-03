#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include <assert.h>
#include <fstream>
#include <vector>
#include <algorithm>
#include "xcsMacros.h"
#include "codeFragment.h"
#include "classifier.h"
#include "env.h"
#include "filter_list.h"

int countNewCFs = 0;
double predictionArray[max_actions]; //prediction array
double sumClfrFitnessInPredictionArray[max_actions]; //The sum of the fitnesses of classifiers that represent each entry in the prediction array.
int gid=0; // global incremental id to  uniquely identify the classifiers for evaluation reuse

void setInitialVariables(Classifier *clfr, double setSize, int time){
    clfr->id = gid++;
    clfr->prediction = predictionIni;
    clfr->predictionError = predictionErrorIni;
    clfr->accuracy = 0.0;
    clfr->fitness = fitnessIni;
    clfr->numerosity = 1;
    clfr->experience = 0;
    clfr->actionSetSize = 1; // chnged to 1 as per paper instead of setSize argument
    clfr->timeStamp = time;
    clfr->specificness = numberOfNonDontcares(clfr->condition);
}

void initializePopulation(ClassifierSet **population, FILE *cfReadingFilePointer){//, FILE *cfWritingFilePointer)

    *population = NULL;

    initializeCFPopulation(cfReadingFilePointer);//,cfWritingFilePointer);
}




int getNumerositySum(ClassifierSet *set){
    int sum = 0;
    for(; set!=NULL; set=set->next)
    {
        sum += set->classifier->numerosity;
    }
    return sum;
}

int getSetSize(ClassifierSet *set){
    int size = 0;
    for(; set!=NULL; set=set->next)
    {
        size++;
    }
    return size;
}

double getAvgFitness(ClassifierSet *set){
    ClassifierSet *setp;
    double fitsum=0.0;
    int setsum=0;

    for(setp=set; setp!=NULL; setp=setp->next)
    {
        fitsum += setp->classifier->fitness;
        setsum += setp->classifier->numerosity;
    }
    return fitsum/(double)setsum;
}

int getNumFitterCFs(ClassifierSet *set, double avgFitness){
    ClassifierSet *setp = set;
    int numCFs = getNumPreviousCFs();
    while(setp!=NULL && setp->classifier->fitness > avgFitness)
    {
        numCFs += setp->classifier->specificness;
        setp = setp->next;
    }
    return numCFs;
}


void remove_classifier(ClassifierList &set, int cl_id) {
    auto previous_it = set.before_begin();
    for(auto it =  set.begin(); it != set.end(); it++){
        if(it->id == cl_id){
            break;
        }else{
            previous_it = it;
        }
    }
    set.erase_after(previous_it);
}

// ####################### match set operations ##########################################

/**
 * Gets the match-set that matches state from pop.
 * If a classifier was deleted, record its address in killset to be
 * able to update former actionsets.
 * The iteration time 'itTime' is used when creating a new classifier
 * due to covering. Covering occurs when not all possible actions are
 * present in the match set. Thus, it is made sure that all actions
 * are present in the match set.
 */
 /*
  * Code review
  * Overall flow is as per algorithm
  */


void *
getMatchSet(ClassifierList &pop, ClassifierList &match_set, float *state, int itTime, int action, int img_id) {
    Classifier *coverClfr= nullptr;
    int popSize=0, pop_numerosity=0, representedActions;

    bool coveredActions[numActions];
    popSize = get_set_size(pop);
    get_matching_classifiers(pop, state, match_set, img_id, true);
    pop_numerosity = get_set_numerosity(match_set);

    representedActions = nrActionsInSet(match_set,coveredActions);

    // TEMP: insert filters every time
    if(popSize < maxPopSize/2) {
        representedActions = 0;
        for(int i=0; i<numActions; i++){
            coveredActions[i] = false;
        }
    }

    while(representedActions < numActions){  // create covering classifiers, if not all actions are covered
        for(int i=0; i<numActions; i++){
            if(!coveredActions[i]){  // make sure that all actions are covered!
                // TEMP: boost covering
                // add large number of classifier in case of filter approach in covering
                int add = 1;
                for(int j=0; j<add; j++){
                    coverClfr = matchingCondAndSpecifiedAct(state, i, pop_numerosity + 1, itTime);
                    // before inserting the new classifier into the population check for subsumption by a generic one
                    if(!subsumeClassifierToSet(*coverClfr, pop, nullptr)) {
                        pop.push_front(*coverClfr);
                    }
                    if(!subsumeClassifierToSet(*coverClfr, match_set, nullptr)) {
                        match_set.push_front(*coverClfr);
                    }
                    free(coverClfr); // this memory was allocated by matchingConditionAndSpecifidAct
                    coverClfr= nullptr;
                    pop_numerosity++;
                    popSize++;
                }
            }
        }

        /* Delete classifier if population is too big and record it in killset */
        while( popSize > maxPopSize )
        {
            /* PL */
            int cl_id = deleteStochClassifier(pop);
            // also remove from match set
            remove_classifier(match_set, cl_id);
            popSize--;
        }
        representedActions = nrActionsInSet(match_set,coveredActions);
    }
}
/**
 * Returns the number of actions in the set and stores which actions are covered in the array coveredActions.
 */
int nrActionsInSet(ClassifierList &set, bool *coveredActions)
{
    int nr=0;
    for(int i=0; i<numActions; i++){
        coveredActions[i] = false;
    }
    for(auto it=set.begin(); it != set.end() && nr < numActions; it++){
        if(!coveredActions[it->action]){
            coveredActions[it->action] = true;
            nr++;
        }
    }
    return nr;
}
float computeDistance(CodeFragment clfrCond[], float cond[]){

  float mid, sum = 0.0;
  for(int i=0; i<clfrCondLength; i++){
            if( isDontcareCF(clfrCond[i]) || evaluateCF(clfrCond[i],cond)==1 )
            {
                sum += 1;
            }
    }
  return sum;
}


bool isConditionMatched(Classifier &cl, float state[], int img_id, bool train)
{
    for(int i=0; i<clfrCondLength; i++)
    {
        //std::cout<<"iscond\n";
        //if( !isDontcareCF(clfrCond[i]) && evaluateCF(clfrCond[i].codeFragment,state)==0 )
        if(!isDontcareCF(cl.condition[i]) && evaluateCF(cl.condition[i], state, cl.id, img_id, train) == 0 )
        {
            return false;
        }
    }
    return true;
}
Classifier* matchingCondAndSpecifiedAct(float state[], int act, int setSize, int time)  //matchingCondAndSpecifiedAct(currentState,i,matchSetNumerositySum+1,time)
{
    Classifier *clfr;
    assert((clfr=( struct Classifier*)calloc(1,sizeof(struct Classifier)))!=NULL); // get memory for the new classifier
    createMatchingCondition(clfr->condition,state);
    clfr->action = act;
    setInitialVariables(clfr,setSize,time);
    return clfr;
}
void createMatchingCondition(CodeFragment cond[], float state[])
{
//    int coun = 0;
    //printf("in.......\n");
    for(int i=0; i<clfrCondLength; i++)
    {
       CodeFragment tempCF;
        do
        {
            //ramped half-and-half, check end for sanity only

            //std::cout<<getNumPreviousCFs();
            tempCF = createNewCF(getNumPreviousCFs()+countNewCFs);
            opType* end = randomProgram(tempCF.codeFragment,irand(2),cfMaxDepth,cfMinDepth);
            validateDepth(tempCF.codeFragment,end); //validate depth
            /////////////
            tempCF = addLeafCF(tempCF, state);
            //std::cout<<"No of eval: "<<coun++<<"\n";
            //printCF(tempCF);

            /////////////
        }
        while( evaluateCF(tempCF,state)!=1 );
//            coun  = 0;
        //while( evaluateCF(tempCF.codeFragment,state)!=1 );
        memmove(&cond[i],&tempCF,sizeof(CodeFragment));
        countNewCFs++;
    }
    //printf("out.......\n");
}

// ######################### prediction array operations ############################################

void getPredictionArray(ClassifierList &match_set)  //determines the prediction array out of the match set ms
{
    for(int i=0; i<numActions; i++)
    {
        predictionArray[i]=0.0;
        sumClfrFitnessInPredictionArray[i]=0.0;
    }
    for(auto it : match_set)
    {
        int actionValue = it.action;
        predictionArray[actionValue]+= it.prediction * it.fitness;
        sumClfrFitnessInPredictionArray[actionValue]+= it.fitness;
    }
    for(int i=0; i<numActions; i++)
    {
        if(sumClfrFitnessInPredictionArray[i]!=0)
        {
            predictionArray[i] /= sumClfrFitnessInPredictionArray[i];
        }
        else
        {
            predictionArray[i]=0;
        }
    }
}
double getBestValue()  //Returns the highest value in the prediction array.
{
    double max = predictionArray[0];
    for(int i=1; i<numActions; i++)
    {
        if(max<predictionArray[i])
        {
            max = predictionArray[i];
        }
    }
    return max;
}
int randomActionWinner()  //Selects an action randomly. The function assures that the chosen action is represented by at least one classifier in the prediction array.
{
    int ret=0;
    do
    {
        ret = irand(numActions);
    }
    while(sumClfrFitnessInPredictionArray[ret]==0);
    return ret;
}
int bestActionWinner()  //Selects the action in the prediction array with the best value.
{
    int ret=irand(numActions);
    for(int i=0; i<numActions; i++)
    {
        if(predictionArray[ret]<predictionArray[i])
        {
            ret=i;
        }
    }
    return ret;
}
int rouletteActionWinner()  //Selects an action in the prediction array by roulette wheel selection.
{
    double bidSum=0.0;
    int i;
    for(i=0; i<numActions; i++)
    {
        bidSum += predictionArray[i];
    }
    bidSum *= drand();
    double bidC=0.0;
    for(i=0; bidC<bidSum; i++)
    {
        bidC += predictionArray[i];
    }
    return i;
}

// ######################## action set operations #########################################

void *getActionSet(int action, ClassifierList &match_set,
                   ClassifierList &action_set)  // constructs an action set out of the match set ms.
{
    for(auto it : match_set){
        if(action == it.action){
            action_set.push_front(it);
        }
    }
}

/**
 * Updates all parameters in the action set.
 * Essentially, reinforcement Learning as well as the fitness evaluation takes place in this set.
 * Moreover, the prediction error and the action set size estimate is updated. Also,
 * action set subsumption takes place if selected. As in the algorithmic description, the fitness is updated
 * after prediction and prediction error. However, in order to be more conservative the prediction error is
 * updated before the prediction.
 * @param maxPrediction The maximum prediction value in the successive prediction array (should be set to zero in single step environments).
 * @param reward The actual resulting reward after the execution of an action.
 */
 /*
  * code review notes
  * The formulas are implemented in a slightly differnt form with same end result. Can be simplified as per the paper.
  */
 void updateActionSet(ClassifierList &action_set, double maxPrediction, double reward, ClassifierList &pop)
{
    double P, setsize=0.0;
    ClassifierSet *setp;

    P = reward + gama*maxPrediction;

    for(auto it : action_set)
    {
        setsize += it.numerosity;
        it.experience++;
    }

    for(auto it : action_set)   // update prediction, prediction error and action set size estimate
    {
        if((double)it.experience < 1.0/beta)
        {
            // !first adjustments! -> simply calculate the average
            it.predictionError = (it.predictionError * ((double)it.experience - 1.0) + absoluteValue(P - it.prediction)) / (double)it.experience;
            it.prediction = (it.prediction * ((double)it.experience - 1.0) + P) / (double)it.experience;
            it.actionSetSize = (it.actionSetSize *((double)(it.experience - 1))+setsize)/(double)it.experience;
        }
        else
        {
            // normal adjustment -> use widrow hoff delta rule
            it.predictionError += beta * (absoluteValue(P - it.prediction) - it.predictionError);
            it.prediction += beta * (P - it.prediction);
            it.actionSetSize += beta * (setsize - it.actionSetSize);
        }
    }
    updateFitness(action_set);
    if(doActSetSubsumption)
    {
        doActionSetSubsumption(action_set, pop);
    }
}

// update the fitnesses of an action set (the previous [A] in multi-step envs or the current [A] in single-step envs.)
/*
 * code review notes
 * correctly implementation
 */
void updateFitness(ClassifierList &action_set)
{
    ClassifierSet *setp;
    double ksum=0.0;


    //First, calculate the accuracies of the classifier and the accuracy sums
    for(auto it : action_set){
        if(it.predictionError <= epsilon_0){
            it.accuracy = 1.0;
        }
        else{
            it.accuracy = alpha * pow(it.predictionError / epsilon_0 , -nu);
        }
        ksum += it.accuracy*(double)it.numerosity;
    }

    //Next, update the fitnesses accordingly
    for(auto it : action_set){
        it.fitness += beta * ( (it.accuracy * it.numerosity) / ksum - it.fitness );
    }
}

// ############################ discovery mechanism #########################################

/**
 * The discovery conmponent with the genetic algorithm
 * note: some classifiers in set could be deleted !
 */
void discoveryComponent(ClassifierList &action_set, ClassifierList &pop, int itTime, float situation[])
{
    ClassifierSet *setp;
    Classifier cl[2], parents[2];
    double fitsum=0.0;
    int i, len, setsum=0, gaitsum=0;
    // if the classifier set is empty, return (due to deletion)
    if(std::distance(action_set.before_begin(), action_set.end()) == 0) return;

    getDiscoversSums(action_set, &fitsum, &setsum, &gaitsum); // get all sums that are needed to do the discovery

    // do not do a GA if the average number of time-steps in the set since the last GA is less or equal than thetaGA
    if( itTime - (double)gaitsum / (double)setsum < theta_GA)
    {
        return;
    }
    setTimeStamps(action_set, itTime); // todo how this timestamp is reflected in population???

    selectTwoClassifiers(cl, parents, action_set, fitsum, setsum); // select two classifiers (tournament selection) and copy them
    // Prediction, prediction error and fitness is only updated if crossover is done instead of always
    // (this is reverted because of slightly decreased performance)
    if(crossover(cl[0], cl[1], situation) || true){
        cl[0].prediction   = (cl[0].prediction + cl[1].prediction) / 2.0;
        cl[0].predictionError = predictionErrorReduction * ( (cl[0].predictionError + cl[1].predictionError) / 2.0 );
        cl[0].fitness = fitnessReduction * ( (cl[0].fitness + cl[1].fitness) / 2.0 );

        cl[1].prediction = cl[0].prediction;
        cl[1].predictionError = cl[0].predictionError;
        cl[1].fitness = cl[0].fitness;
    }
    for(i=0; i<2; i++)  // do mutation
    {
        mutation(cl[i], situation);
    }

    cl[0].specificness = numberOfNonDontcares(cl[0].condition);
    cl[1].specificness = numberOfNonDontcares(cl[1].condition);

    // get the length of the population to check if clasifiers have to be deleted
    len = get_set_numerosity(pop);

    // insert the new two classifiers and delete two if necessary
    insertDiscoveredClassifier(cl, parents, action_set, pop, len, situation);
}
void getDiscoversSums(ClassifierList action_set, double *fitsum, int *setsum, int *gaitsum)  // Calculate all necessary sums in the set for the discovery component.
{
    *fitsum=0.0;
    *setsum=0;
    *gaitsum=0;
    for(auto setp : action_set)
    {
        (*fitsum)+=setp.fitness;
        (*setsum)+=setp.numerosity;
        (*gaitsum) += setp.timeStamp*setp.numerosity;
    }
}
void setTimeStamps(ClassifierList action_set, int itTime)  // Sets the time steps of all classifiers in the set to itTime (because a GA application is occurring in this set!).
{
    for(auto it : action_set)
    {
        it.timeStamp = itTime;
    }
}

void tournament_selection(Classifier *child, Classifier* parent, ClassifierList& set, double setsum)
{
    int first = irand(setsum);
    int second = irand(setsum);
    Classifier* first_classifier = nullptr, *second_classifier = nullptr;
    int i=0;
    for(auto it : set){
        if(i == first){
            first_classifier = &it;
        }
        if(i == second){
            second_classifier = &it;
        }
        i++;
    }
    assert(first_classifier != nullptr);
    assert(second_classifier != nullptr);
    if(first_classifier->fitness > second_classifier->fitness){
        *child = *parent = *first_classifier;
    }else if(first_classifier->fitness < second_classifier->fitness){
        *child = *parent = *second_classifier;
    }else{
        int r = irand(2);
        if(r == 0){
            *child = *parent = *first_classifier;
        }else{
            *child = *parent = *second_classifier;
        }
    }
}

// ########################### selection mechanism ########################################

/**
 * Select two classifiers using the chosen selection mechanism and copy them as offspring.
 */
void selectTwoClassifiers(Classifier cl[], Classifier parents[], ClassifierList &action_set, double fitsum, int setsum)
{
    tournament_selection(&cl[0], &parents[0], action_set, setsum);
    tournament_selection(&cl[1], &parents[1], action_set, setsum);

    for(int i=0; i<2; i++) {
        cl[i].id = gid++;
        cl[i].numerosity = 1;
        cl[i].experience = 0;
    }
}

/**
 * Selects a classifier from 'set' using tournament selection.
 * If 'notMe' is not the NULL pointer and forceDifferentInTournament is set to a value larger 0,
 * this classifier is not selected except if it is the only classifier.
 */
Classifier* selectClassifierUsingTournamentSelection(ClassifierList action_set, int setsum, Classifier *notMe)
{
    ClassifierSet *setp, *winnerSet=NULL;
    Classifier *winner=NULL;
    double fitness=-1.0, value;
    int i, j, *sel, size=0;

    if(notMe!=0)
    {
        if(drand() < forceDifferentInTournament)
        {
            setsum -= notMe->numerosity;
        }
        else
        {
            notMe=0; /* picking the same guy is allowed */
        }
    }
    if(setsum<=0)  /* only one classifier in set */
    {
        return &*action_set.begin();
    }

    if(tournamentSize>1) /* tournament with fixed size */
    {
        assert((sel = (int *)calloc(setsum, sizeof(int)))!=NULL);
        for(i=0; i<tournamentSize; i++)   /* (with replacement) */
        {
            sel[irand(setsum)]=1;
        }
        if(i-tournamentSize != 0 && drand() > i-tournamentSize)
        {
            /* possible probabilistic selection of the last guy */
            sel[irand(setsum)]=1;
        }
        for(auto setp : action_set)
        {
            if(&setp != notMe)
            {
                if(fitness < setp.fitness/setp.numerosity)
                {
                    for(j=0; j<setp.numerosity; j++)
                    {
                        if(sel[i+j])
                        {
                            freeSet(&winnerSet);
                            addClassifierToPointerSet(&setp, &winnerSet);
                            fitness = setp.fitness/setp.numerosity;
                            break; /* go to next classifier since this one is already a winner*/
                        }
                    }
                }
                i += setp.numerosity;
            }
        }
        free(sel);
        assert(winnerSet!=NULL);
        size=1;
    }
    else /* tournament selection with the tournament size approx. equal to tournamentSize*setsum */
    {
        winnerSet=NULL;
        while(winnerSet==NULL)
        {
            size=0;
            for(auto setp : action_set)
            {
                if(&setp != notMe)   /* do not reselect the same classifier -> this only applies if forcedDifferentInTournament is set!*/
                {
                    value = setp.predictionError;
                    if(winnerSet==NULL || (!doGAErrorBasedSelect && fitness - selectTolerance <= setp.fitness/setp.numerosity) || (doGAErrorBasedSelect && fitness + selectTolerance * maxPayoff >= value))
                    {
                        /* if his fitness is worse then do not bother */
                        for(i=0; i<setp.numerosity; i++)
                        {
                            if(drand() < tournamentSize)
                            {
                                /* this classifier is a selection candidate and
                                * his fitness/error is higher/lower or similar to the other classifier */
                                if(winnerSet==NULL)
                                {
                                    /* the first guy in the tournament */
                                    addClassifierToPointerSet(&setp, &winnerSet);
                                    if(doGAErrorBasedSelect)
                                    {
                                        fitness = value;
                                    }
                                    else
                                    {
                                        fitness = setp.fitness/setp.numerosity;
                                    }
                                    size=1;
                                }
                                else
                                {
                                    /* another guy in the tournament */
                                    if( (!doGAErrorBasedSelect && fitness + selectTolerance > setp.fitness/setp.numerosity) ||(doGAErrorBasedSelect && fitness - selectTolerance * maxPayoff < value))
                                    {
                                        /* both classifiers in tournament have a similar fitness/error */
                                        size += addClassifierToPointerSet(&setp, &winnerSet);
                                    }
                                    else
                                    {
                                        /* new classifier in tournament is clearly better */
                                        freeSet(&winnerSet);
                                        winnerSet=NULL;
                                        addClassifierToPointerSet(&setp, &winnerSet);
                                        if(doGAErrorBasedSelect)
                                        {
                                            fitness = value;
                                        }
                                        else
                                        {
                                            fitness = setp.fitness/setp.numerosity;
                                        }
                                        size=1;
                                    }
                                }
                                break; /* go to next classifier since this one is already a winner*/
                            }
                        }
                    }
                }
            }
        }
    }
    /* choose one of the equally best winners at random */
    size = irand(size);
    for(setp=winnerSet; setp!=NULL; setp=setp->next)
    {
        if(size==0)
        {
            break;
        }
        size--;
    }
    assert(setp!=NULL);
    winner = setp->classifier;
    freeSet(&winnerSet);
    return winner;
}

/**
 * Select a classifier for the discovery mechanism using roulette wheel selection
 */
Classifier* selectClassifierUsingRWS(ClassifierSet *set, double fitsum)
{
    ClassifierSet *setp;
    double choicep;

    choicep=drand()*fitsum;
    setp=set;
    fitsum=setp->classifier->fitness;
    while(choicep>fitsum)
    {
        setp=setp->next;
        fitsum+=setp->classifier->fitness;
    }

    return setp->classifier;
}

// ########################## crossover and mutation ########################################

void crossover_filter(Filter& parent1, Filter& parent2)
{
    int point1 = irand(numLeaf);
    int point2 = irand(numLeaf);

    if(point1 > point2){
        int temp = point1;
        point1 = point2;
        point2 = temp;
    }

    float temp_lower, temp_upper;
    for(int i=point1; i<point2; i++){
        temp_lower = parent1.lower_bounds[i];
        temp_upper = parent1.upper_bounds[i];
        parent1.lower_bounds[i] = parent2.lower_bounds[i];
        parent1.upper_bounds[i] = parent2.upper_bounds[i];
        parent2.lower_bounds[i] = temp_lower;
        parent2.upper_bounds[i] = temp_upper;
    }
}


bool crossover(Classifier &cl1, Classifier &cl2, float situation[])
{
    Filter filter1, filter2;
    bool filter1_result = false;
    bool filter2_result = false;

    for (int i = 0; i < clfrCondLength; i++) {
        if(drand() < pX){  // swap CF with pX probability
            CodeFragment temp = cl1.condition[i];
            cl1.condition[i] = cl2.condition[i];
            cl2.condition[i] = temp;
            continue; // if cf crossover done then skip filter crossover
        }
        for(int j=0; j < cl1.condition[i].num_filters && j < cl2.condition[i].num_filters; j++) {
            if(drand() < pX) {  // filter crossover with pX probability
                filter1 = get_filter(cl1.condition[i].filter_id[j]);
                filter2 = get_filter(cl2.condition[i].filter_id[j]);
                // Skip crossover if filters are not of the same size or type
                if(filter1.filter_size != filter2.filter_size ||
                filter1.is_dilated != filter2.is_dilated) continue;
                filter1_result = evaluate_filter(filter1, situation);
                filter2_result = evaluate_filter(filter2, situation);
                for (int tries = 0; tries < 100; tries++) {
                    crossover_filter(filter1, filter2);
                    if (filter1_result == evaluate_filter(filter1, situation)
                        && filter2_result == evaluate_filter(filter2, situation)) {
                        cl1.condition[i].filter_id[j] = add_filter(filter1);
                        cl2.condition[i].filter_id[j] = add_filter(filter2);
                        break;
                    } else {
                        filter1 = get_filter(cl1.condition[i].filter_id[j]);
                        filter2 = get_filter(cl2.condition[i].filter_id[j]);
                    }
                }
            }
        }
    }
    return true;
}


bool crossover_old(Classifier **cl, int crossoverType)  // Determines if crossover is applied and calls then the selected crossover type.
{
    bool changed = false;
    if(drand()<pX)
    {
        changed = true;
        if(crossoverType == 0)
        {
            uniformCrossover(cl);
        }
        else if(crossoverType == 1)
        {
            onePointCrossover(cl);
        }
        else
        {
            twoPointCrossover(cl);
        }
    }
    return changed;
}

void uniformCrossover(Classifier **cl)  // Crosses the two received classifiers using uniform crossover.
{
    CodeFragment help;
    for(int i=0; i<clfrCondLength; i++)
    {
        if(drand() < 0.5)
        {
            memmove(&help,&(cl[0]->condition[i]),sizeof(CodeFragment));
            memmove(&(cl[0]->condition[i]),&(cl[1]->condition[i]),sizeof(CodeFragment));
            memmove(&(cl[1]->condition[i]),&help,sizeof(CodeFragment));
        }
    }
}

void onePointCrossover(Classifier **cl){  // Crosses the two received classifiers using one-point crossover.

    CodeFragment help;
    int sep = irand(clfrCondLength);
    if(sep<0 || sep>=clfrCondLength)
    {
        printf("\ninvalid crossover point\n");
        exit(0);
    }
    for(int i=0; i<=sep; i++)
    {
        memmove(&help,&(cl[0]->condition[i]),sizeof(CodeFragment));
        memmove(&(cl[0]->condition[i]),&(cl[1]->condition[i]),sizeof(CodeFragment));
        memmove(&(cl[1]->condition[i]),&help,sizeof(CodeFragment));
    }

}

void twoPointCrossover(Classifier **cl){  // Crosses the two received classifiers using two-point crossover.

    CodeFragment help;
    int sep1 = irand(clfrCondLength);
    int sep2 = irand(clfrCondLength);
    if(sep1<0 || sep1>=clfrCondLength || sep2<0 || sep2>=clfrCondLength)
    {
        printf("\ninvalid crossover point\n");
        exit(0);
    }
    if(sep1>sep2)
    {
        int help=sep1;
        sep1=sep2;
        sep2=help;
    }
    for(int i=sep1; i<=sep2; i++)
    {
        memmove(&help,&(cl[0]->condition[i]),sizeof(CodeFragment));
        memmove(&(cl[0]->condition[i]),&(cl[1]->condition[i]),sizeof(CodeFragment));
        memmove(&(cl[1]->condition[i]),&help,sizeof(CodeFragment));
    }
}
/**
 * Apply mutation to classifier 'clfr'.
 * If niche mutation is applied, 'state' is considered to constrain mutation.
 * returns if the condition was changed.
 */
/*
bool mutation_old(Classifier *clfr, float state[])
{
    bool changedCond = false, changedAct = false;
    if(mutationType == 0)
    {
        changedCond = applyNicheMutation2(clfr,state);
    }
    else
    {
        changedCond = applyGeneralMutation(clfr,state);
    }
    changedAct = mutateAction(clfr);
    return (changedCond || changedAct);
}
*/

void apply_filter_mutation(Filter& filter, float state[])
{
    float delta = 0;

    for(int i=0; i<filter.filter_size*filter.filter_size; i++){
        // todo: pM should be used only to initiate mutation. Prob of mutating each elle should be less?
        if(drand() < pM) {  // mutation based on mutation probability
            delta = drand()*m;  // how much to mutate
            if(drand() < 0.5){  // revert sign with 50% probability
                delta *= -1;
            }
            if(drand() < 0.5){
                filter.lower_bounds[i] = roundRealValue(fmax(filter.lower_bounds[i] + delta, 0), precisionDigits);
            }else{
                filter.upper_bounds[i] = roundRealValue( fmin( filter.upper_bounds[i] + delta, 1), precisionDigits);
            }
        }
    }
}

bool mutation(Classifier &clfr, float *state)
{
    Filter filter_to_mutate;
    bool previous_evaluation_result = false;

    for(int i=0; i<clfrCondLength; i++){
        // 2 level mutation (CF and filter)
        if(drand() < pM){
           if(mutate_cf(clfr.condition[i])) {
               continue; // if cf is mutated then skip filter mutation
           }
        }
        for(int j=0; j<clfr.condition[i].num_filters; j++) {
            if(drand() >= pM) continue; // mutate only with probability of pM

            // if current filter is not promising and there is promising filter available in filter store
            // then use it with probability p_ol
            if(get_filter(clfr.condition[i].filter_id[j]).fitness < 1 && drand() < p_ol){
                int id = get_promising_filter_id();
                if(id != -1){   // if promising filter found
                   clfr.condition[i].filter_id[j] = id;
                }
            }
            filter_to_mutate = get_filter(clfr.condition[i].filter_id[j]);
            previous_evaluation_result = evaluate_filter(filter_to_mutate, state);
            for(int tries = 0; tries < 100; tries++){
                apply_filter_mutation(filter_to_mutate, state);
                if(previous_evaluation_result == evaluate_filter(filter_to_mutate, state)){
                    clfr.condition[i].filter_id[j] = add_filter(filter_to_mutate);
                    break;
                }else{
                    filter_to_mutate = get_filter(clfr.condition[i].filter_id[j]);
                }
            }
        }
    }
    // mutate action
    mutateAction(clfr);
    return true;
}


/**
 * Mutates the condition of the classifier. If one allele is mutated depends on the constant pM.
 * This mutation is a niche mutation. It assures that the resulting classifier still matches the current situation.
 */
bool applyNicheMutation(Classifier *clfr, float state[])
{
    bool changed = false;

    // Mutation as done paper XOF. It achieves best results among all 3 algos
    for(int i=0; i<clfrCondLength; i++)
    {
        if(drand()<pM)
        {
            changed = true;
            if(drand()<m_0)//(isDontcareCF(clfr->condition[i]))
            {
                CodeFragment tempCF;
                do
                {
                    //ramped half-and-half, check end for sanity only
                    tempCF = createNewCF(getNumPreviousCFs()+countNewCFs);
                    opType* end = randomProgram(tempCF.codeFragment,irand(2),cfMaxDepth,cfMinDepth);
                    validateDepth(tempCF.codeFragment,end); //validate depth
                    ////////
                    tempCF = addLeafCF(tempCF, state);
                    ////////
                }
                while( evaluateCF(tempCF, state)!=1 );
                //while( evaluateCF(tempCF.codeFragment,state)!=1 );
                memmove(&clfr->condition[i],&tempCF,sizeof(CodeFragment));
                countNewCFs++;
            }
        }
    }

   return changed;
}

/*
// copy of original and using first alternative
// Mutation code copied from XCSR code Hassan, although not being used even in that paper
bool applyNicheMutation1(Classifier *clfr, float state[])
{
    bool changed = false;
    // Mutation code copied from XCSR code Hassan, although not being used even in that paper
    for(int i=0; i<clfrCondLength; i++)
    {
        if(drand()<pM)
        {
            changed = true;
            float temp = drand()*m_0;
            if(drand()<m_0)
            {
                temp = -temp;
            }
            if(drand()<m_0)
            {
                for(int k=0; k<numLeaf; k++)
                {
                       clfr->condition[i].leaf[k].lowerBound += temp;
                       if(clfr->condition[i].leaf[k].lowerBound>state[i])
                            clfr->condition[i].leaf[k].lowerBound = state[i] - drand()*m_0;
                }
            }
            else
            {
                for(int k=0; k<numLeaf; k++)
                {
                       clfr->condition[i].leaf[k].upperBound += temp;
                       if(state[i]>clfr->condition[i].leaf[k].upperBound)
                            clfr->condition[i].leaf[k].upperBound = state[i] + drand()*m_0;
                }
            }
        }
    }

   return changed;
}
*/

/*
// copy of original and using 2nd alternative
/////////// original code as per LML paper of Hassan /////////////////
bool applyNicheMutation2(Classifier *clfr, float state[])
{
    bool changed = false;
    // Mutation code copied from XCSR code Hassan, although not being used even in that paper
   /////////// original code as per LML paper of Hassan /////////////////
    for(int i=0; i<clfrCondLength; i++)
    {
        if(drand()<pM)
        {
            changed = true;
            if(isDontcareCF(clfr->condition[i]))
            {
                CodeFragment tempCF;
                do
                {
                    //ramped half-and-half, check end for sanity only
                    tempCF = createNewCF(getNumPreviousCFs()+countNewCFs);
                    opType* end = randomProgram(tempCF.codeFragment,irand(2),cfMaxDepth,cfMinDepth);
                    validateDepth(tempCF.codeFragment,end); //validate depth
                    ////////
                    tempCF = addLeafCF(tempCF, state);
                    ////////
                }
                while( evaluateCF(tempCF, state)!=1 );
                //while( evaluateCF(tempCF.codeFragment,state)!=1 );
                memmove(&clfr->condition[i],&tempCF,sizeof(CodeFragment));
                countNewCFs++;
                clfr->specificness++;
            }
            else
            {
                memmove(&clfr->condition[i],&dontcareCF,sizeof(CodeFragment));//dontcareCF;
                clfr->specificness--;
            }
        }
    }
    return changed;
}
*/

/**
 * Mutates the condition of the classifier. If one allele is mutated depends on the constant pM.
 * This mutation is a general mutation.
 */
/*
bool applyGeneralMutation(Classifier *clfr, float state[])
{
    bool changed = false;
    for(int i=0; i<clfrCondLength; i++)
    {
        if(drand()<pM)
        {
            changed = true;
            if(isDontcareCF(clfr->condition[i]))
            {
                CodeFragment tempCF;
                //ramped half-and-half, check end for sanity only
                tempCF = createNewCF(getNumPreviousCFs()+countNewCFs);
                opType* end = randomProgram(tempCF.codeFragment,irand(2),cfMaxDepth,cfMinDepth);
                validateDepth(tempCF.codeFragment,end); //validate depth
                //////
                tempCF = addLeafCF(tempCF, state);
                //////
                memmove(&clfr->condition[i],&tempCF,sizeof(CodeFragment));
                countNewCFs++;
                clfr->specificness++;
            }
            else
            {
                memmove(&clfr->condition[i],&dontcareCF,sizeof(CodeFragment));//dontcareCF;
                clfr->specificness--;
            }
        }
    }
    return changed;
}
*/

bool mutateAction(Classifier& clfr)  //Mutates the action of the classifier.
{
    bool changed = false;
    if(drand()<pM)
    {
        changed = true;
        int act=0;
        do
        {
            act = irand(numActions);
        }
        while(act==clfr.action);
        clfr.action=act;
    }
    return changed;
}

// ###################### offspring insertion #################################

/**
 * Insert a discovered classifier into the population and respects the population size.
 */
void insertDiscoveredClassifier(Classifier cl[], Classifier parents[], ClassifierList &action_set, ClassifierList &pop,
                                int len, float *state)
{
    Classifier *killedp;
    len+=2;
    if(doGASubsumption)
    {
        subsumeClassifier(cl[0], parents, action_set, pop, state);
        subsumeClassifier(cl[1], parents, action_set, pop, state);
    }
    else
    {
        pop.push_front(cl[0]);
        pop.push_front(cl[1]);
    }

    while(len > maxPopSize)
    {
        len--;
        int cl_id = deleteStochClassifier(pop);

        /* record the deleted classifier to update other sets */
//        if(killedp!=NULL)
//        {
//            addClassifierToPointerSet(killedp,killset);
//            /* update the set */
//            updateSet(action_set, *killset);
//        }
    }
}

// ################################ subsumption deletion #################################

/**
 * Action set subsumption as described in the algorithmic describtion of XCS
 */
void doActionSetSubsumption(ClassifierList &action_set, ClassifierList &pop)
{
    Classifier *subsumer=NULL;

    /* Find the most general subsumer */
    for(auto it : action_set)
    {
        if(isSubsumer(it))
        {
            if(subsumer==NULL || isMoreGeneral(it, *subsumer))
            {
                subsumer = &it;
            }
        }
    }

    /* If a subsumer was found, subsume all classifiers that are more specific. */
    if(subsumer!=NULL)
    {
        auto setp = action_set.begin();
        auto setpl = action_set.begin();
        int cl_id = 0;
        for(; setp != action_set.end(); setp++)
        {
            while(isMoreGeneral(*subsumer, *setp))
            {
                subsumer->numerosity += setp->numerosity;
                remove_classifier(pop, setpl->id);
                if(setp == setpl){
                    setp++;
                    setpl++;
                }else{
                    setp = setpl;
                }
            }
            setpl=setp;
        }
    }
}

/**
 * Tries to subsume the parents.
 */
void subsumeClassifier(Classifier& cl, Classifier parents[], ClassifierList &action_set, ClassifierList &pop, float *state)
{
    int i;
    for(i=0; i<2; i++)
    {
        if(subsumes(parents[i],cl))
        {
            parents[i].numerosity++;
            return;
        }
    }
    // changed from action set subsumption to population subsumption
    if(subsumeClassifierToSet(cl, pop , nullptr))
    {
        return;
    }
    pop.push_front(cl);
}

/**
 * Try to subsume in the specified set.
 */
bool subsumeClassifierToSet(Classifier &cl, ClassifierList &cl_set, ClassifierSet *set)
{
    Classifier *subCl[maxPopSize];
    int numSub=0;

    for(auto & it : cl_set)
    {
        if(subsumes(it,cl))
        {
            subCl[numSub]=&it;
            numSub++;
        }
    }
    /* if there were classifiers found to subsume, then choose randomly one and subsume */
    if(numSub>0)
    {
        numSub = irand(numSub);
        subCl[numSub]->numerosity++;
        return true;
    }
    return false;
}

bool subsumes(Classifier &cl1, Classifier &cl2)  // check if classifier cl1 subsumes cl2
{
    return cl1.action==cl2.action && isSubsumer(cl1) && isMoreGeneral(cl1,cl2);
}

bool isSubsumer(Classifier &cl)
{
    return cl.experience > theta_sub && cl.predictionError <= epsilon_0;
}

// function added by me in new code
// To subsume a CF to more general CF
// filhal randomly selected two CFs, In future will select CFs based on their Fitness
/*
 * code review notes
 * this function recreate a new CF if it sumbsumed by another CF after randomly selecting two CFs
 * it essentially creating a new code fragment if it is already covered by another in the same classifer
 */
void subsumeCFs(Classifier *clfr, float state[])
{
    int tmpindex1 = rand() % clfrCondLength;
    int tmpindex2 = rand() % clfrCondLength;

    if(isGeneralCF(clfr->condition[tmpindex1], clfr->condition[tmpindex2]))
    {
            CodeFragment tempCF;
            do
            {
                //ramped half-and-half, check end for sanity only
                tempCF = createNewCF(getNumPreviousCFs()+countNewCFs);
                opType* end = randomProgram(tempCF.codeFragment,irand(2),cfMaxDepth,cfMinDepth);
                validateDepth(tempCF.codeFragment,end); //validate depth
                ////////
                tempCF = addLeafCF(tempCF, state);
                ////////
            }
            while( evaluateCF(tempCF, state)!=1 );
            //while( evaluateCF(tempCF.codeFragment,state)!=1 );
            memmove(&clfr->condition[tmpindex2],&tempCF,sizeof(CodeFragment));
            //countNewCFs++;
    }
}

bool is_filter_general(const Filter& filter_general, const Filter& filter_to_check)
{
    // only filter of same size and type can be general
    if(filter_general.filter_size != filter_to_check.filter_size ||
    filter_general.is_dilated != filter_to_check.is_dilated) return false;
    for(int i=0; i<filter_general.filter_size*filter_general.filter_size; i++){

            if(filter_to_check.lower_bounds[i] < filter_general.lower_bounds[i]
            || filter_to_check.upper_bounds[i] > filter_general.upper_bounds[i]){
                return false;
            }
    }
    return true;
}

bool is_filter_covered_by_condition(int filter_to_check_id, CodeFragment code_fragments[])
{
    for(int i=0; i<clfrCondLength; i++){
        for(int j=0; j < code_fragments[i].num_filters; j++) {
            if(filter_to_check_id == code_fragments[i].filter_id[j]){
                return true;
            }
        }
    }
    return false;
}

/*
 * This function checks that all the filters of one classifier are present in the second.
 * todo: At this time it does not check the result of the evaluation of the code fragment. It only compares filters.
 */
 bool isMoreGeneral(Classifier &clfr_general, Classifier &clfr_specific)
 {
     for(int i=0; i<clfrCondLength; i++){
         for(int j=0; j < clfr_specific.condition[i].num_filters; j++) {
             if (!is_filter_covered_by_condition(clfr_specific.condition[i].filter_id[j], clfr_general.condition)) {
                 return false;
             }
         }
     }
    return true;
 }

/*
  * Code review notes
  * It appears the isMmoreGeneral does not actually check the the generality but the equality
  * It needs to be updated.
  */
 bool isMoreGeneral_old(Classifier *clfr1, Classifier *clfr2)
{
    /* clfr1 is more general than clfr2 if:
     * Number of dontcares in clfr1 > Number of dontcares in clfr2, and
     * Each non-dontcare in clfr1 is in clfr2.
     */
    if(compareDontcares(clfr1,clfr2) && checkNonDontcares(clfr1->condition,clfr2->condition))
    {
        return true;
    }
    return false;
}
bool compareDontcares(Classifier *clfr1, Classifier *clfr2)  //returns true if "Number of dontcares in clfr1 > Number of dontcares in clfr2".
{
    return (clfr1->specificness < clfr2->specificness)? true : false;
}
bool checkNonDontcares(CodeFragment cond1[], CodeFragment cond2[])  //returns true if "Each non-dontcare in cond1 is in cond2."
{
    for(int i=0; i<clfrCondLength; i++)
    {
        if(!isDontcareCF(cond1[i]) && !isExists(cond1[i],cond2,clfrCondLength))
        {
            return false;
        }
    }
    return true;
}

// ###################### adding classifiers to a set ###################################

/**
 * Adds only the pointers to the pointerset, ensures that no pointer is added twice,
 * returns if the pointer was added
 */
bool addClassifierToPointerSet(Classifier *cl, ClassifierSet **pointerset)
{
    ClassifierSet *setp;
    for(setp=*pointerset; setp!=NULL; setp=setp->next)
    {
        if(setp->classifier == cl)  // classifier is already in Set
        {
            return false;
        }
    }
    // add the classifier, as it is not already in the pointerset
    assert((setp=( struct ClassifierSet*)calloc(1,sizeof(ClassifierSet)))!=NULL);
    setp->classifier=cl;  // todo this has already been deleted ???
    setp->next=*pointerset;
    *pointerset=setp;
    return true;
}

// adds the classifier cl to the population, makes sure that the same classifier does not exist yet.
bool addClassifierToSet(Classifier *cl, ClassifierSet **clSet)
{
    ClassifierSet *setp;
    // Check if classifier exists already. If so, just increase the numerosity and free the space of the new classifier
    for(setp=*clSet; setp!=NULL; setp=setp->next)
    {
        if( equals(setp->classifier,cl) )
        {
            setp->classifier->numerosity++;
            freeClassifier(cl);
            return true;
        }
    }
    // classifier does not exist yet -> add new classifier
    assert((setp=( struct ClassifierSet*)calloc(1,sizeof(struct ClassifierSet)))!=NULL);
    setp->classifier=cl;
    setp->next=*clSet;
    *clSet=setp;
    return false;
}

// adds the new Clasifier cl to the ClassifierSet 'clSet'.
void addNewClassifierToSet(Classifier *cl, ClassifierSet **clSet)
{
    ClassifierSet *setp;
    assert((setp=( struct ClassifierSet *)calloc(1,sizeof(struct ClassifierSet)))!=NULL);
    setp->classifier=cl;
    setp->next=*clSet;
    *clSet=setp;
}
bool equals(Classifier *clfr1, Classifier *clfr2)  //Returns if the two classifiers are identical in condition and action.
{
    if(clfr1->specificness!=clfr2->specificness || clfr1->action != clfr2->action)
    {
        return false;
    }
    if(!checkNonDontcares(clfr1->condition, clfr2->condition))
    {
        return false;
    }
    return true;
}

// ############################## deletion ############################################

/**
 * Deletes one classifier in the population.
 * The classifier that will be deleted is chosen by roulette wheel selection
 * considering the deletion vote. Returns position of the macro-classifier which got decreased by one micro-classifier.
 **/
int deleteStochClassifier(ClassifierList &pop)
{
    double sum=0.0, choicep, meanf=0.0;
    int size=0;

    for(auto it : pop){
        meanf += it.fitness;
        size += it.numerosity;
    }
    meanf/=(double)size;

    /* get the delete proportion, which depends on the average fitness */
    for(auto it : pop){
        sum += getDelProp(&it,meanf);
    }

    /* choose the classifier that will be deleted */
    choicep=drand()*sum;
    /* look for the classifier */
    auto setp = pop.begin();
    auto setpl = pop.begin();
    sum = getDelProp(&*setp,meanf);
    while(sum < choicep)
    {
        setpl=setp;
        setp++;
        sum += getDelProp(&*setp,meanf);
    }

    int cl_id = setp->id; // id of the removed classifier
    /* delete the classifier */
    pop.erase_after(setpl);

    /* return the pointer to the deleted classifier, to be able to update other sets */
    return cl_id;
}
double getDelProp(Classifier *clfr, double meanFitness)  //Returns the vote for deletion of the classifier.
{
    if(clfr->fitness/(double)clfr->numerosity >= delta*meanFitness || clfr->experience < theta_del)
    {
        return (double)(clfr->actionSetSize*clfr->numerosity);
    }
    else
    {
        return (double)clfr->actionSetSize*(double)clfr->numerosity*meanFitness / (clfr->fitness/(double)clfr->numerosity);
    }
}
/**
 * Deletes the classifier setp from the population pop, setpl points to the classifier that is before setp in the list
 */
Classifier* deleteTypeOfClassifier(ClassifierSet *setp, ClassifierSet *setpl, ClassifierSet **pop)
{
    Classifier *killedp=NULL;

    /* setp must point to some classifier! */
    assert(setp!=NULL);

    if(setp->classifier->numerosity>1)
    {
        /* if the numerosity is greater than one -> just decrease it */
        setp->classifier->numerosity--;
    }
    else
    {
        /* if not, delete it and record it in killedp */
        if(setp==setpl)
        {
            *pop=setp->next;
        }
        else
        {
            setpl->next=setp->next;
        }
        killedp=setp->classifier;
        delete setp->classifier;
        delete setp;
    }
    /* return a pointer ot a deleted classifier (NULL if the numerosity was just decreased) */  // todo: what???
    return killedp;
}
/**
 * Check if the classifier pointers that are in killset are in uset - delete the pointers,
 * if they are inside. Note that the classifiers in killset are already deleted, so do not
 * read their values or try to delete them again here!
 */
bool updateSet(ClassifierSet **uset, ClassifierSet *killset)
{
    ClassifierSet *setp,*setpl,*killp,*usetp;
    bool updated = true;

    if(*uset==NULL || killset==NULL)  // if one of the sets is empty -> do not do anything
    {
        return false;
    }
    // check all classifiers in uset
    setp=*uset;
    while(updated && setp!=NULL)
    {
        setp=*uset;
        setpl=*uset;
        updated = false;
        while(setp!=NULL && !updated)
        {
            for(killp=killset; killp!=NULL; killp=killp->next)
            {
                if(killp->classifier == setp->classifier)
                {
                    // if killed classifier found, delete the struct classifier set in uset
                    updated = true;
                    if(setp==setpl)   // first entry in set
                    {
                        usetp=*uset;
                        *uset=usetp->next;
                        free(usetp);
                        break;
                    }
                    else
                    {
                        setpl->next=setp->next;
                        free(setp);
                        setp = *uset;
                        setpl= *uset;
                        break;
                    }
                }
            }
            if(updated)  // check the whole uset again, if one pointer was deleted
            {
                break;
            }
            setpl=setp;
            setp=setp->next;
        }
    }
    return updated;
}

/**
 * Deletes the classifier pointer from the specified set
 * and returns if the pointer was found and deleted.
 */
bool deleteClassifierPointerFromSet(ClassifierSet **set, Classifier *clp)
{
    ClassifierSet *setp, *setpl;
    for(setp=*set, setpl=*set; setp!=NULL; setp=setp->next)
    {
        if(setp->classifier==clp)
        {
            if(setpl==setp)
            {
                *set=(*set)->next;
                free(setp);
            }
            else
            {
                setpl->next=setp->next;
                free(setp);
            }
            return true;
        }
        setpl=setp;
    }
    return false;
}

//############# concrete deletion of a classifier or a whole classifier set ############

/**
 * Frees only the complete ClassifierSet (not the Classifiers itself)!
 */
void freeSet(ClassifierSet **cls)
{
    ClassifierSet *clp;
    while(*cls!=NULL)
    {
        clp=(*cls)->next;
        free(*cls);
        *cls=clp;
    }
}

/**
 * Frees the complete ClassifierSet with the corresponding Classifiers.
 */
void freeClassifierSet(ClassifierSet **cls)
{
    ClassifierSet *clp;
    while(*cls!=NULL)
    {
        freeClassifier((*cls)->classifier);
        clp=(*cls)->next;
        free(*cls);
        *cls=clp;
    }
}

/**
 * Frees one classifier.
 */
void freeClassifier(Classifier *cl)
{
//    // remove filters from the master list
//    for(int i=0; i<cl->condition->num_filters; i++){
//        remove_filter(cl->condition->filter_id[i]);
//    }
    free(cl);
}

// ############################### output operations ####################################

/**
 * print the classifiers in a ClassifierSet
 */
void printClassifierSet(ClassifierSet *head)
{
    for(; head!=NULL; head=head->next)
    {
        printClassifier(head->classifier);
    }
}



/**
 * print the classifier in a ClassifierSet to the file fp
 */
void fprintClassifierSet(FILE *fpClfr, FILE *fpCF, ClassifierSet *head)
{
    print_filter_stats();
    print_filter_evaluation_stats();
    ClassifierSet* set;
    for(set=head; set!=NULL; set=set->next)
    {
        fprintClassifier(fpClfr, set->classifier);
    }
    std::cout << "Global Classifier ID: " << gid << std::endl;
    storeCFs(head,fpCF);
}

/**
 * print a single classifier
 */
void printClassifier(Classifier *clfr)
{
    for(int i=0; i<clfrCondLength; i++)
    {
        printCF(clfr->condition[i]);
        printf("\n");
    }
    printf("Action: %d\n",clfr->action);
    printf("Numerosity: %d ",clfr->numerosity);
    printf("Accuracy: %f ",clfr->accuracy);
    printf("Fitness: %f ",clfr->fitness);
    printf("Prediction Error: %f ",clfr->predictionError);
    printf("Prediction: %f ",clfr->prediction);
    printf("Experience: %d ",clfr->experience);
    printf("specificness: %d\n",clfr->specificness);
}

/**
 * print a single classifier to the file fp
 */
/*
void fprintClassifier(FILE *fp, Classifier *classifier){
	char buf[1000];
	for(int i=0; i<clfrCondLength; i++){
		outprog(classifier->condition[i].codeFragment,cfMaxLength,fp);
		fwrite("\n",strlen("\n"),1,fp);
	}
	sprintf(buf,"Action: %d\n",classifier->action); fwrite(buf,strlen(buf),1,fp);
	sprintf(buf,"Numerosity: %d ",classifier->numerosity); fwrite(buf,strlen(buf),1,fp);
	sprintf(buf,"Accuracy: %f ",classifier->accuracy); fwrite(buf,strlen(buf),1,fp);
	sprintf(buf,"Fitness: %f ",classifier->fitness); fwrite(buf,strlen(buf),1,fp);
	sprintf(buf,"Prediction Error: %f ",classifier->predictionError); fwrite(buf,strlen(buf),1,fp);
	sprintf(buf,"Prediction: %f ",classifier->prediction); fwrite(buf,strlen(buf),1,fp);
	sprintf(buf,"Experience: %d ",classifier->experience); fwrite(buf,strlen(buf),1,fp);
	sprintf(buf,"Specificness: %d\n",classifier->specificness); fwrite(buf,strlen(buf),1,fp);
}
*/

void fprintClassifier(FILE *fp, Classifier *classifier)
{
    char *buf;
    int len;
    for(int i=0; i<clfrCondLength; i++)
    {
        //outprog(classifier->condition[i].codeFragment,cfMaxLength,fp);
        outprog(classifier->condition[i],cfMaxLength,fp);
        fwrite("\n",strlen("\n"),1,fp);
    }

    len = snprintf(NULL,0,"Action: %d\n",classifier->action);
    if(!(buf = (char*)malloc((len + 1) * sizeof(char))))
    {
        printf("\nError in file writing ...\n");
        exit(0);
    }
    len = snprintf(buf,len+1,"Action: %d\n",classifier->action);
    fwrite(buf,strlen(buf),1,fp);
    free(buf);

    len = snprintf(NULL,0,"id: %d ",classifier->id);
    if(!(buf = (char*)malloc((len + 1) * sizeof(char))))
    {
        printf("\nError in file writing ...\n");
        exit(0);
    }
    len = snprintf(buf,len+1,"id: %d ",classifier->id);
    fwrite(buf,strlen(buf),1,fp);
    free(buf);

    len = snprintf(NULL,0,"Numerosity: %d ",classifier->numerosity);
    if(!(buf = (char*)malloc((len + 1) * sizeof(char))))
    {
        printf("\nError in file writing ...\n");
        exit(0);
    }
    len = snprintf(buf,len+1,"Numerosity: %d ",classifier->numerosity);
    fwrite(buf,strlen(buf),1,fp);
    free(buf);

    len = snprintf(NULL,0,"Accuracy: %f ",classifier->accuracy);
    if(!(buf = (char*)malloc((len + 1) * sizeof(char))))
    {
        printf("\nError in file writing ...\n");
        exit(0);
    }
    len = snprintf(buf,len+1,"Accuracy: %f ",classifier->accuracy);
    fwrite(buf,strlen(buf),1,fp);
    free(buf);

    len = snprintf(NULL,0,"Fitness: %f ",classifier->fitness);
    if(!(buf = (char*)malloc((len + 1) * sizeof(char))))
    {
        printf("\nError in file writing ...\n");
        exit(0);
    }
    len = snprintf(buf,len+1,"Fitness: %f ",classifier->fitness);
    fwrite(buf,strlen(buf),1,fp);
    free(buf);

    len = snprintf(NULL,0,"PredictionError: %f ",classifier->predictionError);
    if(!(buf = (char*)malloc((len + 1) * sizeof(char))))
    {
        printf("\nError in file writing ...\n");
        exit(0);
    }
    len = snprintf(buf,len+1,"PredictionError: %f ",classifier->predictionError);
    fwrite(buf,strlen(buf),1,fp);
    free(buf);

    len = snprintf(NULL,0,"Prediction: %f ",classifier->prediction);
    if(!(buf = (char*)malloc((len + 1) * sizeof(char))))
    {
        printf("\nError in file writing ...\n");
        exit(0);
    }
    len = snprintf(buf,len+1,"Prediction: %f ",classifier->prediction);
    fwrite(buf,strlen(buf),1,fp);
    free(buf);

    len = snprintf(NULL,0,"Experience: %d ",classifier->experience);
    if(!(buf = (char*)malloc((len + 1) * sizeof(char))))
    {
        printf("\nError in file writing ...\n");
        exit(0);
    }
    len = snprintf(buf,len+1,"Experience: %d ",classifier->experience);
    fwrite(buf,strlen(buf),1,fp);
    free(buf);

    len = snprintf(NULL,0,"Specificness: %d ",classifier->specificness);
    if(!(buf = (char*)malloc((len + 1) * sizeof(char))))
    {
        printf("\nError in file writing ...\n");
        exit(0);
    }
    len = snprintf(buf,len+1,"Specificness: %d ",classifier->specificness);
    fwrite(buf,strlen(buf),1,fp);
    free(buf);

    len = snprintf(NULL,0,"ActionSetSize: %f ",classifier->actionSetSize);
    if(!(buf = (char*)malloc((len + 1) * sizeof(char))))
    {
        printf("\nError in file writing ...\n");
        exit(0);
    }
    len = snprintf(buf,len+1,"ActionSetSize: %f ",classifier->actionSetSize);
    fwrite(buf,strlen(buf),1,fp);
    free(buf);

    len = snprintf(NULL,0,"TimeStamp: %d\n",classifier->timeStamp);
    if(!(buf = (char*)malloc((len + 1) * sizeof(char))))
    {
        printf("\nError in file writing ...\n");
        exit(0);
    }
    len = snprintf(buf,len+1,"TimeStamp: %d\n",classifier->timeStamp);
    fwrite(buf,strlen(buf),1,fp);
    free(buf);

    fflush(fp);
}

/*###################### sorting the classifier list ###################################*/

/**
 * Sort the classifier set cls in numerosity, prediction, fitness, or error order.
 * type 0 = numerosity order, type 1 = prediction order, type 2 = fitness order, type 3=error order
 */
ClassifierSet* sortClassifierSet(ClassifierSet **cls, int type)
{
    ClassifierSet *clsp, *maxcl, *newcls, *newclshead;
    double max;
    max=0.0;
    assert((newclshead=( struct ClassifierSet *)calloc(1,sizeof(struct ClassifierSet)))!=NULL);
    newcls=newclshead;
    do
    {
        max=-100000.0;
        /* check the classifier set cls for the next maximum -> already inserted classifier are referenced by the NULL pointer */
        for( clsp=*cls, maxcl=NULL; clsp!=NULL; clsp=clsp->next )
        {
            if(clsp->classifier!=NULL && (maxcl==NULL || ((type==0 && clsp->classifier->numerosity>max) || (type==1 && clsp->classifier->prediction>max) || (type==2 && clsp->classifier->fitness/clsp->classifier->numerosity > max) || (type==3 && -1.0*(clsp->classifier->predictionError) > max))))
            {
                if(type==0)
                {
                    max=clsp->classifier->numerosity;
                }
                else if (type==1)
                {
                    max=clsp->classifier->prediction;
                }
                else if (type==2)
                {
                    max=clsp->classifier->fitness/clsp->classifier->numerosity;
                }
                else if(type==3)
                {
                    max=-1.0*(clsp->classifier->predictionError);
                }
                maxcl=clsp;
            }
        }
        if(max>-100000.0)
        {
            assert((newcls->next=( struct ClassifierSet *)calloc(1,sizeof(struct ClassifierSet)))!=NULL);
            newcls=newcls->next;
            newcls->next=NULL;
            newcls->classifier=maxcl->classifier;
            maxcl->classifier=NULL; // do not delete the classifier itself, as it will be present in the new, sorted classifier list
        }
    }
    while(max>-100000.0);

    // set the new ClassifierSet pointer and free the old stuff
    newcls=newclshead->next;
    free(newclshead);
    freeSet(cls);

    return newcls; // return the pointer to the new ClassifierSet
}

void simplifyPopulation(ClassifierSet **population)
{
    ClassifierSet *set, *setp, *killset = NULL;
    int setSize = getSetSize(*population);
    printf("\nPopulation Size: %d\n",setSize);
    //consider only well experienced and accurate classifier rules
    for(set=*population,setp=*population; set!=NULL; set=set->next)
    {
        //printClassifier(set->classifier); printf("\n");
        while(set->classifier->experience <= 1.0/beta || set->classifier->predictionError > epsilon_0)
        {
            //printClassifier(set->classifier);	printf("\n");
            if(setp==set)
            {
                *population=set->next;
                deleteClassifierPointerFromSet(population,set->classifier);
                freeClassifier(set->classifier);
                addClassifierToPointerSet(set->classifier,&killset);
                free(set);
                set=*population;
                setp=*population;
            }
            else
            {
                setp->next=set->next;
                deleteClassifierPointerFromSet(population,set->classifier);
                freeClassifier(set->classifier);
                addClassifierToPointerSet(set->classifier,&killset);
                free(set);
                set=setp;
            }
        }
        setp=set;
    }
    freeSet(&killset);
    setSize = getSetSize(*population);
    printf("Simplified Population Size: %d\n",setSize);
}
/*################################## Utilitiy ##########################################*/

/**
 * Get the absolute value of 'value'.
 */
double absoluteValue(double value)
{
    if(value < 0.0)
    {
        return -1.0*value;
    }
    else
    {
        return value;
    }
}

/*
 * This function is used to update stats of all filters in the filter list
 * Which includes nuemrosity and fitness of the filter
 * If classifier is "promising" then increase the fitness of the filter
 * A promising classifier is one whose error < 10 and experience > 10
 */

void manage_filter_list(ClassifierList &pop){
    ClassifierSet * setp;
    Classifier* cl;

    // reset statistics of all filters before updating
    reset_filter_stats();

    for(auto it : pop){
        for(int i=0; i<clfrCondLength; i++){
            for(int j=0; j<it.condition[i].num_filters; j++){
                Filter& f = get_filter(it.condition[i].filter_id[j]);
                f.numerosity++;
                // if classifier is "promising" then increase the fitness of the filter
                // a promising classifier is one whose error < 10 and experience > 10
                if(it.predictionError < epsilon_0 && it.experience > theta_filter){
                    f.fitness++;
                }

            }
        }
    }
    // remove filters with numerosity=0 from the filter list
    std::forward_list<int> removed_filters;
    remove_unused_filters(removed_filters);
    update_evaluation_cache(removed_filters);
}

int get_set_size(ClassifierList &pop) {
    return std::distance(pop.begin(), pop.end());
}

int get_set_numerosity(ClassifierList &pop) {
    int pop_numerosity = 0;
    std::for_each(pop.begin(), pop.end(),
                  [&pop_numerosity](ClassifierList::value_type& cl)
                  {
                      pop_numerosity+= cl.numerosity;
                  });
    return pop_numerosity;
}

void get_matching_classifiers(ClassifierList &pop, float *state, ClassifierList &match_set, int img_id, bool train) {

    std::for_each(pop.begin(), pop.end(), [&match_set, &state, img_id, train](ClassifierList::value_type cl)
    {
        if(isConditionMatched(cl, state, img_id, train)){
            match_set.push_front(cl);
        }
    });
}
