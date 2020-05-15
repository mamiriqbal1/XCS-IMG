#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include <assert.h>
#include <fstream>
#include <vector>
#include "xcsMacros.h"
#include "codeFragment.h"
#include "classifier.h"
#include "env.h"

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
    clfr->actionSetSize = setSize;
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
 int covering = 0;
ClassifierSet* getMatchSet(ClassifierSet **population, ClassifierSet **killset, float state[], int itTime, int action, int img_id){
    ClassifierSet *mset=NULL, *poppointer;
    Classifier *killedp, *coverClfr;
    int popSize=0, setSize=0, representedActions;

    bool coveredActions[numActions];
    for(poppointer= *population; poppointer!=NULL; poppointer=poppointer->next)
    {
        popSize += poppointer->classifier->numerosity; // calculate the population size
        if(isConditionMatched(poppointer->classifier->condition,state, poppointer->classifier->id, img_id))
        {
            //std::cout<<"Previous\n";
            addNewClassifierToSet(poppointer->classifier, &mset); // add matching classifier to the matchset
            setSize+=poppointer->classifier->numerosity; // calculate size of the match set
        }
    }
    representedActions = nrActionsInSet(mset,coveredActions);

    // TEMP: insert filters every time
    if(popSize < maxPopSize/2) {
        representedActions = 0;
        coveredActions[0] = false;
        coveredActions[1] = false;
    }


    while(representedActions < numActions)  // create covering classifiers, if not all actions are covered
    {
        for(int i=0; i<numActions; i++)
        {
            if(coveredActions[i]==false)  // make sure that all actions are covered!
            {
                // TEMP: boost covering
                // add large number of classifier in case of filter approach in covering
                int add = 1;
                for(int j=0; j<add; j++) {
                    coverClfr = matchingCondAndSpecifiedAct(state, i, setSize + 1, itTime);
                    // before inserting the new classifier into the population check for subsumption by a generic one
                    // todo: setSize and popSize needs to be incremented in case of subsumption?
                    if(!subsumeClassifierToSet(coverClfr,*population)) {
                        addNewClassifierToSet(coverClfr, &mset);
                        setSize++;
                        addNewClassifierToSet(coverClfr, population);
                        popSize++;
                    }
                }
            }
        }

        /* Delete classifier if population is too big and record it in killset */
        while( popSize > maxPopSize )
        {
            /* PL */
            killedp = deleteStochClassifier(population);
            if(killedp!=NULL)
            {
                deleteClassifierPointerFromSet(&mset, killedp);
                addClassifierToPointerSet(killedp, killset);
            }
            popSize--;
        }
        representedActions = nrActionsInSet(mset,coveredActions);
    }
    return mset; // return the match set
}
/**
 * Returns the number of actions in the set and stores which actions are covered in the array coveredActions.
 */
int nrActionsInSet(ClassifierSet *set, bool coveredActions[])
{
    int nr;

    for(int i=0; i<numActions; i++)
    {
        coveredActions[i] = false;
    }

    for(nr=0; nr<numActions && set!=NULL; set=set->next)
    {
        if(coveredActions[set->classifier->action]==false)
        {
            coveredActions[set->classifier->action] = true;
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


bool isConditionMatched(CodeFragment clfrCond[], float state[], int cl_id, int img_id)
{
    for(int i=0; i<clfrCondLength; i++)
    {
        //std::cout<<"iscond\n";
        //if( !isDontcareCF(clfrCond[i]) && evaluateCF(clfrCond[i].codeFragment,state)==0 )
        if( !isDontcareCF(clfrCond[i]) && evaluateCF(clfrCond[i],state, cl_id, img_id)==0 )
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
        if(drand()<P_dontcare)
        {
            memmove(&cond[i],&dontcareCF,sizeof(CodeFragment));
        }
        else
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
    }
    //printf("out.......\n");
}

// ######################### prediction array operations ############################################

void getPredictionArray(ClassifierSet *ms)  //determines the prediction array out of the match set ms
{
    assert(ms!=NULL); // ms should never be NULL (because of covering)
    for(int i=0; i<numActions; i++)
    {
        predictionArray[i]=0.0;
        sumClfrFitnessInPredictionArray[i]=0.0;
    }
    for(; ms!=NULL ; ms=ms->next)
    {
        int actionValue = ms->classifier->action;
        predictionArray[actionValue]+= ms->classifier->prediction*ms->classifier->fitness;
        sumClfrFitnessInPredictionArray[actionValue]+= ms->classifier->fitness;
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

ClassifierSet* getActionSet(int action, ClassifierSet *ms)  // constructs an action set out of the match set ms.
{
    ClassifierSet *aset=NULL;
    for(; ms!=NULL; ms=ms->next)
    {
        if(action == ms->classifier->action)
        {
            addNewClassifierToSet(ms->classifier,&aset);
        }
    }
    return aset;
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
void updateActionSet(ClassifierSet **aset, double maxPrediction, double reward, ClassifierSet **pop, ClassifierSet **killset)
{
    double P, setsize=0.0;
    ClassifierSet *setp;

    P = reward + gama*maxPrediction;

    for(setp=*aset; setp!=NULL; setp=setp->next)
    {
        setsize += setp->classifier->numerosity;
        setp->classifier->experience++;
    }

    for(setp=*aset; setp!=NULL; setp=setp->next)   // update prediction, prediction error and action set size estimate
    {
        if((double)setp->classifier->experience < 1.0/beta)
        {
            // !first adjustments! -> simply calculate the average
            setp->classifier->predictionError = (setp->classifier->predictionError * ((double)setp->classifier->experience - 1.0) + absoluteValue(P - setp->classifier->prediction)) / (double)setp->classifier->experience;
            setp->classifier->prediction = (setp->classifier->prediction * ((double)setp->classifier->experience - 1.0) + P) / (double)setp->classifier->experience;
            setp->classifier->actionSetSize = (setp->classifier->actionSetSize *((double)(setp->classifier->experience - 1))+setsize)/(double)setp->classifier->experience;
        }
        else
        {
            // normal adjustment -> use widrow hoff delta rule
            setp->classifier->predictionError += beta * (absoluteValue(P - setp->classifier->prediction) - setp->classifier->predictionError);
            setp->classifier->prediction += beta * (P - setp->classifier->prediction);
            setp->classifier->actionSetSize += beta * (setsize - setp->classifier->actionSetSize);
        }
    }
    updateFitness(*aset);
    if(doActSetSubsumption)
    {
        doActionSetSubsumption(aset,pop,killset);
    }
}

// update the fitnesses of an action set (the previous [A] in multi-step envs or the current [A] in single-step envs.)
/*
 * code review notes
 * correctly implementation
 */
void updateFitness(ClassifierSet *aset)
{
    ClassifierSet *setp;
    double ksum=0.0;

    if(aset==NULL)  // if the action set got NULL (due to deletion) return
    {
        return;
    }

    //First, calculate the accuracies of the classifier and the accuracy sums
    for(setp=aset; setp!=NULL; setp=setp->next)
    {
        if(setp->classifier->predictionError <= epsilon_0)
        {
            setp->classifier->accuracy = 1.0;
        }
        else
        {
            setp->classifier->accuracy = alpha * pow(setp->classifier->predictionError / epsilon_0 , -nu);
        }
        ksum += setp->classifier->accuracy*(double)setp->classifier->numerosity;
    }

    //Next, update the fitnesses accordingly
    for(setp=aset; setp!=NULL; setp=setp->next)
    {
        setp->classifier->fitness += beta * ( (setp->classifier->accuracy * setp->classifier->numerosity) / ksum - setp->classifier->fitness );
    }
}

// ############################ discovery mechanism #########################################

/**
 * The discovery conmponent with the genetic algorithm
 * note: some classifiers in set could be deleted !
 */
void discoveryComponent(ClassifierSet **set, ClassifierSet **pop, ClassifierSet **killset, int itTime, float situation[])
{
    ClassifierSet *setp;
    Classifier *cl[2], *parents[2];
    double fitsum=0.0;
    int i, len, setsum=0, gaitsum=0;

    if(*set==NULL)  // if the classifier set is empty, return (due to deletion)
    {
        return;
    }

    getDiscoversSums(*set, &fitsum, &setsum, &gaitsum); // get all sums that are needed to do the discovery

    // do not do a GA if the average number of time-steps in the set since the last GA is less or equal than thetaGA
    if( itTime - (double)gaitsum / (double)setsum < theta_GA)
    {
        return;
    }
    setTimeStamps(*set, itTime);

    selectTwoClassifiers(cl, parents, *set, fitsum, setsum); // select two classifiers (tournament selection) and copy them
    if(crossover(cl, situation) || true){
        cl[0]->prediction   = (cl[0]->prediction + cl[1]->prediction) / 2.0;
        cl[0]->predictionError = predictionErrorReduction * ( (cl[0]->predictionError + cl[1]->predictionError) / 2.0 );
        cl[0]->fitness = fitnessReduction * ( (cl[0]->fitness + cl[1]->fitness) / 2.0 );

        cl[1]->prediction = cl[0]->prediction;
        cl[1]->predictionError = cl[0]->predictionError;
        cl[1]->fitness = cl[0]->fitness;
    }
    for(i=0; i<2; i++)  // do mutation
    {
        mutation(cl[i], situation);
    }

    cl[0]->specificness = numberOfNonDontcares(cl[0]->condition);
    cl[1]->specificness = numberOfNonDontcares(cl[1]->condition);

    // get the length of the population to check if clasifiers have to be deleted
    for(len=0, setp=*pop; setp!=NULL; setp=setp->next)
    {
        len += setp->classifier->numerosity;
    }

    // insert the new two classifiers and delete two if necessary
    insertDiscoveredClassifier(cl, parents, set, pop, killset, len, situation);
}
void getDiscoversSums(ClassifierSet *set, double *fitsum, int *setsum, int *gaitsum)  // Calculate all necessary sums in the set for the discovery component.
{
    ClassifierSet *setp;

    *fitsum=0.0;
    *setsum=0;
    *gaitsum=0;
    for(setp=set; setp!=NULL; setp=setp->next)
    {
        (*fitsum)+=setp->classifier->fitness;
        (*setsum)+=setp->classifier->numerosity;
        (*gaitsum) += setp->classifier->timeStamp*setp->classifier->numerosity;
    }
}
void setTimeStamps(ClassifierSet *set, int itTime)  // Sets the time steps of all classifiers in the set to itTime (because a GA application is occurring in this set!).
{
    for( ; set!=NULL; set=set->next)
    {
        set->classifier->timeStamp = itTime;
    }
}

// ########################### selection mechanism ########################################

/**
 * Select two classifiers using the chosen selection mechanism and copy them as offspring.
 */
void selectTwoClassifiers(Classifier **cl, Classifier **parents, ClassifierSet *set, double fitsum, int setsum)
{
    int i; /*,j, length;*/
    Classifier *clp;

    assert(set!=NULL);

    for(i=0; i<2; i++)
    {
        if(tournamentSize == 0)
        {
            clp = selectClassifierUsingRWS(set,fitsum);
        }
        else
        {
            if(i==0)
            {
                clp = selectClassifierUsingTournamentSelection(set, setsum, 0);
            }
            else
            {
                clp = selectClassifierUsingTournamentSelection(set, setsum, parents[0]);
            }
        }

        parents[i]=clp;

        assert((cl[i]=(struct Classifier *)calloc(1,sizeof(struct Classifier)))!=NULL);

        memmove( &(cl[i]->condition),&(clp->condition),sizeof(clp->condition) );

        cl[i]->id = gid++;
        cl[i]->action = clp->action;
        cl[i]->prediction = clp->prediction;
        cl[i]->predictionError = clp->predictionError;
        cl[i]->accuracy = clp->accuracy;
        cl[i]->specificness = clp->specificness;
        cl[i]->fitness = clp->fitness/(double)clp->numerosity;
        cl[i]->numerosity = 1;
        cl[i]->experience = 0;
        cl[i]->actionSetSize = clp->actionSetSize;
        cl[i]->timeStamp = clp->timeStamp;
    }
}

/**
 * Selects a classifier from 'set' using tournament selection.
 * If 'notMe' is not the NULL pointer and forceDifferentInTournament is set to a value larger 0,
 * this classifier is not selected except if it is the only classifier.
 */
Classifier* selectClassifierUsingTournamentSelection(ClassifierSet *set, int setsum, Classifier *notMe)
{
    ClassifierSet *setp, *winnerSet=NULL;
    Classifier *winner=NULL;
    double fitness=-1.0, value;
    int i, j, *sel, size=0;

    assert(set!=NULL);/* there must be at least one classifier in the set */

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
        return set->classifier;
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
        for(setp=set, i=0; setp!=NULL; setp=setp->next)
        {
            if(setp->classifier != notMe)
            {
                if(fitness < setp->classifier->fitness/setp->classifier->numerosity)
                {
                    for(j=0; j<setp->classifier->numerosity; j++)
                    {
                        if(sel[i+j])
                        {
                            freeSet(&winnerSet);
                            addClassifierToPointerSet(setp->classifier, &winnerSet);
                            fitness = setp->classifier->fitness/setp->classifier->numerosity;
                            break; /* go to next classifier since this one is already a winner*/
                        }
                    }
                }
                i += setp->classifier->numerosity;
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
            for(setp=set; setp!=NULL; setp=setp->next)
            {
                if(setp->classifier != notMe)   /* do not reselect the same classifier -> this only applies if forcedDifferentInTournament is set!*/
                {
                    value = setp->classifier->predictionError;
                    if(winnerSet==NULL || (!doGAErrorBasedSelect && fitness - selectTolerance <= setp->classifier->fitness/setp->classifier->numerosity) || (doGAErrorBasedSelect && fitness + selectTolerance * maxPayoff >= value))
                    {
                        /* if his fitness is worse then do not bother */
                        for(i=0; i<setp->classifier->numerosity; i++)
                        {
                            if(drand() < tournamentSize)
                            {
                                /* this classifier is a selection candidate and
                                * his fitness/error is higher/lower or similar to the other classifier */
                                if(winnerSet==NULL)
                                {
                                    /* the first guy in the tournament */
                                    addClassifierToPointerSet(setp->classifier, &winnerSet);
                                    if(doGAErrorBasedSelect)
                                    {
                                        fitness = value;
                                    }
                                    else
                                    {
                                        fitness = setp->classifier->fitness/setp->classifier->numerosity;
                                    }
                                    size=1;
                                }
                                else
                                {
                                    /* another guy in the tournament */
                                    if( (!doGAErrorBasedSelect && fitness + selectTolerance > setp->classifier->fitness/setp->classifier->numerosity) ||(doGAErrorBasedSelect && fitness - selectTolerance * maxPayoff < value))
                                    {
                                        /* both classifiers in tournament have a similar fitness/error */
                                        size += addClassifierToPointerSet(setp->classifier, &winnerSet);
                                    }
                                    else
                                    {
                                        /* new classifier in tournament is clearly better */
                                        freeSet(&winnerSet);
                                        winnerSet=NULL;
                                        addClassifierToPointerSet(setp->classifier, &winnerSet);
                                        if(doGAErrorBasedSelect)
                                        {
                                            fitness = value;
                                        }
                                        else
                                        {
                                            fitness = setp->classifier->fitness/setp->classifier->numerosity;
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

void crossover_filter(Leaf parent1[], Leaf parent2[])
{
    int point1 = irand(numLeaf);
    int point2 = irand(numLeaf);

    if(point1 > point2){
        int temp = point1;
        point1 = point2;
        point2 = temp;
    }

    Leaf temp;
    for(int i=point1; i<point2; i++){
        memmove(&temp, &parent1[i], sizeof(Leaf));
        memmove(&parent1[i], &parent2[i],  sizeof(Leaf));
        memmove(&parent2[i], &temp, sizeof(Leaf));
    }
}


bool crossover(Classifier **cl, float situation[])  // Determines if crossover is applied and calls then the selected crossover type.
{
//    twoPointCrossover(cl);
//    return true;

    Leaf previous1[numLeaf];
    Leaf previous2[numLeaf];
    if(drand() < pX) {
        for (int i = 0; i < clfrCondLength; i++) {
            if (!isDontcareCF(cl[0]->condition[i]) && !isDontcareCF(cl[1]->condition[i])) {
                memmove(&previous1, &cl[0]->condition[i].leaf, sizeof(previous1));
                memmove(&previous2, &cl[1]->condition[i].leaf, sizeof(previous2));
                crossover_filter(cl[0]->condition[i].leaf, cl[1]->condition[i].leaf);
                // revert if the resultant classifiers do not satisfies the situation
                if(!evaluateCF(cl[0]->condition[i], situation) || !evaluateCF(cl[1]->condition[i], situation)){
                    memmove(&cl[0]->condition[i].leaf, &previous1, sizeof(previous1));
                    memmove(&cl[1]->condition[i].leaf, &previous2, sizeof(previous2));
                }
            }
        }
        return true;
    }else{
        return false;
    }
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

void apply_filter_mutation(Leaf filter[], float state[])
{
    float delta = 0;

    for(int i=0; i<numLeaf; i++){
        if(drand() < pM) {  // mutation based on mutation probability
            delta = drand()*m;  // how much to mutate
            if(drand() < 0.5){  // revert sign with 50% probability
                delta *= -1;
            }
            if(drand() < 0.5){
                filter[i].lowerBound = roundRealValue(fmax(filter[i].lowerBound + delta, 0), precisionDigits);
            }else{
                filter[i].upperBound = roundRealValue( fmin( filter[i].upperBound + delta, 1), precisionDigits);
            }
        }
    }


}

bool mutation(Classifier *clfr, float state[])
{
    Leaf previous[numLeaf];

    for(int i=0; i<clfrCondLength; i++){
        memmove(&previous, &clfr->condition[i].leaf,sizeof(previous));
        if(!isDontcareCF(clfr->condition[i]) ) {
            do {
                memmove(&clfr->condition[i].leaf, &previous, sizeof(previous));
                apply_filter_mutation(clfr->condition[i].leaf, state);
            }while(!evaluateCF(clfr->condition[i], state));
        }
    }

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


/**
 * Mutates the condition of the classifier. If one allele is mutated depends on the constant pM.
 * This mutation is a general mutation.
 */
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
bool mutateAction(Classifier *clfr)  //Mutates the action of the classifier.
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
        while(act==clfr->action);
        clfr->action=act;
    }
    return changed;
}

// ###################### offspring insertion #################################

/**
 * Insert a discovered classifier into the population and respects the population size.
 */
void insertDiscoveredClassifier(Classifier **cl, Classifier **parents, ClassifierSet **set, ClassifierSet **pop, ClassifierSet **killset, int len, float state[])
{
    Classifier *killedp;
    len+=2;
    if(doGASubsumption)
    {
        subsumeClassifier(cl[0],parents,*set,pop,state);
        subsumeClassifier(cl[1],parents,*set,pop,state);
    }
    else
    {
        addClassifierToSet(cl[0],pop);
        addClassifierToSet(cl[1],pop);
    }

    while(len > maxPopSize)
    {
        len--;
        killedp=deleteStochClassifier(pop);

        /* record the deleted classifier to update other sets */
        if(killedp!=NULL)
        {
            addClassifierToPointerSet(killedp,killset);
            /* update the set */
            updateSet(set, *killset);
        }
    }
}

// ################################ subsumption deletion #################################

/**
 * Action set subsumption as described in the algorithmic describtion of XCS
 */
void doActionSetSubsumption(ClassifierSet **aset, ClassifierSet **pop, ClassifierSet **killset)
{
    Classifier *subsumer=NULL;
    ClassifierSet *setp, *setpl;

    /* Find the most general subsumer */
    for(setp=*aset; setp!=NULL; setp=setp->next)
    {
        if(isSubsumer(setp->classifier))
        {
            if(subsumer==NULL || isMoreGeneral(setp->classifier, subsumer))
            {
                subsumer = setp->classifier;
            }
        }
    }

    /* If a subsumer was found, subsume all classifiers that are more specific. */
    if(subsumer!=NULL)
    {
        for(setp=*aset, setpl=*aset; setp!=NULL; setp=setp->next)
        {
            while(isMoreGeneral(subsumer, setp->classifier))
            {
                subsumer->numerosity += setp->classifier->numerosity;
                if(setpl==setp)
                {
                    *aset=setp->next;
                    deleteClassifierPointerFromSet(pop,setp->classifier);
                    freeClassifier(setp->classifier);
                    addClassifierToPointerSet(setp->classifier,killset);
                    free(setp);
                    setp=*aset;
                    setpl=*aset;
                }
                else
                {
                    setpl->next=setp->next;
                    deleteClassifierPointerFromSet(pop,setp->classifier);
                    freeClassifier(setp->classifier);
                    addClassifierToPointerSet(setp->classifier,killset);
                    free(setp);
                    setp=setpl;
                }
            }
            setpl=setp;
        }
    }
}

/**
 * Tries to subsume the parents.
 */
void subsumeClassifier(Classifier *cl, Classifier **parents, ClassifierSet *locset, ClassifierSet **pop,  float state[])
{
    int i;
    for(i=0; i<2; i++)
    {
        if(parents[i]!=NULL && subsumes(parents[i],cl))
        {
            parents[i]->numerosity++;
            freeClassifier(cl);
            // code review notes: It appears that code fragments are being subsumed within a classifier
            //Here code for subsume of CFs of Parants
            //subsumeCFs(parents[i], state);
            return;
        }
    }
    // changed from action set subsumption to population subsumption
    if(subsumeClassifierToSet(cl,*pop))
    {
        return;
    }
    addClassifierToSet(cl,pop);
}

/**
 * Try to subsume in the specified set.
 */
bool subsumeClassifierToSet(Classifier *cl, ClassifierSet *set)
{
    ClassifierSet * setp;
    Classifier *subCl[maxPopSize];
    int numSub=0;

    for(setp=set; setp!=NULL; setp=setp->next)
    {
        if(subsumes(setp->classifier,cl))
        {
            subCl[numSub]=setp->classifier;
            numSub++;
        }
    }
    /* if there were classifiers found to subsume, then choose randomly one and subsume */
    if(numSub>0)
    {
        numSub = irand(numSub);
        subCl[numSub]->numerosity++;
        freeClassifier(cl);
        return true;
    }
    return false;
}

bool subsumes(Classifier *cl1, Classifier *cl2)  // check if classifier cl1 subsumes cl2
{
    return cl1->action==cl2->action && isSubsumer(cl1) && isMoreGeneral(cl1,cl2);
}

bool isSubsumer(Classifier *cl)
{
    return cl->experience > theta_sub && cl->predictionError <= epsilon_0;
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

bool is_filter_general(Leaf filter_general[], Leaf filter_to_check[])
{
    int filter_size = (int)sqrt(numLeaf);  // filter size
    for(int i=0; i<filter_size*filter_size; i++){

            if(filter_to_check[i].lowerBound<filter_general[i].lowerBound || filter_to_check[i].upperBound > filter_general[i].upperBound){
                return false;
            }
    }
    return true;
}

bool is_filter_covered_by_condition(Leaf filter_to_check[], CodeFragment code_fragments[])
{
    for(int i=0; i<clfrCondLength; i++){
        if(is_filter_general(code_fragments[i].leaf, filter_to_check)){
           return true;
        }
    }
    return false;
}

/**
 * Returns if the classifier clfr1 is more general than the classifier clfr2. It is made sure that the classifier is indeed more general and
 * not equally general as well as that the more specific classifier is completely included in the more general one (do not specify overlapping regions)
 */

// Writing by amir

 bool isMoreGeneral(Classifier *clfr1, Classifier *clfr2)
 {
     for(int i=0; i<clfrCondLength; i++)
     {
         if(!is_filter_covered_by_condition(clfr2->condition[i].leaf, clfr1->condition)){
             return false;
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
    setp->classifier=cl;
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
Classifier* deleteStochClassifier(ClassifierSet **pop)
{
    ClassifierSet *setp,*setpl;
    Classifier *killedp=NULL;
    double sum=0.0, choicep, meanf=0.0;
    int size=0;

    /* get the sum of the fitness and the numerosity */
    for(setp=*pop; setp!=NULL; setp=setp->next)
    {
        meanf+=setp->classifier->fitness;
        size+=setp->classifier->numerosity;
    }
    meanf/=(double)size;

    /* get the delete proportion, which depends on the average fitness */
    for(setp=*pop; setp!=NULL; setp=setp->next)
    {
        sum += getDelProp(setp->classifier,meanf);
    }

    /* choose the classifier that will be deleted */
    choicep=drand()*sum;

    /* look for the classifier */
    setp=*pop;
    setpl=*pop;
    sum = getDelProp(setp->classifier,meanf);
    while(sum < choicep)
    {
        setpl=setp;
        setp=setp->next;
        sum += getDelProp(setp->classifier,meanf);
    }

    /* delete the classifier */
    killedp=deleteTypeOfClassifier(setp, setpl, pop);

    /* return the pointer to the deleted classifier, to be able to update other sets */
    return killedp;
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
    /* return a pointer ot a deleted classifier (NULL if the numerosity was just decreased) */
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
