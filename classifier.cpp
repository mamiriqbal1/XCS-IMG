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

void setInitialVariables(Classifier &clfr, double setSize, int time){
    clfr.id = gid++;
    clfr.prediction = predictionIni;
    clfr.predictionError = predictionErrorIni;
    clfr.accuracy = 0.0;
    clfr.fitness = fitnessIni;
    clfr.numerosity = 1;
    clfr.experience = 0;
    clfr.actionSetSize = 1; // chnged to 1 as per paper instead of setSize argument
    clfr.timeStamp = time;
    clfr.specificness = numberOfNonDontcares(clfr.condition);
}


double getAvgFitness(ClassifierMap &pop){
    double fitsum=0.0;
    int setsum=0;

    for(auto &item : pop)
    {
        fitsum += item.second.fitness;
        setsum += item.second.numerosity;
    }
    return fitsum/(double)setsum;
}

int getNumFitterCFs(ClassifierMap &pop, double avgFitness){
    int numCFs = getNumPreviousCFs();
    for(auto & item : pop){
        if(item.second.fitness > avgFitness){
            numCFs += item.second.specificness;
        }
    }
    return numCFs;
}


void remove_classifier(ClassifierSet &set, int cl_id) {
    for(auto it = set.ids.begin(); it != set.ids.end(); it++){
        if(*it == cl_id){
            set.ids.erase(it);
            return;
        }
    }
}

int get_set_numerosity(ClassifierSet &set)
{
    int sum = 0;
    for(auto& id : set.ids){
        sum += set.pop[id].numerosity;
    }
    return sum;
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
getMatchSet(ClassifierMap &pop, ClassifierSet &match_set, float *state, int itTime, int action, int img_id) {
    int population_numerosity=0, match_set_numerosity=0, representedActions;

    bool coveredActions[numActions];
    population_numerosity = get_pop_numerosity(pop);
    get_matching_classifiers(pop, state, match_set, img_id, true);
    match_set_numerosity = get_set_numerosity(match_set);

    representedActions = nrActionsInSet(match_set,coveredActions);

    // TEMP: insert filters every time
    if(population_numerosity < maxPopSize / 2) {
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
                    Classifier coverClfr;
                    matchingCondAndSpecifiedAct(coverClfr, state, i, match_set_numerosity + 1,
                                                            itTime);
                    // before inserting the new classifier into the population check for subsumption by a generic one
                    if(!subsumeClassifierToPop(coverClfr, pop)) {
                        pop[coverClfr.id] = coverClfr;
                        match_set.ids.push_back(coverClfr.id);
                    }
                    match_set_numerosity++;
                    population_numerosity++;
                }
            }
        }

        /* Delete classifier if population is too big and record it in killset */
        while(population_numerosity > maxPopSize )
        {
            /* PL */
            int cl_id = deleteStochClassifier(pop);
            // also remove from match set
            remove_classifier(match_set, cl_id);
            population_numerosity--;
        }
        representedActions = nrActionsInSet(match_set,coveredActions);
    }
}
/**
 * Returns the number of actions in the set and stores which actions are covered in the array coveredActions.
 */
int nrActionsInSet(ClassifierSet &match_set, bool *coveredActions)
{
    int nr=0;
    for(int i=0; i<numActions; i++){
        coveredActions[i] = false;
    }
    for(auto & id : match_set.ids){
        if(!coveredActions[match_set.pop[id].action]){
            coveredActions[match_set.pop[id].action] = true;
            nr++;
            if(nr >= numActions) break;
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
void matchingCondAndSpecifiedAct(Classifier &cl, float *state, int act, int setSize,
                                 int time)  //matchingCondAndSpecifiedAct(currentState,i,matchSetNumerositySum+1,time)
{
    createMatchingCondition(cl.condition,state);
    cl.action = act;
    setInitialVariables(cl,setSize,time);
}
void createMatchingCondition(CodeFragment cond[], float state[])
{
    for(int i=0; i<clfrCondLength; i++)
    {
        do
        {
            createNewCF(getNumPreviousCFs() + countNewCFs, cond[i]);
            opType* end = randomProgram(cond[i].codeFragment,irand(2),cfMaxDepth,cfMinDepth);
            validateDepth(cond[i].codeFragment,end); //validate depth
            cond[i] = addLeafCF(cond[i], state);
        }
        while( evaluateCF(cond[i],state)!=1 );
        countNewCFs++;
    }
}

// ######################### prediction array operations ############################################

void getPredictionArray(ClassifierSet &match_set)  //determines the prediction array out of the match set ms
{
    for(int i=0; i<numActions; i++)
    {
        predictionArray[i]=0.0;
        sumClfrFitnessInPredictionArray[i]=0.0;
    }
    for(auto& id : match_set.ids)
    {
        int actionValue = match_set.pop[id].action;
        predictionArray[actionValue]+= match_set.pop[id].prediction * match_set.pop[id].fitness;
        sumClfrFitnessInPredictionArray[actionValue]+= match_set.pop[id].fitness;
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

// ######################## action set operations #########################################

void getActionSet(int action, ClassifierSet &match_set,
                   ClassifierSet &action_set)  // constructs an action set out of the match set ms.
{
    for(auto& id : match_set.ids){
        if(action == match_set.pop[id].action){
            action_set.ids.push_back(id);
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
 void updateActionSet(ClassifierSet &action_set, double maxPrediction, double reward)
{
    double P, setsize=0.0;

    P = reward + gama*maxPrediction;

    for(auto& id : action_set.ids)
    {
        setsize += action_set.pop[id].numerosity;
        action_set.pop[id].experience++;
    }

    for(auto& id : action_set.ids)   // update prediction, prediction error and action set size estimate
    {
        if((double)action_set.pop[id].experience < 1.0 / beta)
        {
            // !first adjustments! -> simply calculate the average
            action_set.pop[id].predictionError = (action_set.pop[id].predictionError * ((double)action_set.pop[id].experience - 1.0) + absoluteValue(P - action_set.pop[id].prediction)) / (double)action_set.pop[id].experience;
            action_set.pop[id].prediction = (action_set.pop[id].prediction * ((double)action_set.pop[id].experience - 1.0) + P) / (double)action_set.pop[id].experience;
            action_set.pop[id].actionSetSize = (action_set.pop[id].actionSetSize * ((double)(action_set.pop[id].experience - 1)) + setsize) / (double)action_set.pop[id].experience;
        }
        else
        {
            // normal adjustment -> use widrow hoff delta rule
            action_set.pop[id].predictionError += beta * (absoluteValue(P - action_set.pop[id].prediction) - action_set.pop[id].predictionError);
            action_set.pop[id].prediction += beta * (P - action_set.pop[id].prediction);
            action_set.pop[id].actionSetSize += beta * (setsize - action_set.pop[id].actionSetSize);
        }
    }
    updateFitness(action_set);
    if(doActSetSubsumption)
    {
        doActionSetSubsumption(action_set);
    }
}

// update the fitnesses of an action set (the previous [A] in multi-step envs or the current [A] in single-step envs.)
/*
 * code review notes
 * correctly implementation
 */
void updateFitness(ClassifierSet &action_set)
{
    double ksum=0.0;


    //First, calculate the accuracies of the classifier and the accuracy sums
    for(auto& id : action_set.ids){
        if(action_set.pop[id].predictionError <= epsilon_0){
            action_set.pop[id].accuracy = 1.0;
        }
        else{
            action_set.pop[id].accuracy = alpha * pow(action_set.pop[id].predictionError / epsilon_0 , -nu);
        }
        ksum += action_set.pop[id].accuracy * (double)action_set.pop[id].numerosity;
    }

    //Next, update the fitnesses accordingly
    for(auto& id : action_set.ids){
        action_set.pop[id].fitness += beta * ((action_set.pop[id].accuracy * action_set.pop[id].numerosity) / ksum - action_set.pop[id].fitness );
    }
}

// ############################ discovery mechanism #########################################

/**
 * The discovery conmponent with the genetic algorithm
 * note: some classifiers in set could be deleted !
 */
void discoveryComponent(ClassifierSet &action_set, ClassifierMap &pop, int itTime, float *situation)
{
    Classifier child[2];
    int parent[2];
    double fitsum=0.0;
    int i, len, setsum=0, gaitsum=0;
    // if the classifier set is empty, return (due to deletion)
    if(action_set.ids.size() == 0) return;

    getDiscoversSums(action_set, &fitsum, &setsum, &gaitsum); // get all sums that are needed to do the discovery

    // do not do a GA if the average number of time-steps in the set since the last GA is less or equal than thetaGA
    if( itTime - (double)gaitsum / (double)setsum < theta_GA)
    {
        return;
    }
    setTimeStamps(action_set, itTime);

    selectTwoClassifiers(child, parent, action_set, fitsum, setsum); // select two classifiers (tournament selection) and copy them
    // Prediction, prediction error and fitness is only updated if crossover is done instead of always
    // (this is reverted because of slightly decreased performance)
    if(crossover(child[0], child[1], situation) || true){
        child[0].prediction   = (child[0].prediction + child[1].prediction) / 2.0;
        child[0].predictionError = predictionErrorReduction * ((child[0].predictionError + child[1].predictionError) / 2.0 );
        child[0].fitness = fitnessReduction * ((child[0].fitness + child[1].fitness) / 2.0 );

        child[1].prediction = child[0].prediction;
        child[1].predictionError = child[0].predictionError;
        child[1].fitness = child[0].fitness;
    }
    for(i=0; i<2; i++)  // do mutation
    {
        mutation(child[i], situation);
    }

    child[0].specificness = numberOfNonDontcares(child[0].condition);
    child[1].specificness = numberOfNonDontcares(child[1].condition);

    // get the length of the population to check if clasifiers have to be deleted
    len = get_pop_numerosity(pop);

    // insert the new two classifiers and delete two if necessary
    insertDiscoveredClassifier(child, parent, pop, len);
}
void getDiscoversSums(ClassifierSet &action_set, double *fitsum, int *setsum, int *gaitsum)  // Calculate all necessary sums in the set for the discovery component.
{
    *fitsum=0.0;
    *setsum=0;
    *gaitsum=0;
    for(auto& id : action_set.ids)
    {
        (*fitsum)+=action_set.pop[id].fitness;
        (*setsum)+=action_set.pop[id].numerosity;
        (*gaitsum) += action_set.pop[id].timeStamp * action_set.pop[id].numerosity;
    }
}
void setTimeStamps(ClassifierSet &action_set, int itTime)  // Sets the time steps of all classifiers in the set to itTime (because a GA application is occurring in this set!).
{
    for(auto& id : action_set.ids)
    {
        action_set.pop[id].timeStamp = itTime;
    }
}

void tournament_selection(Classifier *child, int *parent, ClassifierSet &set, double setsum)
{
    int first = irand(set.ids.size());
    int second = irand(set.ids.size());
    Classifier* first_classifier = nullptr, *second_classifier = nullptr;
    int i=0;
    for(auto& id : set.ids){
        if(i == first){
            first_classifier = &set.pop[id];
        }
        if(i == second){
            second_classifier = &set.pop[id];
        }
        i++;
    }
    assert(first_classifier != nullptr);
    assert(second_classifier != nullptr);
    if(first_classifier->fitness > second_classifier->fitness){
        *child = *first_classifier;
        *parent = first_classifier->id;
    }else if(first_classifier->fitness < second_classifier->fitness){
        *child = *second_classifier;
        *parent = second_classifier->id;
    }else{
        int r = irand(2);
        if(r == 0){
            *child = *first_classifier;
            *parent = first_classifier->id;
        }else{
            *child = *second_classifier;
            *parent = second_classifier->id;
        }
    }
}

// ########################### selection mechanism ########################################

/**
 * Select two classifiers using the chosen selection mechanism and copy them as offspring.
 */
void selectTwoClassifiers(Classifier *cl, int *parents, ClassifierSet &action_set, double fitsum, int setsum)
{
    // todo parents need not be copied externally???
    tournament_selection(&cl[0], &parents[0], action_set, setsum);
    tournament_selection(&cl[1], &parents[1], action_set, setsum);

    for(int i=0; i<2; i++) {
        cl[i].id = gid++;
        cl[i].numerosity = 1;
        cl[i].experience = 0;
    }
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
void insertDiscoveredClassifier(Classifier *child, int *parent, ClassifierMap &pop, int len)
{
    len+=2;
    if(doGASubsumption)
    {
        if(!subsumeClassifier(child[0], pop[parent[0]], pop[parent[1]])){
            pop[child[0].id] = child[0];
        }
        if(!subsumeClassifier(child[1], pop[parent[0]], pop[parent[1]])){
            pop[child[1].id] = child[1];
        }
    }
    else
    {
        pop[child[0].id] = child[0];
        pop[child[1].id] = child[1];
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
void doActionSetSubsumption(ClassifierSet &action_set)
{
    int subsumer = -1;

    /* Find the most general subsumer */
    for(auto& id : action_set.ids)
    {
        if(isSubsumer(action_set.pop[id]))
        {
            if(subsumer== -1 || isMoreGeneral(action_set.pop[id], action_set.pop[subsumer]))
            {
                subsumer = id;
            }
        }
    }

    /* If a subsumer was found, subsume all classifiers that are more specific. */
    if(subsumer!= -1)
    {
        for(auto& id : action_set.ids){
            if(isMoreGeneral(action_set.pop[subsumer], action_set.pop[id]))
                action_set.pop[subsumer].numerosity += action_set.pop[id].numerosity;
                remove_classifier(action_set, id);
                action_set.pop.erase(id);
        }

    }
}

/**
 * Tries to subsume the parents.
 */
bool subsumeClassifier(Classifier &cl, Classifier &p1, Classifier &p2)
{
    int i;
    if(subsumes(p1, cl))
    {
        p1.numerosity++;
        return true;
    }
    if(subsumes(p2, cl))
    {
        p2.numerosity++;
        return true;
    }
    return false;
    // as per algorithm child submsumption in population is not done
    // changed from action set subsumption to population subsumption
//    if(subsumeClassifierToPop(cl, pop))
//    {
//        return;
//    }
//    pop.push_front(cl);
}


/**
 * Try to subsume in the population.
 */
bool subsumeClassifierToPop(Classifier &cl, ClassifierMap &cl_set)
{
    std::vector<int> subsumers;
    subsumers.reserve(maxPopSize);

    for(auto & item : cl_set)
    {
        if(subsumes(item.second, cl))
        {
            subsumers.push_back(item.second.id);
        }
    }
    /* if there were classifiers found to subsume, then choose randomly one and subsume */
    if(subsumers.size() > 0)
    {
        int k = irand(subsumers.size());
        cl_set[subsumers[k]].numerosity++;
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

// ###################### adding classifiers to a set ###################################

// ############################## deletion ############################################

/**
 * Deletes one classifier in the population.
 * The classifier that will be deleted is chosen by roulette wheel selection
 * considering the deletion vote. Returns position of the macro-classifier which got decreased by one micro-classifier.
 **/
int deleteStochClassifier(ClassifierMap &pop)
{
    double vote_sum=0.0, choicep, meanf=0.0;
    int size=0;

    for(auto& item : pop){
        meanf += item.second.fitness;
        size += item.second.numerosity;
    }
    meanf/=(double)size;

    /* get the delete proportion, which depends on the average fitness */
    for(auto& item : pop){
        vote_sum += getDelProp(item.second, meanf);
    }

    /* choose the classifier that will be deleted */
    choicep= drand() * vote_sum;
    /* look for the classifier */
    vote_sum = 0;
    int removed_id = -1;
    for(auto& item : pop){
        vote_sum += getDelProp(item.second, meanf);
        if(vote_sum > choicep){
            removed_id = item.second.id;
            if(item.second.numerosity > 1){
                item.second.numerosity--;
            }else{
                pop.erase(item.second.id);
            }
            return removed_id;
        }
    }
}


double getDelProp(Classifier &clfr, double meanFitness)  //Returns the vote for deletion of the classifier.
{
    if(clfr.fitness/(double)clfr.numerosity >= delta*meanFitness || clfr.experience < theta_del)
    {
        return (double)(clfr.actionSetSize*clfr.numerosity);
    }
    else
    {
        return (double)clfr.actionSetSize*(double)clfr.numerosity*meanFitness / (clfr.fitness/(double)clfr.numerosity);
    }
}

//############# concrete deletion of a classifier or a whole classifier set ############

// ############################### output operations ####################################



/**
 * print the classifier in a delete_ClassifierSet to the file fp
 */
void fprintClassifierSet(FILE *fpClfr, FILE *fpCF, ClassifierMap &pop)
{
    print_filter_stats();
    print_filter_evaluation_stats();
    for(auto& item : pop)
    {
        fprintClassifier(fpClfr, item.second);
    }
    std::cout << "Global Classifier ID: " << gid << std::endl;
    storeCFs(pop, fpCF);
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

void fprintClassifier(FILE *fp, Classifier &classifier)
{
    char *buf;
    int len;
    for(int i=0; i<clfrCondLength; i++)
    {
        //outprog(classifier->condition[i].codeFragment,cfMaxLength,fp);
        outprog(classifier.condition[i],cfMaxLength,fp);
        fwrite("\n",strlen("\n"),1,fp);
    }

    len = snprintf(NULL,0,"Action: %d\n",classifier.action);
    if(!(buf = (char*)malloc((len + 1) * sizeof(char))))
    {
        printf("\nError in file writing ...\n");
        exit(0);
    }
    len = snprintf(buf,len+1,"Action: %d\n",classifier.action);
    fwrite(buf,strlen(buf),1,fp);
    free(buf);

    len = snprintf(NULL,0,"id: %d ",classifier.id);
    if(!(buf = (char*)malloc((len + 1) * sizeof(char))))
    {
        printf("\nError in file writing ...\n");
        exit(0);
    }
    len = snprintf(buf,len+1,"id: %d ",classifier.id);
    fwrite(buf,strlen(buf),1,fp);
    free(buf);

    len = snprintf(NULL,0,"Numerosity: %d ",classifier.numerosity);
    if(!(buf = (char*)malloc((len + 1) * sizeof(char))))
    {
        printf("\nError in file writing ...\n");
        exit(0);
    }
    len = snprintf(buf,len+1,"Numerosity: %d ",classifier.numerosity);
    fwrite(buf,strlen(buf),1,fp);
    free(buf);

    len = snprintf(NULL,0,"Accuracy: %f ",classifier.accuracy);
    if(!(buf = (char*)malloc((len + 1) * sizeof(char))))
    {
        printf("\nError in file writing ...\n");
        exit(0);
    }
    len = snprintf(buf,len+1,"Accuracy: %f ",classifier.accuracy);
    fwrite(buf,strlen(buf),1,fp);
    free(buf);

    len = snprintf(NULL,0,"Fitness: %f ",classifier.fitness);
    if(!(buf = (char*)malloc((len + 1) * sizeof(char))))
    {
        printf("\nError in file writing ...\n");
        exit(0);
    }
    len = snprintf(buf,len+1,"Fitness: %f ",classifier.fitness);
    fwrite(buf,strlen(buf),1,fp);
    free(buf);

    len = snprintf(NULL,0,"PredictionError: %f ",classifier.predictionError);
    if(!(buf = (char*)malloc((len + 1) * sizeof(char))))
    {
        printf("\nError in file writing ...\n");
        exit(0);
    }
    len = snprintf(buf,len+1,"PredictionError: %f ",classifier.predictionError);
    fwrite(buf,strlen(buf),1,fp);
    free(buf);

    len = snprintf(NULL,0,"Prediction: %f ",classifier.prediction);
    if(!(buf = (char*)malloc((len + 1) * sizeof(char))))
    {
        printf("\nError in file writing ...\n");
        exit(0);
    }
    len = snprintf(buf,len+1,"Prediction: %f ",classifier.prediction);
    fwrite(buf,strlen(buf),1,fp);
    free(buf);

    len = snprintf(NULL,0,"Experience: %d ",classifier.experience);
    if(!(buf = (char*)malloc((len + 1) * sizeof(char))))
    {
        printf("\nError in file writing ...\n");
        exit(0);
    }
    len = snprintf(buf,len+1,"Experience: %d ",classifier.experience);
    fwrite(buf,strlen(buf),1,fp);
    free(buf);

    len = snprintf(NULL,0,"Specificness: %d ",classifier.specificness);
    if(!(buf = (char*)malloc((len + 1) * sizeof(char))))
    {
        printf("\nError in file writing ...\n");
        exit(0);
    }
    len = snprintf(buf,len+1,"Specificness: %d ",classifier.specificness);
    fwrite(buf,strlen(buf),1,fp);
    free(buf);

    len = snprintf(NULL,0,"ActionSetSize: %f ",classifier.actionSetSize);
    if(!(buf = (char*)malloc((len + 1) * sizeof(char))))
    {
        printf("\nError in file writing ...\n");
        exit(0);
    }
    len = snprintf(buf,len+1,"ActionSetSize: %f ",classifier.actionSetSize);
    fwrite(buf,strlen(buf),1,fp);
    free(buf);

    len = snprintf(NULL,0,"TimeStamp: %d\n",classifier.timeStamp);
    if(!(buf = (char*)malloc((len + 1) * sizeof(char))))
    {
        printf("\nError in file writing ...\n");
        exit(0);
    }
    len = snprintf(buf,len+1,"TimeStamp: %d\n",classifier.timeStamp);
    fwrite(buf,strlen(buf),1,fp);
    free(buf);

    fflush(fp);
}

/*###################### sorting the classifier list ###################################*/

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

void manage_filter_list(ClassifierMap &pop){
    // reset statistics of all filters before updating
    reset_filter_stats();

    for(auto& item : pop){
        for(int i=0; i<clfrCondLength; i++){
            for(int j=0; j < item.second.condition[i].num_filters; j++){
                Filter& f = get_filter(item.second.condition[i].filter_id[j]);
                f.numerosity++;
                // if classifier is "promising" then increase the fitness of the filter
                // a promising classifier is one whose error < 10 and experience > 10
                if(item.second.predictionError < epsilon_0 && item.second.experience > theta_filter){
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

int get_pop_numerosity(ClassifierMap &pop) {
    int pop_numerosity = 0;
    std::for_each(pop.begin(), pop.end(),
                  [&pop_numerosity](ClassifierMap::value_type& item)
                  {
                      pop_numerosity+= item.second.numerosity;
                  });
    return pop_numerosity;
}

void get_matching_classifiers(ClassifierMap &pop, float *state, ClassifierSet &match_set, int img_id, bool train) {

    std::for_each(pop.begin(), pop.end(), [&match_set, &state, img_id, train](ClassifierMap::value_type& item)
    {
        if(isConditionMatched(item.second, state, img_id, train)){
            match_set.ids.push_back(item.first);
        }
    });
}
