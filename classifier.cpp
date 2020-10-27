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
#include <float.h>
#include <iomanip>
#include <sys/stat.h>

int cf_gid = 0;
double predictionArray[max_actions]; //prediction array
double sumClfrFitnessInPredictionArray[max_actions]; //The sum of the fitnesses of classifiers that represent each entry in the prediction array.
int classifier_gid=0; // global incremental id to  uniquely identify the classifiers for evaluation reuse

void setInitialVariables(Classifier &clfr, double setSize, int time){
    clfr.id = classifier_gid++;
    clfr.prediction = predictionIni;
    clfr.predictionError = predictionErrorIni;
    clfr.accuracy = 0.0;
    clfr.fitness = fitnessIni;
    clfr.numerosity = 1;
    clfr.experience = 0;
    clfr.actionSetSize = 1; // chnged to 1 as per paper instead of setSize argument
    clfr.timeStamp = time;
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
void getMatchSet(ClassifierMap &pop, ClassifierSet &match_set, float *state, int itTime, int action, int img_id) {
    int population_numerosity=0, match_set_numerosity=0, representedActions;

    bool coveredActions[numActions];
    population_numerosity = get_pop_numerosity(pop);
    get_matching_classifiers(pop, state, match_set, img_id, true);
    match_set_numerosity = get_set_numerosity(match_set);

    representedActions = nrActionsInSet(match_set,coveredActions);

//    // TEMP: insert filters every time
//    if(population_numerosity < maxPopSize / 2) {
//        representedActions = 0;
//        for(int i=0; i<numActions; i++){
//            coveredActions[i] = false;
//        }
//    }

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


bool isConditionMatched(Classifier &cl, float state[], int img_id, bool train)
{
    for(int i=0; i < cl.cf.size(); i++)
    {
        //std::cout<<"iscond\n";
        //if( !isDontcareCF(clfrCond[i]) && evaluateCF(clfrCond[i].reverse_polish,state)==0 )
        if(!isDontcareCF(cl.cf[i]) && evaluateCF(cl.cf[i], state, cl.id, img_id, train) == 0 )
        {
            return false;
        }
    }
    return true;
}
void matchingCondAndSpecifiedAct(Classifier &cl, float *state, int act, int setSize, int time)
{
    createMatchingCondition(cl, state);
    cl.action = act;
    setInitialVariables(cl,setSize,time);
}

/*
 * Transfer the filters from kb to running state (i.e. the filter master list)
 */
void transfer_kb_filter(CodeFragment & cf)
{
    for(int i=0; i<cf.num_filters; i++){
        cf.filter_id[i] = add_filter(kb_filter[cf.filter_id[i]]);
    }
}


void createMatchingCondition(Classifier &cl, float *state)
{
    int count = irand(clfrCondMaxLength);
    for(int i=0; i<=count; i++) {
        add_cf(cl, state);
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

    selectTwoClassifiers(child[0], child[1] , parent[0], parent[1], action_set, fitsum, setsum); // select two classifiers (tournament selection) and copy them
    // Prediction, prediction error and fitness is only updated if crossover is done instead of always
    // (this is reverted because of slightly decreased performance)
    crossover(child[0], child[1], situation);
    for(i=0; i<2; i++)  // do mutation
    {
        mutation(child[i], situation);
    }

    child[0].prediction   = (child[0].prediction + child[1].prediction) / 2.0;
    child[0].predictionError = predictionErrorReduction * ((child[0].predictionError + child[1].predictionError) / 2.0 );
    child[0].fitness = fitnessReduction * ((child[0].fitness + child[1].fitness) / 2.0 );

    child[1].prediction = child[0].prediction;
    child[1].predictionError = child[0].predictionError;
    child[1].fitness = child[0].fitness;

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


void tournament_selection(Classifier &child, int &parent, ClassifierSet &set, double setsum)
{
    double best_fitness = -1, prediction_error=0;
    ClassifierList winner_set;

    while(winner_set.empty()) {
        for (auto &id : set.ids) {
            prediction_error = set.pop[id].predictionError;
            if (winner_set.empty() ||
                (!doGAErrorBasedSelect &&
                 best_fitness - selectTolerance <= set.pop[id].fitness / set.pop[id].numerosity) ||
                (doGAErrorBasedSelect && best_fitness + selectTolerance * maxPayoff >= prediction_error)) {
                for (int i = 0; i < set.pop[id].numerosity; i++) {
                    if (drand() < tournamentSize) {
                        if (winner_set.empty()) {
                            winner_set.push_back(id);
                            if (doGAErrorBasedSelect) {
                                best_fitness = prediction_error;
                            } else {
                                best_fitness = set.pop[id].fitness / set.pop[id].numerosity;
                            }
                        } else {
                            /* another guy in the tournament */
                            if ((!doGAErrorBasedSelect &&
                                 best_fitness + selectTolerance > set.pop[id].fitness / set.pop[id].numerosity) ||
                                (doGAErrorBasedSelect &&
                                 best_fitness - selectTolerance * maxPayoff < prediction_error)) {
                                /* both classifiers in tournament have a similar fitness/error */
                                winner_set.push_back(id);
                            } else {
                                /* new classifier in tournament is clearly better */
                                winner_set.clear();
                                winner_set.push_back(id);
                                if (doGAErrorBasedSelect) {
                                    best_fitness = prediction_error;
                                } else {
                                    best_fitness = set.pop[id].fitness / set.pop[id].numerosity;
                                }
                            }
                        }
                        break; /* go to next classifier since this one is already a winner*/
                    }
                }
            }
        }
    }
    /* choose one of the equally best winners at random */
    assert(!winner_set.empty());
    auto random_it = winner_set.begin();
    random_it = std::next(winner_set.begin(), irand(winner_set.size()));
    child = set.pop[*random_it];
    parent = *random_it;
}

void tournament_selection_(Classifier &child, int &parent, ClassifierSet &set, double setsum)
{
    int first_index = irand(set.ids.size());
    int second_index = irand(set.ids.size());
    if(first_index > second_index){
        int temp = first_index;
        first_index = second_index;
        second_index = temp;
    }
    int first=-1, second=-1;
    auto it = set.ids.begin();
    std::next(it, first_index);
    first = *it;
    std::next(it, second_index - first_index);
    second = *it;
    assert(first != -1);
    assert(second != -1);
    int selected = -1;
    if(set.pop[first].fitness > set.pop[second].fitness){
        selected = first;
    }else if(set.pop[first].fitness < set.pop[second].fitness){
        selected = second;
    }else{
        int r = irand(2);
        if(r == 0){
            selected = first;
        }else{
            selected = second;
        }
    }
    child = set.pop[selected];
    parent = selected;
}

// ########################### selection mechanism ########################################

/**
 * Select two classifiers using the chosen selection mechanism and copy them as offspring.
 */
void selectTwoClassifiers(Classifier &child1, Classifier &child2, int &parent1, int &parent2, ClassifierSet &action_set,
                          double fitsum, int setsum)
{
    tournament_selection(child1, parent1, action_set, setsum);
    tournament_selection(child2, parent2, action_set, setsum);

    child1.id = classifier_gid++;
    child1.numerosity = 1;
    child1.experience = 0;
    child1.fitness = child1.fitness / child1.numerosity;
    child2.id = classifier_gid++;
    child2.numerosity = 1;
    child2.experience = 0;
    child2.fitness = child2.fitness / child2.numerosity;

    for(int i=0; i < child1.cf.size(); i++){
        child1.cf[i].cf_id = cf_gid++;
    }
    for(int i=0; i < child2.cf.size(); i++){
        child2.cf[i].cf_id = cf_gid++;
    }

}

// ########################## crossover and mutation ########################################

void union_filter(Filter& parent1, Filter& parent2)
{
    assert(parent1.filter_size == parent2.filter_size);
    assert((parent1.is_dilated == parent2.is_dilated));
    Filter temp = parent1;

    for(int i=0; i<parent1.filter_size; i++){
        temp.lower_bounds[i] = parent1.lower_bounds[i] < parent2.lower_bounds[i] ? parent1.lower_bounds[i] : parent2.lower_bounds[i];
        temp.upper_bounds[i] = parent1.upper_bounds[i] > parent2.upper_bounds[i] ? parent1.upper_bounds[i] : parent2.upper_bounds[i];
    }

    parent1.lower_bounds = parent2.lower_bounds = temp.lower_bounds;
    parent1.upper_bounds = parent2.upper_bounds = temp.upper_bounds;
}

void crossover_filter(Filter& parent1, Filter& parent2)
{
    assert(parent1.filter_size == parent2.filter_size);
    assert((parent1.is_dilated == parent2.is_dilated));
    int point1 = irand(parent1.filter_size*parent1.filter_size);
    int point2 = irand(parent1.filter_size*parent1.filter_size);

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

/*
 * implement one point crossover
 */

bool crossover(Classifier &cl1, Classifier &cl2, float *state) {
    // crossover probability check
    if (drand() >= pX) return false;

    int size = cl1.cf.size() < cl2.cf.size() ? cl1.cf.size() : cl2.cf.size();
    if(size<=1) return false;
    int p = 0;
    do {
        p = irand(size);
    }while(p==0);
    CodeFragmentVector cfv1 = cl1.cf;
    CodeFragmentVector cfv2 = cl2.cf;
    CodeFragmentVector new_cfv1, new_cfv2;
    for(int i=0; i<p; i++){
        new_cfv1.push_back(cfv1[i]);
        new_cfv2.push_back(cfv2[i]);
    }
    for(int i=p; i<cfv2.size(); i++){
        if(!is_cf_covered(cfv2[i], new_cfv1)) {  // ensure that duplicate cf is not inserted
            new_cfv1.push_back(cfv2[i]);
        }
    }
    for(int i=p; i<cfv1.size(); i++){
        if(!is_cf_covered(cfv1[i], new_cfv2)) { // ensure that duplicate cf is not inserted
            new_cfv2.push_back(cfv1[i]);
        }
    }
    cl1.cf = new_cfv1;
    cl2.cf = new_cfv2;
}

bool crossover_old(Classifier &cl1, Classifier &cl2, float *state)
{
    // crossover probability check
    if(drand() >= pX) return false;

    bool success = false;
    int cf_index1 = irand(cl1.cf.size());
    int cf_index2 = irand(cl2.cf.size());
    CodeFragment cf1 = cl1.cf[cf_index1];
    CodeFragment cf2 = cl2.cf[cf_index2];
    // swap code fragments if they do not introduce duplication
    if(!is_cf_covered(cf1, cl2) && !is_cf_covered(cf2, cl1)){
        std::swap(cl1.cf[cf_index1], cl2.cf[cf_index2]);
        success = true;
    }
    // swap filters within code fragments if they do not introduce duplication of cf
    if(!success){
        int filter_index_1 = irand(cf1.num_filters);
        int filter_index_2 = irand(cf2.num_filters);
        // just swap the filters - just as good as crossover of filters
        std::swap(cf1.filter_id[filter_index_1], cf2.filter_id[filter_index_2]);
        if(!is_cf_covered(cf1, cl1) && !is_cf_covered(cf2, cl2)){
            cl1.cf[cf_index1] = cf1;
            cl2.cf[cf_index2] = cf2;
            success = true;
        }
    }
    if(!success){
        cf1 = cl1.cf[cf_index1];
        cf2 = cl2.cf[cf_index2];
        // crossover two filters randomly
        int filter_index_1 = irand(cl1.cf[cf_index1].num_filters);
        int filter_index_2 = irand(cl2.cf[cf_index2].num_filters);
        Filter filter1 = get_filter(cl1.cf[cf_index1].filter_id[filter_index_1]);
        Filter filter2 = get_filter(cl2.cf[cf_index2].filter_id[filter_index_2]);
        // Only do crossover if filters are of the same size or type
        if (filter1.filter_size == filter2.filter_size && filter1.is_dilated == filter2.is_dilated) {
            union_filter(filter1, filter2);
            cf1.filter_id[filter_index_1] = add_filter(filter1);
            cf2.filter_id[filter_index_2] = add_filter(filter2);
            if (evaluateCF(cf1, state) != 1) {
                negate_cf(cf1);
            }
            if (evaluateCF(cf2, state) != 1) {
                negate_cf(cf2);
            }
            if(!is_cf_covered(cf1, cl1) && !is_cf_covered(cf2, cl2)){
                cl1.cf[cf_index1] = cf1;
                cl2.cf[cf_index2] = cf2;
                success = true;
            }
        }
    }
    return success;
}


void apply_filter_mutation(Filter& filter, float state[])
{
    float delta = 0;

    for(int i=0; i<filter.filter_size*filter.filter_size; i++){
        if(drand() < pM_allel) {  // mutation based on mutation probability
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
    // mutation probability check
    if(drand() >= pM) return false;

    // There are 4 alternatives we must try all before returning false
    int alternative = irand(2);
    bool success = false;

    for(int count=0; !success && count < 2 ; count++){

        // alternative - mutate cf operator
        if(!success && (alternative + count) % 2 == 0){
            // acquire cf_index again since cf list might have changed
            int cf_index = irand(clfr.cf.size());
            CodeFragment cf = clfr.cf[cf_index];
            cf = clfr.cf[cf_index];
            if(mutate_cf(cf, state)) {
                if(!is_cf_covered(cf, clfr)){
                    clfr.cf[cf_index] = cf;
                    success = true;
                }
            }
        }
        // alternative - add new filter
        if(!success && (alternative + count) % 2 == 1){
            // acquire cf_index again since cf list might have changed
            int cf_index = irand(clfr.cf.size());
            CodeFragment cf = clfr.cf[cf_index];
            int filter_index = irand(cf.num_filters);
            cf.filter_id[filter_index] = get_new_filter(state);
            if (evaluateCF(cf, state) != 1) {
                negate_cf(cf);
            }
            if(!is_cf_covered(cf, clfr)){
                clfr.cf[cf_index] = cf;
                success = true;
            }
        }
    }

    // mutate action
    bool action_success = false; // mutateAction(clfr);
    return success || action_success;
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
        if(!subsumeClassifier(child[0], pop[parent[0]], pop[parent[1]], pop)){
            pop[child[0].id] = child[0];
        }
        if(!subsumeClassifier(child[1], pop[parent[0]], pop[parent[1]], pop)){
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
bool subsumeClassifier(Classifier &cl, Classifier &p1, Classifier &p2, ClassifierMap &pop)
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
    // as per algorithm child submsumption in population is not done
    // changed from action set subsumption to population subsumption
    if(subsumeClassifierToPop(cl, pop)){
        return true;
    }
    return false;
}


/**
 * Try to subsume in the population.
 */
bool subsumeClassifierToPop(Classifier &cl, ClassifierMap &cl_set)
{
    std::list<int> subsumers;

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
        auto random_it = subsumers.begin();
        random_it = std::next(random_it, irand(subsumers.size()));
        cl_set[*random_it].numerosity++;
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

bool is_filter_covered_by_condition(int filter_to_check_id, Classifier &cl)
{
    for(int i=0; i < cl.cf.size(); i++){
        for(int j=0; j < cl.cf[i].num_filters; j++) {
            if(filter_to_check_id == cl.cf[i].filter_id[j]){
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

// changed it to compare cfs for equality
bool isMoreGeneral(Classifier &clfr_general, Classifier &clfr_specific)
{
    bool more_general = true;
    for(int i=0; i < clfr_specific.cf.size(); i++){
        if(!is_cf_covered(clfr_specific.cf[i], clfr_general)){
            more_general = false;
            break;
        }
    }
    return more_general;
}



bool isMoreGeneral_old(Classifier &clfr_general, Classifier &clfr_specific)
 {
     for(int i=0; i < clfr_specific.cf.size(); i++){
         for(int j=0; j < clfr_specific.cf[i].num_filters; j++) {
             if (!is_filter_covered_by_condition(clfr_specific.cf[i].filter_id[j], clfr_general)) {
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


void print_population_stats(ClassifierMap &pop, std::ofstream &output_stats_file)
{
    int size = pop.size();
    //std::cout<<"\n--- Population Stats ---\n";
    //std::cout << "Global Classifier ID: " << classifier_gid << std::endl;
    //std::cout << "Population set size: " << size << std::endl;
    //std::cout << "Population numerosity size: " << get_pop_numerosity(pop) << std::endl;

    output_stats_file<<"\n--- Population Stats ---\n";
    output_stats_file << "Global Classifier ID: " << classifier_gid << std::endl;
    output_stats_file<< "Population set size: " << size << std::endl;
    output_stats_file<< "Population numerosity size: " << get_pop_numerosity(pop) << std::endl;
    int n_total = 0, n_min = INT16_MAX, n_max = -1;
    float f_total = 0, f_min = FLT_MAX, f_max = -1;
    std::for_each(pop.begin(), pop.end(),
                  [&n_total, &n_min, &n_max, &f_total, &f_min, &f_max]
                          (const ClassifierMap::value_type & item)
                  {
                      n_total+= item.second.numerosity;
                      if(n_min > item.second.numerosity) n_min = item.second.numerosity;
                      if(n_max < item.second.numerosity) n_max = item.second.numerosity;
                      f_total+= item.second.fitness;
                      if(f_min > item.second.fitness) f_min = item.second.fitness;
                      if(f_max < item.second.fitness) f_max = item.second.fitness;
                  });
    //std::cout<<"avg numerosity: "<<n_total/(float)size<<" , max numerosity: "<<n_max<<" , min numerosity: "<<n_min<<std::endl;
    //std::cout<<"avg fitness: "<<f_total/(float)size<<" , max fitness: "<<f_max<<" , min fitness: "<<f_min<<std::endl;
    //std::cout<<"--- Population Stats ---\n\n";

    output_stats_file<<"avg numerosity: "<<n_total/(float)size<<" , max numerosity: "<<n_max<<" , min numerosity: "<<n_min<<std::endl;
    output_stats_file<<"avg fitness: "<<f_total/(float)size<<" , max fitness: "<<f_max<<" , min fitness: "<<f_min<<std::endl;
    output_stats_file<<"--- Population Stats ---\n\n";
}

/**
 * This function saves the classifier population and outputs various stats.
 * This function also saves promising code fragments and filters for reuse by the subsequent experiments
 */
void save_experiment_results(ClassifierMap &pop, std::string path_postfix)
{
    std::string output_full_path = output_path + path_postfix;
    mkdir(output_full_path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    std::ofstream output_classifier_file;
    output_classifier_file.open(output_full_path + output_classifier_file_name);
    if(!output_classifier_file.is_open()){
        std::cout << "Could not open output classifier file";
        exit(1);
    }
    std::ofstream output_code_fragment_file;
    output_code_fragment_file.open(output_full_path + output_code_fragment_file_name);
    if(!output_code_fragment_file.is_open()){
        std::cout << "Could not open output code fragment file";
        exit(1);
    }
    std::ofstream output_promising_code_fragment_file;
    output_promising_code_fragment_file.open(output_full_path + output_promising_code_fragment_file_name);
    if(!output_promising_code_fragment_file.is_open()){
        std::cout << "Could not open output promising code fragment file";
        exit(1);
    }
    std::ofstream output_filter_file;
    output_filter_file.open(output_full_path + output_filter_file_name);
    if(!output_filter_file.is_open()){
        std::cout << "Could not open output code filter file";
        exit(1);
    }
    std::ofstream output_promising_filter_file;
    output_promising_filter_file.open(output_full_path + output_promising_filter_file_name);
    if(!output_promising_filter_file.is_open()){
        std::cout << "Could not open output code filter file";
        exit(1);
    }
    std::ofstream output_stats_file;
    output_stats_file.open(output_full_path + output_stats_file_name);
    if(!output_stats_file.is_open()){
        std::cout << "Could not open output stats file";
        exit(1);
    }
    print_population_stats(pop, output_stats_file);
    print_filter_stats(output_stats_file);
    print_filter_evaluation_stats(output_stats_file);
    for(auto& item : pop)
    {
        fprintClassifier(item.second, output_classifier_file, output_code_fragment_file,
                output_promising_code_fragment_file, output_promising_filter_file);
    }
    output_filters(output_filter_file, output_promising_filter_file);
    //storeCFs(pop, fpCF);
    output_classifier_file.close();
    output_code_fragment_file.close();
    output_promising_code_fragment_file.close();
    output_filter_file.close();
    output_promising_filter_file.close();
    output_stats_file.close();
}

/**
 * print a single classifier to the file fp
 */
/*
void fprintClassifier(FILE *fp, Classifier *classifier){
	char buf[1000];
	for(int i=0; i<clfrCondMaxLength; i++){
		outprog(classifier->code_fragment[i].reverse_polish,cfMaxLength,fp);
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



void fprintClassifier(Classifier &classifier, std::ofstream &output_classifier_file,
                      std::ofstream &output_code_fragment_file, std::ofstream &output_promising_code_fragment_file,
                      std::ofstream &output_promising_filter_file)
{
    output_classifier_file << "id ";
    output_classifier_file.width(5);
    output_classifier_file << classifier.id;
    output_classifier_file << " num ";
    output_classifier_file << classifier.numerosity;
    output_classifier_file << " exp ";
    output_classifier_file.width(5);
    output_classifier_file << classifier.experience;
    output_classifier_file << " num_cf ";
    output_classifier_file << classifier.cf.size();
    output_classifier_file << " fitness ";
    output_classifier_file.width(11);
    output_classifier_file << classifier.fitness;
    output_classifier_file << " accuracy ";
    output_classifier_file.width(11);
    output_classifier_file << classifier.accuracy;
    output_classifier_file << " prediction ";
    output_classifier_file.width(11);
    output_classifier_file << classifier.prediction;
    output_classifier_file << " error ";
    output_classifier_file.width(11);
    output_classifier_file << classifier.predictionError;
    output_classifier_file << " action_set_size ";
    output_classifier_file.width(7);
    output_classifier_file << classifier.actionSetSize;
    output_classifier_file << " time_stamp ";
    output_classifier_file.width(5);
    output_classifier_file << classifier.timeStamp;
    output_classifier_file << " action ";
    output_classifier_file << classifier.action << std::endl;

    for(auto & cf : classifier.cf)
    {
        output_classifier_file << cf.cf_id << " ";
        output_code_fragment_to_file(cf, output_code_fragment_file);
        // output promising code fragments separately
        if(is_promising_classifier(classifier)){
            output_code_fragment_to_file(cf, output_promising_code_fragment_file);
        }
    }
    output_classifier_file << std::endl;
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


inline bool is_promising_classifier(Classifier& cl)
{
    cl.predictionError < epsilon_0 && cl.experience > theta_promising;

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
        for(int i=0; i < item.second.cf.size(); i++){
            for(int j=0; j < item.second.cf[i].num_filters; j++){
                Filter& f = get_filter(item.second.cf[i].filter_id[j]);
                f.numerosity++;
                // if classifier is "promising" then increase the fitness of the filter
                // a promising classifier is one whose error < 10 and experience > 10
                if(is_promising_classifier(item.second)){
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
