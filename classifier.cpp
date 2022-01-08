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
#include <stack>
#include "cf_list.h"

CharMatrix evaluation_cache;

ClassifierVector population; //(maxPopSize + 10);
int classifier_gid=0; // global incremental id to  uniquely identify the classifiers for evaluation reuse

DoubleVector predictionArray; // [max_actions]; //prediction array
DoubleVector sumClfrFitnessInPredictionArray; //[max_actions]; //The sum of the fitnesses of classifiers that represent each entry in the prediction array.

std::vector<int> cl_gid_vector;
std::stack<int, std::vector<int>> classifier_gid_stack(cl_gid_vector);
//std::stack<int, std::vector<int>> classifier_gid_stack;


int get_next_cl_gid()
{
    int result = -1;
    if(classifier_gid_stack.size()>0){
        int val = classifier_gid_stack.top();
        classifier_gid_stack.pop();
        result = val;
    }else{
        // only grow the vector when a new element is needed
        population.resize(classifier_gid+1);
        result = classifier_gid++;
    }
    // reset evaluation cache
//    evaluation_cache[result].assign(trainNumInstances, UNKNOWN);
    std::fill(evaluation_cache[result].begin(), evaluation_cache[result].end(), UNKNOWN);
    return result;
}


void initialize_parameters()
{
    DoubleVector prediction_array(numActions);
    predictionArray = prediction_array;
    DoubleVector sum(numActions);
    sumClfrFitnessInPredictionArray = sum;
    // if class map not initialized through config
    if(class_map.empty()) {
        for (int i = 0; i < numActions; i++) {
            class_map[i] = i;
        }
    }
}


void initialize_population(int size)
{
    population.reserve(size);
    // create a matrix with rows=population size and columns = number of images
    evaluation_cache = CharMatrix(size, CharVector (trainNumInstances, UNKNOWN));
}

void setInitialVariables(Classifier &clfr, double setSize, int time){
//    clfr.id = get_next_cl_gid();  // it will be set just before adding to population
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
        sum += population[id].numerosity;
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
 void getMatchSet(ClassifierSet &match_set, float *state, int itTime, int action, int img_id) {
    int population_numerosity=0, match_set_numerosity=0, representedActions;

    bool coveredActions[numActions];
    population_numerosity = get_pop_size(true);
     get_matching_classifiers(state, match_set, img_id, true);
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
                     if(!subsumeClassifierToPop(coverClfr)) {
                         coverClfr.id = get_next_cl_gid();
                         population[coverClfr.id] = coverClfr;
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
            int cl_id = deleteStochClassifier(population);
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
        if(!coveredActions[population[id].action]){
            coveredActions[population[id].action] = true;
            nr++;
            if(nr >= numActions) break;
        }
    }
    return nr;
}


bool isConditionMatched(Classifier &cl, float state[], int img_id, bool train)
{
    bool result = false;
    if(train && evaluation_cache[cl.id][img_id] == NOT_MATCHED) result = false;
    else if(train && evaluation_cache[cl.id][img_id] == MATCHED)  result = true;
    else {
        bool matched = true;
        for (int i = 0; i < clfrCondMaxLength && matched; i++) {
            if (cl.cf_ids[i] != -1 && evaluate_cf_slide(get_cf(cl.cf_ids[i]), state, cl.id, img_id, train) == 0) {
                matched = false;
            }
        }
        if(matched)  result = true;
        else  result = false;

        if(train){
            if(matched)  evaluation_cache[cl.id][img_id] = MATCHED;
            else  evaluation_cache[cl.id][img_id] = NOT_MATCHED;
        }
    }
    return result;
}
void matchingCondAndSpecifiedAct(Classifier &cl, float *state, int act, int setSize, int time)
{
    createMatchingCondition(cl, state);
    cl.action = act;
    setInitialVariables(cl,setSize,time);
}

//todo: create_new_cf should be done to temporary location and should be added to main cf_list only when cl is added to pop

/*
 * Default cf has id == -1 that represent don't care
 */
void createMatchingCondition(Classifier &cl, float *state)
{
    bool cf_added = false; // to ensure that at least one cf is added
    do {
        for (int i = 0; i < clfrCondMaxLength; i++) {
            if (drand() >= P_dontcare) {
                cf_added = true;
                CodeFragment new_cf;
                cl.cf_ids[i] = create_new_cf(state);
            }
        }
    }while(!cf_added);
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
        int actionValue = population[id].action;
        predictionArray[actionValue]+= population[id].prediction * population[id].fitness;
        sumClfrFitnessInPredictionArray[actionValue]+= population[id].fitness;
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
        if(action == population[id].action){
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
  * code is as per paper
  */
 void updateActionSet(ClassifierSet &action_set, double maxPrediction, double reward)
{
    double P, action_set_numerosity=0.0;

    P = reward + gama*maxPrediction;

    for(auto& id : action_set.ids)
    {
        action_set_numerosity += population[id].numerosity;
        population[id].experience++;
    }

    for(auto& id : action_set.ids)   // update prediction, prediction error and action set size estimate
    {
        if((double)population[id].experience < 1.0 / beta)
        {
            // !first adjustments! -> simply calculate the average
            population[id].predictionError += (absoluteValue(P - population[id].prediction) - population[id].predictionError) / (double)population[id].experience;
            population[id].prediction += (P - population[id].prediction) / (double)population[id].experience;
            population[id].actionSetSize += (action_set_numerosity - population[id].actionSetSize) / (double)population[id].experience;
        }
        else
        {
            // normal adjustment -> use widrow hoff delta rule
            population[id].predictionError += beta * (absoluteValue(P - population[id].prediction) - population[id].predictionError);
            population[id].prediction += beta * (P - population[id].prediction);
            population[id].actionSetSize += beta * (action_set_numerosity - population[id].actionSetSize);
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
        if(population[id].predictionError <= epsilon_0){
            population[id].accuracy = 1.0;
        }
        else{
            population[id].accuracy = alpha * pow(population[id].predictionError / epsilon_0 , -nu);
        }
        ksum += population[id].accuracy * (double)population[id].numerosity;
    }

    //Next, update the fitnesses accordingly
    for(auto& id : action_set.ids){
        population[id].fitness += beta * ((population[id].accuracy * population[id].numerosity) / ksum - population[id].fitness );
    }
}

// ############################ discovery mechanism #########################################

/**
 * The discovery conmponent with the genetic algorithm
 * note: some classifiers in set could be deleted !
 */
void discoveryComponent(ClassifierSet &action_set, int itTime, float *situation, int action)
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

//    add_new_classifiers_to_population(situation, action, itTime);
//    return;

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
    len = get_pop_size(true);

    // insert the new two classifiers and delete two if necessary
    insertDiscoveredClassifier(child, parent, action_set, len);
}
void getDiscoversSums(ClassifierSet &action_set, double *fitsum, int *setsum, int *gaitsum)  // Calculate all necessary sums in the set for the discovery component.
{
    *fitsum=0.0;
    *setsum=0;
    *gaitsum=0;
    for(auto& id : action_set.ids)
    {
        (*fitsum)+=population[id].fitness;
        (*setsum)+=population[id].numerosity;
        (*gaitsum) += population[id].timeStamp * population[id].numerosity;
    }
}
void setTimeStamps(ClassifierSet &action_set, int itTime)  // Sets the time steps of all classifiers in the set to itTime (because a GA application is occurring in this set!).
{
    for(auto& id : action_set.ids)
    {
        population[id].timeStamp = itTime;
    }
}


void tournament_selection(Classifier &child, int &parent, ClassifierSet &set, double setsum)
{
    double best_fitness = -1, prediction_error=0;
    ClassifierIDVector winner_set;

    while(winner_set.empty()) {
        for (auto &id : set.ids) {
            prediction_error = population[id].predictionError;
            if (winner_set.empty() ||
                (!doGAErrorBasedSelect &&
                 best_fitness - selectTolerance <= population[id].fitness / population[id].numerosity) ||
                (doGAErrorBasedSelect && best_fitness + selectTolerance * maxPayoff >= prediction_error)) {
                for (int i = 0; i < population[id].numerosity; i++) {
                    if (drand() < tournamentSize) {
                        if (winner_set.empty()) {
                            winner_set.push_back(id);
                            if (doGAErrorBasedSelect) {
                                best_fitness = prediction_error;
                            } else {
                                best_fitness = population[id].fitness / population[id].numerosity;
                            }
                        } else {
                            /* another guy in the tournament */
                            if ((!doGAErrorBasedSelect &&
                                 best_fitness + selectTolerance > population[id].fitness / population[id].numerosity) ||
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
                                    best_fitness = population[id].fitness / population[id].numerosity;
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
    child = population[*random_it];
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
    if(population[first].fitness > population[second].fitness){
        selected = first;
    }else if(population[first].fitness < population[second].fitness){
        selected = second;
    }else{
        int r = irand(2);
        if(r == 0){
            selected = first;
        }else{
            selected = second;
        }
    }
    child = population[selected];
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

//    child1.id = get_next_cl_gid();
    child1.numerosity = 1;
    child1.experience = 0;
    child1.fitness = child1.fitness / child1.numerosity;
//    child2.id = get_next_cl_gid();
    child2.numerosity = 1;
    child2.experience = 0;
    child2.fitness = child2.fitness / child2.numerosity;
}


/*
 * implement two point crossover
 */

bool crossover(Classifier &cl1, Classifier &cl2, float *state) {
    // crossover probability check
    if (drand() >= pX) return false;

    int size = clfrCondMaxLength;
    int p1 = irand(size);
    int p2 = irand(size);
    if(p1 > p2){
        std::swap(p1,p2);
    }

    for(int i=p1; i<p2; i++){
        std::swap(cl1.cf_ids[i], cl2.cf_ids[i]);
    }

    return true;
}


/*
 * Sync with original code. Toggle one code fragment
 */
bool mutation(Classifier &clfr, float *state)
{
    bool changed = false;
    for(int i=0; i<clfrCondMaxLength; i++){
        if(drand() < pM){
            changed = true;
            if(clfr.cf_ids[i] != -1){
                clfr.cf_ids[i] = -1; // set as don't care
            }else{
                CodeFragment new_cf;
                clfr.cf_ids[i] = create_new_cf(state);
            }
        }
    }

    return changed;
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
void insertDiscoveredClassifier(Classifier *child, int *parent, ClassifierSet &action_set, int len)
{
    len+=2;
    if(doGASubsumption)
    {
        if(!subsumeClassifier(child[0], population[parent[0]], population[parent[1]], action_set)){
            child[0].id = get_next_cl_gid();
            population[child[0].id] = child[0];
        }
        if(!subsumeClassifier(child[1], population[parent[0]], population[parent[1]], action_set)){
            child[1].id = get_next_cl_gid();
            population[child[1].id] = child[1];
        }
    }
    else
    {
        child[0].id = get_next_cl_gid();
        child[1].id = get_next_cl_gid();
        population[child[0].id] = child[0];
        population[child[1].id] = child[1];
    }

    while(len > maxPopSize)
    {
        len--;
        int cl_id = deleteStochClassifier(population);
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
        if(isSubsumer(population[id]))
        {
            if(subsumer== -1 || isMoreGeneral(population[id], population[subsumer]))
            {
                subsumer = id;
            }
        }
    }

    /* If a subsumer was found, subsume all classifiers that are more specific. */
    if(subsumer!= -1)
    {
        for(auto& id : action_set.ids){
            if(isMoreGeneral(population[subsumer], population[id]))
                population[subsumer].numerosity += population[id].numerosity;
                remove_classifier(action_set, id);
                population[id].id = -1;
        }

    }
}

/**
 * Tries to subsume the parents.
 */
bool subsumeClassifier(Classifier &cl, Classifier &p1, Classifier &p2, ClassifierSet &action_set)
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
    if(subsumeClassifierToSet(cl, action_set)){
        return true;
    }
    return false;
}


/**
 * Try to subsume in the set.
 */
bool subsumeClassifierToSet(Classifier &cl, ClassifierSet &cl_set)
{
    std::list<int> subsumers;

    for(auto & id : cl_set.ids)
    {
        if(subsumes(population[id], cl))
        {
            subsumers.push_back(id);
        }
    }
    /* if there were classifiers found to subsume, then choose randomly one and subsume */
    if(subsumers.size() > 0)
    {
        auto random_it = subsumers.begin();
        random_it = std::next(random_it, irand(subsumers.size()));
        population[*random_it].numerosity++;
        return true;
    }
    return false;
}

int count_classifier_cfs(const Classifier &cl)
{
    int count = 0;
    for(int id : cl.cf_ids){
        if(id != -1) count++;
    }
    return count;
}

void add_classifier_cfs_to_list(Classifier &cl)
{
    for(int cf_id : cl.cf_ids){
        if(cf_id != -1) {
            add_cf_to_list(get_cf(cf_id));
        }
    }
}

void remove_classifier_cfs_from_list(Classifier &cl)
{
    for(int cf_id : cl.cf_ids){
        if(cf_id != -1) {
            remove_cf_from_list(cf_id);
        }
    }
}

/**
 * Try to subsume in the population.
 */
bool subsumeClassifierToPop(Classifier &cl)
{
    std::list<int> subsumers;

    for(auto & item : population)
    {
        if(item.id == -1) continue; // skip empty slots in the array
        if(subsumes(item, cl))
        {
            subsumers.push_back(item.id);
        }
    }
    /* if there were classifiers found to subsume, then choose randomly one and subsume */
    if(subsumers.size() > 0)
    {
        auto random_it = subsumers.begin();
        random_it = std::next(random_it, irand(subsumers.size()));
        population[*random_it].numerosity++;
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


/*
 * This function checks that all the filters of one classifier are present in the second.
 */

// changed it to compare cfs for equality
bool isMoreGeneral(Classifier &clfr_general, Classifier &clfr_specific)
{
    bool more_general = true;
    for(int i=0; i < clfrCondMaxLength; i++){
        if(clfr_general.cf_ids[i] != -1 && !is_cf_covered(get_cf(clfr_general.cf_ids[i]), clfr_specific)){
            more_general = false;
            break;
        }
    }
    return more_general;
}


// ###################### adding classifiers to a set ###################################

// ############################## deletion ############################################

/**
 * Deletes one classifier in the population.
 * The classifier that will be deleted is chosen by roulette wheel selection
 * considering the deletion vote. Returns position of the macro-classifier which got decreased by one micro-classifier.
 **/
int deleteStochClassifier(ClassifierVector &pop)
{
    double vote_sum=0.0, choicep, meanf=0.0;
    int size=0;

    for(auto& item : pop){
        if(item.id == -1) continue; // skip empty slots in the array
        meanf += item.fitness;
        size += item.numerosity;
    }
    meanf/=(double)size;

    /* get the delete proportion, which depends on the average fitness */
    for(auto& item : pop){
        if(item.id == -1) continue; // skip empty slots in the array
        vote_sum += getDelProp(item, meanf);
    }

    /* choose the classifier that will be deleted */
    choicep= drand() * vote_sum;
    /* look for the classifier */
    vote_sum = 0;
    int removed_id = -1;
    for(auto& item : pop){
        if(item.id == -1) continue; // skip empty slots in the array
        vote_sum += getDelProp(item, meanf);
        if(vote_sum > choicep){
            removed_id = item.id;
            if(item.numerosity > 1){
                item.numerosity--;
            }else{
                classifier_gid_stack.push(removed_id);
                pop[item.id].id = -1;
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


void print_population_stats(std::ofstream &output_stats_file)
{
    int size = 0; //get_pop_size(pop, false);

    output_stats_file<<"\n--- Population Stats ---\n";
    output_stats_file << "Global Classifier ID: " << classifier_gid << std::endl;
    int n_total = 0, n_min = INT16_MAX, n_max = -1;
    int cf_count = 0;
    float f_total = 0, f_min = FLT_MAX, f_max = -1;
    std::for_each(population.begin(), population.end(),
                  [&size, &cf_count, &n_total, &n_min, &n_max, &f_total, &f_min, &f_max]
                          (const ClassifierVector::value_type & item)
                  {
                      if(item.id == -1) return; // skip empty slots in the array
                      size++;
                      cf_count += count_classifier_cfs(item);
                      n_total+= item.numerosity;
                      if(n_min > item.numerosity) n_min = item.numerosity;
                      if(n_max < item.numerosity) n_max = item.numerosity;
                      f_total+= item.fitness;
                      if(f_min > item.fitness) f_min = item.fitness;
                      if(f_max < item.fitness) f_max = item.fitness;
                  });

    output_stats_file<< "Population set size: " << size << std::endl;
    output_stats_file << "Population numerosity size: " << n_total << std::endl;
    output_stats_file<< "Avg cf count: " << cf_count/(float)size << std::endl;
    output_stats_file<<"avg numerosity: "<<n_total/(float)size<<" , max numerosity: "<<n_max<<" , min numerosity: "<<n_min<<std::endl;
    output_stats_file<<"avg fitness: "<<f_total/(float)size<<" , max fitness: "<<f_max<<" , min fitness: "<<f_min<<std::endl;
    output_stats_file<<"--- Population Stats ---\n\n";
}

/**
 * This function saves the classifier population and outputs various stats.
 * This function also saves promising code fragments and filters for reuse by the subsequent experiments
 */
void save_experiment_results(std::string path_postfix)
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
    std::ofstream output_parameter_file;
    output_parameter_file.open(output_full_path + output_parameter_file_name);
    if(!output_parameter_file.is_open()){
        std::cout << "Could not open output parameter file";
        exit(1);
    }
    output_parameter_file<<"pM "<<pM<<std::endl;
    print_population_stats(output_stats_file);
    print_code_fragment_stats(output_stats_file);
    print_filter_stats(output_stats_file);
    print_filter_evaluation_stats(output_stats_file);
    write_classifier_header(output_classifier_file);
    for(auto& item : population)
    {
        if(item.id == -1) continue; // skip empty slots in the array
        fprintClassifier(item, output_classifier_file);
    }
//    output_filters(output_filter_file, output_promising_filter_file);
    output_cf_list(output_code_fragment_file, output_promising_code_fragment_file);
    //storeCFs(pop, fpCF);
    output_classifier_file.close();
    output_code_fragment_file.close();
    output_promising_code_fragment_file.close();
    output_filter_file.close();
    output_promising_filter_file.close();
    output_stats_file.close();
    output_parameter_file.close();
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


void write_classifier_header(std::ofstream &output_classifier_file)
{
    output_classifier_file << "id ";
    output_classifier_file << "numerosity ";
    output_classifier_file << "experience ";
    output_classifier_file << "cf_count ";
    output_classifier_file << "fitness ";
    output_classifier_file << "accuracy ";
    output_classifier_file << "prediction ";
    output_classifier_file << "error ";
    output_classifier_file << "action_set_size ";
    output_classifier_file << "time_stamp ";
    output_classifier_file << "action ";
    output_classifier_file << "cfs... ";
    output_classifier_file << std::endl;

}

void fprintClassifier(Classifier &classifier, std::ofstream &output_classifier_file)
{
    output_classifier_file << std::fixed;
    output_classifier_file << classifier.id << " ";
    output_classifier_file << classifier.numerosity << " ";
    output_classifier_file << classifier.experience << " ";
    output_classifier_file << count_classifier_cfs(classifier) << " ";
    output_classifier_file << classifier.fitness << " ";
    output_classifier_file << classifier.accuracy << " ";
    output_classifier_file << classifier.prediction << " ";
    output_classifier_file << classifier.predictionError << " ";
    output_classifier_file << classifier.actionSetSize << " ";
    output_classifier_file << classifier.timeStamp << " ";
    output_classifier_file << classifier.action << " ";

    for(auto & id : classifier.cf_ids)
    {
        output_classifier_file << id << " ";

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
    return cl.predictionError < epsilon_0 && cl.experience > theta_promising;

}


/*
 * This function is used to update stats of all filters in the filter list
 * Which includes nuemrosity and fitness of the filter
 * If classifier is "promising" then increase the fitness of the filter
 * A promising classifier is one whose error < 10 and experience > 10
 */

void manage_filter_and_cf_list() {
    // reset statistics of all filters before updating
    reset_filter_stats();
    reset_cf_stats();

    for(auto& item : population){
        if(item.id == -1) continue; // skip empty slots in the array
        bool promising = is_promising_classifier(item);
        for(int i=0; i < clfrCondMaxLength; i++){
            if(item.cf_ids[i] != -1) {
                CodeFragment & cf = get_cf(item.cf_ids[i]);
                cf.numerosity++;
                if(promising){
                    cf.fitness++;
                }
                for (int j = 0; j < cf.num_filters; j++) {
                    Filter &f = get_filter(cf.filter_ids[j]);
                    f.numerosity++;
                    // if classifier is "promising" then increase the fitness of the filter
                    // a promising classifier is one whose error < 10 and experience > 10
                    if (promising) {
                        f.fitness++;
                    }

                }
            }
        }
    }
    // remove filters with numerosity=0 from the filter list
    std::forward_list<int> removed_filters;
    remove_unused_filters(removed_filters);
    update_evaluation_cache(removed_filters);
    prepare_promising_filter_list();
    remove_unused_cf();
    prepare_promising_cf_list();
}

int get_pop_size(bool numerosity) {
    int pop_numerosity = 0;
    int pop_size = 0;
    std::for_each(population.begin(), population.end(),
                  [&pop_numerosity, &pop_size](ClassifierVector::value_type& item)
                  {
                      if(item.id == -1) return; // skip empty slots in the array
                      pop_numerosity+= item.numerosity;
                      pop_size += 1;
                  });
    if(numerosity) return pop_numerosity;
    else return pop_size;
}

void get_matching_classifiers(float *state, ClassifierSet &match_set, int img_id, bool train) {
    std::for_each(population.begin(), population.end(), [&match_set, &state, img_id, train](ClassifierVector::value_type& item)
    {
        if(item.id == -1) return; // skip empty slots in the array
        if(isConditionMatched(item, state, img_id, train)){
            match_set.ids.push_back(item.id);
        }
    });
}


void load_classifier(std::string classifier_file_name)
{
    int loaded_cl_gid = -1;
    std::string line;
    std::ifstream cl_file(classifier_file_name);
    if (!cl_file.is_open()) {
        std::string error("Error opening input file: ");
        error.append(classifier_file_name).append(", could not load data!");
        throw std::runtime_error(error);
    }

    // skip header line
    getline(cl_file, line);
    while(getline(cl_file, line)) {
        // load classifier
        Classifier cl;
        int num_cf=-1;
        std::stringstream line1(line);
        line1>>cl.id;
        line1>>cl.numerosity;
        line1>>cl.experience;
        line1>>num_cf;
        line1>>cl.fitness;
        line1>>cl.accuracy;
        line1>>cl.prediction;
        line1>>cl.predictionError;
        line1>>cl.actionSetSize;
        line1>>cl.timeStamp;
        line1>>cl.action;

        for(int i=0; i<clfrCondMaxLength; i++){
            line1>>cl.cf_ids[i];
        }
        population.resize(cl.id + 1);
        population[cl.id] = cl;
        if(loaded_cl_gid < cl.id){
            loaded_cl_gid = cl.id;
        }
    }
    classifier_gid = 1 + loaded_cl_gid;
    // populate stack with available slots till classifier_gid
    for(int i=0; i<classifier_gid; i++){
        if(population[i].id == -1) classifier_gid_stack.push(i);
    }
}


// add new classifiers to population
void add_new_classifiers_to_population(float* state, int action, int itTime)
{
    Classifier coverClfr;
    matchingCondAndSpecifiedAct(coverClfr, state, action, 0, itTime);
    if(!subsumeClassifierToPop(coverClfr)) {
        coverClfr.id = get_next_cl_gid();
        population[coverClfr.id] = coverClfr;
    }
    int len = get_pop_size(true);
    while(len > maxPopSize)
    {
        len--;
        int cl_id = deleteStochClassifier(population);
    }
}