#include <iostream>
#include <string.h>
#include <string>
#include <assert.h>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <math.h>
#include "xcsMacros.h"
#include "codeFragment.h"
#include "env.h"
#include "filter_list.h"
#include <unordered_map>
#include <utility>
#include <algorithm>
#include <stack>
//using namespace std;

char globalBuf[1000];
int numPreviousCFs = 0;     //When to set
int startingPreviousCFID = 0;   //what is this value

CodeFragment *previousCFPopulation;

// <image_id, location_where_matched> if location >= 0 then matched otherwise not matched.
typedef std::unordered_map<int, int> ImageEvaluationMap;
typedef std::unordered_map<int, ImageEvaluationMap> FilterEvaluationMap;  // <filter_id, <img_id, bool>>
FilterEvaluationMap evaluation_map;
FilterEvaluationMap evaluation_validation_map;
unsigned long map_hits = 0;


void print_filter_evaluation_stats(std::ofstream &output_stats_file) {
    //std::cout<<"--- Filter Evaluation Stats ---\n";

    output_stats_file<<"--- Filter Evaluation Stats ---\n";
    int size = std::distance(evaluation_map.begin(), evaluation_map.end());
    int total = 0, positive = 0, min = INT16_MAX, max = -1;
    std::for_each(evaluation_map.begin(), evaluation_map.end(),
                  [&total, &positive, &min, &max](const FilterEvaluationMap::value_type & item)
                  {
                      int size = item.second.size();
                      int yes = std::count_if(item.second.begin(), item.second.end(),
                              [](const ImageEvaluationMap::value_type & item2)
                              {
                                return item2.second>=0;
                              });
                      total += size;
                      positive += yes;
                      if(min > yes) min = yes;
                      if(max < yes) max = yes;
                 });
    //std::cout<<"map hits: "<<map_hits<<std::endl;
    //std::cout<<"total evaluations recorded: "<<total<<" , total evaluated filters: "<<size<<std::endl;
    //std::cout<<"avg positive: "<<positive/(float)size<<" , max positive: "<<max<<" , min positive: "<<min<<std::endl;
    //std::cout<<"--- Filter Evaluation Stats ---\n\n";

    output_stats_file<<"map hits: "<<map_hits<<std::endl;
    output_stats_file<<"total evaluations recorded: "<<total<<" , total evaluated filters: "<<size<<std::endl;
    output_stats_file<<"avg positive: "<<positive/(float)size<<" , max positive: "<<max<<" , min positive: "<<min<<std::endl;
    output_stats_file<<"--- Filter Evaluation Stats ---\n\n";
}


opType str_to_opt(std::string str)
{
    if(!str.compare("o")) return OPNOP;
    if(!str.compare("&")) return OPAND;
    if(!str.compare("|")) return OPOR;
    if(!str.compare("d")) return OPNAND;
    if(!str.compare("r")) return OPNOR;
    if(!str.compare("~")) return OPNOT;
    if(!str.compare("^")) return OPXOR;
    if(!str.compare("n")) return OPXNOR;
    std::string error("Error: wrong operator ");
    error.append(str).append(" exiting...");
    throw std::runtime_error(error);
}

int getNumPreviousCFs()
{
    return numPreviousCFs;
}

bool isDontcareCF(CodeFragment &cf)
{
    return (cf.cf_id == -1) ? true : false;
    //return (cf.reverse_polish[0] == OPUNITY) ? true : false;
}

int validateDepth(opType *cf)
{
    //display 'cf' for debugging
    /*
    int i=0;
    char* temp = NULL;
    while(cf[i]!=OPNOP){
    	temp = opchar(cf[i++]);
    	printf("%s ",temp);
    }
    */
    opType* end = nullptr;
    for(int i=0; i<cfMaxLength; i++){
        if(cf[i] == OPNOP){
            end = &cf[i];
            break;
        }
    }
    const opType* start = cf;
    opType* p = end -1;
    int a =1;
    int depth = -1;
    DepthMax(start,&p,a,depth);
    //printf("Depth: %d\n",depth);
    //assert(depth<=cfMaxDepth);
    return depth; // return depth of the cf
}

void DepthMax(const opType* const end,opType** prog, int& argstogo, int& depth)
{
    if(*prog<end || argstogo<=0) return;
    const int a = getNumberOfArguments(*prog[0]);
    argstogo += a-1;
    *prog = (*prog-1);
    depth++;
    const int d0 = depth;
    for(int i=0; i<a; i++)
    {
        int d1 = d0;
        DepthMax(end,prog,argstogo,d1);
        depth = (depth>d1)? depth : d1; //max(depth,d1);
    }
}

void initializeNewCF(int id, CodeFragment &cf)
{
    for(int i=0; i<cfMaxLength; i++)
    {
        cf.reverse_polish[i] = OPNOP;
    }
    for(int i=0; i < cfMaxLeaf; i++){
        cf.filter_id[i] = -1;
    }
    cf.cf_id = id;
}

// create without regard to state
void create_new_filter_from_input_random(Filter& filter, float *state)
{
    // randomly selects a position in the image to create filter bounds
    //int filter_size = (int)sqrt(cfMaxLeaf);  // filter size
    int new_filter_size = filter_sizes[irand(num_filter_sizes)]; // select a filter size randomly
    bool is_dilated = false;
    if(allow_dilated_filters){
        is_dilated = irand(2) != 0;
    }

    filter.filter_size = new_filter_size;
    filter.lower_bounds.reserve(filter.filter_size*filter.filter_size);
    filter.lower_bounds.assign(filter.filter_size*filter.filter_size, -1);
    filter.upper_bounds.reserve(filter.filter_size*filter.filter_size);
    filter.upper_bounds.assign(filter.filter_size*filter.filter_size, -1);
    filter.is_dilated = is_dilated;
    for(int i=0; i<new_filter_size*new_filter_size; i++){
        float lower = drand();
        float upper = drand();
        if(lower > upper) std::swap(lower,upper);

        filter.lower_bounds[i] = roundRealValue(fmax(lower, 0), precisionDigits);
        filter.upper_bounds[i] = roundRealValue(fmin(upper, 1),precisionDigits);
    }
}
// new function for setting filter bounds
void create_new_filter_from_input(Filter& filter, float *state)
{
    // randomly selects a position in the image to create filter bounds
    //int filter_size = (int)sqrt(cfMaxLeaf);  // filter size
    int new_filter_size = filter_sizes[irand(num_filter_sizes)]; // select a filter size randomly
    bool is_dilated = false;
    if(allow_dilated_filters){
        is_dilated = irand(2) != 0;
    }
    int effective_filter_size = new_filter_size;
    if(is_dilated){
        effective_filter_size = new_filter_size + new_filter_size -1;
    }
    int step = is_dilated ? 2 : 1;  // this will be used to map normal coordinates to dilated coordinates
    float pixel_values[new_filter_size*new_filter_size];
    float sum = 0;
    int x_position = -1;
    int y_position = -1;
    do{
        sum = 0;
        int index = 0;
        x_position = irand(image_width-effective_filter_size);
        y_position = irand(image_height-effective_filter_size);
        for(int y=y_position; y<y_position+effective_filter_size; y+=step){
            for(int x=x_position; x<x_position+effective_filter_size; x+=step){
                pixel_values[index] = state[y*image_width+x];
                sum += pixel_values[index];
                index++;
            }
        }
    }while(false && sum <= 0.1); // get to some interesting area in the image. All blanks will be ignored.
    filter.x = x_position;
    filter.y = y_position;

    filter.filter_size = new_filter_size;
    filter.lower_bounds.reserve(filter.filter_size*filter.filter_size);
    filter.lower_bounds.assign(filter.filter_size*filter.filter_size, -1);
    filter.upper_bounds.reserve(filter.filter_size*filter.filter_size);
    filter.upper_bounds.assign(filter.filter_size*filter.filter_size, -1);
    filter.is_dilated = is_dilated;
    float delta = drand();
    for(int i=0; i<new_filter_size*new_filter_size; i++){
        filter.lower_bounds[i] = roundRealValue(fmax(pixel_values[i] - delta, 0), precisionDigits);
        filter.upper_bounds[i] = roundRealValue(fmin(pixel_values[i] + delta, 1),precisionDigits);
    }
}

// new function that randomly selects a position on filter and create a matching filter at that position
void addLeafCF(CodeFragment &cf, float *state){

    // count number of leaves
    int count = 0;
    for(int i=0; i<cfMaxLength; i++){
        const opType opcode = cf.reverse_polish[i];
        if(opcode == OPNOP)
        {
            break;
        }
        if(0<=opcode && opcode<condLength)  //code_fragment bit
        {
            count++;
        }
    }
    cf.num_filters = count;
    //cf.filter = new Filter[count];

    int leafNum = 0;
    for(int i=0; i<cfMaxLength; i++)
    {
        const opType opcode = cf.reverse_polish[i];
        if(opcode == OPNOP)
        {
            break;
        }
        if(0<=opcode && opcode<condLength)  //code_fragment bit
        {
            cf.reverse_polish[i] = leafNum;
            cf.filter_id[leafNum] = get_new_filter(state);
            leafNum++;
        }
    }
    assert(cf.num_filters == leafNum);
}

/*
 * Return -1 if not matched otherwise return the location of match
 */
int evaluate_filter_actual(const Filter& filter, float *state)
{
    int step = 1; // this will be used to map normal coordinates to dilated coordinates
    int effective_filter_size = filter.filter_size;
    if(filter.is_dilated){
        step = 2;
        effective_filter_size = filter.filter_size + filter.filter_size -1;
    }
    bool match_failed = false; // flag that controls if the next position to be evaluated when current does not match
    int i = filter.y;
    int j = filter.x;
    for(int k=0; k<filter.filter_size && !match_failed; k++){  // k is filter y coordinate
        for(int l=0; l<filter.filter_size && !match_failed; l++){  // l is filter x coordinate
            if(state[i*image_width+j + k*step*image_width+l*step] < filter.lower_bounds[k*filter.filter_size+l]
               || state[i*image_width+j + k*step*image_width+l*step] > filter.upper_bounds[k*filter.filter_size+l]){
                match_failed = true;
            }
        }
    }
    if(!match_failed){
        // return the actual position where the filter matched
        return i*image_width+j;
        //return true;
    }
    return -1;
    //return false;
}

int evaluate_filter_actual_old(const Filter& filter, float *state)
{
    int step = 1; // this will be used to map normal coordinates to dilated coordinates
    int effective_filter_size = filter.filter_size;
    if(filter.is_dilated){
        step = 2;
        effective_filter_size = filter.filter_size + filter.filter_size -1;
    }
    bool match_failed = false; // flag that controls if the next position to be evaluated when current does not match
    for(int i=0; i<image_height - effective_filter_size; i++){  // i is image y coordinate
        for(int j=0; j<image_width - effective_filter_size; j++){  // j is image x coordinate
            match_failed = false;
            for(int k=0; k<filter.filter_size && !match_failed; k++){  // k is filter y coordinate
                for(int l=0; l<filter.filter_size && !match_failed; l++){  // l is filter x coordinate
                    if(state[i*image_width+j + k*step*image_width+l*step] < filter.lower_bounds[k*filter.filter_size+l]
                    || state[i*image_width+j + k*step*image_width+l*step] > filter.upper_bounds[k*filter.filter_size+l]){
                        match_failed = true;
                    }
                }
            }
            if(!match_failed){
                // return the actual position where the filter matched
                return i*image_width+j;
                //return true;
            }
        }
    }
    return -1;
    //return false;
}

void update_evaluation_cache(std::forward_list<int>& removed_filters){
    std::forward_list<int>::iterator it;
    for(it = removed_filters.begin(); it != removed_filters.end(); it++){
        evaluation_map.erase(*it);
        evaluation_validation_map.erase((*it));
    }
}

bool evaluate_filter(const Filter& filter, float state[], int cl_id, int img_id, bool train)
{
//    // if cl_id or img_id is -1 then do not check evaluation map otherwise check for prior results
//    if(img_id >=0){
//        // return prior result if it is found
//        if(train){
//            ImageEvaluationMap& inner_map = evaluation_map[filter.id];
//            if(inner_map.count(img_id) > 0){  // the image evaluation exist - unordered map always return  1
//                map_hits++;
//                return inner_map[img_id]>=0;
//            }
//        }else {
//            ImageEvaluationMap &inner_map = evaluation_validation_map[filter.id];
//            if (inner_map.count(img_id) > 0) {  // the image evaluation exist - unordered map always return  1
//                map_hits++;
//                return inner_map[img_id]>=0;
//            }
//        }
//    }

    int evaluation = evaluate_filter_actual(filter, state);

//    // set hasmap entry for re-using evaluation
//    if(img_id >=0) {
//        if(train){
//            ImageEvaluationMap& inner_map = evaluation_map[filter.id];
//            inner_map[img_id] = evaluation;
//        }else{
//            ImageEvaluationMap& inner_map = evaluation_validation_map[filter.id];
//            inner_map[img_id] = evaluation;
//        }
//    }
    return evaluation>=0;
}

/*
 * This function will save data that can be used to visualize classifiers and filters for an image
 * that was predicted correctly or incorrectly
 */
void save_visualization_data(ClassifierSet &match_set, int img_id, std::ofstream &output_visualization_file) {

    for(auto & id:match_set.ids){
        output_visualization_file<<id<<" "<<match_set.pop[id].action<<" ";
    }
    output_visualization_file<<std::endl;
//    for(auto & item:evaluation_validation_map){
//        // save filterid, location, size, isdilated
//        output_visualization_file<<item.first<<" "<<item.second[img_id]<< " "
//        <<get_filter(item.first).filter_size<<" "<<get_filter(item.first).is_dilated<<" ";
//    }
//    output_visualization_file<<std::endl;
}

opType negate_operator(opType op){
    switch(op){
        case OPAND:
            return OPNAND;
        case OPNAND:
            return OPAND;
        case OPOR:
            return OPNOR;
        case OPNOR:
            return OPOR;
        case OPXOR:
            return OPXNOR;
        case OPXNOR:
            return OPXOR;
        default:
            assert(false);
    }// end switch
}


/*
 * Replace the operator with its negative.
 * If the depth is zero then add a NOT operator
 * If the top most operator is OPNOT then depth must be 1 so remove it
 */
bool negate_cf(CodeFragment &cf){
    if(cfMaxDepth == 0) return false;
    bool success = false;
    int depth = validateDepth(cf.reverse_polish.data());
    if(depth == 0){ // add a NOT operator
       cf.reverse_polish[1] = OPNOT;
       cf.reverse_polish[2] = OPNOP;
       return true;
    }
    for(int i=0; i<cfMaxLength; i++){
        if(cf.reverse_polish[i] == OPNOP){ // reached end
            if(cf.reverse_polish[i-1] == OPNOT){ // it must be depth 1 cf
                assert(depth == 1 && cf.num_filters == 1);
                cf.reverse_polish[i-1] = OPNOP;
                success = true;
            }else if(i < cfMaxLength){  // negate the top operator
                cf.reverse_polish[i-1] = negate_operator(cf.reverse_polish[i-1]);
                success = true;
            }
            break;
        }
    }
    return success;
}

/*
 * Get a new filter either from promising filters or create a new from the state
 */
int get_new_filter(float *state) {
    int id = -1;
    if(drand() < p_promising_filter){
        id = get_promising_filter_id();
    }
    if(id == -1) {
        Filter new_filter;
        create_new_filter_from_input(new_filter, state);
        id = add_filter(new_filter);
    }
    return id;
}

bool is_full(CodeFragment& cf){
    int i=0;
    for(i=0; i<cfMaxLength; i++){
        if(cf.reverse_polish[i] == OPNOP){
            break;
        }
    }
    if(i >= cfMaxLength - 2) {
        return true;
    }else {
        return false;
    }
}

/*
 * Removes an operator from the cf.
 * Do nothing in case of depth 0 or depth 1 with NOT operator
 */
bool remove_operator(CodeFragment& cf, float* state){
    int depth = validateDepth(cf.reverse_polish.data());
    if(depth == 0) return false;

    if(cf.reverse_polish[1] == OPNOT){ // NOT will only be added to a zero depth cf via negation
        return false;
    }
    std::vector<int> temp_reverse_polish;
    temp_reverse_polish.reserve(cfMaxLength);
    temp_reverse_polish.assign(cfMaxLength, OPNOP);
    temp_reverse_polish = cf.reverse_polish;
    for(int i=0; i<cfMaxLength; i++){
        if(cf.reverse_polish[i] < 0){ // its an operator
            // remove leaf from leaves list
            int leave_index = cf.reverse_polish[i-1];
            cf.filter_id[leave_index] = -1;
            for(int j=leave_index+1; j<cfMaxLeaf; j++){
                cf.filter_id[j-1] = cf.filter_id[j];
            }
            cf.filter_id[cfMaxLeaf-1] = -1;
            // now adjust indexes of all leaves which were greater than leave_index
            for(int k=0; k<cfMaxLength; k++){
                if(cf.reverse_polish[k] > leave_index){
                    cf.reverse_polish[k]--;
                    temp_reverse_polish[k]--;
                }
            }

            // shift left all the contents to remove the operator and one operand
            std::copy(
                    std::next(temp_reverse_polish.begin(),i+1),
                    temp_reverse_polish.end(),
                    std::next(cf.reverse_polish.begin(),i-1)
            );
            // now reset last two  slots
            auto it = cf.reverse_polish.end();
            it--;
            *it = OPNOP;
            it--;
            *it = OPNOP;
            cf.num_filters--;
            if (evaluateCF(cf, state) != 1){
                negate_cf(cf);
                assert(evaluateCF(cf, state) == 1);
            }
            return true;
        }
    }
    assert(false); // should not reach here
}




/*
 * Shrink cf by removing an operator from it
 */
bool shrink_cf(CodeFragment &cf, float* state){
    return remove_operator(cf, state);
}


/*
 * Adds a new operator to the cf from the operator list. The operator list does not include OPNOT operator.
 * OPNOT will only be added as a result of negation of zero depth cf.
 */
bool add_operator(CodeFragment& cf, float* state){
    int depth = validateDepth(cf.reverse_polish.data());
    if(is_full(cf)) return false;
//    if(depth >= cfMaxDepth) return false;

    if(cf.reverse_polish[1] == OPNOT){ // NOT will only be added to a zero depth cf via negation
        assert(depth == 1 && cf.num_filters == 1);
        negate_cf(cf); // remove NOT operator
    }
    opType selected_operator = functionCodes[irand(totalFunctions)];
    std::vector<int> temp_reverse_polish;
    temp_reverse_polish.reserve(cfMaxLength);
    temp_reverse_polish.assign(cfMaxLength, OPNOP);
    temp_reverse_polish = cf.reverse_polish;
    //std::copy(cf.reverse_polish.begin(), cf.reverse_polish.end(), temp_reverse_polish.begin());
    int new_filter_id = get_new_filter(state);
    for(int i=0; i<cfMaxLength; i++){
        if(cf.reverse_polish[i] >= 0){ // its an operand
            // shift right all the contents to make room for one operand and one operator
            std::copy(
                    std::next(temp_reverse_polish.begin(),i+1),
                    std::prev(temp_reverse_polish.end(),2),
                    std::next(cf.reverse_polish.begin(),i+3)
                    );
            cf.reverse_polish[i+1] = cf.num_filters;
            cf.filter_id[cf.num_filters] = new_filter_id;
            cf.num_filters++;
            cf.reverse_polish[i+2] = selected_operator;
            if(validateDepth(cf.reverse_polish.data()) <= cfMaxDepth) {
                // now evaluate cf before returning
                if (evaluateCF(cf, state) != 1){
                    negate_cf(cf);
                    assert(evaluateCF(cf, state) == 1);
                }
                return true;
            }else{  // reset and try the next operand
                cf.reverse_polish = temp_reverse_polish;
                //std::copy(std::begin(temp_reverse_polish), std::end(temp_reverse_polish), std::begin(cf.reverse_polish));
                cf.num_filters--;
                cf.filter_id[cf.num_filters] = -1;
            }
        }
    }
    assert(false); // should not reach here
}




/*
 * Grow cf by adding a new operator to it respecting the maximum depth
 */
bool grow_cf(CodeFragment &cf, float* state){
    return add_operator(cf, state);
}

bool remove_cf(Classifier &cl, float* state){
   if(cl.cf.size() <= 1) return false;

   int cf_index = irand(cl.cf.size());
   cl.cf.erase(cl.cf.begin()+cf_index);
   return true;
}


bool add_cf(Classifier &cl, float* state){
    int num_cf = cl.cf.size();
    if(num_cf >= clfrCondMaxLength) return false;

    // add new cf to the classifier
    CodeFragment temp;
    initializeNewCF(cf_gid, temp);
    temp.cf_id = -1;
    if (use_kb && drand() < 0.5) {
        CodeFragment received_cf = get_kb_code_fragment(state);;
        if (received_cf.cf_id != -1) {
            received_cf.cf_id = cf_gid;
            temp = received_cf;
            // add the filters from kb to master filter list
            transfer_kb_filter(temp);
        }
    }
    if (temp.cf_id == -1) { // if cf not received from kb
        temp.cf_id = cf_gid;
        // create a cf of depth zero to start with
        opType *end = randomProgram(temp.reverse_polish.data(), 0, cfMaxDepth, 0);
        addLeafCF(temp, state);
    }
    if (evaluateCF(temp, state) != 1) {
        negate_cf(temp);
    }

    cl.cf.push_back(temp);
    cf_gid++;
    return true;
}


/*
 * Mutate CF by randomly changing an operator
 */
bool mutate_cf(CodeFragment &cf, float *state) {
    int functions_index[cfMaxLength];  // collect indices of functions for possible mutation
    int functions_i = 0;
    for(int i=0; /*i<cfMaxLength*/; i++){
        const opType opcode = cf.reverse_polish[i];
        if(opcode == OPNOP){
            break;
        }else if(getNumberOfArguments(opcode) > 1){  // only select binary functions for now
           functions_index[functions_i++] = i;
        }
    }
    if(functions_i == 0) return false;
    // number of functions in the CF are functions_i
    int k = irand(functions_i); // select a function index at random for mutation
    opType new_function=0;
    do{
        new_function = randomFunction(); // there is a possiblity of selecting the same function
        if(getNumberOfArguments(new_function) == getNumberOfArguments(cf.reverse_polish[functions_index[k]])){
            cf.reverse_polish[functions_index[k]] = new_function;
            break;
        }
    }while(true);
    if (evaluateCF(cf, state) != 1){
        negate_cf(cf);
        assert(evaluateCF(cf, state) == 1);
    }
    return true;
}

bool is_cf_equal(CodeFragment& cf1, CodeFragment& cf2)
{
    if(cf1.num_filters != cf2.num_filters ||
       cf1.reverse_polish != cf2.reverse_polish ||
       cf1.filter_id != cf2.filter_id) {
        return false;
    }else {
        return true;
    }
}


bool is_cf_covered(CodeFragment& cf, CodeFragmentVector & cfv)
{
    bool covered = false;
    for(int i=0; i<cfv.size(); i++){
        if(is_cf_equal(cf, cfv[i])){
            covered = true;
            break;
        }
    }
    return covered;
}

bool is_cf_covered(CodeFragment& cf, Classifier& cl)
{
    return is_cf_covered(cf, cl.cf);

}

int evaluateCF(CodeFragment &cf, float *state, int cl_id, int img_id, bool train){
    int stack[cfMaxStack];
    stack[0] = 0;
    int SP = 0;
    int tmptmpfval1 = 0;
    for(int i=0; /*i<cfMaxLength*/; i++)
    {
        const opType opcode = cf.reverse_polish[i];
        if(opcode == OPNOP)
        {
            break;
        }
        if(isPreviousLevelsCode(opcode))  //CF from any previous level
        {

            int valueOfCF = evaluateCF(previousCFPopulation[opcode - condLength], state, train);
            stack[SP++] = valueOfCF;
        }
        else if(0<=opcode && opcode<condLength)  //code_fragment bit
        {

            //if(cf.leaf[opcode].lowerBound<=state[cf.leaf[opcode].featureNumber] && state[cf.leaf[opcode].featureNumber]<=cf.leaf[opcode].upperBound)
            if(evaluate_filter(get_filter(cf.filter_id[opcode]), state, cl_id, img_id, train))
            {
                stack[SP++] = 1;   //changed
            }
            else
            {
                stack[SP++] = 0;   //changed
            }

        }
        else if(opcode == OPNOT)
        {
            const int sp = stack[--SP];
            stack[SP++] = (!sp)?1:0;
        }
        else
        {
            const int sp2 = stack[--SP];
            const int sp1 = stack[--SP];
            switch(opcode)
            {
            case OPAND:
                stack[SP++] = (sp1&&sp2)?1:0;
                break;
            case OPOR:
                stack[SP++] = (sp1||sp2)?1:0;
                break;
            case OPNAND:
                stack[SP++] = (sp1&&sp2)?0:1;
                break;
            case OPNOR:
                stack[SP++] = (sp1||sp2)?0:1;
                break;
            case OPXOR:
                stack[SP++] = (sp1!=sp2)?1:0;
                break;
            case OPXNOR:
                stack[SP++] = (sp1!=sp2)?0:1;
                break;
            }//end switch
        }
    }
    int value = stack[--SP];
    //std::cout<<"SP: "<<SP<<"\n";
    assert(SP==0);

    return value;
}

/*
int evaluateCF(CodeFragment cf, float state[], int cl_id, int img_id){
    // if cl_id or img_id is -1 then do not check evaluation  map otherwise check for prior results
    if(cl_id >=0 and img_id >=0){
        // return prior result if it is found
        if(evaluation_map.find(pair(cl_id, img_id)) != evaluation_map.end()){  // if found
            return evaluation_map[pair(cl_id, img_id)];
        }
    }
    // set featureNumber appropriately and then call evaluateCF_old
    int size = (int)sqrt(cfMaxLeaf);  // filter size
    int index[size*size];
    bool done = false;
    int return_value = 0;

    for(int i=0; i<image_height - size && !done; i++){
        for(int j=0; j<image_width - size && !done; j++){
           for(int k=0; k<size; k++){
               for(int l=0; l<size; l++){
                   index[k*size+l] = i*image_width+j + k*image_width+l;
                   cf.leaf[k*size+l].featureNumber = i*image_width+j + k*image_width+l;
               }
           }
           if(evaluateCF_old(cf, state) == 1){
               done = true;
               return_value = 1;
           }else{
               continue;
           }
        }
    }
    // reset after evaluation
    for(int m=0; m<size*size; m++){
        cf.leaf[m].featureNumber = 0;
    }
    // set hasmap entry for re-using evaluation across epochs
    if(cl_id >=0 and img_id >=0) {
        evaluation_map[pair(cl_id, img_id)] = return_value;
    }
    return return_value;
}
*/

bool isPreviousLevelsCode(const opType code){
    return false;
}

// ######################################## tree operations #################################################

inline int getNumberOfArguments(const opType code){
    switch(code)
    {
    case OPAND:
    case OPOR:
    case OPNAND:
    case OPNOR:
    case OPXOR:
    case OPXNOR:
        return 2;
    case OPNOT:
        return 1;
    default:
        return 0;
    }//end switch code
}

inline opType randomLeaf(){
    opType leaf = OPNOP;
    leaf = 0; // feature number is no more used irand(condLength);
    return leaf;
}

inline opType randomFunction(){
    return functionCodes[irand(totalFunctions)];
}
//generate reverse polish
opType* randomProgram(opType prog[], const int isfull, const int maxDepth, const int minDepth){
    //if reached max depth or probabilistically reached mindepth
    if( maxDepth<=0 || ( (!isfull) && minDepth<=0 && irand(2) ) )
    {
        *prog = randomLeaf();
        return prog+1;
    }
    else
    {
        opType* pc = prog;
        opType newFunction;
        newFunction = randomFunction();
        const int numArgs = getNumberOfArguments(newFunction);
        for(int i=0; i<numArgs; i++)
        {
            pc = randomProgram(pc,isfull,maxDepth-1,minDepth-1);
        }
        *pc = newFunction;
        return pc+1;
    }
}//end randomProgram
// --------------------------------- Start Display Functions -----------------------------

//end leafname

//end leafname


//end opchar


inline std::string op_to_str(opType code)
{
    switch(code)
    {
        case OPAND:
            return std::string("&");
        case OPOR:
            return std::string("|");
        case OPNAND:
            return std::string("d");
        case OPNOR:
            return std::string("r");
        case OPNOT:
            return std::string("~");
        case OPNOP:
            return std::string("o");
        case OPXOR:
            return std::string("^");
        case OPXNOR:
            return std::string("n");
        default:
            return std::string("[" + std::to_string(code) + "!!]");
    }//end switch code

}

void output_code_fragment_to_file(CodeFragment &cf, std::ofstream &output_code_fragment_file)
{
    output_code_fragment_file << cf.cf_id << " ";
    std::string str;
    opType code = 0;
    for(int i=0; i<cfMaxLength; i++){
        code = cf.reverse_polish[i];
        if(code == OPNOP){
            break;
        }else if(0<=code && code < cfMaxLeaf){  // if code is zero then it is a filter_id index
            output_code_fragment_file << "D" << cf.filter_id[code] << " "; // print filter id
        }else if (code>=condLength && isPreviousLevelsCode(code)){
            // output previous when implemented
        }else{
            // output code str
            output_code_fragment_file << op_to_str(code) << " ";
        }
    }
    output_code_fragment_file << std::endl;
}


/*
 * Return a matching filter randomly.
 * Returns filder with id -1 if not found in multiple tries
 */

Filter get_kb_filter(float* state)
{
    Filter f;
    auto random_it = kb_filter.begin();
    bool matched = false;
    int tries=0;
    do{
        tries++;
        random_it = std::next(kb_filter.begin(), irand(kb_filter.size()));
        // ignore the state for now.
        //matched = evaluate_filter(random_it->second, state);
        matched = true;
    }while(!matched && tries < 100);
    if(matched){
        f = random_it->second;
    }else{
        f.id = -1;
    }
    return f;
}


CodeFragment get_kb_code_fragment(float* state)
{
    CodeFragment cf;
    auto random_it = kb_cf.begin();
    bool matched = false;
    int tries = 0;
    do{
        tries++;
        random_it = std::next(kb_cf.begin(), irand(kb_cf.size()));
        // ignore the state for now
        // transfer all filters to the list before evaluation
        //matched = evaluateCF(random_it->second, state);
        matched = true;
    }while(!matched && tries < 100);
    if(matched){
        cf = random_it->second;
    }else{
        cf.cf_id = -1;
    }
    return cf;
}

