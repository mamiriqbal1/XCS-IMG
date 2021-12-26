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
#include "cf_list.h"
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
typedef std::unordered_map<int, ImageEvaluationMap> FilterEvaluationMap;  // <filter_ids, <img_id, bool>>
FilterEvaluationMap evaluation_map;
FilterEvaluationMap evaluation_validation_map;
unsigned long map_hits = 0;

/*
 * define constant patterns for vertical, horizontal and diagonal edges
 * eventually these should go to some configuration with more generic rotated edges
 * VERTICAL, HORIZONTAL, DIAGONAL_45, DIAGONAL_135 = 0, 1, 2, 3
 */
int EDGE_PATTERNS[4][6][6] =
        {
        {
                {0,0,0,1,1,1},
                {0,0,0,1,1,1},
                {0,0,0,1,1,1},
                {0,0,0,1,1,1},
                {0,0,0,1,1,1},
                {0,0,0,1,1,1},
        },
        {
                {0,0,0,0,0,0},
                {0,0,0,0,0,0},
                {0,0,0,0,0,0},
                {1,1,1,1,1,1},
                {1,1,1,1,1,1},
                {1,1,1,1,1,1},
        },
        {
                {0,0,0,0,0,1},
                {0,0,0,0,1,1},
                {0,0,0,1,1,1},
                {0,0,0,1,1,1},
                {0,0,1,1,1,1},
                {0,1,1,1,1,1},
        },
        {
                {0,1,1,1,1,1},
                {0,0,1,1,1,1},
                {0,0,0,1,1,1},
                {0,0,0,1,1,1},
                {0,0,0,0,1,1},
                {0,0,0,0,0,1},
        }
        };


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
    cf.cf_id = id;
    cf.matching_threshold = filter_matching_threshold;
}


/*
 * Set bb as per parameter and initialize matrix accordingly
 */
void set_cf_bounding_box(CodeFragment &cf, BoundingBox bb)
{
    cf.bb = bb;
    FloatMatrix pattern(cf.bb.size_y, FloatVector (cf.bb.size_x, NOT_INITIALIZED));
    IntMatrix mask(cf.bb.size_y, IntVector (cf.bb.size_x, CELL_LIGHT));
}




bool detect_new_feature(CodeFragment &cf, float *state)
{
    int tries = 0;
    bool found = false;
    do{
        float metric[4] = {0};
        float max = 0;
        int max_i = 0;
        initialize_cf_bounding_box(cf);
        for(int i=0; i<4; i++) {
            cf.edge_pattern_id = i;
            if(evaluate_cf_bb(cf, state, &metric[i])){
                found = true;
            }
            if(metric[i] > max) {
                max = metric[i];
                max_i = i;
            }
        }
        if(found){
            cf.edge_pattern_id = max_i;
            return true;
        }
        tries++;
    }while(tries < 1000);
    assert(false); // it must detect a feature
    return false;
}


/*
 * Initialize bb randomly
 */
void initialize_cf_bounding_box(CodeFragment &cf)
{
    BoundingBox bb;
    bb.size_x = cf_min_bounding_box_size + irand(cf_max_bounding_box_size - cf_min_bounding_box_size+1);
    bb.size_y = cf_min_bounding_box_size + irand(cf_max_bounding_box_size - cf_min_bounding_box_size+1);
    bb.x = irand(image_width - bb.size_x+1);
    bb.y = irand(image_height - bb.size_y+1);
    set_cf_bounding_box(cf, bb);
}


// translate bounding box coordinates to state coordinates
inline int translate(BoundingBox& bb, int x, int y)
{
    return (bb.y+y)*image_width + bb.x +x;
}


// new function for setting filter bounds
void create_new_filter_from_input(Filter &filter, float *state, BoundingBox bb, Position &relative_position)
{
    filter.size_x = irand(bb.size_x)+1;
    filter.size_y = irand(bb.size_y)+1;
    filter.values.reserve(filter.size_x*filter.size_y);
    filter.values.assign(filter.size_x*filter.size_y, -1);

    relative_position = generate_relative_position(filter, bb);
    int state_x = bb.x + relative_position.x;
    int state_y = bb.y + relative_position.y;

    int i = 0;
    for(int y=state_y; y<state_y + filter.size_y; y++){
        for(int x=state_x; x<state_x + filter.size_x; x++){
           filter.values[i] = state[y*image_width+x];
           i++;
        }
    }
}

// new function that randomly selects a position on filter and create a matching filter at that position
void addLeafCF(CodeFragment &cf, float *state, BoundingBox bb) {

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
            Position p;
            cf.filter_ids[leafNum] = get_new_filter(state, bb, p);
            cf.filter_positions[leafNum] = p;
            leafNum++;
        }
    }
    assert(cf.num_filters == leafNum);
}


void update_evaluation_cache(std::forward_list<int>& removed_filters){
    std::forward_list<int>::iterator it;
    for(it = removed_filters.begin(); it != removed_filters.end(); it++){
        evaluation_map.erase(*it);
        evaluation_validation_map.erase((*it));
    }
}



/*
 * This function will save data that can be used to visualize classifiers and filters for an image
 * that was predicted correctly or incorrectly
 */
void save_visualization_data(ClassifierSet &action_set, int img_id, std::ofstream &output_visualization_file)
{
    for(auto & id:action_set.ids){
        output_visualization_file << id << " ";
    }
    output_visualization_file<<std::endl;
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
int get_new_filter(float *state, BoundingBox bb, Position &relative_position) {
    int id = -1;
//    if(use_kb && drand() < p_kb){
//        Filter kb_filter = get_kb_filter(state);
//        if(kb_filter.id != -1) {
//            relative_position = generate_relative_position(kb_filter, bb);
//            id = add_filter(kb_filter);
//        }
//    }
//    if(false && id == -1 && drand() < p_promising){
//        id = get_promising_filter_id();
//        if(id != -1) {
//            Filter temp_filter = get_filter(id);
//            relative_position = generate_relative_position(temp_filter, bb);
//            id = add_filter(temp_filter);
//        }
//    }
    if(id == -1) {
        Filter new_filter;
        Position p;
        create_new_filter_from_input(new_filter, state, bb, p);
        relative_position = p;
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
//bool remove_operator(CodeFragment& cf, float* state){
//    int depth = validateDepth(cf.reverse_polish.data());
//    if(depth == 0) return false;
//
//    if(cf.reverse_polish[1] == OPNOT){ // NOT will only be added to a zero depth cf via negation
//        return false;
//    }
//    std::vector<int> temp_reverse_polish;
//    temp_reverse_polish.reserve(cfMaxLength);
//    temp_reverse_polish.assign(cfMaxLength, OPNOP);
//    temp_reverse_polish = cf.reverse_polish;
//    for(int i=0; i<cfMaxLength; i++){
//        if(cf.reverse_polish[i] < 0){ // its an operator
//            // remove leaf from leaves list
//            int leave_index = cf.reverse_polish[i-1];
//            cf.filter_ids[leave_index] = -1;
//            for(int j=leave_index+1; j<cfMaxLeaf; j++){
//                cf.filter_ids[j - 1] = cf.filter_ids[j];
//            }
//            cf.filter_ids[cfMaxLeaf - 1] = -1;
//            // now adjust indexes of all leaves which were greater than leave_index
//            for(int k=0; k<cfMaxLength; k++){
//                if(cf.reverse_polish[k] > leave_index){
//                    cf.reverse_polish[k]--;
//                    temp_reverse_polish[k]--;
//                }
//            }
//
//            // shift left all the contents to remove the operator and one operand
//            std::copy(
//                    std::next(temp_reverse_polish.begin(),i+1),
//                    temp_reverse_polish.end(),
//                    std::next(cf.reverse_polish.begin(),i-1)
//            );
//            // now reset last two  slots
//            auto it = cf.reverse_polish.end();
//            it--;
//            *it = OPNOP;
//            it--;
//            *it = OPNOP;
//            cf.num_filters--;
//            if (evaluateCF(cf, state) != 1){
//                negate_cf(cf);
//                assert(evaluateCF(cf, state) == 1);
//            }
//            return true;
//        }
//    }
//    assert(false); // should not reach here
//}




/*
 * Shrink cf by removing an operator from it
 */
//bool shrink_cf(CodeFragment &cf, float* state){
//    return remove_operator(cf, state);
//}


/*
 * Adds a new operator to the cf from the operator list. The operator list does not include OPNOT operator.
 * OPNOT will only be added as a result of negation of zero depth cf.
 */
//bool add_operator(CodeFragment& cf, float* state){
//    int depth = validateDepth(cf.reverse_polish.data());
//    if(is_full(cf)) return false;
////    if(depth >= cfMaxDepth) return false;
//
//    if(cf.reverse_polish[1] == OPNOT){ // NOT will only be added to a zero depth cf via negation
//        assert(depth == 1 && cf.num_filters == 1);
//        negate_cf(cf); // remove NOT operator
//    }
//    opType selected_operator = functionCodes[irand(totalFunctions)];
//    std::vector<int> temp_reverse_polish;
//    temp_reverse_polish.reserve(cfMaxLength);
//    temp_reverse_polish.assign(cfMaxLength, OPNOP);
//    temp_reverse_polish = cf.reverse_polish;
//    //std::copy(cf.reverse_polish.begin(), cf.reverse_polish.end(), temp_reverse_polish.begin());
//    int new_filter_id = get_new_filter(state, BoundingBox(), <#initializer#>);
//    for(int i=0; i<cfMaxLength; i++){
//        if(cf.reverse_polish[i] >= 0){ // its an operand
//            // shift right all the contents to make room for one operand and one operator
//            std::copy(
//                    std::next(temp_reverse_polish.begin(),i+1),
//                    std::prev(temp_reverse_polish.end(),2),
//                    std::next(cf.reverse_polish.begin(),i+3)
//                    );
//            cf.reverse_polish[i+1] = cf.num_filters;
//            cf.filter_ids[cf.num_filters] = new_filter_id;
//            cf.num_filters++;
//            cf.reverse_polish[i+2] = selected_operator;
//            if(validateDepth(cf.reverse_polish.data()) <= cfMaxDepth) {
//                // now evaluate cf before returning
//                if (evaluateCF(cf, state) != 1){
//                    negate_cf(cf);
//                    assert(evaluateCF(cf, state) == 1);
//                }
//                return true;
//            }else{  // reset and try the next operand
//                cf.reverse_polish = temp_reverse_polish;
//                //std::copy(std::begin(temp_reverse_polish), std::end(temp_reverse_polish), std::begin(cf.reverse_polish));
//                cf.num_filters--;
//                cf.filter_ids[cf.num_filters] = -1;
//            }
//        }
//    }
//    assert(false); // should not reach here
//}


int create_new_cf(float *state) {
    int new_cf_id = -1;
//    initializeNewCF(-1, new_cf);
    if (use_kb && drand() < p_kb) {
        CodeFragment received_cf = get_kb_code_fragment(state);
        if (received_cf.cf_id != -1) {
            if(evaluate_cf_slide(received_cf, state)){  // todo slide now depends on a parameter in the evaluate_cf_slide which is currently zero
                received_cf.cf_id = get_next_cf_gid();
                received_cf.numerosity = 1;
                received_cf.fitness = 0;
                // add the filters from kb to master filter list
                transfer_filters_from_kb_cf(received_cf);
                add_cf_to_list(received_cf);
                new_cf_id = received_cf.cf_id;
            }
        }
    }
    if(new_cf_id == -1 && drand() < p_promising){
        int id = get_promising_cf_id();
        if(id != -1) {
            CodeFragment& cf = get_cf(id);
            if(evaluate_cf_slide(cf, state)) {     // only use promising cf if it evaluates to true
                new_cf_id = cf.cf_id;
            }
//            else{  // copy and negate cf
//                CodeFragment temp = cf;
//                temp.cf_id = get_next_cf_gid();
//                temp.numerosity = 1;
//                temp.fitness = 0;
//                negate_cf(temp);
//                add_cf_to_list(temp);
//            }
        }
    }
    if (new_cf_id == -1) { // if cf not received from kb
        CodeFragment new_cf;
        initializeNewCF(-1, new_cf);
        detect_new_feature(new_cf, state);
//        opType *end = randomProgram(new_cf.reverse_polish.data(), 0, cfMaxDepth, cfMinDepth);
//        addLeafCF(new_cf, state, new_cf.bb);
//        if (evaluateCF(new_cf, state) != 1) {
//            negate_cf(new_cf);
//        }
        new_cf.cf_id = get_next_cf_gid();
        add_cf_to_list(new_cf);
        new_cf_id = new_cf.cf_id;
    }
    return new_cf_id;
}


/*
 * Mutate CF by changing the mask
 */
bool mutate_cf(CodeFragment &cf, float *state) {
    float change = drand() / 100;
    int sign = irand(2) == 0 ? 1 : -1;
    cf.matching_threshold += (change * sign);
    cf.matching_threshold = std::fmax(cf.matching_threshold, 0);
    cf.matching_threshold = std::fmin(cf.matching_threshold, 1);
    return true;
}


bool is_cf_equal(CodeFragment& cf1, CodeFragment& cf2)
{
    if(cf1.cf_id == cf2.cf_id) return true;

    if(cf1.bb.x == cf2.bb.x&&
        cf1.bb.y == cf2.bb.y &&
        cf1.bb.size_x == cf2.bb.size_x &&
        cf1.bb.size_y == cf2.bb.size_y &&
        cf1.edge_pattern_id == cf2.edge_pattern_id){
        return true;
    }else{
        return false;
    }
}


bool is_cf_covered(CodeFragment& cf, Classifier& cl)
{
    if(cf.cf_id == -1) return true;

    bool covered = false;
    for(int i=0; i<clfrCondMaxLength; i++){
        if(cl.cf_ids[i] != -1) {
            if (is_cf_equal(cf, get_cf(cl.cf_ids[i]))) {
                covered = true;
                break;
            }
        }
    }
    return covered;
}


int evaluate_cf_slide(CodeFragment &cf, float *state, int cl_id, int img_id, bool train)
{
    // evaluate cf in the +/- 4 (8x8 region)
    int region_shift = 0; // INT16_MAX; // 0;
    // save cf coordinates
    int cf_x = cf.bb.x;
    int cf_y = cf.bb.y;

    // define region of search
    int y_start = std::max(cf_y - region_shift, 0);
    int x_start = std::max(cf_x - region_shift, 0);
    int y_end = std::min(cf_y + region_shift, image_height - cf.bb.size_y);
    int x_end = std::min(cf_x + region_shift ,image_width - cf.bb.size_x);

    bool matched = false;
    int evaluation = -1;

    for(int y=y_start; y<=y_end && !matched; y++){
        for(int x=x_start; x<=x_end && !matched; x++){
            cf.bb.y = y;
            cf.bb.x = x;
            if(evaluateCF(cf, state, cl_id, img_id, train)){  // keep x,y coordinates for successful evaluation and return true
                matched = true;
            }
        }
    }
    if(matched) evaluation = 1;
    else evaluation = 0;
    // restore coordinates
    cf.bb.y = cf_y;
    cf.bb.x = cf_x;
    return evaluation;
}


bool evaluate_cf_bb(CodeFragment &cf, float *state, float *metric)
{
    float region_sum_1 = 0;
    float region_sum_2 = 0;
    int region_count_1 = 0;
    int region_count_2 = 0;

    for(int y=0; y<cf.bb.size_y; y++){
        for(int x=0; x<cf.bb.size_x; x++){
            int region = EDGE_PATTERNS[cf.edge_pattern_id][y][x];
            float val = state[translate(cf.bb, x, y)];
            if(region == 0){
                region_sum_1 += val;
                region_count_1++;
            }else{
                region_sum_2 += val;
                region_count_2++;
            }
        }
    }
    assert(region_count_1 == region_count_2);
    float region_mean_1 = region_sum_1 / region_count_1;
    float region_mean_2 = region_sum_2 / region_count_2;
    float diff = std::abs(region_mean_1 - region_mean_2);
    if(metric != nullptr){
        *metric = diff;
    }

    if(diff > cf.matching_threshold){
        return true;
    }else{
        return false;
    }
}


int evaluateCF(CodeFragment &cf, float *state, int cl_id, int img_id, bool train) {

    if(evaluate_cf_bb(cf, state, nullptr)) return 1;
    else return 0;


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
            Position p;
            p.x = cf.bb.x + cf.filter_positions[opcode].x;
            p.y = cf.bb.y + cf.filter_positions[opcode].y;
            int result = evaluate_filter(get_filter(cf.filter_ids[opcode]), state, p, cl_id, img_id, train);

            if(result != -1)  // if matched
            {
                stack[SP++] = 1;   //changed
            }
            else // if not matched
            {
                stack[SP++] = 0;   //changed
            }

        }
        else if(opcode == OPNOT)
        {
            const int sp = stack[--SP];
            stack[SP++] = (!sp)?1:0;
            // no change in the filter stack because NOT retains all filters regardless its bool value
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
    output_code_fragment_file << cf.cf_id << " " << cf.numerosity << " " << cf.fitness << " "
    << cf.bb.x << " " << cf.bb.y << " " << cf.bb.size_x << " " << cf.bb.size_y << " "
    << cf.matching_threshold << " " << cf.edge_pattern_id << std::endl;
    return;
}


/*
 * Transfer the filters from kb to running state (i.e. the filter master list)
 */
void transfer_filters_from_kb_cf(CodeFragment & kb_cf)
{
    for(int i=0; i < kb_cf.num_filters; i++){
        kb_cf.filter_ids[i] = transfer_kb_filter(kb_cf.filter_ids[i]);
    }
}