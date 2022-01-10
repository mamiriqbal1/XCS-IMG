#include <iostream>
#include <string>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include "xcsMacros.h"
#include "codeFragment.h"
#include "filter_list.h"
#include "cf_list.h"
#include <unordered_map>
#include <stack>
//using namespace std;

char globalBuf[1000];
int numPreviousCFs = 0;     //When to set
int startingPreviousCFID = 0;   //what is this value

CodeFragment *previousCFPopulation;

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
    cf.pattern = FloatMatrix(cf.bb.size_y, FloatVector (cf.bb.size_x, NOT_INITIALIZED));
    cf.mask = IntMatrix(cf.bb.size_y, IntVector (cf.bb.size_x, ENABLED));
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


bool set_cf_pattern_and_mask(CodeFragment &cf, float* state)
{
    bool result = false;
    float max = 0;
    float min = 1;
    // initialize the pattern and mask
    // pattern gets all values from the state
    for(int y=0; y<cf.bb.size_y; y++){
        for(int x=0; x<cf.bb.size_x; x++){
            cf.pattern[y][x] = state[translate(cf.bb, x, y)];
            if(cf.pattern[y][x] > max) max = cf.pattern[y][x];
            if(cf.pattern[y][x] < min) min = cf.pattern[y][x];
        }
    }
    // see if region is not all homogeneous
    if(max - min > 0.1) result = true;

    // mask enables 25% to 75% pixels in the pattern
    double enabled_pixels_percentage = (drand() / 2) + 0.5;
    for(int y=0; y<cf.bb.size_y; y++){
        for(int x=0; x<cf.bb.size_x; x++){
            if(true || drand() < enabled_pixels_percentage){
                cf.mask[y][x] = ENABLED;
            }else{
                cf.mask[y][x] = DISABLED;
            }
        }
    }

    return result;
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
        bool found_interesting_region = false;
        do {
            initialize_cf_bounding_box(new_cf);
            found_interesting_region = set_cf_pattern_and_mask(new_cf, state);
        }while(!found_interesting_region);
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
        cf1.pattern == cf2.pattern &&
        cf1.mask == cf2.mask){
        return true;
    }else{
        return false;
    }

//    if(cf1.num_filters != cf2.num_filters ||
//       cf1.reverse_polish != cf2.reverse_polish ||
//       cf1.bb.x != cf2.bb.x ||
//       cf1.bb.y != cf2.bb.y ||
//       cf1.bb.size != cf2.bb.size ||
////       cf1.filter_positions != cf2.filter_positions ||
//       cf1.filter_ids != cf2.filter_ids) {
//        return false;
//    }else {
//        for(int i=0; i<cf1.num_filters; i++){
//            if(cf1.filter_positions[i].x != cf2.filter_positions[i].x ||
//                cf1.filter_positions[i].y != cf2.filter_positions[i].y){
//                return false;
//            }
//        }
//        return true;
//    }
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
            if(evaluateCF(cf, state, 0, 0, cl_id, img_id, train)){  // keep x,y coordinates for successful evaluation and return true
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


bool evaluate_cf_bb(CodeFragment &cf, float *state, int shift_x, int shift_y)
{
    double distance = 0;
    int mask_size = 0;

    for(int y=0; y<cf.bb.size_y; y++){
        for(int x=0; x<cf.bb.size_x; x++){
            if(cf.mask[y][x] == ENABLED) {
                mask_size += 1;
                BoundingBox bb = cf.bb;
                bb.y = bb.y + shift_y;
                bb.x = bb.x + shift_x;
                distance += std::abs(cf.pattern[y][x] - state[translate(bb, x, y)]);
            }
        }
    }

    if(distance/(mask_size) < cf.matching_threshold){  // average distance per pixel
        return true;
    }else{
        return false;
    }
}


int evaluateCF(CodeFragment &cf, float *state, int shift_x, int shift_y, int cl_id, int img_id, bool train) {
    if(evaluate_cf_bb(cf, state, shift_x, shift_y)) return 1;
    else return 0;
}


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
    << cf.bb.x << " " << cf.bb.y << " " << cf.bb.size_x << " " << cf.bb.size_y << " " << cf.matching_threshold << " ";
    for(int y=0; y<cf.bb.size_y; y++){
        for(int x=0; x<cf.bb.size_x; x++){
            output_code_fragment_file << cf.pattern[y][x] <<  " ";
        }
    }
    output_code_fragment_file << std::endl;
    return;

    std::string str;
    opType code = 0;
    for(int i=0; i<cfMaxLength; i++){
        code = cf.reverse_polish[i];
        if(code == OPNOP){
            break;
        }else if(0<=code && code < cfMaxLeaf){  // if code is zero then it is a filter_ids index
            output_code_fragment_file << "D" << cf.filter_ids[code] << " " << cf.filter_positions[code].x << " " << cf.filter_positions[code].y << " "; // print filter id
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
 * Transfer the filters from kb to running state (i.e. the filter master list)
 */
void transfer_filters_from_kb_cf(CodeFragment & kb_cf)
{
    for(int i=0; i < kb_cf.num_filters; i++){
        kb_cf.filter_ids[i] = transfer_kb_filter(kb_cf.filter_ids[i]);
    }
}