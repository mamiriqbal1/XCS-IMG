#include <iostream>
#include <string>
#include <assert.h>
#include <stdlib.h>
#include "xcsMacros.h"
#include "codeFragment.h"
#include "cf_list.h"
#include <stack>

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
    cf.mask = IntMatrix(cf.bb.size_y, IntVector (cf.bb.size_x, DISABLED));
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


bool set_cf_pattern(CodeFragment &cf, float* state)
{
    float INTERESTING_THRESHOLD = 0.05;
    float EDGE_THRESHOLD = 0.5;
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
    if(max - min <= INTERESTING_THRESHOLD) return false; // result = true;

    for(int y=0; y<cf.bb.size_y - 1; y++){
        for(int x=0; x<cf.bb.size_x - 1; x++){
            if(std::abs(cf.pattern[y][x] - cf.pattern[y][x+1]) > EDGE_THRESHOLD){
                cf.mask[y][x] = ENABLED;
                cf.mask[y][x+1] = ENABLED;
            }
        }
    }

    for(int x=0; x<cf.bb.size_x - 1; x++){
        for(int y=0; y<cf.bb.size_y - 1; y++){
            if(std::abs(cf.pattern[y][x] - cf.pattern[y+1][x]) > EDGE_THRESHOLD){
                cf.mask[y][x] = ENABLED;
                cf.mask[y+1][x] = ENABLED;
            }
        }
    }

    return true;
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
        }
    }
    if (new_cf_id == -1) { // if cf not received from kb
        CodeFragment new_cf;
        initializeNewCF(-1, new_cf);
        bool found_interesting_region = false;
        do {
            initialize_cf_bounding_box(new_cf);
            found_interesting_region = set_cf_pattern(new_cf, state);
        }while(!found_interesting_region);

        new_cf.cf_id = get_next_cf_gid();
        add_cf_to_list(new_cf);
        new_cf_id = new_cf.cf_id;
    }
    return new_cf_id;
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
    bool matched = false;
    // evaluate cf with +/- 1 shift flexibility
    int cf_shift_region = 1;
    for(int cf_shift_y = cf_shift_region*-1; cf_shift_y <= cf_shift_region && !matched ; cf_shift_y++){
        int new_shift_y = std::min(std::max(shift_y + cf_shift_y, 0), image_height - cf.bb.size_y);
        for(int cf_shift_x = cf_shift_region*-1; cf_shift_x <= cf_shift_region && !matched ; cf_shift_x++) {
            int new_shift_x = std::min(std::max(shift_x + cf_shift_x, 0), image_width - cf.bb.size_x);
            matched = evaluate_cf_bb(cf, state, new_shift_x, new_shift_y);
        }
    }
    if(matched) return 1;
    else return 0;
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
}


