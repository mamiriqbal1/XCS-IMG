#ifndef XCS_IMG_CF_H
#define XCS_IMG_CF_H

#include <fstream>
#include "configuration.h"
#include "classifier.h"

void initializeNewCF(int id, CodeFragment &cf);
int evaluateCF(CodeFragment &cf, float *state, int shift_x, int shift_y, int cl_id = -1, int img_id = -1,
               bool train = true);
int evaluate_cf_slide(CodeFragment &cf, float *state, int cl_id = -1, int img_id = -1, bool train = true);
int create_new_cf(float *state);
void output_code_fragment_to_file(CodeFragment &cf, std::ofstream &output_code_fragment_file);
void save_visualization_data(ClassifierSet &action_set, int img_id, std::ofstream &output_visualization_file);
bool is_cf_covered(CodeFragment& cf, Classifier& cl);
void set_cf_bounding_box(CodeFragment &cf, BoundingBox bb);
void initialize_cf_bounding_box(CodeFragment &cf);
int translate(BoundingBox& bb, int x, int y);
bool set_cf_pattern(CodeFragment &cf, float* state);
bool evaluate_cf_bb(CodeFragment &cf, float *state, int shift_x, int shift_y);



#endif //XCS_IMG_CF_H
