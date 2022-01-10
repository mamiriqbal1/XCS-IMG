#ifndef XCS_IMG_CF_H
#define XCS_IMG_CF_H

#include <fstream>
#include "configuration.h"
#include "classifier.h"

void addLeafCF(CodeFragment &cf, float *state, BoundingBox bb);

int getNumPreviousCFs();
bool isDontcareCF(CodeFragment &cf);

int validateDepth(opType *cf);
void initializeNewCF(int id, CodeFragment &cf);

int evaluateCF(CodeFragment &cf, float *state, int shift_x, int shift_y, int cl_id = -1, int img_id = -1,
               bool train = true);
int evaluate_cf_slide(CodeFragment &cf, float *state, int cl_id = -1, int img_id = -1, bool train = true);
bool isPreviousLevelsCode(const opType code);

int getNumberOfArguments(const opType opcode);

opType randomLeaf();

opType randomFunction();
opType* randomProgram(opType prog[], const int isfull, const int maxDepth, const int minDepth);

void DepthMax(const opType* const end,opType** prog, int& argstogo, int& depth);


bool mutate_cf(CodeFragment &cf, float *state);

bool shrink_cf(CodeFragment &cf, float* state);
bool negate_cf(CodeFragment &cf);
int get_new_filter(float *state, BoundingBox bb, Position &relative_position);
int create_new_cf(float *state);

void output_code_fragment_to_file(CodeFragment &cf, std::ofstream &output_code_fragment_file);
opType str_to_opt(std::string str);
void save_visualization_data(ClassifierSet &action_set, int img_id, std::ofstream &output_visualization_file);
bool is_cf_covered(CodeFragment& cf, Classifier& cl);
void transfer_filters_from_kb_cf(CodeFragment & kb_cf);
void set_cf_bounding_box(CodeFragment &cf, BoundingBox bb);
void initialize_cf_bounding_box(CodeFragment &cf);
int translate(BoundingBox& bb, int x, int y);
bool set_cf_pattern_and_mask(CodeFragment &cf, float* state);
bool evaluate_cf_bb(CodeFragment &cf, float *state, int shift_x, int shift_y);



#endif //XCS_IMG_CF_H
