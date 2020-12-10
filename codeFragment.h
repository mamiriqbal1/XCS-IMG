#include <fstream>
#include "configuration.h"
#include "classifier.h"

void addLeafCF(CodeFragment &cf, float *state);

int getNumPreviousCFs();
bool isDontcareCF(CodeFragment &cf);

int validateDepth(opType *cf);
void initializeNewCF(int id, CodeFragment &cf);

int evaluateCF(CodeFragment &cf, float *state, int cl_id=-1, int img_id=-1, bool train=true);
bool isPreviousLevelsCode(const opType code);

int getNumberOfArguments(const opType opcode);

opType randomLeaf();

opType randomFunction();
opType* randomProgram(opType prog[], const int isfull, const int maxDepth, const int minDepth);

void DepthMax(const opType* const end,opType** prog, int& argstogo, int& depth);


bool evaluate_filter(const Filter& filter, float state[], int cl_id=-1, int img_id=-1, bool train=true);
void update_evaluation_cache(std::forward_list<int>& removed_filters);
void print_filter_evaluation_stats(std::ofstream &output_stats_file);
bool mutate_cf(CodeFragment &cf, float *state);

bool shrink_cf(CodeFragment &cf, float* state);
bool negate_cf(CodeFragment &cf);
int get_new_filter(float *state);
bool create_new_cf(CodeFragment &cf, float* state);

void output_code_fragment_to_file(CodeFragment &cf, std::ofstream &output_code_fragment_file);
Filter get_kb_filter(float* state);
opType str_to_opt(std::string str);
CodeFragment get_kb_code_fragment(float* state);
void save_visualization_data(ClassifierSet &action_set, int img_id, std::ofstream &output_visualization_file);
bool is_cf_covered(CodeFragment& cf, Classifier& cl);
