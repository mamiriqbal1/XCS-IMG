//
// Created by amir on 5/25/20.
//

#ifndef RCFC_KB_FILTER_LIST_H
#define RCFC_KB_FILTER_LIST_H
#include "configuration.h"
#include <stdexcept>

int add_filter(Filter filter_to_add);
Filter& get_filter(int filter_id);
void reset_filter_stats();
void remove_unused_filters(std::forward_list<int>& removed_filters);
void print_filter_stats(std::ofstream &output_stats_file);
int get_promising_filter_id();
void output_filters(std::ofstream &output_filter_file, std::ofstream &output_promising_filter_file);
void initialize_filter_list(int size);
void output_cf_list(std::ofstream &output_code_fragment_file, std::ofstream &output_promising_code_fragment_file);
void load_filter_for_resume(std::string filter_file_name);
void load_filter_for_kb(std::string filter_file_name);
Filter get_kb_filter(float* state);
int evaluate_filter(const Filter& filter, float *state, int cl_id=-1, int img_id=-1, bool train=true);
int evaluate_filter_actual(const Filter& filter, float *state);
int evaluate_filter_actual_slide(Filter& filter, float *state);
void prepare_promising_filter_list();


#endif //RCFC_KB_FILTER_LIST_H
