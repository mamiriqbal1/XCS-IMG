//
// Created by amir on 11/26/20.
//

#ifndef RCFC_KB_CF_LIST_H
#define RCFC_KB_CF_LIST_H

#include "configuration.h"

int get_next_cf_gid();
void initialize_cf_list(int size);
CodeFragment& get_cf(int id);
void add_cf_to_list(CodeFragment & cf);
void remove_cf_from_list(int id);
void print_code_fragment_stats(std::ofstream &output_stats_file);
void load_cf_for_resume(std::string cf_file_name);
void update_cf_list_parameters(ClassifierVector pop);
void reset_cf_stats();
void prepare_promising_cf_list();
int get_promising_cf_id();
void remove_unused_cf();
void load_cf_for_kb(std::string cf_file_name);
CodeFragment get_kb_code_fragment(float* state);
int count_total_cf();
#endif //RCFC_KB_CF_LIST_H