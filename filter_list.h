//
// Created by amir on 5/25/20.
//

#ifndef RCFC_KB_FILTER_LIST_H
#define RCFC_KB_FILTER_LIST_H
#include "configuration.h"
#include <stdexcept>

extern FilterList master_filter_list; // The main filter list that is maintained

int add_filter(Filter filter_to_add);
Filter& get_filter(int filter_id);
void reset_filter_stats();
void remove_unused_filters(std::forward_list<int>& removed_filters);
void print_filter_stats(std::ofstream &output_stats_file);
int get_promising_filter_id();
void output_filters(std::ofstream &output_filter_file, std::ofstream &output_promising_filter_file);
void initialize_filter_list(int size);
#endif //RCFC_KB_FILTER_LIST_H
