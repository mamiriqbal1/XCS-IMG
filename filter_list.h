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
void print_filter_stats();
int get_promising_filter_id();
void output_filter_to_file(std::ofstream& output_filter_file);

#endif //RCFC_KB_FILTER_LIST_H
