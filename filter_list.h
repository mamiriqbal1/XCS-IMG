//
// Created by amir on 5/25/20.
//

#ifndef RCFC_KB_FILTER_LIST_H
#define RCFC_KB_FILTER_LIST_H
#include "configuration.h"
#include <stdexcept>


int add_filter(Filter filter_to_add);
const Filter& get_filter(int filter_id);
int push_front(Filter filter_to_push);
void pop_front();

#endif //RCFC_KB_FILTER_LIST_H
