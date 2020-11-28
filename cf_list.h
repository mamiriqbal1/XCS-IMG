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

#endif //RCFC_KB_CF_LIST_H
