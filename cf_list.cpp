//
// Created by amir on 11/26/20.
//
#include <stack>
#include <assert.h>
#include "cf_list.h"
#include "configuration.h"


int cf_gid = 0;
CodeFragmentVector main_cf_list;

std::vector<int> cf_gid_vector;
std::stack<int, std::vector<int>> cf_gid_stack(cf_gid_vector);

int get_next_cf_gid()
{
    if(cf_gid_stack.size()>0){
        int val = cf_gid_stack.top();
        cf_gid_stack.pop();
        return val;
    }else{
        return cf_gid++;
    }
}


void initialize_cf_list(int size)
{
    main_cf_list.reserve(size);
    CodeFragment temp;
    main_cf_list.assign(size, temp);
}

CodeFragment& get_cf(int id)
{
    return main_cf_list[id];
}


void add_cf_to_list(CodeFragment & cf)
{
    assert(cf.cf_id != -1);
    if(main_cf_list[cf.cf_id].cf_id != -1){
        main_cf_list[cf.cf_id].numerosity++;
    }else{
        cf.numerosity = 1;
        main_cf_list[cf.cf_id] = cf;
    }
}

void remove_cf_from_list(int id)
{
    if(main_cf_list[id].cf_id != -1) {
        main_cf_list[id].numerosity--;
        if (main_cf_list[id].numerosity == 0) {
            main_cf_list[id].cf_id = -1;
        }
    }else{
        assert(false);
    }
}