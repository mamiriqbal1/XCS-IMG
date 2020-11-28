//
// Created by amir on 11/26/20.
//
#include <stack>
#include <assert.h>
#include "cf_list.h"
#include "configuration.h"
#include <iostream>
#include <fstream>


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
        // only grow the vector when a new element is needed
        main_cf_list.resize(cf_gid + 1);
        return cf_gid++;
    }
}


void initialize_cf_list(int size)
{
    main_cf_list.reserve(size);
//    CodeFragment temp;
//    main_cf_list.assign(size, temp);
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
            cf_gid_stack.push(id);
        }
    }else{
        assert(false);
    }
}


/*
 * Print stats about the code fragment list
 */
void print_code_fragment_stats(std::ofstream &output_stats_file) {
    output_stats_file<<"\n--- Code Fragment Stats ---\n";
    output_stats_file<<"cf_gid: "<<cf_gid<<std::endl;
    int n_total = 0, n_min = INT16_MAX, n_max = -1;
    int f_total = 0, f_min = INT16_MAX, f_max = -1;
    int size = 0;
    std::for_each(main_cf_list.begin(), main_cf_list.end(),
                  [&n_total, &n_min, &n_max, &f_total, &f_min, &f_max, &size]
                          (const CodeFragmentVector::value_type & item)
                  {
                      if(item.cf_id == -1) return; // skip empty slots in the array
                      size++;
                      n_total+= item.numerosity;
                      if(n_min > item.numerosity) n_min = item.numerosity;
                      if(n_max < item.numerosity) n_max = item.numerosity;
                      f_total+= item.num_filters;
                      if(f_min > item.num_filters) f_min = item.num_filters;
                      if(f_max < item.num_filters) f_max = item.num_filters;
                  });
    output_stats_file<<"cf list size: "<<size<<std::endl;
    output_stats_file<<"avg numerosity: "<<n_total/(float)size<<" , max numerosity: "<<n_max<<" , min numerosity: "<<n_min<<std::endl;
    output_stats_file<<"avg num filters: "<<f_total/(float)size<<" , max num filters: "<<f_max<<" , min num filters: "<<f_min<<std::endl;
    output_stats_file<<"--- Code Fragment Stats ---\n\n";
}
