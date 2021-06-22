//
// Created by amir on 11/26/20.
//
#include <stack>
#include <assert.h>
#include "cf_list.h"
#include "configuration.h"
#include "codeFragment.h"
#include <iostream>
#include <fstream>
#include <sstream>


int cf_gid = 0;
CodeFragmentVector main_cf_list;

std::vector<int> cf_gid_vector;
std::stack<int, std::vector<int>> cf_gid_stack(cf_gid_vector);

// lists of promising filter and cf ids updated periodically
std::vector<int> promising_cf_ids;


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
        if (main_cf_list[id].numerosity <= 0) {
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


void output_cf_list(std::ofstream &output_code_fragment_file, std::ofstream &output_promising_code_fragment_file) {
    for(CodeFragment & item : main_cf_list){
        if(item.cf_id == -1) continue; // skip empty slots in the array
        output_code_fragment_to_file(item, output_code_fragment_file);
        //output promising code fragments separately
//        if (is_promising_classifier(classifier)) {
//            output_code_fragment_to_file(id, output_promising_code_fragment_file);
//        }

    }
}

void load_code_fragment(std::string cf_file_name)
{
    int loaded_cf_gid = -1;
    std::string line;
    std::ifstream cf_file(cf_file_name);
    if (!cf_file.is_open()) {
        std::string error("Error opening input file: ");
        error.append(cf_file_name).append(", could not load data!");
        throw std::runtime_error(error);
    }

    while(getline(cf_file, line)) {
        // load cf
        CodeFragment cf;
        int id=0;
        std::stringstream line1(line);
        line1>>id;
        initializeNewCF(id, cf);
        int val=0;
        line1>>val;
        cf.bb.x = val;
        line1>>val;
        cf.bb.y = val;
        line1>>val;
        cf.bb.size = val;

        int index = 0, leaf_index = 0;
        while(!line1.eof()){
            std::string token;
            line1>>token;
            // last token is "" that needs to be handled
            if(token.empty()) break;
            if(token.substr(0,1) == "D"){ // this is filter id
                int filter_id = std::stoi(token.substr(1));
                cf.filter_ids[leaf_index] = filter_id;
                int val = 0;
                line1 >> val;
                cf.filter_positions[leaf_index].x = val;
                line1 >> val;
                cf.filter_positions[leaf_index].y = val;
                cf.filter_ids[leaf_index] = filter_id;
                cf.reverse_polish[index] = leaf_index;
                leaf_index++;
            }else{ // this is operator
                cf.reverse_polish[index] = str_to_opt(token);
            }
            index++;
        }
        cf.reverse_polish[index] = OPNOP; // terminate the reverse polish
        cf.num_filters = leaf_index;
        main_cf_list.resize(cf.cf_id+1);
        main_cf_list[cf.cf_id] = cf;
        if(cf.cf_id > loaded_cf_gid){
            loaded_cf_gid = cf.cf_id;
        }
    }
    cf_gid = 1 + loaded_cf_gid;
    // populate stack with available slots
    for(int i=0; i<main_cf_list.size(); i++){
        if(main_cf_list[i].cf_id == -1) cf_gid_stack.push(i);
    }
}


void update_cf_list_parameters(ClassifierVector pop)
{
    for(Classifier& cl: pop){
        if(cl.id == -1) continue; // skip empty slot
        for(int cf_id : cl.cf_ids){
            if(cf_id == -1) continue; // skip empty slot
            main_cf_list[cf_id].numerosity++;
        }
    }
}

void reset_cf_stats()
{
    for(CodeFragment& cf: main_cf_list){
        if(cf.cf_id == -1) continue; // skip empty slot
        cf.numerosity = 0;
        cf.fitness = 0;
    }
    promising_cf_ids.clear();
}


void remove_unused_cf()
{
    for(CodeFragment& cf: main_cf_list){
        if(cf.cf_id == -1) continue; // skip empty slot
        if(cf.numerosity <= 0){
            remove_cf_from_list(cf.cf_id);
        }
    }
}

void prepare_promising_cf_list()
{
    for(CodeFragment& cf: main_cf_list){
        if(cf.cf_id == -1) continue; // skip empty slot
        if(cf.fitness > 0) promising_cf_ids.push_back(cf.cf_id);
    }

}

int get_promising_cf_id()
{
    int size = promising_cf_ids.size();
    if(size == 0) return -1;
    int index = irand(promising_cf_ids.size());
    // select two filters for tournament
    int selected[2];
    selected[0] = promising_cf_ids[irand(size)];
    selected[1] = promising_cf_ids[irand(size)];
    if(main_cf_list[selected[0]].fitness > main_cf_list[selected[1]].fitness){
        return selected[0];
    }else if (main_cf_list[selected[0]].fitness < main_cf_list[selected[1]].fitness){
        return selected[1];
    }else{
        return selected[irand(2)];  // return randomly if both have equal fitness
    }
}
