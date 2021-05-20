//
// Created by amir on 5/25/20.
//

#include "filter_list.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <numeric>
#include <vector>
#include <float.h>
#include <stack>
#include <sstream>

//FilterList master_filter_list(maxPopSize*clfrCondMaxLength*cfMaxLeaf, 0); // The main filter list that is maintained
FilterList master_filter_list; // The main filter list that is maintained
FilterList kb_filter_list; // The knowledge base filter list

std::vector<int> filter_gid_vector;
std::stack<int, std::vector<int>> filter_gid_stack(filter_gid_vector);

// lists of promising filter and cf ids updated periodically
std::vector<int> promising_filter_ids;

int get_next_filter_gid()
{
    if(filter_gid_stack.size()>0){
        int val = filter_gid_stack.top();
        filter_gid_stack.pop();
        return val;
    }else{
        // only grow the vector when a new element is needed
        master_filter_list.filters.resize(master_filter_list.gid+1);
        return master_filter_list.gid++;
    }
}

void initialize_filter_list(int size)
{
    master_filter_list.filters.reserve(size);
}



/*
 * Inserts a filter into the list.
 * First it checks if the filter is subsumed by another in the list.
 * Only a filter which is part of at least on "promising" classifier (i.e. with fitness > 0) can subsume the new filter
 * If yes then return the id of more general filter.
 * If no then return the id of new filter.
 */
int add_filter_plain(Filter filter_to_add){
    filter_to_add.id = get_next_filter_gid();
    filter_to_add.numerosity = 1;
    master_filter_list.filters[filter_to_add.id] = filter_to_add;
    return filter_to_add.id;
}


int add_filter(Filter filter_to_add){
    // use find_if from <algorithm> to find a more general filter.
    auto general_filter_iterator = std::find_if(master_filter_list.filters.begin(), master_filter_list.filters.end(),
            [&filter_to_add](const FilterMap::value_type & filter_item) -> bool
            {
                if(filter_item.id == -1) return false; // skip empty slots in the array
                if(filter_item.fitness < 1) return false; // only a promising filter can subsume;
                // only filter of same size and type and position can subsume
                if(filter_item.x != filter_to_add.x || filter_item.y != filter_to_add.y ||
                   filter_item.filter_size != filter_to_add.filter_size ||
                   filter_item.is_dilated != filter_to_add.is_dilated) return false;
                for(int i=0; i<filter_to_add.filter_size*filter_to_add.filter_size; i++){
                    if(filter_item.lower_bounds[i] > filter_to_add.lower_bounds[i]
                       || filter_item.upper_bounds[i] < filter_to_add.upper_bounds[i]){
                        return false;
                    }
                }
                return true;
            });
    if(general_filter_iterator == master_filter_list.filters.end()){  // did not find any general filter
        filter_to_add.id = get_next_filter_gid();
        filter_to_add.numerosity = 1;
        master_filter_list.filters[filter_to_add.id] = filter_to_add;
        return filter_to_add.id;
    }else{
        general_filter_iterator->numerosity++;
        return general_filter_iterator->id;
    }
}



/*
 * Returns reference to an object of Filter against the ID.
 * The copied filter can be used for crossover or mutation
 * Throws std::runtime_error in case id is not found
 */
Filter& get_filter(int filter_id){
    return master_filter_list.filters.at(filter_id);
}


/*
 * This function resets the numerosity and fitness of all filters
 */
void reset_filter_stats(){
    std::for_each(master_filter_list.filters.begin(), master_filter_list.filters.end(),
                            [](FilterMap::value_type & filter_item)
                            {
                                filter_item.numerosity=0;
                                filter_item.fitness=0;
                            });
    promising_filter_ids.clear();
}


void prepare_promising_filter_list()
{
    for(auto it=master_filter_list.filters.begin(); it!=master_filter_list.filters.end();++it){
        if(it->id == -1) continue; // skip empty slots in the array
        if(it->fitness > 0){
            promising_filter_ids.push_back(it->id);
        }
    }
}


/*
 * This function removes any unused filters
 */
void remove_unused_filters(std::forward_list<int>& removed_filters){
    for(auto it=master_filter_list.filters.begin(); it!=master_filter_list.filters.end();++it){
        if(it->id == -1) continue; // skip empty slots in the array
        if(it->numerosity <= 0){
            removed_filters.push_front(it->id);
            filter_gid_stack.push(it->id);
            master_filter_list.filters[it->id].id = -1;  // remove unused filter by making it empty slot
        }
    }
}


/*
 * Print stats about the filter list
 */
void print_filter_stats(std::ofstream &output_stats_file) {
    output_stats_file<<"\n--- Filter Stats ---\n";
    output_stats_file<<"filter_gid: "<<master_filter_list.gid<<std::endl;
    int size = 0;
    int n_total = 0, n_min = INT16_MAX, n_max = -1;
    float f_total = 0, f_min = FLT_MAX, f_max = -1;
    int promising_filters = 0;
    int filter_sizes_count[100] = {0};  // ASSUME MAX SIZE OF FILTER TO BE 100 :)
    int num_dilated = 0;
    std::for_each(master_filter_list.filters.begin(), master_filter_list.filters.end(),
            [&size, &n_total, &n_min, &n_max, &f_total, &f_min, &f_max, &promising_filters, &filter_sizes_count, &num_dilated]
            (const FilterMap::value_type & filter_item)
            {
                if(filter_item.id == -1) return; // skip empty slots in the array
                size++;
                n_total+= filter_item.numerosity;
                if(n_min > filter_item.numerosity) n_min = filter_item.numerosity;
                if(n_max < filter_item.numerosity) n_max = filter_item.numerosity;
                f_total+= filter_item.fitness;
                if(f_min > filter_item.fitness) f_min = filter_item.fitness;
                if(f_max < filter_item.fitness) f_max = filter_item.fitness;
                if(filter_item.fitness>0) promising_filters++;
                filter_sizes_count[filter_item.filter_size]++;
                if(filter_item.is_dilated) num_dilated++;
            });

    output_stats_file<<"filter list size: "<<size<<std::endl;
    output_stats_file<<"avg numerosity: "<<n_total/(float)size<<" , max numerosity: "<<n_max<<" , min numerosity: "<<n_min<<std::endl;
    output_stats_file<<"avg fitness: "<<f_total/(float)size<<" , max fitness: "<<f_max<<" , min fitness: "<<f_min<<std::endl;
    output_stats_file<<"promising filters: "<<promising_filters<<std::endl;
    for(int i=0; i<100; i++){
        if(filter_sizes_count[i] > 0){
            //std::cout<<"filter size "<<i<<" count: "<<filter_sizes_count[i]<<std::endl;
            output_stats_file<<"filter size "<<i<<" count: "<<filter_sizes_count[i]<<std::endl;
        }
    }
    //std::cout<<"dilated filters: "<<num_dilated<<std::endl;
    //std::cout<<"--- Filter Stats ---\n\n";

    output_stats_file<<"dilated filters: "<<num_dilated<<std::endl;
    output_stats_file<<"--- Filter Stats ---\n\n";
}

/*
 * Returns a promising filter using binary tournament selection to be used in mutation or covering
 * Tournament size = 2
 * Returns -1 if no promising filter found
 */
int get_promising_filter_id(){
    // select two filters for tournament
    int size = promising_filter_ids.size();
    if(size == 0) return -1;
    int selected[2];
    selected[0] = promising_filter_ids[irand(size)];
    selected[1] = promising_filter_ids[irand(size)];
    if(master_filter_list.filters[selected[0]].fitness > master_filter_list.filters[selected[1]].fitness){
        return selected[0];
    }else if (master_filter_list.filters[selected[0]].fitness < master_filter_list.filters[selected[1]].fitness){
        return selected[1];
    }else{
        return selected[irand(2)];  // return randomly if both have equal fitness
    }
}

void output_bounds_to_file(Filter& filter, std::ofstream &output_filter_file) {
    output_filter_file << "lower_bound ";
    for(int i=0; i<filter.filter_size*filter.filter_size; i++){
        output_filter_file.width(4);
        output_filter_file << filter.lower_bounds[i] << " ";
    }
    output_filter_file << "\nupper_bound ";
    for(int i=0; i<filter.filter_size*filter.filter_size; i++){
        output_filter_file.width(4);
        output_filter_file << filter.upper_bounds[i] << " ";
    }
}

void output_filter_to_file(Filter& filter, std::ofstream &output_filter_file) {
    output_filter_file << "id ";
    output_filter_file.width(5);
    output_filter_file << filter.id;
    output_filter_file << " x ";
    output_filter_file.width(5);
    output_filter_file << filter.x;
    output_filter_file << " y ";
    output_filter_file.width(5);
    output_filter_file << filter.y;
    output_filter_file << " size ";
    output_filter_file << filter.filter_size;
    output_filter_file << " dilated ";
    output_filter_file << filter.is_dilated;
    output_filter_file << " fitness ";
    output_filter_file.width(11);
    output_filter_file << filter.fitness;
    output_filter_file << " num ";
    output_filter_file << filter.numerosity << std::endl;
    output_bounds_to_file(filter, output_filter_file);
    output_filter_file << std::endl;
}

void output_filters(std::ofstream &output_filter_file, std::ofstream &output_promising_filter_file) {
    for(auto & item : master_filter_list.filters){
        if(item.id == -1) continue; // skip empty slots in the array
        output_filter_to_file(item, output_filter_file);
        if(item.fitness > 0) {
            output_filter_to_file(item, output_promising_filter_file);
        }
    }
}

void extract_filter_attributes(std::string& line, Filter& f) {
    std::string str;
    std::stringstream line1(line);
    line1 >> str;
    line1 >> f.id;
    line1 >> str;
    line1 >> f.x;
    line1 >> str;
    line1 >> f.y;
    line1 >> str;
    line1 >> f.filter_size;
    line1 >> str;
    line1 >> f.is_dilated;
    line1 >> str;
    line1 >> f.fitness;
    line1 >> str;
    line1 >> f.numerosity;
}


void extract_filter_lb(std::string& line, Filter& f) {
    std::string str;
    std::stringstream line2(line);
    line2 >> str;
    f.lower_bounds.reserve(f.filter_size * f.filter_size);
    f.lower_bounds.assign(f.filter_size * f.filter_size, -1);
    f.upper_bounds.reserve(f.filter_size * f.filter_size);
    f.upper_bounds.assign(f.filter_size * f.filter_size, -1);
    for (int i = 0; i < f.filter_size * f.filter_size; i++) {
        line2 >> f.lower_bounds[i];
    }
}


void extract_filter_ub(std::string& line, Filter& f) {
    std::string str;
    std::stringstream line3(line);
    line3 >> str;
    for(int i=0; i<f.filter_size*f.filter_size; i++){
        line3 >> f.upper_bounds[i];
    }
}

void load_filter_for_resume(std::string filter_file_name)
{
    int loaded_gid = -1;
    std::string line;
    std::ifstream filter_file(filter_file_name);
    if (!filter_file.is_open()) {
        std::string error("Error opening input file: ");
        error.append(filter_file_name).append(", could not load data!");
        throw std::runtime_error(error);
    }

    while(getline(filter_file, line)){
        Filter f;
        extract_filter_attributes(line, f);
        getline(filter_file, line);
        extract_filter_lb(line, f);
        getline(filter_file, line);
        extract_filter_ub(line, f);

        master_filter_list.filters.resize(f.id + 1);
        master_filter_list.filters[f.id] = f;
        if(f.id > loaded_gid){
            loaded_gid = f.id;
        }
    }
    master_filter_list.gid = 1 + loaded_gid;
    // populate stack with available slots
    for(int i=0; i<master_filter_list.filters.size(); i++){
        if(master_filter_list.filters[i].id == -1) filter_gid_stack.push(i);
    }
}
/*
 * Only load filter that are sufficiently specific.
 * Avoid don't care and too wide filters
 */
bool should_load_filter(Filter& f)
{
    float thresh_hold = 0.5;
    for(int i=0; i<f.filter_size*f.filter_size; i++){
        if(f.upper_bounds[i] - f.lower_bounds[i] <= thresh_hold){
            return true;
        }
    }
    return false;
}

void load_filter_for_kb(std::string filter_file_name)
{
    int filter_id_seq = 0;
    std::string line;
    std::ifstream filter_file(filter_file_name);
    if (!filter_file.is_open()) {
        std::string error("Error opening input file: ");
        error.append(filter_file_name).append(", could not load data!");
        throw std::runtime_error(error);
    }

    while(getline(filter_file, line)){
        Filter f;
        extract_filter_attributes(line, f);
        getline(filter_file, line);
        extract_filter_lb(line, f);
        getline(filter_file, line);
        extract_filter_ub(line, f);

        // check if the filter is specific enough to be loaded
        if(should_load_filter(f)) {
            // we set the filter ids sequentially
            f.id = filter_id_seq++;
            kb_filter_list.filters.resize(f.id + 1);
            kb_filter_list.filters[f.id] = f;
        }
    }
    kb_filter_list.gid = filter_id_seq;
}


/*
 * Return a matching filter randomly.
 * Returns filter with id -1 if could not match with the image at any location
 * Otherwise return filter with correctly matching x,y coordinates
 */

Filter get_kb_filter(float* state)
{
    Filter f;
    int random_index = irand(kb_filter_list.gid);
    f = kb_filter_list.filters[random_index];
    int result = evaluate_filter_actual_slide(f, state);
    if(result == -1){
        f.id = -1;
    }

    return f;
}


/*
 * Return -1 if not matched otherwise return the location of match
 */
int evaluate_filter_actual(const Filter& filter, float *state)
{
    int step = 1; // this will be used to map normal coordinates to dilated coordinates
    int effective_filter_size = filter.filter_size;
    if(filter.is_dilated){
        step = 2;
        effective_filter_size = filter.filter_size + filter.filter_size -1;
    }
    bool match_failed = false; // flag that controls if the next position to be evaluated when current does not match
    int i = filter.y;
    int j = filter.x;
    int y = 0, x = 0;
    for(int k=0; k<filter.filter_size && !match_failed; k++){  // k is filter y coordinate
        for(int l=0; l<filter.filter_size && !match_failed; l++){  // l is filter x coordinate
            if(state[i*image_width+j + k*step*image_width+l*step] < filter.lower_bounds[k*filter.filter_size+l]
               || state[i*image_width+j + k*step*image_width+l*step] > filter.upper_bounds[k*filter.filter_size+l]){
                match_failed = true;
                y = k;
                x = l;
            }
        }
    }
    if(match_failed){
        // return the actual position of filter that did not match
        return y*filter.filter_size+x;
    }else {  // if matched return -1
        return -1;
    }
}

/*
 * Returns -1 if filter could not be matched to any location
 * Otherwise return position
 * The filter in argument is always modified to the last location checked
 * If the filter matches, the filer has correct x,y coordinates when function ends
 */

int evaluate_filter_actual_slide(Filter& filter, float *state)
{
    int step = 1; // this will be used to map normal coordinates to dilated coordinates
    int effective_filter_size = filter.filter_size;
    if(filter.is_dilated){
        step = 2;
        effective_filter_size = filter.filter_size + filter.filter_size -1;
    }
    std::vector<int> result_list;
    for(int i=0; i<image_height - effective_filter_size; i++){  // i is image y coordinate
        for(int j=0; j<image_width - effective_filter_size; j++){  // j is image x coordinate
            filter.x = i;
            filter.y = j;
            int result = evaluate_filter_actual(filter, state);
            if(result != -1){
                result_list.push_back(result);
            }
        }
    }
    if(result_list.size() > 0){
        int index = irand(result_list.size());
        filter.x = result_list[index] % image_width;
        filter.y = result_list[index] / image_width;
        return result_list[index];
    }else {
        return -1;
    }
}


int evaluate_filter(const Filter& filter, float *state, int cl_id, int img_id, bool train)
{
//    // if cl_id or img_id is -1 then do not check evaluation map otherwise check for prior results
//    if(img_id >=0){
//        // return prior result if it is found
//        if(train){
//            ImageEvaluationMap& inner_map = evaluation_map[filter.id];
//            if(inner_map.count(img_id) > 0){  // the image evaluation exist - unordered map always return  1
//                map_hits++;
//                return inner_map[img_id]>=0;
//            }
//        }else {
//            ImageEvaluationMap &inner_map = evaluation_validation_map[filter.id];
//            if (inner_map.count(img_id) > 0) {  // the image evaluation exist - unordered map always return  1
//                map_hits++;
//                return inner_map[img_id]>=0;
//            }
//        }
//    }

    int evaluation = evaluate_filter_actual(filter, state);

//    // set hasmap entry for re-using evaluation
//    if(img_id >=0) {
//        if(train){
//            ImageEvaluationMap& inner_map = evaluation_map[filter.id];
//            inner_map[img_id] = evaluation;
//        }else{
//            ImageEvaluationMap& inner_map = evaluation_validation_map[filter.id];
//            inner_map[img_id] = evaluation;
//        }
//    }
    return evaluation;
}


