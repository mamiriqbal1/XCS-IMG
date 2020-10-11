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

FilterList master_filter_list; // The main filter list that is maintained


/*
 * Inserts a filter into the list.
 * First it checks if the filter is subsumed by another in the list.
 * Only a filter which is part of at least on "promising" classifier (i.e. with fitness > 0) can subsume the new filter
 * If yes then return the id of more general filter.
 * If no then return the id of new filter.
 */
int add_filter(Filter filter_to_add){
    // use find_if from <algorithm> to find a more general filter.
    auto general_filter_iterator = std::find_if(master_filter_list.filters.begin(), master_filter_list.filters.end(),
            [&filter_to_add](const FilterMap::value_type & filter_item) -> bool
            {
                if(filter_item.second.fitness < 1) return false; // only a promising filter can subsume;
                // only filter of same size and type and position can subsume
                if(filter_item.second.filter_size != filter_to_add.filter_size ||
                        filter_item.second.is_dilated != filter_to_add.is_dilated ||
                        filter_item.second.x != filter_to_add.x ||
                        filter_item.second.y != filter_to_add.y) return false;
                for(int i=0; i<filter_to_add.filter_size*filter_to_add.filter_size; i++){
                    if(filter_item.second.lower_bounds[i] > filter_to_add.lower_bounds[i]
                       || filter_item.second.upper_bounds[i] < filter_to_add.upper_bounds[i]){
                        return false;
                    }
                }
                return true;
            });
    if(general_filter_iterator == master_filter_list.filters.end()){  // did not find any general filter
        filter_to_add.id = master_filter_list.gid++;
        filter_to_add.numerosity = 1;
        master_filter_list.filters[filter_to_add.id] = filter_to_add;
        return filter_to_add.id;
    }else{
        general_filter_iterator->second.numerosity++;
        return general_filter_iterator->second.id;
    }
}



/*
 * Returns reference to an object of Filter against the ID.
 * The copied filter can be used for crossover or mutation
 * Throws std::runtime_error in case id is not found
 */
Filter& get_filter(int filter_id){
    if(master_filter_list.filters.count(filter_id)){
        return master_filter_list.filters.at(filter_id);
    }else{
        std::cout<<std::endl<<filter_id<<std::endl;
        //print_filter_stats(<#initializer#>);
        // print stack trace by causing a segmentation fault that will be handled by our handler
        int *foo = (int *) -1; // make a bad pointer
        printf("%d\n", *foo);       // causes segfault
    }
}


/*
 * This function resets the numerosity and fitness of all filters
 */
void reset_filter_stats(){
    std::for_each(master_filter_list.filters.begin(), master_filter_list.filters.end(),
                            [](FilterMap::value_type & filter_item)
                            {
                                filter_item.second.numerosity=0;
                                filter_item.second.fitness=0;
                            });
}

/*
 * This function removes any unused filters
 */
void remove_unused_filters(std::forward_list<int>& removed_filters){
    for(auto it=master_filter_list.filters.begin(); it!=master_filter_list.filters.end();it++){
        if(it->second.numerosity == 0){
            removed_filters.push_front(it->second.id);
        }
    }
    for(auto it=removed_filters.begin(); it!=removed_filters.end();it++){
        master_filter_list.filters.erase(*it);
    }
    //print_filter_stats();
}


/*
 * Print stats about the filter list
 */
void print_filter_stats(std::ofstream &output_stats_file) {
    //std::cout<<"\n--- Filter Stats ---\n";
    //std::cout<<"classifier_gid: "<<master_filter_list.gid<<std::endl;

    output_stats_file<<"\n--- Filter Stats ---\n";
    output_stats_file<<"filter_gid: "<<master_filter_list.gid<<std::endl;
    int size = master_filter_list.filters.size();
    //std::cout<<"filter list size: "<<size<<std::endl;

    output_stats_file<<"filter list size: "<<size<<std::endl;
    int n_total = 0, n_min = INT16_MAX, n_max = -1;
    float f_total = 0, f_min = FLT_MAX, f_max = -1;
    int promising_filters = 0;
    int filter_sizes_count[100] = {0};  // ASSUME MAX SIZE OF FILTER TO BE 100 :)
    int num_dilated = 0;
    std::for_each(master_filter_list.filters.begin(), master_filter_list.filters.end(),
            [&n_total, &n_min, &n_max, &f_total, &f_min, &f_max, &promising_filters, &filter_sizes_count, &num_dilated]
            (const FilterMap::value_type & filter_item)
            {
                n_total+= filter_item.second.numerosity;
                if(n_min > filter_item.second.numerosity) n_min = filter_item.second.numerosity;
                if(n_max < filter_item.second.numerosity) n_max = filter_item.second.numerosity;
                f_total+= filter_item.second.fitness;
                if(f_min > filter_item.second.fitness) f_min = filter_item.second.fitness;
                if(f_max < filter_item.second.fitness) f_max = filter_item.second.fitness;
                if(filter_item.second.fitness>0) promising_filters++;
                filter_sizes_count[filter_item.second.filter_size]++;
                if(filter_item.second.is_dilated) num_dilated++;
            });
    //std::cout<<"avg numerosity: "<<n_total/(float)size<<" , max numerosity: "<<n_max<<" , min numerosity: "<<n_min<<std::endl;
    //std::cout<<"avg fitness: "<<f_total/(float)size<<" , max fitness: "<<f_max<<" , min fitness: "<<f_min<<std::endl;
    //std::cout<<"promising filters: "<<promising_filters<<std::endl;

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
    // create a vector of promising filters
    std::vector<int> promising_filters;
    for(auto it=master_filter_list.filters.begin(); it!=master_filter_list.filters.end();it++){
        if(it->second.fitness > 0){
            promising_filters.push_back(it->second.id);
        }
    }
    // select two filters for tournament
    int size = promising_filters.size();
    if(size == 0) return -1;
    int selected[2];
    selected[0] = promising_filters[irand(size)];
    selected[1] = promising_filters[irand(size)];
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
        output_filter_to_file(item.second, output_filter_file);
        if(item.second.fitness > 0) {
            output_filter_to_file(item.second, output_promising_filter_file);
        }
    }
}
