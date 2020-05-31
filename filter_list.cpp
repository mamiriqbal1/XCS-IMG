//
// Created by amir on 5/25/20.
//

#include "filter_list.h"
#include <algorithm>
#include <iostream>
#include <numeric>
#include <vector>

FilterList master_filter_list; // The main filter list that is maintained

/*
 * This functions is used to push a filter at the front of the list.
 * This function does not check if filter is subsumed by existing filters
 * It can be used for crossover and mutation to push filters temporarily and then pop if validation fails
 * It also returns ID of the filter so that all other functions with filter ID can work just like for normal filters
 */
int push_front(Filter filter_to_push){
    filter_to_push.id = master_filter_list.gid++;
    master_filter_list.filters.push_front(filter_to_push);
    return filter_to_push.id;
}

/*
 * This functions is used to pop the last pushed filter
 * It can be used for crossover and mutation to pop filters if validation fails
 */
void pop_front(){
    master_filter_list.filters.pop_front();
}

/*
 * Inserts a filter into the list.
 * First it checks if the filter is subsumed by another in the list.
 * Only a filter which is part of at least on "promising" classifier (i.e. with fitness > 0) can subsume the new filter
 * If yes then return the id of more general filter.
 * If no then return the id of new filter.
 */
int add_filter(Filter filter_to_add){
    // use find_if from <algorithm> to find a more general filter.
    std::forward_list<Filter>::iterator general_filter_iterator;
    general_filter_iterator = std::find_if(master_filter_list.filters.begin(), master_filter_list.filters.end(),
            [&filter_to_add](Filter& filter_from_list) -> bool
            {
                if(filter_from_list.fitness == 0) return false; // only a fitter filter can subsume;
                for(int i=0; i<filter_size*filter_size; i++){
                    if(filter_from_list.lower_bounds[i] > filter_to_add.lower_bounds[i]
                       || filter_from_list.upper_bounds[i] < filter_to_add.upper_bounds[i]){
                        return false;
                    }
                }
                return true;
            });
    if(general_filter_iterator == master_filter_list.filters.end()){  // did not find any general filter
        filter_to_add.id = master_filter_list.gid++;
        filter_to_add.numerosity = 1;
        master_filter_list.filters.push_front(filter_to_add);
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
    std::forward_list<Filter>::iterator it;
    it = std::find_if(master_filter_list.filters.begin(), master_filter_list.filters.end(),
            [filter_id](Filter& filter_from_list) -> bool
            {
                return filter_from_list.id == filter_id;
            });
    if(it == master_filter_list.filters.end()){
        throw std::runtime_error("Filter ID not found in the list: " + std::to_string(filter_id));
    }else{
        return *it;
    }
}

/*
 * This function first checks the numerosity of the filter and decreases it if it is more than 1 otherwise removes
 * the filter.
 */
void remove_filter(int filter_id){
    std::forward_list<Filter>::iterator iterator, previous_iterator;
    for(iterator = master_filter_list.filters.before_begin();
    iterator != master_filter_list.filters.end(); ){
       previous_iterator = iterator++;
       if(iterator->id == filter_id){
           if(iterator->numerosity > 1){
               iterator->numerosity--;
           }else {
               master_filter_list.filters.erase_after(previous_iterator);
           }
           break;
       }
    }
}

/*
 * This function resets the numerosity and fitness of all filters
 */
void reset_filter_stats(){
    std::for_each(master_filter_list.filters.begin(), master_filter_list.filters.end(),
                            [](Filter& filter_from_list)
                            {
                                filter_from_list.numerosity=0;
                                filter_from_list.fitness=0;
                            });
}

/*
 * This function removes any unused filters
 */
void remove_unused_filters(std::forward_list<int>& removed_filters){
    master_filter_list.filters.remove_if(
            [&removed_filters](Filter& filter_from_list) -> bool
            {
                if(filter_from_list.numerosity == 0){
                    removed_filters.push_front(filter_from_list.id);
                    return true;
                }else{
                    return false;
                }

            });

}

/*
 * Print stats about the filter list
 */
void print_filter_stats(){
    std::cout<<"\n--- Filter Stats ---\n";
    std::cout<<"gid: "<<master_filter_list.gid<<std::endl;
    int size = std::distance(master_filter_list.filters.begin(), master_filter_list.filters.end());
    std::cout<<"filter list size: "<<size<<std::endl;
    int n_total = 0, n_min = INT16_MAX, n_max = -1;
    int f_total = 0, f_min = INT16_MAX, f_max = -1;
    std::for_each(master_filter_list.filters.begin(), master_filter_list.filters.end(),
            [&n_total, &n_min, &n_max, &f_total, &f_min, &f_max](Filter& filter_from_list)
            {
                n_total+= filter_from_list.numerosity;
                if(n_min > filter_from_list.numerosity) n_min = filter_from_list.numerosity;
                if(n_max < filter_from_list.numerosity) n_max = filter_from_list.numerosity;
                f_total+= filter_from_list.fitness;
                if(f_min > filter_from_list.fitness) f_min = filter_from_list.fitness;
                if(f_max < filter_from_list.fitness) f_max = filter_from_list.fitness;

            });
    std::cout<<"avg numerosity: "<<n_total/(float)size<<" max numerosity: "<<n_max<<" min numerosity: "<<n_min<<std::endl;
    std::cout<<"avg fitness: "<<f_total/(float)size<<" max fitness: "<<f_max<<" min fitness: "<<f_min<<std::endl;
    std::cout<<"--- Filter Stats ---\n\n";
}