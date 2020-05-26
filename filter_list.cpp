//
// Created by amir on 5/25/20.
//

#include "filter_list.h"
#include <algorithm>

FilterList filter_list; // The main filter list that is maintained
Filter prior_filter; // The prior filter that will be checked in the predicate function




/*
 * This function checks if the filter provided in the argument by the container is more general
 * Than an already specified filter.
 *
 */
bool is_already_covered(Filter& filter_from_list){
    for(int i=0; i<filter_size*filter_size; i++){
        if(filter_from_list.lower_bounds[i] > prior_filter.lower_bounds[i]
        || filter_from_list.upper_bounds[i] < prior_filter.upper_bounds[i]){
            return false;
        }
    }
    return true;
}

/*
 * Inserts a filter into the list.
 * First it checks if the filter is subsumed by another in the list.
 * If yes then return the id of more general filter.
 * If no then return the id of new filter.
 */
int add_filter(Filter& filter_to_add){
    // use find_if from <algorithm> to find a more general filter.
    prior_filter = filter_to_add;  // set prior filter to be used by the predicate function
    std::forward_list<Filter>::iterator general_filter_iterator;
    general_filter_iterator = std::find_if(filter_list.filters.begin(), filter_list.filters.end(), is_already_covered);
    if(general_filter_iterator == filter_list.filters.end()){  // did not find any general filter
        filter_to_add.id = filter_list.gid++;
        filter_list.filters.push_front(filter_to_add);
        return filter_to_add.id;
    }else{
        return general_filter_iterator->id;
    }
}