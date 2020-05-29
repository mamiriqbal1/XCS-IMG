//
// Created by amir on 5/25/20.
//

#include "filter_list.h"
#include <algorithm>

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
 * If yes then return the id of more general filter.
 * If no then return the id of new filter.
 */
int add_filter(Filter filter_to_add){
    // use find_if from <algorithm> to find a more general filter.
    std::forward_list<Filter>::iterator general_filter_iterator;
    general_filter_iterator = std::find_if(master_filter_list.filters.begin(), master_filter_list.filters.end(),
            [&filter_to_add](Filter& filter_from_list) -> bool
            {
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
        master_filter_list.filters.push_front(filter_to_add);
        return filter_to_add.id;
    }else{
        return general_filter_iterator->id;
    }
}



/*
 * Returns an object of Filter against the ID.
 * The copied filter can be used for crossover or mutation
 * Throws std::runtime_error in case id is not found
 */
const Filter& get_filter(int filter_id){
    std::forward_list<Filter>::iterator iterator;
    iterator = std::find_if(master_filter_list.filters.begin(), master_filter_list.filters.end(),
            [filter_id](Filter& filter_from_list) -> bool
            {
                return filter_from_list.id == filter_id;
            });
    if(iterator == master_filter_list.filters.end()){  // did not find any general filter
        throw std::runtime_error("Filter ID not found in the list");
    }else{
        return *iterator;
    }
}


