#include <iostream>
#include <string.h>
#include <string>
#include <assert.h>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <math.h>
#include "xcsMacros.h"
#include "codeFragment.h"
#include "env.h"
#include "filter_list.h"
#include <unordered_map>
#include <utility>
#include <algorithm>
//using namespace std;

char globalBuf[1000];
int numPreviousCFs = 0;     //When to set
int startingPreviousCFID = 0;   //what is this value

CodeFragment *previousCFPopulation;

typedef std::unordered_map<int, bool> ImageMap;
typedef std::unordered_map<int, ImageMap> FilterMap;
FilterMap evaluation_map;
FilterMap evaluation_validation_map;
int map_hits = 0;


void print_filter_evaluation_stats(std::ofstream &output_stats_file) {
    std::cout<<"--- Filter Evaluation Stats ---\n";

    output_stats_file<<"--- Filter Evaluation Stats ---\n";
    int size = std::distance(evaluation_map.begin(), evaluation_map.end());
    int total = 0, positive = 0, min = INT16_MAX, max = -1;
    std::for_each(evaluation_map.begin(), evaluation_map.end(),
                  [&total, &positive, &min, &max](const FilterMap::value_type & item)
                  {
                      int size = item.second.size();
                      int yes = std::count_if(item.second.begin(), item.second.end(),
                              [](const ImageMap::value_type & item2)
                              {
                                return item2.second;
                              });
                      total += size;
                      positive += yes;
                      if(min > yes) min = yes;
                      if(max < yes) max = yes;
                 });
    std::cout<<"map hits: "<<map_hits<<std::endl;
    std::cout<<"total evaluations recorded: "<<total<<" , total evaluated filters: "<<size<<std::endl;
    std::cout<<"avg positive: "<<positive/(float)size<<" , max positive: "<<max<<" , min positive: "<<min<<std::endl;
    std::cout<<"--- Filter Evaluation Stats ---\n\n";

    output_stats_file<<"map hits: "<<map_hits<<std::endl;
    output_stats_file<<"total evaluations recorded: "<<total<<" , total evaluated filters: "<<size<<std::endl;
    output_stats_file<<"avg positive: "<<positive/(float)size<<" , max positive: "<<max<<" , min positive: "<<min<<std::endl;
    output_stats_file<<"--- Filter Evaluation Stats ---\n\n";
}

void initializeCFPopulation(FILE *cfReadingFilePointer)//, FILE *code_fragment_file)
{
    if(use_kb)
    {
        printf("\nLoading Previous Population! Please wait ....\n");
       getPreviousCFPopulation(cfReadingFilePointer);
    }
}

void getPreviousCFPopulation(FILE *cfReadingFilePointer)
{
    //read from file
    //std::cout<<"read prev\n";
    //int cntr = 0;
   int lineSize = 1000; // large enough to cover spaces in the CF
    char previousCFContents[lineSize];
    CodeFragment previousCF;
    int numReadPreviousCFs = 0;
    if( fgets(previousCFContents,lineSize,cfReadingFilePointer) !=NULL )  //first line is the number of CFs
    {
        numPreviousCFs = atoi(previousCFContents); // maximum number of stored CFs
        //printf("\n%d\n",numPreviousCFs);

        previousCFPopulation = new CodeFragment[numPreviousCFs];
        //previousCFPopulation = (CodeFragment*)malloc(numPreviousCFs*sizeof(CodeFragment));
    }
    else
    {
        printf("\nError in reading previous CFs.\n");
        exit(0);
    }
    if( fgets(previousCFContents,lineSize,cfReadingFilePointer) !=NULL )  //second line is the starting ID number of previous CFs
    {
        startingPreviousCFID = atoi(previousCFContents);
        //std::cout<<previousCFContents;

    }
    else
    {
        printf("\nError in reading previous CFs.\n");
        exit(0);
    }

    while(fgets(previousCFContents, lineSize, cfReadingFilePointer) != NULL)   //get each line from the file
    {

        //std::cout<<previouwcsCFContents<<std::endl;
        // strip trailing '\n' if it exists
        int len = strlen(previousCFContents)-1;
        //printf("\n%d\n",len);
        if(previousCFContents[len] == '\n')
        {
            previousCFContents[len] = 0;
        }
        //printf("\n%s",previousCFContents);

        char *pch = NULL;
        opType c = 0;
        int i = 0;
        int lfNum = 0;
        createNewCF(numReadPreviousCFs + condLength, previousCF); // start CFs ID numbers from condLength to avoid confusion with CFs and code_fragment bits
        //tokenize the previousCFContents
        //int counting = 0;

        pch = strtok(previousCFContents,",");
        while (pch != NULL)
        {
            //printf("\n%d\n",*pch);
            //std::cout<<counting++<<"\n";
            if(pch[0] == 'D')  //D0, D1, D2, etc
            {
                //int length = strlen(pch);
                //std::cout<<"\n"<<length<<"\n";
                //std::cout<<pch;
                // todo: loading of previous CF to be adapted
                //previousCF.leaf[lfNum] = leafNode(pch);
                previousCF.reverse_polish[i++] = lfNum;
                lfNum++;
                //std::cout<<"Pch with D: "<<pch<<std::endl;

                //std::cout<<"\n"<<pch<<"\n";
                //getchar();
            }
            else
            {
                //printf("\n%d\n",*pch);
                c = getOpType(pch);
                previousCF.reverse_polish[i++] = c;

                //std::cout<<"C: "<<c<<std::endl;


            }
            //printf("\n%d\n",c);
            //std::cout<<pch;

            if(c == OPNOP)
            {
                //std::cout<<"C: "<<c<<std::endl;
                if(!isExists(previousCF,previousCFPopulation,numReadPreviousCFs)) //insert the code fragment
                {
                    //std::cout<<"is exit false: "<<std::endl;

                    //printf("\n%s",previousCFContents);
                    //printCF(previousCF);
                    memmove(&previousCFPopulation[numReadPreviousCFs],&previousCF,sizeof(CodeFragment));
                    numReadPreviousCFs++;
                    //printf("\n%d\n",numReadPreviousCFs);
                    //std::cout<<"C: "<<c<<std::endl;
                }
                break;
            }
            //printf("\n%d\n",*pch);
            //pch = strtok (previousCFContents,",");

            pch = strtok (NULL,",");
            //std::cout<<"\n"<<pch<<"\n";
            //getchar();


        } //end while pch
        //std::cout<<"Line No: "<<cntr++<<std::endl;
        //getchar();
    } //end while fgets
    //printf("\n%s",previousCFContents);
    numPreviousCFs = numReadPreviousCFs; //number of distint CFs
    printf("\nPrevious CFs loaded: %d\n",numPreviousCFs);
    fclose(cfReadingFilePointer); //close
    //exit(0);


}

opType getOpType(char str[])
{
    opType ret;
    if(strcmp(str,"o")==0) return OPNOP;
    if(strcmp(str,"&")==0) return OPAND;
    if(strcmp(str,"|")==0) return OPOR;
    if(strcmp(str,"d")==0) return OPNAND;
    if(strcmp(str,"r")==0) return OPNOR;
    if(strcmp(str,"~")==0) return OPNOT;
/*
    if(str[0] == 'D')  //D0, D1, D2, etc
    {
        ret = leafOpCode(atoi(str+1));
    }*/
    else if(str[0] == 'C')  //CF from a previous level
    {
        ret = leafOpCode( atoi(str+2) + (condLength - startingPreviousCFID) );
        //std::cout<<"CF: "<<ret<<std::endl;
    }
    else
    {
        printf("\nError in getOpType ... \n");
        exit(0);
    }
    return ret;
}

bool isExists(CodeFragment &newCF, CodeFragment *cfPopulation, int numCFs)
{

    for(int i=0; i<numCFs; i++)
    {

        if(equalTwoCFs(newCF,cfPopulation[i]))
        {
            //std::cout<<"i"<<std::endl;
            return true;
        }
    }
    return false;
}

bool equalTwoCFs(CodeFragment &cf1, CodeFragment &cf2)
{

    for(int i=0; i<cfMaxLength; i++)
    {
        if(cf1.reverse_polish[i] == cf2.reverse_polish[i])
        {
            if(0<=cf1.reverse_polish[i] && cf1.reverse_polish[i] < numLeaf)
            {
                if(cf1.filter_id[cf1.reverse_polish[i]] != cf2.filter_id[cf2.reverse_polish[i]])
                    return false;
            }
        }else
        {
            return false;
        }
    }

    return true;
}

int getNumPreviousCFs()
{
    return numPreviousCFs;
}

bool isDontcareCF(CodeFragment &cf)
{
    return (cf.cf_id == -1) ? true : false;
    //return (cf.reverse_polish[0] == OPUNITY) ? true : false;
}

int numberOfNonDontcares(CodeFragment cf[])  //returns the number of specific CFs in cf
{
    int count=0;
    for(int i=0; i<clfrCondLength; i++)
    {
        if(!isDontcareCF(cf[i]))
        {
            count++;
        }
    }
    return count;
}

void validateDepth(opType* cf, opType* end)
{
    //display 'cf' for debugging
    /*
    int i=0;
    char* temp = NULL;
    while(cf[i]!=OPNOP){
    	temp = opchar(cf[i++]);
    	printf("%s ",temp);
    }
    */
    const opType* start = cf;
    opType* p = end -1;
    int a =1;
    int depth = -1;
    DepthMax(start,&p,a,depth);
    //printf("Depth: %d\n",depth);
    assert(depth<=cfMaxDepth);
}

void DepthMax(const opType* const end,opType** prog, int& argstogo, int& depth)
{
    if(*prog<end || argstogo<=0) return;
    const int a = getNumberOfArguments(*prog[0]);
    argstogo += a-1;
    *prog = (*prog-1);
    depth++;
    const int d0 = depth;
    for(int i=0; i<a; i++)
    {
        int d1 = d0;
        DepthMax(end,prog,argstogo,d1);
        depth = (depth>d1)? depth : d1; //max(depth,d1);
    }
}

void createNewCF(int id, CodeFragment &cf)
{
    for(int i=0; i<cfMaxLength; i++)
    {
        cf.reverse_polish[i] = OPNOP;
    }
    cf.cf_id = id;
}

void storeCFs(ClassifierMap &pop, FILE *cfWritingFilePointer)
{
    double avgFitness = getAvgFitness(pop);
    int numFitterCFs = getNumFitterCFs(pop, avgFitness);
    int firstCFID = 0;

    char *buf;
    int len;

    // number of CFs
    //sprintf(buf,"%d\n",numFitterCFs); fwrite(buf,strlen(buf),1,code_fragment_file); // number of CFs
    len = snprintf(NULL,0,"%d\n",numFitterCFs);

    if(!(buf = (char*)malloc((len + 1) * sizeof(char))))
    {
        printf("\nError in file writing ...\n");
        exit(0);
    }
    len = snprintf(buf,len+1,"%d\n",numFitterCFs);

    fwrite(buf,strlen(buf),1,cfWritingFilePointer);
    free(buf);

    if(use_kb)
    {
        firstCFID = previousCFPopulation[0].cf_id;
    }
    // cf_id of first CF
    //sprintf(buf,"%d\n",firstCFID); fwrite(buf,strlen(buf),1,code_fragment_file); // cf_id of first CF
    len = snprintf(NULL,0,"%d\n",firstCFID);
    if(!(buf = (char*)malloc((len + 1) * sizeof(char))))
    {
        printf("\nError in file writing ...\n");
        exit(0);
    }


    len = snprintf(buf,len+1,"%d\n",firstCFID);

    fwrite(buf,strlen(buf),1,cfWritingFilePointer);
    free(buf);
  //  std::cout<<"numPrevCFs: "<<numPreviousCFs<<"\n";
//    getchar();

    for(int i=0; i<numPreviousCFs; i++)  //first store CFs from previous level problems
    {
        outprog(previousCFPopulation[i], cfWritingFilePointer);

        //sprintf(buf," ---------> %d",previousCFPopulation[i].cf_id); fwrite(buf,strlen(buf),1,code_fragment_file);
        len = snprintf(NULL,0," ---------> %d",previousCFPopulation[i].cf_id);
        if(!(buf = (char*)malloc((len + 1) * sizeof(char))))
        {
            printf("\nError in file writing ...\n");
            exit(0);
        }
        len = snprintf(buf,len+1," ---------> %d",previousCFPopulation[i].cf_id);

        fwrite(buf,strlen(buf),1,cfWritingFilePointer);
        free(buf);

        fwrite("\n",strlen("\n"),1,cfWritingFilePointer);

        fflush(cfWritingFilePointer);
    }


    for(auto & item : pop)  //CFs from the current problem
    {
        if(item.second.fitness <= avgFitness)  //store CFs from classifiers with fitness > average fitness of the classifiers population
        {
            break;
        }
        for(int i=0; i<clfrCondLength; i++)
        {
        if(!isDontcareCF(item.second.code_fragment[i]))
            {
                //outprog(set->classifier->code_fragment[i].reverse_polish,cfMaxLength,code_fragment_file);
                outprog(item.second.code_fragment[i], cfWritingFilePointer);
               //sprintf(buf," ---------> %d",set->classifier->code_fragment[i].cf_id); fwrite(buf,strlen(buf),1,code_fragment_file);
                len = snprintf(NULL,0," ---------> %d",item.second.code_fragment[i].cf_id);
                if(!(buf = (char*)malloc((len + 1) * sizeof(char))))
                {
                    printf("\nError in file writing ...\n");
                    exit(0);
                }
                len = snprintf(buf,len+1," ---------> %d",item.second.code_fragment[i].cf_id);

                fwrite(buf,strlen(buf),1,cfWritingFilePointer);
                free(buf);

                fwrite("\n",strlen("\n"),1,cfWritingFilePointer);
                fflush(cfWritingFilePointer);
            }
        }
    }

    fwrite("\n\n",strlen("\n\n"),1,cfWritingFilePointer);
    fflush(cfWritingFilePointer);
    delete[] previousCFPopulation;
}

// new function for setting filter bounds
void create_new_filter_from_input(Filter& filter, float *state)
{
    // randomly selects a position in the image to create filter bounds
    //int filter_size = (int)sqrt(numLeaf);  // filter size
    int new_filter_size = filter_sizes[irand(num_filter_sizes)]; // select a filter size randomly
    bool is_dilated = false;
    if(allow_dilated_filters){
        is_dilated = irand(2) != 0;
    }
    int effective_filter_size = new_filter_size;
    if(is_dilated){
        effective_filter_size = new_filter_size + new_filter_size -1;
    }
    int step = is_dilated ? 2 : 1;  // this will be used to map normal coordinates to dilated coordinates
    float pixel_values[new_filter_size*new_filter_size];
    float sum = 0;
    do{
        sum = 0;
        int index = 0;
        int x_position = irand(image_width-effective_filter_size);
        int y_position = irand(image_height-effective_filter_size);
        for(int y=y_position; y<y_position+effective_filter_size; y+=step){
            for(int x=x_position; x<x_position+effective_filter_size; x+=step){
                pixel_values[index] = state[y*image_width+x];
                sum += pixel_values[index];
                index++;
            }
        }
    }while(sum <= 0.1); // get to some interesting area in the image. All blanks will be ignored.

    for(int i=0; i<new_filter_size*new_filter_size; i++){
        float delta = drand();
        filter.lower_bounds[i] = roundRealValue(fmax(pixel_values[i] - delta, 0), precisionDigits);
        filter.upper_bounds[i] = roundRealValue(fmin(pixel_values[i] + delta, 1),precisionDigits);
    }
    filter.filter_size = new_filter_size;
    filter.is_dilated = is_dilated;
}

// new function that randomly selects a position on filter and create a matching filter at that position
CodeFragment addLeafCF(CodeFragment &cf, float *state){

    // count number of leaves
    int count = 0;
    for(int i=0; i<cfMaxLength; i++){
        const opType opcode = cf.reverse_polish[i];
        if(opcode == OPNOP)
        {
            break;
        }
        if(0<=opcode && opcode<condLength)  //code_fragment bit
        {
            count++;
        }
    }
    cf.num_filters = count;
    //cf.filter = new Filter[count];

    int leafNum = 0;
    for(int i=0; i<cfMaxLength; i++)
    {
        const opType opcode = cf.reverse_polish[i];
        if(opcode == OPNOP)
        {
            break;
        }
        if(0<=opcode && opcode<condLength)  //code_fragment bit
        {
            cf.reverse_polish[i] = leafNum;
            // With probability p_ol if there is promising filter available in filter store then use it
            int id = -1;
            if(drand() < p_ol){
                id = get_promising_filter_id();
            }
            if(id != -1){
               cf.filter_id[leafNum] = id;
            }else {
                Filter new_filter;
                create_new_filter_from_input(new_filter, state);
                cf.filter_id[leafNum] = add_filter(new_filter);
            }
            leafNum++;
        }
    }
    return cf;
}


bool evaluate_filter_actual(const Filter& filter, float state[])
{
    int step = filter.is_dilated ? 2 : 1;  // this will be used to map normal coordinates to dilated coordinates
    int effective_filter_size = filter.filter_size;
    if(filter.is_dilated){
        effective_filter_size = filter.filter_size + filter.filter_size -1;
    }
    bool match_failed = false; // flag that controls if the next position to be evaluated when current does not match
    int k = 0;
    int l = 0;
    for(int i=0; i<image_height - effective_filter_size; i++){  // i is image y coordinate
        for(int j=0; j<image_width - effective_filter_size; j++){  // j is image x coordinate
            match_failed = false;
            for(k=0; k<filter.filter_size && !match_failed; k+=step){  // k is filter y coordinate
                for(l=0; l<filter.filter_size && !match_failed; l+=step){  // l is filter x coordinate
                    if(state[i*image_width+j + k*image_width+l] < filter.lower_bounds[k*filter.filter_size+l]
                    || state[i*image_width+j + k*image_width+l] > filter.upper_bounds[k*filter.filter_size+l]){
                        match_failed = true;
                    }
                }
            }
            if(!match_failed){
                return true;
            }
        }
    }
    return false;
}

void update_evaluation_cache(std::forward_list<int>& removed_filters){
    std::forward_list<int>::iterator it;
    for(it = removed_filters.begin(); it != removed_filters.end(); it++){
        evaluation_map.erase(*it);
        evaluation_validation_map.erase((*it));
    }
}

bool evaluate_filter(const Filter& filter, float state[], int cl_id, int img_id, bool train)
{
    // if cl_id or img_id is -1 then do not check evaluation map otherwise check for prior results
    if(img_id >=0){
        // return prior result if it is found
        if(train){
            ImageMap& inner_map = evaluation_map[filter.id];
            if(inner_map.count(img_id) > 0){  // the image evaluation exist - unordered map always return  1
                map_hits++;
                return inner_map[img_id];
            }
        }else {
            ImageMap &inner_map = evaluation_validation_map[filter.id];
            if (inner_map.count(img_id) > 0) {  // the image evaluation exist - unordered map always return  1
                map_hits++;
                return inner_map[img_id];
            }
        }
    }

    bool evaluation = evaluate_filter_actual(filter, state);

    // set hasmap entry for re-using evaluation
    if(img_id >=0) {
        if(train){
            ImageMap& inner_map = evaluation_map[filter.id];
            inner_map[img_id] = evaluation;
        }else{
            ImageMap& inner_map = evaluation_validation_map[filter.id];
            inner_map[img_id] = evaluation;
        }
    }
    return evaluation;
}

bool mutate_cf(CodeFragment &cf){
    int functions_index[cfMaxLength];  // collect indices of functions for possible mutation
    int functions_i = 0;
    for(int i=0; /*i<cfMaxLength*/; i++){
        const opType opcode = cf.reverse_polish[i];
        if(opcode == OPNOP){
            break;
        }else if(getNumberOfArguments(opcode) > 1){  // only select binary functions for now
           functions_index[functions_i++] = i;
        }
    }
    if(functions_i == 0) return false;
    // number of functions in the CF are functions_i
    int k = irand(functions_i); // select a function index at random for mutation
    opType new_function=0;
    do{
        new_function = randomFunction(); // there is a possiblity of selecting the same function
        if(getNumberOfArguments(new_function) == getNumberOfArguments(cf.reverse_polish[functions_index[k]])){
            cf.reverse_polish[functions_index[k]] = new_function;
            break;
        }
    }while(true);
    return true;
}


int evaluateCF(CodeFragment &cf, float *state, int cl_id, int img_id, bool train){
    int stack[cfMaxStack];
    stack[0] = 0;
    int SP = 0;
    int tmptmpfval1 = 0;
    for(int i=0; /*i<cfMaxLength*/; i++)
    {
        const opType opcode = cf.reverse_polish[i];
        if(opcode == OPNOP)
        {
            break;
        }
        if(isPreviousLevelsCode(opcode))  //CF from any previous level
        {

            int valueOfCF = evaluateCF(previousCFPopulation[opcode - condLength], state, train);
            stack[SP++] = valueOfCF;
        }
        else if(0<=opcode && opcode<condLength)  //code_fragment bit
        {

            //if(cf.leaf[opcode].lowerBound<=state[cf.leaf[opcode].featureNumber] && state[cf.leaf[opcode].featureNumber]<=cf.leaf[opcode].upperBound)
            if(evaluate_filter(get_filter(cf.filter_id[opcode]), state, cl_id, img_id, train))
            {
                stack[SP++] = 1;   //changed
            }
            else
            {
                stack[SP++] = 0;   //changed
            }

        }
        else if(opcode == OPNOT)
        {
            const int sp = stack[--SP];
            stack[SP++] = (!sp)?1:0;
        }
        else
        {
            const int sp2 = stack[--SP];
            const int sp1 = stack[--SP];
            switch(opcode)
            {
            case OPAND:
                stack[SP++] = (sp1&&sp2)?1:0;
                break;
            case OPOR:
                stack[SP++] = (sp1||sp2)?1:0;
                break;
            case OPNAND:
                stack[SP++] = (sp1&&sp2)?0:1;
                break;
            case OPNOR:
                stack[SP++] = (sp1||sp2)?0:1;
                break;
            }//end switch
        }
    }
    int value = stack[--SP];
    //std::cout<<"SP: "<<SP<<"\n";
    assert(SP==0);

    return value;
}

/*
int evaluateCF(CodeFragment cf, float state[], int cl_id, int img_id){
    // if cl_id or img_id is -1 then do not check evaluation  map otherwise check for prior results
    if(cl_id >=0 and img_id >=0){
        // return prior result if it is found
        if(evaluation_map.find(pair(cl_id, img_id)) != evaluation_map.end()){  // if found
            return evaluation_map[pair(cl_id, img_id)];
        }
    }
    // set featureNumber appropriately and then call evaluateCF_old
    int size = (int)sqrt(numLeaf);  // filter size
    int index[size*size];
    bool done = false;
    int return_value = 0;

    for(int i=0; i<image_height - size && !done; i++){
        for(int j=0; j<image_width - size && !done; j++){
           for(int k=0; k<size; k++){
               for(int l=0; l<size; l++){
                   index[k*size+l] = i*image_width+j + k*image_width+l;
                   cf.leaf[k*size+l].featureNumber = i*image_width+j + k*image_width+l;
               }
           }
           if(evaluateCF_old(cf, state) == 1){
               done = true;
               return_value = 1;
           }else{
               continue;
           }
        }
    }
    // reset after evaluation
    for(int m=0; m<size*size; m++){
        cf.leaf[m].featureNumber = 0;
    }
    // set hasmap entry for re-using evaluation across epochs
    if(cl_id >=0 and img_id >=0) {
        evaluation_map[pair(cl_id, img_id)] = return_value;
    }
    return return_value;
}
*/

bool isPreviousLevelsCode(const opType code){
    if(use_kb)
    {
        return (previousCFPopulation[0].cf_id <= code && code < (previousCFPopulation[0].cf_id + numPreviousCFs) ) ? true : false;
    }
    return false;
}

// ######################################## tree operations #################################################

inline int getNumberOfArguments(const opType code){
    switch(code)
    {
    case OPAND:
    case OPOR:
    case OPNAND:
    case OPNOR:
        return 2;
    case OPNOT:
        return 1;
    default:
        return 0;
    }//end switch code
}

inline opType leafOpCode(const int r){
    return r;
}

inline opType randomLeaf(){
    opType leaf = OPNOP;
    if(numPreviousCFs==0 || !use_kb)
    {
        leaf = 0; // feature number is no more used irand(condLength);
        //printf("leaf_1 %d\n",leaf);
        return leaf;
    }
    double p = drand();
    if(p < 0.5)
    {
        leaf =  0; // feature number is no more irand(condLength);
        //printf("leaf_2 %d\n",leaf);
        return leaf;
    }
    int n = irand(numPreviousCFs);
    leaf = previousCFPopulation[n].cf_id;
    //printf("n: %d leaf_3 %d\n",n,leaf);
    return leaf;
}

int validLeaf(const opType opcode){
    if( 0<=opcode && opcode<condLength )
    {
        return opcode;
    }
    if(isPreviousLevelsCode(opcode))
    {
        return opcode;
    }
    return -1;
}

inline opType randomFunction(){
    return functionCodes[irand(totalFunctions)];
}
//generate reverse polish
opType* randomProgram(opType prog[], const int isfull, const int maxDepth, const int minDepth){
    //if reached max depth or probabilistically reached mindepth
    if( maxDepth<=0 || ( (!isfull) && minDepth<=0 && irand(2) ) )
    {
        *prog = randomLeaf();
        return prog+1;
    }
    else
    {
        opType* pc = prog;
        opType newFunction;
        newFunction = randomFunction();
        const int numArgs = getNumberOfArguments(newFunction);
        for(int i=0; i<numArgs; i++)
        {
            pc = randomProgram(pc,isfull,maxDepth-1,minDepth-1);
        }
        *pc = newFunction;
        return pc+1;
    }
}//end randomProgram
// --------------------------------- Start Display Functions -----------------------------

char* leafname(const opType code){
    if(0<=code && code<condLength)
    {
        sprintf(globalBuf,"D%d ",validLeaf(code));
    }
    else if(isPreviousLevelsCode(code))
    {
        sprintf(globalBuf,"CF_%d ",validLeaf(code));
    }
    else
    {
        printf("\nERROR! invlalid leaf name\n");
        exit(0);
    }
    return globalBuf;
}//end leafname

char* leafname(const Leaf leaf){
    opType code;
    code = leaf.featureNumber;
//    std::cout<<"leaf No: "<<leaf.featureNumber<<"\n";
    if(0<=code && code<condLength)
    {
        sprintf(globalBuf,"D%d ",validLeaf(code));
    }
    else if(isPreviousLevelsCode(code))
    {
        sprintf(globalBuf,"CF_%d ",validLeaf(code));
    }
    else
    {
        printf("\nERROR! invlalid leaf name\n");
        exit(0);
    }
    return globalBuf;
}//end leafname


char* opchar(const opType code){
    switch(code)
    {
    case OPAND:
        return (char*)"&,";
    case OPOR:
        return (char*)"|,";
    case OPNAND:
        return (char*)"d,";
    case OPNOR:
        return (char*)"r,";
    case OPNOT:
        return (char*)"~,";
    case OPNOP:
        return (char*)"o,";
    //case OPUNITY:
    //	return (char*)"1 ";
    default:
        sprintf(globalBuf,"[%d!!]",code);
        return globalBuf;
    }//end switch code
}//end opchar

void outprog_bin(const opType* prog, int size){
    for(int j = 0; j<size; j++)
        printf("%d ", prog[j]);
}

char* leafInterval(const Leaf leaf){
    opType opcode = leaf.featureNumber;
    std::string leafInv, str, str1,str2,str3;
    //std::stringstream ss;
    if( 0<=opcode && opcode<condLength )
    {
        std::stringstream ss;
        str = "D";
        ss << opcode;
        str1 = ss.str();
        std::stringstream ss1;
        ss1 << leaf.lowerBound;
        str2 = ss1.str();
        std::stringstream ss2;
        ss2 << leaf.upperBound;
        str3 = ss2.str();
        leafInv = str+""+str1+" "+str2+" "+str3+",";//+"]";

        char *cstr = new char[leafInv.length() + 1];
        strcpy(cstr, leafInv.c_str());
        //return (char*)leafInv;

        return cstr;
        //delete [] cstr;
    }
    /*if(isPreviousLevelsCode(opcode))
    {
        str = "CF";
        std::stringstream ssP;
        //str = "D";
        ssP << opcode;
        str1 = ssP.str();
        leafInv = str+""+str1+",";//+"]";

        char *cstr = new char[leafInv.length() + 1];
        strcpy(cstr, leafInv.c_str());
        //return (char*)leafInv;

        return cstr;
        //ss << opcode;
        //return (char*)opcode;
    }*/
    return (char*)"-1";
}


void outprog(CodeFragment &cf, FILE *fp) {
    opType code;

    for(int j = 0; j < cfMaxLength; j++)
    {
        char* temp = NULL;
        code = cf.reverse_polish[j];
        if(0<=code && code<numLeaf)
        {
            // todo: printing for filter to be implemented
            temp = new char[strlen("filter to be printed") + 1];
            strcpy(temp, "filter to be printed");

            //temp = "D"+leafDes;
            //sprintf(temp,"D%d ",leafDes);

            //sprintf(globalBuf,"D%d ",validLeaf(code));
            //temp = leafInterval(cf.leaf);
            delete temp;
        }
        else if (code>=condLength&&isPreviousLevelsCode(code))
        {
            std::string str = "CF";
            std::stringstream ssP;
            //str = "D";
            ssP << code;
            std::string str1 = ssP.str();
            std::string leafInv = str+""+str1+",";//+"]";

            temp = new char[leafInv.length() + 1];
            strcpy(temp, leafInv.c_str());

        }else
        {
            temp = opchar(code);
        }

        /*if(validLeaf(code)>=0)
        {
            temp = leafname(code);
        }*/

        //printf("%s ",temp);
        fwrite(temp,strlen(temp),1,fp);
        if(cf.reverse_polish[j] == OPNOP)
            break; //reduce size of output
    }
    //printf("\n");
    //fwrite("\n",strlen("\n"),1,fp);
}

inline std::string op_to_str(opType code)
{
    switch(code)
    {
        case OPAND:
            return std::string("&");
        case OPOR:
            return std::string("|");
        case OPNAND:
            return std::string("d");
        case OPNOR:
            return std::string("r");
        case OPNOT:
            return std::string("~");
        case OPNOP:
            return std::string("o");
        default:
            return std::string("[" + std::to_string(code) + "!!]");
    }//end switch code

}

void output_code_fragment_to_file(CodeFragment &cf, std::ofstream &output_code_fragment_file)
{
    output_code_fragment_file << cf.cf_id << " ";
    std::string str;
    opType code = 0;
    for(int i=0; i<cfMaxLength; i++){
        code = cf.reverse_polish[i];
        if(code == OPNOP){
            break;
        }else if(0<=code && code<numLeaf){  // if code is zero then it is a filter_id index
            output_code_fragment_file << "D" << cf.filter_id[code] << " "; // print filter id
        }else if (code>=condLength && isPreviousLevelsCode(code)){
            // output previous when implemented
        }else{
            // output code str
            output_code_fragment_file << op_to_str(code) << " ";
        }
    }
    output_code_fragment_file << std::endl;
}
