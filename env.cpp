#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <cstring>
#include <assert.h>
#include "xcsMacros.h"
#include "configuration.h"
#include "env.h"
#include "codeFragment.h"
#include "filter_list.h"


CodeFragmentMap kb_cf;
FilterMap kb_filter;
//DataSource allData[1];
float lowerLimit[condLength];
float upperLimit[condLength];

void initializeInput(DataSource inputArray[],int numofRows){
    int featureIndex,docIndex;
    for(docIndex=0;docIndex<numofRows;docIndex++)
    {
        inputArray[docIndex].state = new float[condLength];
        for(featureIndex=0;featureIndex<condLength;featureIndex++)
        {
            inputArray[docIndex].state[featureIndex] = 0.0;
        }
        inputArray[docIndex].action = 0;
    }
}

float roundRealValue(float val, int num){
    float p = (float)pow(10.0,num);
    val = val * p;
    float tmp = roundf(val);
    return tmp/p;
}

float ScaleRange(float Value,float FromMinValue, float FromMaxValue, float ToMinValue,float ToMaxValue){
    if ((FromMaxValue - FromMinValue) + ToMinValue == 0){
        return 0;
    }else {
        return (Value - FromMinValue) * (ToMaxValue - ToMinValue) / (FromMaxValue - FromMinValue) + ToMinValue;
    }
}

void updateRange(DataSource data[],int totalRows){
    for(int docNo = 0; docNo<totalRows;docNo++)
        {
            for(int featureNo=0;featureNo<condLength;featureNo++)
            {
                data[docNo].state[featureNo] = ScaleRange(data[docNo].state[featureNo],lowerLimit[featureNo],upperLimit[featureNo],0.0,1.0);
                data[docNo].state[featureNo] = roundRealValue(data[docNo].state[featureNo],precisionDigits);
            }
        }
}

// This is my improved re-written function for laoding data
int LoadDataFromFile(DataSource data[], const char inputFile[])
{
    std::string line;
    std::ifstream infile(inputFile);

    if(infile.is_open()){
        while(getline(infile, line)){
            std::stringstream line_stream(line);

        }
    }else{
        std::string error("Error opening input file: ");
        error.append(inputFile).append(", could not load data!");
        throw std::runtime_error(error);
    }
}


// my function
void loadDataFromFile(DataSource data[], const char inputFile[], const int numInstances){

    int featureIndex = 0;
    float tfIDF = 0.0;
    int docIndex = 0;

    std::string line;
    std::ifstream infile(inputFile);


    if (infile.is_open())
    {
        while (!infile.eof())
        {
            if(docIndex<numInstances)
            {
                featureIndex = 0;
                std::string line;
                getline(infile,line);
                std::stringstream singleValue(line);
                //while(singleValue>>featureIndex)
                while(featureIndex <= condLength)
                {
                    singleValue>>tfIDF;
                    if(featureIndex==condLength)
                        //if(featureIndex==5000)
                        data[docIndex].action = int(tfIDF);
                    else
                    {
                        //tfIDF = roundRealValue(tfIDF,precisionDigits);
                        data[docIndex].state[featureIndex] = tfIDF;

                        if(tfIDF>upperLimit[featureIndex])
                            upperLimit[featureIndex] = tfIDF;
                        if(tfIDF<lowerLimit[featureIndex])
                            lowerLimit[featureIndex] = tfIDF;
                    }
                    featureIndex++;
                }
                // std::cout<<std::endl;
                docIndex++;
            }else{
                break;
            }
        }
    }else{
        std::string error("Error opening input file: ");
        error.append(inputFile).append(", could not load data!");
        throw std::runtime_error(error);
    }
    infile.close();
    //updateRange(inputArray,totalNumInstances);
}

int load_filter(std::string filter_file_name, FilterMap& filters)
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
        getline(filter_file, line);
        std::stringstream line2(line);
        line2 >> str;
        f.lower_bounds.reserve(f.filter_size*f.filter_size);
        f.lower_bounds.assign(f.filter_size*f.filter_size, -1);
        f.upper_bounds.reserve(f.filter_size*f.filter_size);
        f.upper_bounds.assign(f.filter_size*f.filter_size, -1);
        for(int i=0; i<f.filter_size*f.filter_size; i++){
            line2 >> f.lower_bounds[i];
        }
        getline(filter_file, line);
        std::stringstream line3(line);
        line3 >> str;
        for(int i=0; i<f.filter_size*f.filter_size; i++){
            line3 >> f.upper_bounds[i];
        }
        filters[f.id] = f;
        if(f.id > loaded_gid){
            loaded_gid = f.id;
        }
    }
    return loaded_gid;
}

int load_code_fragment(std::string cf_file_name, CodeFragmentMap& code_fragments)
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
        int index = 0, leaf_index = 0;
        while(!line1.eof()){
            std::string token;
            line1>>token;
            // last token is "" that needs to be handled
            if(token.empty()) break;
            if(token.substr(0,1) == "D"){ // this is filter id
                int filter_id = std::stoi(token.substr(1));
                cf.filter_id[leaf_index] = filter_id;
                cf.reverse_polish[index] = leaf_index;
                leaf_index++;
            }else{ // this is operator
                cf.reverse_polish[index] = str_to_opt(token);
            }
            index++;
        }
        cf.reverse_polish[index] = OPNOP; // terminate the reverse polish
        cf.num_filters = leaf_index;
        code_fragments[cf.cf_id] = cf;
        if(cf.cf_id > loaded_cf_gid){
            loaded_cf_gid = cf.cf_id;
        }
    }
    return loaded_cf_gid;
}

void load_kb(std::string kb_cf_file_name, std::string kb_filter_file_name) {

    load_filter(kb_filter_file_name, kb_filter);
    load_code_fragment(kb_cf_file_name, kb_cf);
}

int load_classifier(std::string classifier_file_name, ClassifierMap& pop, CodeFragmentMap& code_fragments)
{
    int loaded_cl_gid = -1;
    std::string line;
    std::ifstream cl_file(classifier_file_name);
    if (!cl_file.is_open()) {
        std::string error("Error opening input file: ");
        error.append(classifier_file_name).append(", could not load data!");
        throw std::runtime_error(error);
    }

    while(getline(cl_file, line)) {
        // load classifier
        Classifier cl;
        int num_cf=-1;
        std::string str;
        std::stringstream line1(line);
        line1>>str;
        line1>>cl.id;
        line1>>str;
        line1>>cl.numerosity;
        line1>>str;
        line1>>cl.experience;
        line1>>str;
        line1>>num_cf;
        line1>>str;
        line1>>cl.fitness;
        line1>>str;
        line1>>cl.accuracy;
        line1>>str;
        line1>>cl.prediction;
        line1>>str;
        line1>>cl.predictionError;
        line1>>str;
        line1>>cl.actionSetSize;
        line1>>str;
        line1>>cl.timeStamp;
        line1>>str;
        line1>>cl.action;

        getline(cl_file, line);
        std::stringstream line2(line);
        for(int i=0; i<num_cf; i++){
            int cf_id = -1;
            line2>>cf_id;
            cl.cf.push_back(code_fragments[cf_id]);

        }
        pop[cl.id] = cl;
        if(loaded_cl_gid < cl.id){
            loaded_cl_gid = cl.id;
        }
    }
    return loaded_cl_gid;
}