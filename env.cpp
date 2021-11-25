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
#include "cf_list.h"


FloatVector lowerLimit; //[condLength];
FloatVector upperLimit; //[condLength];

void initialize_env()
{
    FloatVector l(condLength);
    lowerLimit = l;
    FloatVector  u(condLength);
    upperLimit = u;

}


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
    // do nothing
    return val;
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

void normalize_image(DataSource *data, int totalRows){
    for(int docNo = 0; docNo<totalRows;docNo++)
        {
            for(int featureNo=0;featureNo<condLength;featureNo++)
            {
                data[docNo].state[featureNo] = data[docNo].state[featureNo] / IMAGE_MAX_VALUE;
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
                    if(featureIndex==condLength) {
                        data[docIndex].action = class_map[int(tfIDF)];
                    }else
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
    //normalize_image(inputArray,totalNumInstances);
}



void load_kb(std::string kb_cf_file_name, std::string kb_filter_file_name) {

    load_filter_for_kb(kb_filter_file_name);
    load_cf_for_kb(kb_cf_file_name);
}

