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
//FilterMap kb_filter(1024, 0);
FilterMap2 kb_filter;
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



void load_kb(std::string kb_cf_file_name, std::string kb_filter_file_name) {

    //load_filter(kb_filter_file_name, kb_filter);
//    load_code_fragment(kb_cf_file_name, kb_cf);
}

