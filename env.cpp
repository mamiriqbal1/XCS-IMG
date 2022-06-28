#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include "configuration.h"
#include "env.h"
#include "cf_list.h"
#include "codeFragment.h"


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

void extract_features_edge(DataSource *data, int totalRows){
    FloatMatrix temp_img = FloatMatrix(image_height, FloatVector (image_width, NOT_INITIALIZED));
    BoundingBox full_img_bb;
    full_img_bb.size_y = image_height;
    full_img_bb.size_x = image_width;
    full_img_bb.y = 0;
    full_img_bb.x = 0;
    for(int img_id = 0; img_id<totalRows;img_id++){
        // make a copy in temp
        for(int y=0; y<image_height; y++) {
            for(int x=0; x<image_width; x++) {
                temp_img[y][x] = data[img_id].state[translate(full_img_bb, x, y)];
            }
        }
        // mark all pixels as 0 initially
        for(int y=0; y<full_img_bb.size_y - 1; y++){
            for(int x=0; x<full_img_bb.size_x - 1; x++){
                data[img_id].state[translate(full_img_bb, x, y)] = 0;
            }
        }
        // mark edge pixels with 1
        for(int y=0; y<full_img_bb.size_y - 1; y++){
            for(int x=0; x<full_img_bb.size_x - 1; x++){
                if(std::abs(temp_img[y][x] - temp_img[y][x+1]) > EDGE_THRESHOLD){
                    data[img_id].state[translate(full_img_bb, x, y)] = 1;
                }
            }
        }
        for(int x=0; x<full_img_bb.size_x - 1; x++){
            for(int y=0; y<full_img_bb.size_y - 1; y++){
                if(std::abs(temp_img[y][x] - temp_img[y+1][x]) > EDGE_THRESHOLD){
                    data[img_id].state[translate(full_img_bb, x, y)] = 1;
                }
            }
        }
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

    load_cf_for_kb(kb_cf_file_name);
}

