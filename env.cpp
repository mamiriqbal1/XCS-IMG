#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <cstring>
#include "xcsMacros.h"
#include "configuration.h"
#include "env.h"

//DataSource allData[1];
float lowerLimit[condLength];
float upperLimit[condLength];

void initializeLimit(){
    for (int i = 0; i<condLength; i++)
    {
        lowerLimit[i] = 0.0;
        upperLimit[i] = 0.0;
    }
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
    }
    infile.close();
    //updateRange(inputArray,totalNumInstances);
}


void loadData(DataSource inputArray[]){
//    int featureIndex = 0;
//    float tfIDF = 0.0;
//    int docIndex = 0;
//    int class1 = 0;
//    int class2 = 1;
//    int class3 = 2;
//    int class4 = 3;
//    int class5 = 4;
//    /*int class6 = 5;
//    int class7 = 6;
//    int class8 = 7;
//    int class9 = 8;
//    int class10 = 9;
//    /*int class11 = 10;
//    int class12 = 11;
//    int class13 = 12;
//    int class14 = 13;
//    int class15 = 14;
//    int class16 = 15;
//    int class17 = 16;
//    int class18 = 17;
//    int class19 = 18;
//    int class20 = 19;
//    int class21 = 20;
//    int class22 = 21;
//    int class23 = 22;*/
//
//    std::string line;
//    __glibcxx_assert(false); // this function is not to be used
//    std::ifstream infile("not needed");
//
//     if (infile.is_open())
//    {
//      while (!infile.eof())
//        {
//            if(docIndex<totalNumInstances)
//            {
//                featureIndex = 0;
//                getline(infile,line);
//                std::stringstream singleValue(line);
//                //while(singleValue>>featureIndex)
//                while(featureIndex <= condLength)
//                {
//                    singleValue>>tfIDF;
//                     if(featureIndex==condLength)
//                    //if(featureIndex==5000)
//                        //inputArray[docIndex].action = int(tfIDF);
//                        {
//                            if(docIndex < class1Instances)
//                                inputArray[docIndex].action = class1;
//                            else if((docIndex >= class1Instances) && (docIndex < (class1Instances+class2Instances)))
//                                inputArray[docIndex].action = class2;
//                            else if((docIndex >= class1Instances+class2Instances) && (docIndex < (class1Instances+class2Instances+class3Instances)))
//                                inputArray[docIndex].action = class3;
//                            else if((docIndex >= class1Instances+class2Instances+class3Instances) && (docIndex < (class1Instances+class2Instances+class3Instances+class4Instances)))
//                                inputArray[docIndex].action = class4;
//                            else if((docIndex >= class1Instances+class2Instances+class3Instances+class4Instances) && (docIndex < (class1Instances+class2Instances+class3Instances+class4Instances+class5Instances)))
//                                inputArray[docIndex].action = class5;
//                            /*else if((docIndex >= class1Instances+class2Instances+class3Instances+class4Instances+class5Instances) && (docIndex < (class1Instances+class2Instances+class3Instances+class4Instances+class5Instances+class6Instances)))
//                                inputArray[docIndex].action = class6;
//                            else if((docIndex >= class1Instances+class2Instances+class3Instances+class4Instances+class5Instances+class6Instances) && (docIndex < (class1Instances+class2Instances+class3Instances+class4Instances+class5Instances+class6Instances+class7Instances)))
//                                inputArray[docIndex].action = class7;
//                            else if((docIndex >= class1Instances+class2Instances+class3Instances+class4Instances+class5Instances+class6Instances+class7Instances) && (docIndex < (class1Instances+class2Instances+class3Instances+class4Instances+class5Instances+class6Instances+class7Instances+class8Instances)))
//                                inputArray[docIndex].action = class8;
//                            else if((docIndex >= class1Instances+class2Instances+class3Instances+class4Instances+class5Instances+class6Instances+class7Instances+class8Instances) && (docIndex < (class1Instances+class2Instances+class3Instances+class4Instances+class5Instances+class6Instances+class7Instances+class8Instances+class9Instances)))
//                                inputArray[docIndex].action = class9;
//                            else if((docIndex >= class1Instances+class2Instances+class3Instances+class4Instances+class5Instances+class6Instances+class7Instances+class8Instances+class9Instances) && (docIndex < (class1Instances+class2Instances+class3Instances+class4Instances+class5Instances+class6Instances+class7Instances+class8Instances+class9Instances+class10Instances)))
//                                inputArray[docIndex].action = class10;
//                            /*else if((docIndex >= class1Instances+class2Instances+class3Instances+class4Instances+class5Instances+class6Instances+class7Instances+class8Instances+class9Instances+class10Instances) && (docIndex < (class1Instances+class2Instances+class3Instances+class4Instances+class5Instances+class6Instances+class7Instances+class8Instances+class9Instances+class10Instances+class11Instances)))
//                                inputArray[docIndex].action = class11;
//                            else if((docIndex >= class1Instances+class2Instances+class3Instances+class4Instances+class5Instances+class6Instances+class7Instances+class8Instances+class9Instances+class10Instances+class11Instances) && (docIndex < (class1Instances+class2Instances+class3Instances+class4Instances+class5Instances+class6Instances+class7Instances+class8Instances+class9Instances+class10Instances+class11Instances+class12Instances)))
//                                inputArray[docIndex].action = class12;
//                            else if((docIndex >= class1Instances+class2Instances+class3Instances+class4Instances+class5Instances+class6Instances+class7Instances+class8Instances+class9Instances+class10Instances+class11Instances+class12Instances) && (docIndex < (class1Instances+class2Instances+class3Instances+class4Instances+class5Instances+class6Instances+class7Instances+class8Instances+class9Instances+class10Instances+class11Instances+class12Instances+class13Instances)))
//                                inputArray[docIndex].action = class13;
//                            else if((docIndex >= class1Instances+class2Instances+class3Instances+class4Instances+class5Instances+class6Instances+class7Instances+class8Instances+class9Instances+class10Instances+class11Instances+class12Instances+class13Instances) && (docIndex < (class1Instances+class2Instances+class3Instances+class4Instances+class5Instances+class6Instances+class7Instances+class8Instances+class9Instances+class10Instances+class11Instances+class12Instances+class13Instances+class14Instances)))
//                                inputArray[docIndex].action = class14;
//                            else if((docIndex >= class1Instances+class2Instances+class3Instances+class4Instances+class5Instances+class6Instances+class7Instances+class8Instances+class9Instances+class10Instances+class11Instances+class12Instances+class13Instances+class14Instances) && (docIndex < (class1Instances+class2Instances+class3Instances+class4Instances+class5Instances+class6Instances+class7Instances+class8Instances+class9Instances+class10Instances+class11Instances+class12Instances+class13Instances+class14Instances+class15Instances)))
//                                inputArray[docIndex].action = class15;
//                            else if((docIndex >= class1Instances+class2Instances+class3Instances+class4Instances+class5Instances+class6Instances+class7Instances+class8Instances+class9Instances+class10Instances+class11Instances+class12Instances+class13Instances+class14Instances+class15Instances) && (docIndex < (class1Instances+class2Instances+class3Instances+class4Instances+class5Instances+class6Instances+class7Instances+class8Instances+class9Instances+class10Instances+class11Instances+class12Instances+class13Instances+class14Instances+class15Instances+class16Instances)))
//                                inputArray[docIndex].action = class16;
//                            else if((docIndex >= class1Instances+class2Instances+class3Instances+class4Instances+class5Instances+class6Instances+class7Instances+class8Instances+class9Instances+class10Instances+class11Instances+class12Instances+class13Instances+class14Instances+class15Instances+class16Instances) && (docIndex < (class1Instances+class2Instances+class3Instances+class4Instances+class5Instances+class6Instances+class7Instances+class8Instances+class9Instances+class10Instances+class11Instances+class12Instances+class13Instances+class14Instances+class15Instances+class16Instances+class17Instances)))
//                                inputArray[docIndex].action = class17;
//                            else if((docIndex >= class1Instances+class2Instances+class3Instances+class4Instances+class5Instances+class6Instances+class7Instances+class8Instances+class9Instances+class10Instances+class11Instances+class12Instances+class13Instances+class14Instances+class15Instances+class16Instances+class17Instances) && (docIndex < (class1Instances+class2Instances+class3Instances+class4Instances+class5Instances+class6Instances+class7Instances+class8Instances+class9Instances+class10Instances+class11Instances+class12Instances+class13Instances+class14Instances+class15Instances+class16Instances+class17Instances+class18Instances)))
//                                inputArray[docIndex].action = class18;
//                            else if((docIndex >= class1Instances+class2Instances+class3Instances+class4Instances+class5Instances+class6Instances+class7Instances+class8Instances+class9Instances+class10Instances+class11Instances+class12Instances+class13Instances+class14Instances+class15Instances+class16Instances+class17Instances+class18Instances) && (docIndex < (class1Instances+class2Instances+class3Instances+class4Instances+class5Instances+class6Instances+class7Instances+class8Instances+class9Instances+class10Instances+class11Instances+class12Instances+class13Instances+class14Instances+class15Instances+class16Instances+class17Instances+class18Instances+class19Instances)))
//                                inputArray[docIndex].action = class19;
//                            else if((docIndex >= class1Instances+class2Instances+class3Instances+class4Instances+class5Instances+class6Instances+class7Instances+class8Instances+class9Instances+class10Instances+class11Instances+class12Instances+class13Instances+class14Instances+class15Instances+class16Instances+class17Instances+class18Instances+class19Instances) && (docIndex < (class1Instances+class2Instances+class3Instances+class4Instances+class5Instances+class6Instances+class7Instances+class8Instances+class9Instances+class10Instances+class11Instances+class12Instances+class13Instances+class14Instances+class15Instances+class16Instances+class17Instances+class18Instances+class19Instances+class20Instances)))
//                                inputArray[docIndex].action = class20;
//                            else if((docIndex >= class1Instances+class2Instances+class3Instances+class4Instances+class5Instances+class6Instances+class7Instances+class8Instances+class9Instances+class10Instances+class11Instances+class12Instances+class13Instances+class14Instances+class15Instances+class16Instances+class17Instances+class18Instances+class19Instances+class20Instances) && (docIndex < (class1Instances+class2Instances+class3Instances+class4Instances+class5Instances+class6Instances+class7Instances+class8Instances+class9Instances+class10Instances+class11Instances+class12Instances+class13Instances+class14Instances+class15Instances+class16Instances+class17Instances+class18Instances+class19Instances+class20Instances+class21Instances)))
//                                inputArray[docIndex].action = class21;
//                            else if((docIndex >= class1Instances+class2Instances+class3Instances+class4Instances+class5Instances+class6Instances+class7Instances+class8Instances+class9Instances+class10Instances+class11Instances+class12Instances+class13Instances+class14Instances+class15Instances+class16Instances+class17Instances+class18Instances+class19Instances+class20Instances+class21Instances) && (docIndex < (class1Instances+class2Instances+class3Instances+class4Instances+class5Instances+class6Instances+class7Instances+class8Instances+class9Instances+class10Instances+class11Instances+class12Instances+class13Instances+class14Instances+class15Instances+class16Instances+class17Instances+class18Instances+class19Instances+class20Instances+class21Instances+class22Instances)))
//                                inputArray[docIndex].action = class22;
//                            else if((docIndex >= class1Instances+class2Instances+class3Instances+class4Instances+class5Instances+class6Instances+class7Instances+class8Instances+class9Instances+class10Instances+class11Instances+class12Instances+class13Instances+class14Instances+class15Instances+class16Instances+class17Instances+class18Instances+class19Instances+class20Instances+class21Instances+class22Instances) && (docIndex < (class1Instances+class2Instances+class3Instances+class4Instances+class5Instances+class6Instances+class7Instances+class8Instances+class9Instances+class10Instances+class11Instances+class12Instances+class13Instances+class14Instances+class15Instances+class16Instances+class17Instances+class18Instances+class19Instances+class20Instances+class21Instances+class22Instances+class23Instances)))
//                                inputArray[docIndex].action = class23;*/
//                        }
//                    else
//                    {
//                        //tfIDF = roundRealValue(tfIDF,precisionDigits);
//                        inputArray[docIndex].state[featureIndex] = tfIDF;
//                        //std::cout<<featureIndex<<":"<<inputArray[docIndex].state[featureIndex]<<" ";
//                        if(tfIDF>upperLimit[featureIndex])
//                            upperLimit[featureIndex] = tfIDF;
//                        if(tfIDF<lowerLimit[featureIndex])
//                            lowerLimit[featureIndex] = tfIDF;
//                     }
//                     featureIndex++;
//                }
//                docIndex++;
//                std::cout<<"instance: "<<docIndex<<std::endl;
//            }
//        }
//    }
//    infile.close();
//    //updateRange(inputArray,totalNumInstances);
}
/*
void loadDataFile(){

    initializeLimit();
    //  initializeInput(allData,totalNumInstances);
    loadData(allData);

    if(Testing)
    {/*
        splitData();
        writeData(trainData,trainNumInstances,trainingFile);
        writeData(testData,testNumInstances,testingFile);
        updateRange(trainData,trainNumInstances);
        updateRange(testData,testNumInstances); /
    }
    else
    {
            updateRange(allData,totalNumInstances);
    }

}
*/
void splitData(DataSource allData[],DataSource trainData[], DataSource testData[]){ //Split totaldata into training data and test data
//	int indexsUsed[totalNumInstances];
//	//int pos = 0, neg = 0;
//	int class1train = 0, class2train = 0, class3train = 0, class4train = 0, class5train = 0;
//	int class6train = 0, class7train = 0, class8train = 0, class9train = 0, class10train = 0;
//	int class11train = 0, class12train = 0, class13train = 0, class14train = 0, class15train = 0;
//	int class16train = 0, class17train = 0, class18train = 0, class19train = 0, class20train = 0, class21train = 0, class22train = 0, class23train = 0;;
//	for(int i=0; i<totalNumInstances; i++){//initialize to -1 to avoid 0 problem
//		indexsUsed[i]=-1;
//	}
//	for(int i=0; i<trainNumInstances; i++)
//	{//make training dataset
//        std::cout<<"in split data : i:"<<i<<std::endl; // added by me
//		int index = irand(totalNumInstances);
//		if(isIntegerInArray(index, indexsUsed, i))
//		{
//			i--;
//			continue;
//		}
//		else
//        {
//            if(allData[index].action == 0)
//            {
//                if(class1train < class1TrainNumInstance)
//                {
//                    indexsUsed[i]=index;
//                    trainData[i] = allData[index];
//                    class1train++;
//                }
//                else
//                {
//                    i--;
//                    continue;
//                }
//            }
//            else if(allData[index].action == 1)
//            {
//                if(class2train < class2TrainNumInstance)
//                {
//                    indexsUsed[i]=index;
//                    trainData[i] = allData[index];
//                    class2train++;
//                }
//                else
//                {
//                    i--;
//                    continue;
//                }
//            }
//            else if(allData[index].action == 2)
//            {
//                if(class3train < class3TrainNumInstance)
//                {
//                    indexsUsed[i]=index;
//                    trainData[i] = allData[index];
//                    class3train++;
//                }
//                else
//                {
//                    i--;
//                    continue;
//                }
//            }
//            else if(allData[index].action == 3)
//            {
//                if(class4train < class4TrainNumInstance)
//                {
//                    indexsUsed[i]=index;
//                    trainData[i] = allData[index];
//                    class4train++;
//                }
//                else
//                {
//                    i--;
//                    continue;
//                }
//            }
//            else if(allData[index].action == 4)
//            {
//                if(class5train < class5TrainNumInstance)
//                {
//                    indexsUsed[i]=index;
//                    trainData[i] = allData[index];
//                    class5train++;
//                }
//                else
//                {
//                    i--;
//                    continue;
//                }
//            }
//            /*else if(allData[index].action == 5)
//            {
//                if(class6train < class6TrainNumInstance)
//                {
//                    indexsUsed[i]=index;
//                    trainData[i] = allData[index];
//                    class6train++;
//                }
//                else
//                {
//                    i--;
//                    continue;
//                }
//            }
//            else if(allData[index].action == 6)
//            {
//                if(class7train < class7TrainNumInstance)
//                {
//                    indexsUsed[i]=index;
//                    trainData[i] = allData[index];
//                    class7train++;
//                }
//                else
//                {
//                    i--;
//                    continue;
//                }
//            }
//            else if(allData[index].action == 7)
//            {
//                if(class8train < class8TrainNumInstance)
//                {
//                    indexsUsed[i]=index;
//                    trainData[i] = allData[index];
//                    class8train++;
//                }
//                else
//                {
//                    i--;
//                    continue;
//                }
//            }
//            else if(allData[index].action == 8)
//            {
//                if(class9train < class9TrainNumInstance)
//                {
//                    indexsUsed[i]=index;
//                    trainData[i] = allData[index];
//                    class9train++;
//                }
//                else
//                {
//                    i--;
//                    continue;
//                }
//            }
//            else if(allData[index].action == 9)
//            {
//                if(class10train < class10TrainNumInstance)
//                {
//                    indexsUsed[i]=index;
//                    trainData[i] = allData[index];
//                    class10train++;
//                }
//                else
//                {
//                    i--;
//                    continue;
//                }
//            }
//            /*else if(allData[index].action == 10)
//            {
//                if(class11train < class11TrainNumInstance)
//                {
//                    indexsUsed[i]=index;
//                    trainData[i] = allData[index];
//                    class11train++;
//                }
//                else
//                {
//                    i--;
//                    continue;
//                }
//            }
//            else if(allData[index].action == 11)
//            {
//                if(class12train < class12TrainNumInstance)
//                {
//                    indexsUsed[i]=index;
//                    trainData[i] = allData[index];
//                    class12train++;
//                }
//                else
//                {
//                    i--;
//                    continue;
//                }
//            }
//            else if(allData[index].action == 12)
//            {
//                if(class13train < class13TrainNumInstance)
//                {
//                    indexsUsed[i]=index;
//                    trainData[i] = allData[index];
//                    class13train++;
//                }
//                else
//                {
//                    i--;
//                    continue;
//                }
//            }
//            else if(allData[index].action == 13)
//            {
//                if(class14train < class14TrainNumInstance)
//                {
//                    indexsUsed[i]=index;
//                    trainData[i] = allData[index];
//                    class14train++;
//                }
//                else
//                {
//                    i--;
//                    continue;
//                }
//            }
//            else if(allData[index].action == 14)
//            {
//                if(class15train < class15TrainNumInstance)
//                {
//                    indexsUsed[i]=index;
//                    trainData[i] = allData[index];
//                    class15train++;
//                }
//                else
//                {
//                    i--;
//                    continue;
//                }
//            }
//            else if(allData[index].action == 15)
//            {
//                if(class16train < class16TrainNumInstance)
//                {
//                    indexsUsed[i]=index;
//                    trainData[i] = allData[index];
//                    class16train++;
//                }
//                else
//                {
//                    i--;
//                    continue;
//                }
//            }
//            else if(allData[index].action == 16)
//            {
//                if(class17train < class17TrainNumInstance)
//                {
//                    indexsUsed[i]=index;
//                    trainData[i] = allData[index];
//                    class17train++;
//                }
//                else
//                {
//                    i--;
//                    continue;
//                }
//            }
//            else if(allData[index].action == 17)
//            {
//                if(class18train < class18TrainNumInstance)
//                {
//                    indexsUsed[i]=index;
//                    trainData[i] = allData[index];
//                    class18train++;
//                }
//                else
//                {
//                    i--;
//                    continue;
//                }
//            }
//            else if(allData[index].action == 18)
//            {
//                if(class19train < class19TrainNumInstance)
//                {
//                    indexsUsed[i]=index;
//                    trainData[i] = allData[index];
//                    class19train++;
//                }
//                else
//                {
//                    i--;
//                    continue;
//                }
//            }
//            else if(allData[index].action == 19)
//            {
//                if(class20train < class20TrainNumInstance)
//                {
//                    indexsUsed[i]=index;
//                    trainData[i] = allData[index];
//                    class20train++;
//                }
//                else
//                {
//                    i--;
//                    continue;
//                }
//            }
//            else if(allData[index].action == 20)
//            {
//                if(class21train < class21TrainNumInstance)
//                {
//                    indexsUsed[i]=index;
//                    trainData[i] = allData[index];
//                    class21train++;
//                }
//                else
//                {
//                    i--;
//                    continue;
//                }
//            }
//            else if(allData[index].action == 21)
//            {
//                if(class22train < class22TrainNumInstance)
//                {
//                    indexsUsed[i]=index;
//                    trainData[i] = allData[index];
//                    class22train++;
//                }
//                else
//                {
//                    i--;
//                    continue;
//                }
//            }
//            else if(allData[index].action == 22)
//            {
//                if(class23train < class23TrainNumInstance)
//                {
//                    indexsUsed[i]=index;
//                    trainData[i] = allData[index];
//                    class23train++;
//                }
//                else
//                {
//                    i--;
//                    continue;
//                }
//            }*/
//        }
//	}
//	std::cout<<"In Split data after first loop:";
//	//make test dataset from allData-trainData
//	int num=0;
//	for(int i=0; i<totalNumInstances; i++){
//		if(!isIntegerInArray(i, indexsUsed, totalNumInstances)){
//			testData[i-num] = allData[i];
//		}
//		else{
//			num++;
//		}
//	}
//	//std::cout<<"pos:"<<pos<<"---neg:"<<neg<<std::endl;
//	std::cout<<"\n" <<"class1train:"<<class1train<<"---class2train:"<<class2train<<"---class3train:"<<class3train<<"---class4train:"<<class4train<<"---class5train:"<<class5train<<"\n";
//	std::cout<<"class6train:"<<class6train<<"---class7train:"<<class7train<<"---class8train:"<<class8train<<"---class9train:"<<class9train<<"---class10train:"<<class10train<<"\n";
//	//std::cout<<"class11train:"<<class11train<<"---class12train:"<<class12train<<"---class13train:"<<class13train<<"---class14train:"<<class14train<<"---class15train:"<<class15train<<"\n";
//	//std::cout<<"class16train:"<<class16train<<"---class17train:"<<class17train<<"---class18train:"<<class18train<<"---class19train:"<<class19train<<"---class20train:"<<class20train<<"\n";
//	//std::cout<<"class21train:"<<class21train<<"---class22train:"<<class22train<<"---class23train:"<<class23train<<"\n";
}

bool isIntegerInArray(int integer, int array[totalNumInstances], int highestPosition){
	for(int i=0; i<highestPosition; i++){
		if(array[i]==integer)
			return true;
	}
	return false;
}


/*
DataSource resetState(DataSource inputArray[]){ // generates a new random problem instance.
    int index=0;
    DataSource state;
    if(Testing)
    {
        index = irand(trainNumInstances);
        state = trainData[index];/
    }
    else
    {
        index = irand(totalNumInstances);
        std::cout<<index<<" ";
        std::memcpy(state.state, inputArray[index].state, condLength* sizeof(float));
        state.action = inputArray[index].action;
         //std::memcpy(state, inputArray[index], sizeof state);
        //state = inputArray[index];

    }
    return state;
}
*/
