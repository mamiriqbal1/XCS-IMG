#include <stdlib.h>
#include <sys/types.h>
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <assert.h>
#include <math.h>
#include <cstring>
#include <vector>
#include "xcsMacros.h"
#include "time.h"
#include "configuration.h"
#include "classifier.h"
#include "env.h"
#include "xcs.h"
//using namespace std;

struct distanceInputClassifier{
  int posClassifier;
  float distance;
};

double executeAction(int action, int stateAction, bool &wasCorrect){  // Executes the action and determines the reward.

    int ret=0;
   // std::cout<<"--- "<<action<<" = " <<stateAction<<"--- ";

    if(action == stateAction)
    {
       // std::cout<<"-- "<<action << " = "<<stateAction<<"--";
        wasCorrect=true;
        ret = maxPayoff;
    }
    else
    {
        wasCorrect=false;
        ret = 0;
    }
    return (double)ret;
}

void header(){
    time_t t;
    char hostName[] = "act-PC"; //[MAXHOSTNAMELEN];
    time(&t);
    //gethostname(hostName,sizeof(hostName));
    printf("%s\n Seed: %ld %s",hostName,getSeed(),ctime(&t));
}//end header

void Exit(FILE *fp){
    printf("\nShutting down...\n");
    elapsedTime();
    header();
    printf("Elapsed Time: ");
    elapsed(fp);
    printf(" Secs\n");
}//end exit

/*##########---- Output ----##########*/
/**
 * Writes the performance of the XCS to the specified file.
 * The function writes time performance systemError actualPopulationSizeInMacroClassifiers.
 * Performance and system error are averaged over the last fifty trials.
 *
 * @param performance The performance in the last fifty exploration trials.
 * @param sysError The system error in the last fifty exploration trials.
 * @param exploreProbC The number of exploration trials executed so far.
 */
void writePerformance(ClassifierSet *population, int performance[], double sysError[], int exploreProbC){
    char buf[1000];
    double perf=0.0, serr=0.0;
    int setSize;
    for(int i=0; i<testFrequency; i++)
    {
        perf+=performance[i];
        serr+=sysError[i];
    }
    perf/=testFrequency;
    serr/=testFrequency;
    setSize = getSetSize(population);

    snprintf(buf,strlen(buf),"%d ",exploreProbC);

    sprintf(buf,"%d ",exploreProbC);
    fwrite(buf,strlen(buf),1,filePerformance);
    sprintf(buf,"%f ",perf);
    fwrite(buf,strlen(buf),1,filePerformance);
    sprintf(buf,"%f ",serr);
    fwrite(buf,strlen(buf),1,filePerformance);
    sprintf(buf,"%d ",setSize);
    fwrite(buf,strlen(buf),1,filePerformance);
    fwrite("\n",strlen("\n"),1,filePerformance);

      //int numerositySum = getNumerositySum(population); printf("%d %f %f %d %d\n",exploreProbC,perf,serr,setSize,numerositySum);
    //printf("%d %f %f %d\n",exploreProbC,perf,serr,setSize);
}

/*******************************Write Test Performance*************************************/
void writeTestPerformance(ClassifierSet *population, int performance[], double sysError[], int exploreProbC){
    char buf[100];
    double perf=0.0, serr=0.0;

    for(int i=0; i<testFrequency; i++)
    {
        perf+=performance[i];
        serr+=sysError[i];
    }

    perf/=testFrequency;
    serr/=testFrequency;
   // std::cout<<"testing";

    int setSize = getSetSize(population);

    sprintf(buf,"%d ",exploreProbC);
    fwrite(buf,strlen(buf),1,testPerformance);
    sprintf(buf,"%f ",perf);
    fwrite(buf,strlen(buf),1,testPerformance);
    sprintf(buf,"%f ",serr);
    fwrite(buf,strlen(buf),1,testPerformance);
    sprintf(buf,"%d ",setSize);
    fwrite(buf,strlen(buf),1,testPerformance);
    fwrite("\n",strlen("\n"),1,testPerformance);

    //int numerositySum = getNumerositySum(population); printf("%d %f %f %d %d\n",exploreProbC,perf,serr,setSize,numerositySum);
    //printf("%d %f %f %d\n",exploreProbC,perf,serr,setSize);
}


/**************************** Single Step Experiments ***************************/
void doOneSingleStepProblemExplore(ClassifierSet **population, DataSource *object, int counter){ // Executes one main learning loop for a single step problem.

    bool wasCorrect = false;
    ClassifierSet *mset, *aset, *killset=NULL;

    mset = getMatchSet(population,&killset,object->state,counter);
    freeSet(&killset);
    //cout<<"test2";
    //getchar();
    getPredictionArray(mset);
    int actionWinner = randomActionWinner();
    aset = getActionSet(actionWinner,mset);
    double reward = executeAction(actionWinner,object->action,wasCorrect);

    updateActionSet(&aset,0.0,reward,population,&killset);
    freeSet(&killset);

    discoveryComponent(&aset,population,&killset,counter,object->state);
    freeSet(&killset);

    freeSet(&mset);
    freeSet(&aset);
}

void doOneSingleStepProblemExploit(ClassifierSet **population, DataSource *object, int counter, int correct[], double sysError[]){  //Executes one main performance evaluation loop for a single step problem.

    bool wasCorrect = false;
    ClassifierSet *mset, *killset=NULL;

    mset = getMatchSet(population,&killset,object->state,counter);
    freeSet(&killset);

    getPredictionArray(mset);
    int actionWinner = bestActionWinner();
    double reward = executeAction(actionWinner,object->action,wasCorrect);

    //std::cout<<"-- "<<wasCorrect<<wasCorrect<<wasCorrect<<wasCorrect<<wasCorrect<<" --- ";
    if(wasCorrect)
    {
        correct[counter%testFrequency]=1;
    }
    else
    {
        correct[counter%testFrequency]=0;
    }
    sysError[counter%testFrequency] = absoluteValue(reward - getBestValue());

    freeSet(&mset);
}


 void doOneSingleStepExperiment(ClassifierSet **population){  //Executes one single-step experiment monitoring the performance.

    int explore=0;
    int correct[testFrequency];
    double sysError[testFrequency];

    DataSource *state = NULL;
    //int counter = 0;
    //int index = 0;

    for(int exploreProbC=0; exploreProbC <= maxProblems; exploreProbC+=explore)
    {
        explore = (explore+1)%2;
        // state = inputArray[irand(totalNumInstances)];
        //index = ;
        //state = &inputArray[irand(totalNumInstances)];
        if(Testing)
        {
            state = &trainingData[irand(trainNumInstances)];
        }
        else
        {
            state = &input[irand(totalNumInstances)];
        }

        //std::cout<<index<<" ";
        /*for (int i=0;i<condLength;i++)
        {
            if(state->state[i] != 0.0)
                std::cout<<state->state[i]<<" ";
        }
        std::cout<<std::endl;
        getchar();*/

        if(explore==1)
        {
            doOneSingleStepProblemExplore(population,state,exploreProbC);
        }
        else
        {
            doOneSingleStepProblemExploit(population,state,exploreProbC, correct, sysError);
        }
        if(exploreProbC%testFrequency==0 && explore==0 && exploreProbC>0)
        {
            writePerformance(*population,correct,sysError,exploreProbC);
        }
        if((exploreProbC%10000) == 0)
            std::cout<<"ExploreprobC: "<<exploreProbC<<std::endl;
    }
//    delete state;
}

void sortAll(distanceInputClassifier arrayToBeSorted[], int size){ // Bubble Sort Function for Ascending Order
  int i, j;
  bool flag = true;    // set flag to true to start first pass
  distanceInputClassifier temp;          // holding variable
  for(i = 1; i<=size && flag==true; i++){
    flag = false;
    for(j=0; j < size-1; j++){
      if(arrayToBeSorted[j+1].distance < arrayToBeSorted[j].distance){
	temp = arrayToBeSorted[j];  // swap elements
        arrayToBeSorted[j] = arrayToBeSorted[j+1];
        arrayToBeSorted[j+1] = temp;
        flag = true;  // indicates that a swap occurred.
      }
    }
  }
  return;   //arrays are passed to functions by address; nothing is returned
}

void sortK(distanceInputClassifier arrayToBeSorted[], int size){ // Bubble Sort Function for Ascending Order
  int i, j;
  bool flag = true;    // set flag to true to start first pass
  distanceInputClassifier temp;          // holding variable
  for(i = 1; i<=size && flag==true; i++){
    flag = false;
    for(j=0; j < size-1; j++){
      if(arrayToBeSorted[j+1].posClassifier < arrayToBeSorted[j].posClassifier){
	temp = arrayToBeSorted[j];  // swap elements
        arrayToBeSorted[j] = arrayToBeSorted[j+1];
        arrayToBeSorted[j+1] = temp;
        flag = true;  // indicates that a swap occurred.
      }
    }
  }
  return;   //arrays are passed to functions by address; nothing is returned
}

std::string NumberToString(int num){
    std::stringstream ss;
    ss << num;
    return ss.str();
}

void doOneSingleStepTest(ClassifierSet *population){
	bool wasCorrect = false;
	int correctCounter = 0;
	int correct[testNumInstances];
	double sysError[testFrequency];
	DataSource *testState = NULL;
	int tmpcorrectcounter = 0;
    int tmpnotmatched = 0;
	int TP = 0, TN = 0, FP = 0, FN = 0, TNeu = 0, FNeu = 0;
    int class1TP=0, class2TP=0, class3TP=0, class4TP=0, class5TP=0, class6TP=0, class7TP=0, class8TP=0, class9TP=0, class10TP=0;
    int class11TP=0, class12TP=0, class13TP=0, class14TP=0, class15TP=0, class16TP=0, class17TP=0, class18TP=0, class19TP=0, class20TP=0, class21TP=0, class22TP=0, class23TP=0;
	int class1FP=0, class2FP=0, class3FP=0, class4FP=0, class5FP=0, class6FP=0, class7FP=0, class8FP=0, class9FP=0, class10FP=0;
	int class11FP=0, class12FP=0, class13FP=0, class14FP=0, class15FP=0, class16FP=0, class17FP=0, class18FP=0, class19FP=0, class20FP=0, class21FP=0, class22FP=0, class23FP=0;

    //loadDataFile(false);
    int cc = 0;

	for(int t=0; t<testNumInstances; t++){
		ClassifierSet *mset=NULL, *poppointer;
		bool isMatched = false;
		//resetStateTesting(testState,t);
		testState = &testingData[t];
    /*
		for(int i =0;i<condLength;i++)
            if(testState->state[i]!=0.0)
                std::cout<<i<<":"<<testState->state[i]<<" ";
            std::cout<<"=="<<testState->action;
            std::cout<<std::endl;
            getchar();*/

        for(poppointer= population; poppointer!=NULL; poppointer=poppointer->next)
        {
            if(isConditionMatched(poppointer->classifier->condition,testState->state))
            {
                isMatched = true;tmpcorrectcounter++;
                addNewClassifierToSet(poppointer->classifier, &mset); // add matching classifier to the matchset
            }
        }
        if(isMatched == false){
            cc++;
            tmpnotmatched++;
            int popSize = getSetSize(population);
			distanceInputClassifier distanceArray[popSize];
			int i = 0;
			int k = popSize*tournamentSize;
			assert(k > 0); // make sure that k > 0

			for(poppointer= population; poppointer!=NULL; poppointer=poppointer->next){
			  distanceArray[i].posClassifier = i;
			  distanceArray[i].distance = computeDistance(poppointer->classifier->condition,testState->state);
			  i++;
			}
			sortAll(distanceArray,popSize);
			sortK(distanceArray,k);
			int mK = 0, mT = 0;
			for(poppointer=population; poppointer!=NULL&&mK<k; poppointer=poppointer->next,mT++){
			  if(distanceArray[mK].posClassifier == mT){
			    addNewClassifierToSet(poppointer->classifier, &mset); // add classifier to the matchset
			    mK++;
			  }
		   	}
		}
		getPredictionArray(mset);
        int actionWinner = bestActionWinner();
        double reward = executeAction(actionWinner,testState->action,wasCorrect);
        sysError[t%testFrequency] = absoluteValue(reward - getBestValue());

		//float error = computeError(computedSaliencyMap,state.output);
		if(wasCorrect)
        {
            /////  For TWO Classes ///////////////////////
            /*if(actionWinner==1)
            {
                TP++;
            }else if(actionWinner==0)
            {
                TN++;
            }else
            {
                TNeu++;
            }*/
            ///// For Multiple Classes/////// ////////////////
                if(actionWinner==0)
                {
                    class1TP++;
                    //TN++;
                }
                else if(actionWinner==1)
                {
                    class2TP++;
                //TP++;
                }
                else if(actionWinner==2)
                {
                    class3TP++;
                //TP++;
                }
                else if(actionWinner==3)
                {
                    class4TP++;
                //TP++;
                }
                else if(actionWinner==4)
                {
                    class5TP++;
                //TP++;
                }
                /*else if(actionWinner==5)
                {
                    class6TP++;
                    //TN++;
                }
                else if(actionWinner==6)
                {
                    class7TP++;
                //TP++;
                }
                else if(actionWinner==7)
                {
                    class8TP++;
                //TP++;
                }
                else if(actionWinner==8)
                {
                    class9TP++;
                //TP++;
                }
                else if(actionWinner==9)
                {
                    class10TP++;
                //TP++;
                }
                /*else if(actionWinner==10)
                {
                    class11TP++;
                    //TN++;
                }
                else if(actionWinner==11)
                {
                    class12TP++;
                //TP++;
                }
                else if(actionWinner==12)
                {
                    class13TP++;
                //TP++;
                }
                else if(actionWinner==13)
                {
                    class14TP++;
                //TP++;
                }
                else if(actionWinner==14)
                {
                    class15TP++;
                //TP++;
                }
                else if(actionWinner==15)
                {
                    class16TP++;
                    //TN++;
                }
                else if(actionWinner==16)
                {
                    class17TP++;
                //TP++;
                }
                else if(actionWinner==17)
                {
                    class18TP++;
                //TP++;
                }
                else if(actionWinner==18)
                {
                    class19TP++;
                //TP++;
                }
                else if(actionWinner==19)
                {
                    class20TP++;
                //TP++;
                }
                 else if(actionWinner==20)
                {
                    class21TP++;
                //TP++;
                }
                 else if(actionWinner==21)
                {
                    class22TP++;
                //TP++;
                }
                 else if(actionWinner==22)
                {
                    class23TP++;
                //TP++;
                }*/
                else
                {
                TNeu++;
                }

                correct[t%testFrequency]=1;
                correctCounter++;
        }
        else
        {
            //// For TWO CLASSES ///////////////////////
            /*if(actionWinner==1)
            {
                FP++;
            }else if(actionWinner==0)
            {
                FN++;
            }
            else
            {
                FNeu++;
            }*/
            //// For Multiple Classes ///////////////////////
                if(actionWinner==0)
                {
                    class1FP++;
                    //TN++;
                }
                else if(actionWinner==1)
                {
                    class2FP++;
                //TP++;
                }
                else if(actionWinner==2)
                {
                    class3FP++;
                //TP++;
                }
                else if(actionWinner==3)
                {
                    class4FP++;
                //TP++;
                }
                else if(actionWinner==4)
                {
                    class5FP++;
                //TP++;
                }
                /*else if(actionWinner==5)
                {
                    class6FP++;
                    //TN++;
                }
                else if(actionWinner==6)
                {
                    class7FP++;
                //TP++;
                }
                else if(actionWinner==7)
                {
                    class8FP++;
                //TP++;
                }
                else if(actionWinner==8)
                {
                    class9FP++;
                //TP++;
                }
                else if(actionWinner==9)
                {
                    class10FP++;
                //TP++;
                }
                /*else if(actionWinner==10)
                {
                    class11FP++;
                    //TN++;
                }
                else if(actionWinner==11)
                {
                    class12FP++;
                //TP++;
                }
                else if(actionWinner==12)
                {
                    class13FP++;
                //TP++;
                }
                else if(actionWinner==13)
                {
                    class14FP++;
                //TP++;
                }
                else if(actionWinner==14)
                {
                    class15FP++;
                //TP++;
                }
                else if(actionWinner==15)
                {
                    class16FP++;
                    //TN++;
                }
                else if(actionWinner==16)
                {
                    class17FP++;
                //TP++;
                }
                else if(actionWinner==17)
                {
                    class18FP++;
                //TP++;
                }
                else if(actionWinner==18)
                {
                    class19FP++;
                //TP++;
                }
                else if(actionWinner==19)
                {
                    class20FP++;
                //TP++;
                }
                else if(actionWinner==20)
                {
                    class21FP++;
                //TP++;
                }
                else if(actionWinner==21)
                {
                    class22FP++;
                //TP++;
                }
                else if(actionWinner==22)
                {
                    class23FP++;
                //TP++;
                }*/

                correct[t%testFrequency]=0;
        }
        if((t+1) % testFrequency == 0 && t>0)
        {
//            writeTestPerformance(population,correct, sysError, t+1);
        }
		freeSet(&mset);
    }
	std::ofstream resFile;
    resFile.open(resultFile,std::ios::app);
    std::string wLine = "";
    std::string wLine2 = "";

	std::cout<<"Number of KNNs = "<<cc<<std::endl;
	std::cout<<"Number of correct Instances = "<<correctCounter<<std::endl;
	std::cout<<"TP: "<<TP<<"--TN: "<<TN<<"--FP: "<<FP<<"--FN: "<<FN;
    // For two classes
	//wLine = NumberToString(correctCounter) + " " + NumberToString(tmpnotmatched) + " " + NumberToString(TP) + " " + NumberToString(TN)  + " " + NumberToString(FP)  + " " + NumberToString(FN) + "\n";
	//resFile<<wLine;
	// for multiple classes
	wLine = NumberToString(correctCounter) + " " + NumberToString(cc) + " " + NumberToString(class1TP) + " " + NumberToString(class1FP)  + " " + NumberToString(class2TP) + " " + NumberToString(class2FP) +" " + NumberToString(class3TP) + " " + NumberToString(class3FP)+ " " +NumberToString(class4TP) + " " + NumberToString(class4FP)+ " " + NumberToString(class5TP) + " " + NumberToString(class5FP);
	//" " + NumberToString(class6TP) + " " + NumberToString(class6FP)  + " " + NumberToString(class7TP) + " " + NumberToString(class7FP) +" " + NumberToString(class8TP) + " " + NumberToString(class8FP)+ " " +NumberToString(class9TP) + " " + NumberToString(class9FP)+ " " + NumberToString(class10TP) + " " + NumberToString(class10FP);
	//+ " " + NumberToString(class11TP) + " " + NumberToString(class11FP)  + " " + NumberToString(class12TP) + " " + NumberToString(class12FP) +" " + NumberToString(class13TP) + " " + NumberToString(class13FP)+ " " +NumberToString(class14TP) + " " + NumberToString(class14FP)+ " " + NumberToString(class15TP) + " " + NumberToString(class15FP)+
	//+ " " + NumberToString(class16TP) + " " + NumberToString(class16FP)  + " " + NumberToString(class17TP) + " " + NumberToString(class17FP) +" " + NumberToString(class18TP) + " " + NumberToString(class18FP)+ " " +NumberToString(class19TP) + " " + NumberToString(class19FP)+ " " + NumberToString(class20TP) + " " + NumberToString(class20FP)+ " " + NumberToString(class21TP) + " " + NumberToString(class21FP) + " " + NumberToString(class22TP) + " " + NumberToString(class22FP) +" " + NumberToString(class23TP) + " " + NumberToString(class23FP) + "\n";
	resFile<<wLine;
	resFile.close();
}

void loadDataFile(DataSource inputArray[]){

    loadData(inputArray);
    if (Testing)
    {
		std::cout<<"split Data"<<std::endl;
        trainingData = new DataSource[trainNumInstances];
        testingData = new DataSource[testNumInstances];
        splitData(inputArray,trainingData,testingData);
        updateRange(trainingData,trainNumInstances);
        updateRange(testingData,testNumInstances);
        //writeData(trainingData,trainNumInstances,trainingFile);
        //writeData(testingData,testNumInstances,testingFile);

    }
    else
    {
        updateRange(inputArray,totalNumInstances);
    }
}


//void startXCS(int inputFile[][condLength+1])
void startXCS(){
    startTimer();
    ClassifierSet *population;
    initializePopulation(&population,cfReadingFilePointer);//,cfWritingFilePointer);
    printf("\nLoading Input! Please wait ....\n");
    //inp = new ds[totalNumInstances];
    //loadData(inp);

    input = new DataSource[totalNumInstances];
    initializeInput(input,totalNumInstances);
    //loadData(input);
    //loadDataFile();
    loadDataFile(input);

    printf("\nIt is in progress! Please wait ....\n");

    doOneSingleStepExperiment(&population);
    //simplifyPopulation(&population);

    if(Testing){
            printf("\nTesting\n");
            //loadDataFile(false); //load test data
            doOneSingleStepTest(population); // testing
    }

    population = sortClassifierSet(&population,2); // sort according to 'fitness'
    fprintClassifierSet(fileClassifierPopulation,cfWritingFilePointer,population);
    freeClassifierSet(&population); // free population for this experiment
    delete []input;

}//end startXCS

int main(int argc, char **argv){

    //loadFileString();

    char path[10] = "";
    char curLevel[2] = "";
    sprintf(curLevel, "%d", currentProblemLevel);
    char prevLevel[2] = "";
    sprintf(prevLevel, "%d", currentProblemLevel-1);

    for(int j=0; j<run; j++)
    {
        sprintf(path, "%d", j);
        //mkdir(path, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        mkdir(path);
        setSeed(seeds[j]);

        char outputFile[200];
        sprintf(outputFile,"%s%s",path,outputFileName);
        filePerformance=fopen(outputFile, "w"); //open outputFile in writing mode
        if (filePerformance == NULL)
        {
            printf("Error in opening a file.. %s", outputFile);
            exit(1);
        }

        char outputFile2[200];
        sprintf(outputFile2,"%s%s",path,featureFileName);
        cfWritingFilePointer=fopen(outputFile2, "w"); //open outputFile2 in writing mode
        if(cfWritingFilePointer == NULL)
        {
            printf("Error in opening a file.. %s", outputFile2);
            exit(1);
        }

        char outputFile3[200];
        sprintf(outputFile3,"%s%s",path,ruleFileName);
        fileClassifierPopulation=fopen(outputFile3, "w"); //open outputFile3 in writing mode
        if (fileClassifierPopulation == NULL)
        {
            printf("Error in opening a file.. %s", outputFile3);
            exit(1);
        }
        char testFile[100];
        sprintf(testFile,"%s%s",path,"/test_4626sts_10k_100k_200CF_33_testing_1k_1034.txt");
        testPerformance = fopen(testFile, "w"); //open outputFile1 in writing mode
        if (testPerformance == NULL)
        {
            printf("Error in opening a file.. %s", testFile);
            exit(1);
        }
        if (currentProblemLevel>1)
        {
            char inputFile[200];
            sprintf(inputFile,"%s%s",path,"/Previous_features.txt");
            cfReadingFilePointer=fopen(inputFile, "r"); //open inputFile in reading mode
            if(cfReadingFilePointer == NULL)
            {
                printf("Error in opening a file.. %s", inputFile);
                exit(1);
            }
        }


  //      printf("%d have %d\n", omp_get_thread_num(), j);
        printf("%dText Classification : %d\n",condLength,j);
        startXCS();
        printf("Done\n");
        Exit(filePerformance);

    }

    fclose(filePerformance);//close outputFile
    fclose(cfWritingFilePointer);//close
    fclose(fileClassifierPopulation);//close

    exit(0);
}//end main
