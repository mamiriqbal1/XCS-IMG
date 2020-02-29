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
//using namespace std;

char globalBuf[1000];
int numPreviousCFs = 0;     //When to set
int startingPreviousCFID = 0;   //what is this value

CodeFragment *previousCFPopulation;

void initializeCFPopulation(FILE *cfReadingFilePointer)//, FILE *cfWritingFilePointer)
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
        previousCF = createNewCF(numReadPreviousCFs+condLength); // start CFs ID numbers from condLength to avoid confusion with CFs and condition bits
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
                previousCF.leaf[lfNum] = leafNode(pch);
                previousCF.codeFragment[i++] = lfNum;
                lfNum++;
                //std::cout<<"Pch with D: "<<pch<<std::endl;

                //std::cout<<"\n"<<pch<<"\n";
                //getchar();
            }
            else
            {
                //printf("\n%d\n",*pch);
                c = getOpType(pch);
                previousCF.codeFragment[i++] = c;

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
Leaf leafNode(std::string str)
{
    Leaf leaf;
    //char *cstr = new char[strlen(str) + 1];
    //strcpy(cstr, str);
    int length = str.length();
    std::size_t found = str.find(" ");

    std::string strFN = str.substr(1,found);
    std::stringstream convertInt(strFN);
    convertInt>> leaf.featureNumber;

    std::size_t foundLB = str.find(" ",found+1,1);
    std::string strLB = str.substr(found+1,foundLB);
    std::stringstream convertFloat(strLB);
    convertFloat>> leaf.lowerBound;

    std::string strUB = str.substr(foundLB+1,length-foundLB);
    std::stringstream convertFloatU(strUB);
    convertFloatU>> leaf.upperBound;

    //std::cout<<"leaf No: "<<leaf.featureNumber<<" LB: "<<leaf.lowerBound<<" UB: "<<leaf.upperBound;
    return leaf;
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

bool isExists(CodeFragment newCF, CodeFragment cfPopulation[], int numCFs)
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
bool equalTwoCFs(CodeFragment cf1, CodeFragment cf2)
{

    for(int i=0; i<cfMaxLength; i++)
    {
        if(cf1.codeFragment[i] == cf2.codeFragment[i])
        {
            if(0<=cf1.codeFragment[i] && cf1.codeFragment[i]<numLeaf)
            {
            //std::cout<<"e"<<std::endl;
                if(!equalTwoLeaf(cf1.leaf[cf1.codeFragment[i]], cf2.leaf[cf2.codeFragment[i]]))
                    return false;
            }
        }else
        {
            return false;
        }
    }

    return true;
}
bool equalTwoLeaf(Leaf lf1, Leaf lf2)
{
    if(lf1.featureNumber == lf2.featureNumber)
    {
        if(lf1.lowerBound==lf2.lowerBound && lf1.upperBound==lf2.upperBound)
        {
            return true;
        }
        else
            return false;
    }
    else
        return false;
}

// function added by me in new code
bool isGeneralCF(CodeFragment cf1, CodeFragment cf2)
{

    for(int i=0; i<cfMaxLength; i++)
    {
        if(cf1.codeFragment[i] == cf2.codeFragment[i])
        {
            if(0<=cf1.codeFragment[i] && cf1.codeFragment[i]<numLeaf)
            {
            //std::cout<<"e"<<std::endl;
                if(!isMoreGeneralLeaf(cf1.leaf[cf1.codeFragment[i]], cf2.leaf[cf2.codeFragment[i]]))
                    return false;
            }
        }
        else
        {
            return false;
        }
    }

    return true;
}

// function added by me in new code
bool isMoreGeneralLeaf(Leaf lf1, Leaf lf2)
{
    if(lf1.featureNumber == lf2.featureNumber)
    {
        if( (lf1.lowerBound<lf2.lowerBound && lf1.upperBound>=lf2.upperBound) || (lf1.lowerBound<=lf2.lowerBound && lf1.upperBound>lf2.upperBound) )
        //if(lf1.lowerBound==lf2.lowerBound && lf1.upperBound==lf2.upperBound)
        {
            return true;
        }
        else
            return false;
    }
    else
        return false;
}

/*
bool equalTwoCFs(CodeFragment cf1, CodeFragment cf2)
{
    for(int i=0; i<cfMaxLength; i++)
    {
        if(cf1.codeFragment[i] != cf2.codeFragment[i])
        {
            return false;
        }
    }
    return true;
}
*/
int getNumPreviousCFs()
{
    return numPreviousCFs;
}

bool isDontcareCF(CodeFragment cf)
{
    return (cf.cfID == -1) ? true : false;
    //return (cf.codeFragment[0] == OPUNITY) ? true : false;
}

int numberOfNonDontcares(CodeFragment cond[])  //returns the number of specific CFs in cond
{
    int count=0;
    for(int i=0; i<clfrCondLength; i++)
    {
        if(!isDontcareCF(cond[i]))
        {
            count++;
        }
    }
    return count;
}

void printCF(CodeFragment cf)
{
    printf("Printing CF ...... ");
    char* temp = NULL;
    for(int j = 0; j<cfMaxLength; j++)
    {
        if (cf.codeFragment[j]>=0 && cf.codeFragment[j]<numLeaf)
        {
            temp = leafname(cf.leaf[cf.codeFragment[j]]);
            printf("%s ",temp);
        }
        else
        {
            temp = opchar(cf.codeFragment[j]);
            printf("%s ",temp);
        }
        if(cf.codeFragment[j]==OPNOP)
                break; //reduce size of output
    }
    printf(" ------> %d",cf.cfID);
    printf("\n");
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

CodeFragment createNewCF(int id)
{
    CodeFragment newCF;
    for(int i=0; i<cfMaxLength; i++)
    {
        newCF.codeFragment[i] = OPNOP;
    }
    newCF.cfID = id;
    return newCF;
}

void storeCFs(ClassifierSet *population, FILE *cfWritingFilePointer)
{
    double avgFitness = getAvgFitness(population);
    int numFitterCFs = getNumFitterCFs(population,avgFitness);
    int firstCFID = 0;

    char *buf;
    int len;

    // number of CFs
    //sprintf(buf,"%d\n",numFitterCFs); fwrite(buf,strlen(buf),1,cfWritingFilePointer); // number of CFs
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
        firstCFID = previousCFPopulation[0].cfID;
    }
    // cfID of first CF
    //sprintf(buf,"%d\n",firstCFID); fwrite(buf,strlen(buf),1,cfWritingFilePointer); // cfID of first CF
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
        outprog(previousCFPopulation[i],cfMaxLength, cfWritingFilePointer);

        //sprintf(buf," ---------> %d",previousCFPopulation[i].cfID); fwrite(buf,strlen(buf),1,cfWritingFilePointer);
        len = snprintf(NULL,0," ---------> %d",previousCFPopulation[i].cfID);
        if(!(buf = (char*)malloc((len + 1) * sizeof(char))))
        {
            printf("\nError in file writing ...\n");
            exit(0);
        }
        len = snprintf(buf,len+1," ---------> %d",previousCFPopulation[i].cfID);

        fwrite(buf,strlen(buf),1,cfWritingFilePointer);
        free(buf);

        fwrite("\n",strlen("\n"),1,cfWritingFilePointer);

        fflush(cfWritingFilePointer);
    }


    for(ClassifierSet* set=population; set!=NULL; set=set->next)  //CFs from the current problem
    {
        if(set->classifier->fitness <= avgFitness)  //store CFs from classifiers with fitness > average fitness of the classifiers population
        {
            break;
        }
        for(int i=0; i<clfrCondLength; i++)
        {
        if(!isDontcareCF(set->classifier->condition[i]))
            {
                //outprog(set->classifier->condition[i].codeFragment,cfMaxLength,cfWritingFilePointer);
                outprog(set->classifier->condition[i],cfMaxLength,cfWritingFilePointer);
               //sprintf(buf," ---------> %d",set->classifier->condition[i].cfID); fwrite(buf,strlen(buf),1,cfWritingFilePointer);
                len = snprintf(NULL,0," ---------> %d",set->classifier->condition[i].cfID);
                if(!(buf = (char*)malloc((len + 1) * sizeof(char))))
                {
                    printf("\nError in file writing ...\n");
                    exit(0);
                }
                len = snprintf(buf,len+1," ---------> %d",set->classifier->condition[i].cfID);

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
CodeFragment addLeafCF(CodeFragment cf){
    int leafNum = 0;
    for(int i=0; i<cfMaxLength; i++)
    {
        const opType opcode = cf.codeFragment[i];
        //printf("%d ",opcode);
        if(opcode == OPNOP)
        {
            break;
        }
        if(0<=opcode && opcode<condLength)  //condition bit
        {
            cf.leaf[leafNum].featureNumber = opcode;
            float lower = drand();
            float upper = drand();
            if(lower<=upper)
            {
                cf.leaf[leafNum].lowerBound = roundRealValue(lower,precisionDigits);
                cf.leaf[leafNum].upperBound = roundRealValue(upper,precisionDigits);
            }
            else
            {
                cf.leaf[leafNum].lowerBound = roundRealValue(upper,precisionDigits);
                cf.leaf[leafNum].upperBound = roundRealValue(lower,precisionDigits);
            }
            cf.codeFragment[i] = leafNum;
            leafNum++;
        }/*
        if(opcode>=condLength)
        {
            cf.leaf[leafNum].featureNumber = opcode;
            cf.codeFragment[i] = leafNum;
            leafNum++;
        }*/

    }
return cf;

}

int evaluateCF(CodeFragment cf, float state[]){
    //printCF(cf);
    int stack[cfMaxStack];
    stack[0] = 0;
    int SP = 0;
    int tmptmpfval1 = 0;
    //int x = 0;
    for(int i=0; /*i<cfMaxLength*/; i++)
    {
        const opType opcode = cf.codeFragment[i];
        //std::cout<<"opcode no: "<<opcode<<"\n";
        //printf("%d ",opcode);
        if(opcode == OPNOP)
        {
            break;
        }
        if(isPreviousLevelsCode(opcode))  //CF from any previous level
        {

            int valueOfCF = evaluateCF(previousCFPopulation[opcode - condLength], state);
            stack[SP++] = valueOfCF;
        }
        else if(0<=opcode && opcode<numLeaf)  //condition bit
        {

            if(cf.leaf[opcode].lowerBound<=state[cf.leaf[opcode].featureNumber] && state[cf.leaf[opcode].featureNumber]<=cf.leaf[opcode].upperBound)
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

bool isPreviousLevelsCode(const opType code){
    if(use_kb)
    {
        return ( previousCFPopulation[0].cfID<=code && code<(previousCFPopulation[0].cfID + numPreviousCFs) ) ? true : false;
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
    if(numPreviousCFs==0 || use_kb)
    {
        leaf = irand(condLength);
        //printf("leaf_1 %d\n",leaf);
        return leaf;
    }
    double p = drand();
    if(p < 0.5)
    {
        leaf = irand(condLength);
        //printf("leaf_2 %d\n",leaf);
        return leaf;
    }
    int n = irand(numPreviousCFs);
    leaf = previousCFPopulation[n].cfID;
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
opType* randomProgram(opType* prog,const int isfull,const int maxDepth, const int minDepth){
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

//void outprog(const opType* const prog, int size, FILE *fp){
void outprog(CodeFragment prog, int size, FILE *fp){
    //printf("\nDisplaying Program...\n");
    //char* temp = NULL;
    opType code;

    for(int j = 0; j<size; j++)
    //for(int j = 0; j<cfMaxLength; j++)
    {
        char* temp = NULL;
        code = prog.codeFragment[j];
        if(0<=code && code<numLeaf)
        {
            temp = leafInterval(prog.leaf[code]);

            //temp = "D"+leafDes;
            //sprintf(temp,"D%d ",leafDes);

            //sprintf(globalBuf,"D%d ",validLeaf(code));
            //temp = leafInterval(prog.leaf);
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
        if(prog.codeFragment[j]==OPNOP)
            break; //reduce size of output
    }
    //printf("\n");
    //fwrite("\n",strlen("\n"),1,fp);
}
// ------------------------- End Display Functions ----------------------------------
