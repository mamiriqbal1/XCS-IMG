#include <fstream>
#include "configuration.h"
#include "classifier.h"

//const CodeFragment dontcareCF = {{0,0,OPNOT,OPOR,OPNOP,OPNOP,OPNOP,OPNOP},{0,0,0,0},-1}; //dontcareCF: it is just for completeness, not evaluated at all. It's output value is always 1.
//const CodeFragment dontcareCF = {{OPUNITY,OPNOP,OPNOP,OPNOP,OPNOP,OPNOP,OPNOP,OPNOP},-1}; //dontcareCF: it is just for completeness, not evaluated at all. It's output value is always 1.
CodeFragment addLeafCF(CodeFragment cf, float state[]);
bool equalTwoLeaf(Leaf lf1, Leaf lf2);

void initializeCFPopulation(FILE *cfReadingFilePointer);//, FILE *cfWritingFilePointer);
void getPreviousCFPopulation(FILE *cfReadingFilePointer);
opType getOpType(char str[]);
bool isExists(CodeFragment newCF, CodeFragment cfPopulation[], int numCFs);
bool equalTwoCFs(CodeFragment cf1, CodeFragment cf2);
int getNumPreviousCFs();
bool isDontcareCF(CodeFragment cf);
int numberOfNonDontcares(CodeFragment cond[]);
void printCF(CodeFragment cf);
void validateDepth(opType* cf, opType* end);
CodeFragment createNewCF(int id);
void storeCFs(ClassifierSet* population, FILE *cfWritingFilePointer);

bool isGeneralCF(CodeFragment cf1, CodeFragment cf2);
bool isMoreGeneralLeaf(Leaf lf1, Leaf lf2);

//int evaluateCF(opType cf[], int state[]);
int evaluateCF(CodeFragment cf, float state[], int cl_id=-1, int img_id=-1, bool train=true);
bool isPreviousLevelsCode(const opType code);

int getNumberOfArguments(const opType opcode);
opType leafOpCode(const int r);
opType randomLeaf();
int validLeaf(const opType opcode);
opType randomFunction();
opType* randomProgram(opType* prog,const int isfull,const int maxDepth, const int minDepth);
char* leafname(const opType code);

void DepthMax(const opType* const end,opType** prog, int& argstogo, int& depth);
///////////////////
//char* opchar(const opType code);
//void outprog_bin(const opType* prog, int size);
//void outprog(const opType* const prog, int size, FILE *fp);
Leaf leafNode(std::string str);
char* leafname(const Leaf leaf);
char* leafInterval(const Leaf leaf);
char* opchar(opType code);
void outprog_bin(const opType* prog, int size);
void outprog(CodeFragment prog,int size, FILE *fp);
/////////////


bool evaluate_filter(const Filter& filter, float state[], int cl_id=-1, int img_id=-1, bool train=true);
void update_evaluation_cache(std::forward_list<int>& removed_filters);
void print_filter_evaluation_stats();