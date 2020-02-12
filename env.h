


float ScaleRange(float Value,float FromMinValue, float FromMaxValue, float ToMinValue,float ToMaxValue);
//void updateRange(DataSource data[]);
void updateRange(DataSource data[],int totalRows);
void initializeInput(DataSource inputArray[],int numofRows);
//DataSource resetState(DataSource inputArray[]);
void resetStateTesting(DataSource &state, int index);
float roundRealValue(float val, int num);
void loadData(DataSource data[]);
//void loadData(ds data[]);
void initializeLimit();
void loadDataFile();
void loadFileString();
void splitData(DataSource totalData[],DataSource trainData[], DataSource testData[]);
bool isIntegerInArray(int integer, int array[], int highestPosition);
