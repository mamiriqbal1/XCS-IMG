#ifndef XCS_IMG_ENV_H
#define XCS_IMG_ENV_H

float ScaleRange(float Value,float FromMinValue, float FromMaxValue, float ToMinValue,float ToMaxValue);
void normalize_image(DataSource *data, int totalRows);
void initializeInput(DataSource inputArray[],int numofRows);
void loadDataFromFile(DataSource data[], const char inputFile[], const int numInstances);
float roundRealValue(float val, int num);
void load_kb(std::string kb_cf_file_name, std::string kb_filter_file_name);
void initialize_env();

#endif //XCS_IMG_ENV_H
