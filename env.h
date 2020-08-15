float ScaleRange(float Value,float FromMinValue, float FromMaxValue, float ToMinValue,float ToMaxValue);
void updateRange(DataSource data[],int totalRows);
void initializeInput(DataSource inputArray[],int numofRows);
void loadDataFromFile(DataSource data[], const char inputFile[], const int numInstances);
float roundRealValue(float val, int num);
void load_kb(std::string kb_cf_file_name, std::string kb_filter_file_name);
int load_filter(std::string filter_file_name, FilterMap& filters);
int load_code_fragment(std::string cf_file_name, CodeFragmentMap& code_fragments);
int load_classifier(std::string classifier_file_name, ClassifierMap& pop, CodeFragmentMap& code_fragments);
