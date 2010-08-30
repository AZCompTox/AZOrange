#include <string>
#include <vector>
#include <list>

#include "ml.h"
#include "cv.h"
#include "cxcore.h"

using namespace std;


class randomForest
{
public:

    randomForest();

    bool train(string trainFile, const char* modelFileName, string nVarsStr, string nExamplesStr, \
          string maxDepthStr, string minSampleStr, string useSurrogatesStr, string getVarVarianceStr, \
          string nActVarsStr, string nTreesStr, string forestAccStr, string termCritStr, string stratifyStr);
          
    double predict(string example, char token);

    bool save(const char* modelFileName);

    bool load(const char* modelFileName);

protected:	

    string orngFileToCvMat(string orngFilePath, CvMat** cvResp, CvMat** cvData, CvMat** var_type);
    void testVarTypeVector(vector<int> varTypeVector, bool isClassifier);
    void testCvMat(CvMat* cvFileData);
    CvMat* getTrainData(vector<double> DataVector, int NdataRows);
    CvMat* getResponse(list<string> RespList, bool isClassifier, bool reMapClasses);
    CvMat* getVarType(vector<int> varTypeVector, bool isClassifier);
    void printClassAccuracy(CvMat* Data, CvMat* Resp, CvRTrees* RTmodel);
    float GetMatData(const CvMat* Mat, int row, int col);
    int Tokenize(const string* InputStr, string** pSplitedStr ,char Token='\t');

private:

    CvRTrees RTmodel;          //Random Trees instance

};	
