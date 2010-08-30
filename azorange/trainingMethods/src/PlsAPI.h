/** API for using the PLearn PLS Algorithm.
*
*#include "PlsAPI.h"<BR>
*-llib
*
*This class contains the methods for interfacing with PLearn PLS algorithm and
*methods for wraping Orange standard formated files into PLearn standard
*formated files.
*
*/
#ifndef PlsAPI_h
#define PlsAPI_h

#include "PLS.h"

using namespace std;
using namespace PLearn;

class PlsAPI
{

public:

    /** Default constructor.
    */
    PlsAPI();
    
    /** Destructor.
    */
    ~PlsAPI();
    
    /** Set a PLS Parameter.
    *
    *@param ParameterName Name of the name of the parameter to set.
    *  Available parameter names are:
    *    k             -Number of components (factors) to compute.
    *                   Defaullt = 1
    *    method        -PLS method to use. Can be 'pls1','kernel'
    *                   or 'simpls' 
    *                   Defaullt = kernel
    *    precision     -Precision to compute the eigenvectors
    *                   Defaullt = 1e-6.
    *    v              -Set verbose flag. '0' or '1'
    *                   Defaullt = 0.
    *    debug         -set debug flag. '0' or '1'
    *                   Defaullt = 0.
    *                   Note: On debug, the scratch dirs are not deleted!
    *    sDir          -set scratch directory to sDir.
    *                   Defaullt = "." .
    *
    *@param ParameterValue Value of the value to assig to the parameter specified.
    *
    *@return true if parameter was acceped, false if parameter set failed.
    */
    bool SetParameter(string ParameterName, string ParameterValue);

    /** Get a PLS Parameter.
    *
    *@param ParameterName Name of the name of the parameter to read
    *  Available parameter names are:
    *    components    -Number of components (factors) to compute.
    *                   Defaullt = 1
    *    method        -PLS method to use. Can be 'pls1','kernel'
    *                   or 'simpls' 
    *                   Defaullt = kernel
    *    precision     -Precision to compute the eigenvectors
    *                   Defaullt = 1e-6.
    *    v              -Get verbose flag. '0' or '1'
    *                   Defaullt = 0.
    *    debug         -Get debug flag. '0' or '1'
    *                   Defaullt = 0.
    *    sDir          -Get actual scratchdir
    *                   Defaullt = "." .
    *
    *@return The value of the parameter specified if exists, or empty string 
    *if parameter does not exists.
    */
    string GetParameter(string ParameterName);

    /** Train the PLS algorithm.
    *
    *@param TrainFile Orange formated file for training the PLS algorithm.
    *
    *@return True if the train was successful or false if it failed or the
    *file didn't exist. 
    */
    bool Train(string TrainFile);

    /** Run the PLS algorithm.
    *
    *@param InputVector The input vector for the PLS algorithm, the format is
    *a string with separated values Ex:
    *       "2 2.3 1 5 23 0.12"  
    * NOTE: Missing Values in variables will be set to zero.
    *
    *@param Token The char that is used to separate the values in the input
    *vector, default is ' ' (space)
    *
    *@return The output vector predicted from the PLS algorithm as a string 
    *ex:
    *  "POS"  or  "2.32" 
    *If the InputVector was in wrong format or the Prediction has generated error,
    *it returns a string strating with "ERROR" and the followed by the error 
    *description, it can also return "nan", "inf" or "-inf" in value field.
    */
    string Run(string InputVector, char Token = ' ');
    
    /** Load a PLS model from a file.
    *
    *@param PLSModelDirPath The path (directory) with the PLS model data.
    *
    *@return True if the PLS model was loaded successfully or false if the
    *model was not loaded correctly. 
    */
    bool LoadPLSModel(string PLSModelDirPath);

    /** Save the PLS model to a file.
    *
    *@param PLSModelDirPath The path (a directory: the directory specified 
    *will be created) to save the PLS model data.
    *
    *@return True if the PLS model was saved successfully or false if the
    *model was not saved.
    */
    bool SavePLSModel(string PLSModelDirPath);

    /** Get the Class Var Name.
    *
    *@return The name of the class variable.
    */
    string GetClassVarName();

    /** Get the variables names.
    *
    *@return A tab separated string with the var names if they are set,
    *the descrete variables are preceded by the char '#'
    */
    string GetVarNames(){return mVarNames;};

    /** Get the number of examples in train data.
    *
    *@return The number of examples in train data.
    */
    int GetTrainExamplesNumber();
    
    /** Get the number of descriptors in train data.
    *
    *@return The number of descriptors in train data.
    */
    int GetTrainDescriptorsNumber();    

    /** Get the number of outputs.
    *
    *@return The number of outputs to train the data.
    */
    int GetTrainOutputsNumber();

    /** Get the number of symbols in a variable.
    *
    *@param The variable column number.
    *
    *@return The number of symbols found on the variable of column 'col'.
    */
    int GetVarSymbNum(int col);


    /** Get the processor time used to train the PLS
    *
    *@return The time in seconds used by the processor
    *to train the algorithm, If the PLS was loaded instead 
    *of trained, this function will return the time spent 
    *reading the model; Load data and format conversions are 
    *excluded; The clock() function is used.
    */
    double GetCPUTrainTime(){return mCPUTrainTime;};
    
    /** Get the overall time spent to train the PLS
    *
    *@return The time in seconds passed from the start of 
    *the last train until the train of the PLS was complete
    *If the PLS was loaded instead of trained, this 
    *function will return the last time spent reading the 
    *model;  The clock() function is used.
    */
    double GetWallTrainTime(){return mWallTrainTime;};
    
    /** Get the overall time spent to process the train data files
    *
    *@return The time in seconds spent to process the last train file 
    *(parsing and reading) so that it becames a format accepted 
    *by PLearn; time() function is used.
    */
    double GetWallFileProcessTime(){return mWallFileProcessTime;};

    /** Check if the PLS algorithm is already trained.
    *
    *@return True if the PLS algorithm is trained and false if it is 
    *not trained or the parameters have changeed since last train.
    */
    bool IsTrained();

protected:

    /** Convert file from Orange format to PLearn Format.
    *
    *@param orngFilePath native Orange formated file.
    *
    *@param PLrnFilrPath optional path for output ".amat" fprmated file. The
    *default is ./tmp.amat.
    *
    *@return The path for the output .amat formated file or string "ERROR" if
    *the file could no be converted or didn't exist.
    */
    string FileOrngToPLearn(string orngFilePath,string PLrnFilePath="./tmp.amat");

    /** Splits a string by the token specified.
    *                   
    *@param InputStr The input string to split by the Token.
    *                       
    *@param pSplitStr A pointer to an array of strings to place the splited strings.
    *                       
    *@param Token The token to use for spliting the input string, default token is set to TAB. 
    *
    *@return returns the number of strings found between tokens.
    */
    int Tokenize(const string* InputStr, string** pSplitedStr ,char Token='\t');

    /** Converts an integer to a string.
    *                   
    *@param Number The number to convert to string.
    *
    *@return The string of the input number.
    */
    string NumToString(int Number);

    /** Converts a real to a string.
    *                   
    *@param Number The number to convert to string.
    *
    *@return The string of the input number.
    */
    string NumToString(real Number);

    /** Check if a specific file exists.
    *
    *@param Path The path of the file to test if exists.
    *
    *@return True if the file exists or false otherwise.
    */
    bool IsFile(string Path);

private:
    bool mVerbose;          //flag for verbosity on stdout
    bool mDEBUG;            //flag for DEBUG proposes
    bool mTrained;          //Flag indicating the status of the PLS predictor
    string mSDir;             //Directory to use as scratchDir
    string mVarNames;       //The var names (when data has been loaded)
    int mOutSymbNum;        //The number of symbols in output Variable
    VMat mTrainMatrix;      //The matrix with the data to train the PLS
    PLearn::PLS Predictor;  //The PLearn PLS instance    
    double mCPUTrainTime;    //CPU time used on last PLS training
    double mWallTrainTime;   //Wall time spent on last PLS Training
    double mWallFileProcessTime; //Wall time spent on processing the  
                                 //input Data for Training
    time_t rawtime;         //variable used as aux. for time stamps
};

#endif // _PlsAPI_h_

// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:encoding=utf-8:textwidth=79 :
