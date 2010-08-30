#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <list>
#include <set>

#include "ml.h"
#include "cv.h"
#include "randomForest.h"

using namespace std;

int main( int argc, char** argv )
// * Program to train an OpenCV random forest model and save the model to disk in modelFileName.
{
   
    // ***********************************************************************************
    // *                    Handle the input arguments                                   *
    // ***********************************************************************************
    
    if (argc != 13)
    {
        cout<<"\nERROR!"<<endl;
        cout<<"Wrong number of input arguments!"<<endl;
        cout<<"No random forest model built!"<<endl;
        return 1;
    }
    else
    { 
        // Location of train data in Orange format
        string trainFile = argv[1];                      

        // Location in which to save the computed RF model
        const char* modelFileName = argv[2];              

        string nVars = argv[3];                  
        string nExamples = argv[4];                  
        string maxDepth = argv[5];                  
        string minSample = argv[6];                    
        string useSurrogates = argv[7];           
        string getVarVariance = argv[8];               
        string nActVars = argv[9];                     
        string nTrees = argv[10];                       
        string forestAcc = argv[11];                 
        string termCrit = argv[12];                    


        // *****************************************************************************************
        // *                         Train the model and save to modelFileName                     *
        // *****************************************************************************************
        randomForest RFmodel;
        RFmodel.train(trainFile, modelFileName, nVars, nExamples, maxDepth, minSample, useSurrogates,\
                     getVarVariance, nActVars, nTrees, forestAcc, termCrit, "false");
        return 0;	
    }
}


randomForest::randomForest()
{
}	


bool randomForest::train(string trainFile, const char* modelFileName, string nVarsStr, string nExamplesStr, \
                 string maxDepthStr, string minSampleStr, string useSurrogatesStr, string getVarVarianceStr, \
                 string nActVarsStr, string nTreesStr, string forestAccStr, string termCritStr, string stratifyStr)
//
//  Function which trains an OpenCV random forest model from a data set on disk in Orange format. 
//  The data set location is given by trainFile. 
//  The trained model is written to modelFileName. Several model parameters are given as input arguments (see def below). 
//
{

    bool isTrained = false;

    //**************************************************************************************************//
    //       Convert the input arguments from string (from python) to the type required by OpenCV. 
    //**************************************************************************************************//

    
    // The number of variables in the data set (used to check the reading of data from disk)
    int nVars = atoi(nVarsStr.c_str());

    // The number of examples in the data set
    int nExamples = atoi(nExamplesStr.c_str());

    // Max number of nodes along one branch. Default 5.
    int maxDepth = atoi(maxDepthStr.c_str());

    // Minimum number of samples for which to continue branching. Default 10.
    int minSample = atoi(minSampleStr.c_str());

    // Missing value handling. Default false.
    bool useSurrogates = false;
    if (useSurrogatesStr == "true")
       useSurrogates = true; 

    // Calculate variable variance. Default false.
    bool getVarVariance = false;
    if (getVarVarianceStr == "true")
       getVarVariance = true; 

    // Number of variables considered at each node. 0 => default sqrt(Nvariables)
    int nActVars = atoi(nActVarsStr.c_str());

    // The max number of trees to build - used if termCrit = 0. Default 50.
    int nTrees = atoi(nTreesStr.c_str());

    // Desired out of bag error - used if termCrit = 1. Default 0.1.
    double forestAcc = atof(forestAccStr.c_str());

    // Determines which termination criteria to use. 0 => nTrees, 1 => forestAcc.
    int termCrit = atoi(termCritStr.c_str());

    // Stratified data sampling, 1 => stratify sampling
    bool stratify = false;
    if (stratifyStr == "true")
       stratify = true; 
	
    
    //**************************************************************************************************//
    //                      Check for irrational input arguments
    //**************************************************************************************************//
    if (minSample >= nExamples)
    {
        cout<<"ERROR! Invalid minSample: "<<minSample<<endl;
        cout<<"minSample must be smaller than the number of examples."<<endl;
        cout<<"The number of examples is: "<<nExamples<<endl;
        if (nExamples > 10)          
        {
            cout<<"minSample assigned to default value: 10"<<endl;
            minSample = 10;
        }
        else
        {
            cout<<"Too few examples!!"<<endl;
            cout<<"Terminating"<<endl;
            cout<<"No random forest model built"<<endl;
            return isTrained;
        }
    }
    
    if (nActVars > nVars)
    {
        cout<<"ERROR! Invalid nActVars: "<<nActVars<<endl;
        cout<<"nActVars must be smaller than or equal to the number of variables."<<endl;
        cout<<"The number of variables is: "<<nVars<<endl;
        cout<<"nActVars assigned to default value: sqrt(nVars)"<<endl;
        nActVars = 0;
    }
        
	
    //**************************************************************************************************//
    // Read the data from disk into CvMat objects. Orange format of exampleTable with classes is assumed. 
    //**************************************************************************************************//

    CvMat* cvResp = NULL;      // Response vector 
    CvMat* cvData = NULL;      // Train data
    CvMat* var_type = NULL;    // Variable type vector 
    //CvMat* cvRespTest = NULL;      // Response vector test data
    //CvMat* cvDataTest = NULL;      // Test data 
    string rc = orngFileToCvMat(trainFile, &cvResp, &cvData, &var_type);


    //**************************************************************************************************//
    //                 Check the CvMat data read from disk
    //**************************************************************************************************//
    if (rc != "SUCCESS")
    {
        cout<<"ERROR! Could not get training data."<<endl;
        cout<<"Terminating. No random forest model built."<<endl;
        return isTrained;
    }
        
    // Check the number of variables in the training set read from disk
    if (cvData->cols != nVars)
    {
        cout<<"ERROR! Wrong number of variables in the training data."<<endl;
        cout<<"Terminating. No random forest model built."<<endl;
        return isTrained;
    }
 
    // Check the number of examples in the training set read from disk
    if (cvData->rows != nExamples)
    {
        cout<<"ERROR! Wrong number of examples in the training data."<<endl;
        cout<<"Terminating. No random forest model built."<<endl;
        return isTrained;
    }

    // Check the lenght of the response vector
    if (cvResp->rows != nExamples)
    {
        cout<<"ERROR! Wrong number of responses in the training data."<<endl;
        cout<<"Terminating. No random forest model built."<<endl;
        return isTrained;
    }

    // Check the lenght of the var_type vector
    if (var_type->rows != (nVars+1))
    {
        cout<<"ERROR! Wrong number of elements in the var_type vector."<<endl;
        cout<<"Terminating. No random forest model built."<<endl;
        return isTrained;
    }


    //**************************************************************************************************//
    //                                 Set model parameters
    //**************************************************************************************************//
    float priors[2] = {1,1};
    int cvTermCrit;
    if (termCrit == 1)
        cvTermCrit = CV_TERMCRIT_EPS;
    else
        cvTermCrit = CV_TERMCRIT_ITER; 
    CvRTParams Params( maxDepth,                //Max depth 
                       minSample,               //Min sample count 
                       0,                       //Regression accuracy 
                       useSurrogates,           //Use Surrogates 
                       15,                      //Max Categories 
                       priors,                  //Priors  
                       getVarVariance,          //Calculate Variables variance 
                       nActVars,                //Nactive Vars 
                       nTrees,                  //Max Tree Count 
                       forestAcc,               //Accuracy to archive 
                       cvTermCrit               //Criteria to stop algorithm : CV_TERMCRIT_ITER  or   CV_TERMCRIT_EPS(Rate determining)
                       );
                                                	
    //**************************************************************************************************//
    //                                  Train
    //**************************************************************************************************//
    try
    {
        //cout<<"stratify in randomForest.cc "<<stratify<<endl;
        RTmodel.train(cvData, CV_ROW_SAMPLE, cvResp, NULL, NULL, var_type, NULL, Params, stratify);
        RTmodel.save(modelFileName);
        isTrained = true;
    }
    catch (...)
    {
        cout<<"Training of random forest model failed."<<endl;
    }

    // Predict
    //printClassAccuracy(cvData, cvResp, &RTmodel);
    //printClassAccuracy(cvDataTest, cvRespTest, &RTmodel);
    
    ////scPA
    cvReleaseMat(&cvResp);
    cvReleaseMat(&cvData);
    cvReleaseMat(&var_type);
    ////ecPA

    return isTrained;
}	


bool randomForest::load(const char* modelFileName)
{
    bool isLoaded;

    try
    {
        RTmodel.load(modelFileName);
        isLoaded = true;
    }
    catch (...)
    {
        isLoaded = false;
        cout<<"Loading of random forest model failed."<<endl;
    }
    return isLoaded;
   
}


double randomForest::predict(string example, char token)
// Predicts one example given as a token separated string. 
{
    string* InputValues;           //Array to store the splited input values

    // Create a CvMat object from example which is a tab sep string holding the values of an example to be predicted
    int InVecSize=Tokenize(&example, &InputValues , token); 

    CvMat* predData = cvCreateMat(InVecSize, 1, CV_32FC1);   // Allocate memory for a CvMat object
    double value; 
    for(int idx=0; idx<InVecSize; idx++)
    {
        value = atof((InputValues[idx]).c_str());
        cvmSet(predData, idx, 0, value);
    }
    //cout<<"CvMat object to be predicted!"<<endl;
    //testCvMat(predData);

    double prediction = RTmodel.predict(predData);
    ////scPA
    cvReleaseMat(&predData);
    delete[] InputValues;
    ////ecPA
    return prediction;

}


bool randomForest::save(const char* modelFileName)
{
    bool isSaved;

    try
    {
        cout<<"Saving to file: "<<modelFileName<<endl;
        RTmodel.save(modelFileName);
        isSaved = true;
    }
    catch (...)
    {
        isSaved = false;
        cout<<"Saving of random forest model failed"<<endl;
    }
    return isSaved;

}


string randomForest::orngFileToCvMat(string orngFilePath, CvMat** cvResp, CvMat** cvData, CvMat** var_type)
{
    string line;
    int Cols=0;
    int Rows=0;
    int actPos=-1;
    int nxPos=0;
    char token='\t';
    int inFileCols=0;
    int dataCols=0;
    string classElement="";
    string elementName="";
    string element="";
    //Class attribute specifications
    int classCol=-1;
    string className="";
    //Ignored Rows mask
    string ignoredRows=""; 
    
    string* firstThreeOrngRows[3];

    ifstream orngFile(orngFilePath.c_str());
    //ofstream PLrnFile(PLrnFilePath.c_str());

    // Create an stl contrainers to store the data set as stl strings
    vector<double> DataVector;     	
    vector<int> varTypeVector;
    list<string> RespList;       	
    
    // Count the number of data rows    	
    int NdataRows = 0;   	
    
    // Categorical or numeric response column
    bool isClassifier;
        	
    if (orngFile.is_open())
    {
        //.amat commants line
        while (! orngFile.eof() )     // Rows is incremented
        {
            getline (orngFile,line);
            if(line.find(token)!=string::npos)
            {                                           
		    if (Rows<3)
            {
		        inFileCols=Tokenize(&line,&(firstThreeOrngRows[Rows]));
		    }
            else if(Rows==3)//Compute Header
            {
			    if(inFileCols<2)
			    {	
				delete[] firstThreeOrngRows[0];			
				delete[] firstThreeOrngRows[1];			
				delete[] firstThreeOrngRows[2];			
			        return "ERROR";
			    }
			    dataCols=inFileCols;
			    ignoredRows=string(inFileCols,'0');
			    for(int i=0;i<inFileCols;i++)
			    {
			    	//find ignored rows and mark them on mask
			        if(((firstThreeOrngRows[2][i]).compare("i")==0) || 
					    ((firstThreeOrngRows[2][i]).compare("m")==0) || 
					    ((firstThreeOrngRows[2][i]).compare("meta")==0) || 
					    ((firstThreeOrngRows[2][i]).compare("ignore")==0))
			        {
				        ignoredRows[i]='1';
				        dataCols--;
			        }
			       //find class column
			       if((((firstThreeOrngRows[2][i]).compare("c")==0) || ((firstThreeOrngRows[2][i]).compare("class")==0)))
				   classCol=i;
			    }    // end for loop inFileCols
        			
        	    //if class column not found, choose the last column as class
		        if((classCol<0) || (ignoredRows[classCol]=='1'))
				classCol=ignoredRows.rfind('0');
			         
			    //header of columns names
			    for(int i=0;i<inFileCols;i++)
                {
				   if(ignoredRows[i]=='0')  // Variable names
				   {
				      if( ((firstThreeOrngRows[1][i]).compare("c")==0) || ((firstThreeOrngRows[1][i]).compare("continuous")==0) )
                                elementName=firstThreeOrngRows[0][i];
                      else
                            elementName="#"+firstThreeOrngRows[0][i];
                      //elementName+=(firstThreeOrngRows[0][i]);
                      for(unsigned int j=0;j<elementName.length();j++)
                             elementName[j] = elementName[j]==' '?'_':elementName[j];                                       
				      if(i==classCol)
				           className+=elementName;
				    }  // end if ignoredRows
				         	
				    // Assign the varType vector (iterates over all cols here). Remove class column and ignored columns
				    if ((i != classCol) && (ignoredRows[i] != '1'))
				    {
				  	     if( (((firstThreeOrngRows[1][i]).compare("c")==0) || ((firstThreeOrngRows[1][i]).compare("continuous")==0)))
				   	     {
				  	         varTypeVector.push_back(CV_VAR_NUMERICAL); 
				   	     }
				   	     else 
				   	         varTypeVector.push_back(CV_VAR_CATEGORICAL);
				    }    
				    else
				    {
				      	if( (((firstThreeOrngRows[1][i]).compare("c")==0) || ((firstThreeOrngRows[1][i]).compare("continuous")==0)))
				       	    isClassifier = false;
				       	else
				      	    isClassifier = true;    
				    }	
				              
			      }   // End for inFileCols
		        }  // end else Rows == 3
		    if(Rows>=3)//Compute Data Rows
		    {
			    Cols=0;
			    classElement="";
             	actPos=-1;
              	nxPos=line.find(token,0);
               	while((unsigned int)nxPos!=string::npos)//read all row elements but the last
               	{
                    element =line.substr(actPos+1,nxPos-actPos-1);
                    element = element==""?"0":element;
			        if(Cols==classCol)
			        {
			            classElement=element;
			            RespList.push_back(element);				            
			        }    
			        else if(ignoredRows[Cols]=='0')
			        {
                        DataVector.push_back(atof(element.c_str()));
			        }
			        //cout<<"DataVector "<<DataVector[0]<<endl;
                  	actPos=nxPos;
                   	nxPos=line.find(token,actPos+1);
                   	Cols++;
            	}  // end while loop reading all row elements  
                if((unsigned int)actPos<(line.length()-1))//read the last column element
                {   
                    element =line.substr(actPos+1,nxPos-actPos-1);
                    element = element==""?"0":element;
                    if(Cols==classCol)
                    {
                       classElement = element;
                       RespList.push_back(element);
                    }   
			        else if(ignoredRows[Cols]=='0')
			        {
                 	   DataVector.push_back(atof(element.c_str()));
                    }   
                  	Cols++;
                   	NdataRows++;
                  }  // end if reading last element
		       }  // end if Compute Data Rows   
               Rows++;
               } // end  if(line.find(token)!=string::npos)
          } // end while loop
            
        // Collect the varType vector in a CvMat object
        *var_type = getVarType(varTypeVector, isClassifier);
        //testVarTypeVector(varTypeVector, isClassifier);
            
        // Collect the training data in a CvMat object
        *cvData = getTrainData(DataVector, NdataRows);
            
        // Collect the response in a CvMat object
        bool reMapClasses = false;  // Remap in the python layer because of back transformation
        *cvResp = getResponse(RespList, isClassifier, reMapClasses);
                
	    orngFile.close();
    } // orngFile open
    else 
    {
	cout << "Unable to open Orange file for read!"<<endl;
	delete[] firstThreeOrngRows[0];
	delete[] firstThreeOrngRows[1];
	delete[] firstThreeOrngRows[2];
    	return "ERROR";
    }  // orngFile open
    delete[] firstThreeOrngRows[0];
    delete[] firstThreeOrngRows[1];
    delete[] firstThreeOrngRows[2]; 
    return "SUCCESS";
}


CvMat* randomForest::getVarType(vector<int> varTypeVector, bool isClassifier)
/*
 * Create the var_type vector specifying the type of each variable
 * The values of varTypeVector is not used at the moment. All variables are assumed to be either num or cat. 
 * Below is what Intel used in the code they sent. 
 */
{
	int length = varTypeVector.size();
    CvMat* var_type = cvCreateMat(length + 1, 1, CV_8U );
    cvSet(var_type, cvScalarAll(CV_VAR_NUMERICAL));
    if (isClassifier)
        cvSetReal1D(var_type, length, CV_VAR_CATEGORICAL);
    else
        cvSetReal1D(var_type, length, CV_VAR_NUMERICAL);

    return var_type;
}	

CvMat* randomForest::getResponse(list<string> RespList, bool isClassifier, bool reMapClasses)
/*
 * Return a CvMat object holding the responses. If the response vector is categorical and reMapClasses is true
 * the values are mapped to values in the range 0,...,N where N is the number of classes.  
 */
{
	int NdataRows = RespList.size();
    CvMat* cvResp = cvCreateMat(NdataRows, 1, CV_32FC1); 
    map<string, double> classMap;
        
    // Create a map between the strings of class values and doubles in the range 0-N    
	if ((isClassifier) && (reMapClasses))
	{
	    // Get a sorted list of unique values (the different class values) - a set container has this functionality
	    set<string> RespSet;
	    list<string>::iterator i;
        for(i=RespList.begin(); i != RespList.end(); ++i)
        {
            RespSet.insert(*i);	
        }	

        // Create a map between the unique values and consecutive numbers (double because of CvMat type req.)
        double classIdx = 0;
        set<string>::iterator ii;
        for(ii=RespSet.begin(); ii != RespSet.end(); ++ii)
        {
            classMap[*ii] = classIdx;
            classIdx++; 	
        }	
        
	}   // End isClassifier   
    
    // Assign the list to a CvMat object and remap the responses with classMap if classifier 
    int rowNr = 0;
    int colNr = 0;
    double value;
    list<string>::iterator i;
    for(i=RespList.begin(); i != RespList.end(); ++i)
    {
        // Use the classMap to get the value of a classifer
        if (isClassifier && reMapClasses)
            value = classMap[*i];
    	// Transform to double for regression
        else
            value = atof((*i).c_str());
 
        cvmSet(cvResp, rowNr, colNr, value);
        rowNr++;
    }		
		
	return cvResp;
}


CvMat* randomForest::getTrainData(vector<double> DataVector, int NdataRows)
/*
 * Stores the values of DataVector in a CvMat object. 
 * DataVector is a stl vector container where the data is stored row wise. 
 * NdataRows determines the partitioning into rows and columns in CvMat.
 * Response is true if the response vector has categorical values in which case 
 * the response values will be strings. 
 */
{

    int lenDataVector = DataVector.size();
    int NdataCols = (int)lenDataVector/NdataRows;

	// Allocate memory for the CvMat object
    CvMat* cvFileData = cvCreateMat(NdataRows, NdataCols, CV_32FC1);
	
	// Assign values to the elements of CvMat
	for (int idx = 0; idx < lenDataVector; idx++)
	{
		int rowNr = (int)floor(double(idx)/double(NdataCols));
		int colNr = idx - (rowNr*NdataCols);
	    double value = DataVector[idx];
	    cvmSet(cvFileData, rowNr, colNr, value);
	}	 
	
	return cvFileData;
}	

void randomForest::testCvMat(CvMat* cvFileData)
{
	// Assess the value of cvFileData
    int Ncols = cvFileData->cols;
    int Nrows = cvFileData->rows;
    cout<<"cols "<<Ncols<<endl;
    cout<<"rows "<<Nrows<<endl;
    int oldRowIdx = 0; 
    for (int idx = 0; idx < Ncols*Nrows; idx++)
	{
	    int rowNr = (int)floor(double(idx)/double(Ncols));
	    if (rowNr != oldRowIdx)
	    {
	        oldRowIdx = rowNr;
	        cout<<endl;
	    }	 
		int colNr = idx - (rowNr*Ncols);
		cout<<cvmGet(cvFileData, rowNr, colNr)<<" ";
    }	 
    cout<<endl;
}	

void randomForest::testVarTypeVector(vector<int> varTypeVector, bool isClassifier)
{
    // Print varTypeVector
    int lenDataVector = varTypeVector.size();
    cout<<"In testVarTypeVector "<<endl;
    if (isClassifier)
        cout<<"This is a classifier "<<endl;
    else    
        cout<<"This is not a classifier "<<endl;
    cout<<"Length "<<lenDataVector<<endl;
    //for (int idx = 0; idx < lenDataVector; idx++)
    //    cout<<varTypeVector[idx]<<" ";
    cout<<endl;    
}

void randomForest::printClassAccuracy(CvMat* Data, CvMat* Resp, CvRTrees* RTmodel)
{

    double test_hr = 0;
    double r;
    double rateSum;
    double prediction;

    int nTestSamples=Data->rows;
    cout<<"Nrows "<<nTestSamples<<endl;
    cout<<"nVars "<<Data->cols;
    for( int i = 0; i < nTestSamples; i++ )
    {
        CvMat _sample;
        cvGetRow( Data, &_sample, i );
        r = RTmodel->predict( &_sample );
        prediction = r;
        //cout<<"Predicted value "<<r<<endl;
        r = (int)r - (int)GetMatData(Resp, i, 0) ==0 ? 1 : 0;           //For Classifiers
        cout<<prediction<<" "<<GetMatData(Resp, i, 0)<<" "<<r<<endl;
        //cout<<"Actural value "<<GetMatData(Resp, i, 0)<<endl;
        //cout<<"Expression "<<((int)r - (int)GetMatData(Resp, i, 0) ==0 ? 1 : 0)<<endl;;               //For Classifiers
        //out<<"r "<<r<<endl;
        test_hr += r;
        //cout<<"Acc sum "<<test_hr<<endl;
    }
    rateSum = test_hr / (double)(nTestSamples);
    cout<<"Number of correctly classified samples "<<test_hr<<endl;
    cout<<"Number of test examples "<<nTestSamples<<endl;
    cout<<"Classification accuracy "<<rateSum<<endl;
    //printf("%d Trees %.1f%% \n",RTmodel.get_tree_count(),(test_hr / (double)(spPerFold))*100.);


    //For Classifiers
    printf( "Mean recognition rate: %.1f%%\n",rateSum*100. );

    cout<<"----------------------------------------------------------------"<<endl;

}

float randomForest::GetMatData(const CvMat* Mat, int row, int col)//zero based indexes
{
    if((row>=Mat->rows) || (col>=Mat->cols))
    {
       cout<<"ERROR: GetMatData: Matrix index out of limits!"<< endl;
       return 0;
    }
        return cvmGet( Mat, row, col );
}

int randomForest::Tokenize(const string* InputStr, string** pSplitedStr ,char Token)
{                                               
        int Cols=0;                             
        int actPos=-1;      
        try{                
                int nxPos=InputStr->find(Token,0);      
                while((unsigned int)nxPos!=string::npos)//Count all row elements but the last
                {                                       
                        actPos=nxPos;
                        nxPos=InputStr->find(Token,actPos+1);
                        Cols++;
                }
                Cols++;                         //Count the last element
                                
                *pSplitedStr=new string[Cols];
                actPos=-1;  
                nxPos=InputStr->find(Token,0);
                Cols=0;     
                while((unsigned int)nxPos!=string::npos)//Get all row elements but the last
        {               
                        (*pSplitedStr)[Cols]=string(InputStr->substr(actPos+1,nxPos-actPos-1).c_str());
            actPos=nxPos;   
            nxPos=InputStr->find(Token,actPos+1);
                        Cols++;
        }                       
        if((unsigned int)actPos<(InputStr->length()-1))//Get the last column element
        {
                (*pSplitedStr)[Cols]=string(InputStr->substr(actPos+1).c_str());
        }
                else
                        (*pSplitedStr)[Cols]=string("");
                Cols++;
        }
        catch(...)
        {
                Cols=0;
                //cout<<"ERROR:Probably bad formated Orange File!"<<endl;
        }
        return Cols;
}
