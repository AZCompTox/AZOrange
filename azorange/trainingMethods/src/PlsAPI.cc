
/*
#include <plearn/io/pl_log.h>
#include <plearn/math/TMat_maths.h>    //!< For dist.
#include <plearn/math/pl_erf.h>
#include <plearn/math/plapack.h>
#include <plearn/vmat/ShiftAndRescaleVMatrix.h>
#include <plearn/vmat/SubVMatrix.h>
#include <plearn/vmat/VMat_linalg.h>
 */
#include <plearn/db/getDataSet.h>
#include <time.h>
#include <plearn/io/load_and_save.h>
#include <iostream>
#include <string>
#include <stdio.h>
#include <fstream>
#include <sys/stat.h>
#include <map>
#include "PlsAPI.h"

using namespace std;
using namespace PLearn;


string& strip(string& context)
{
    string::size_type pos1 = context.find_first_not_of(' ');
    string::size_type pos2;
    if (pos1>0)
      context.erase(0,pos1);
    pos2 = context.find_last_not_of(' ');
    if (pos2 < (context.size()-1))
      context.erase(pos2+1,context.size()-pos2-1);
    return context;
}


string& replaceAll(string& context, const string& from, const string& to)
{
    size_t lookHere = 0;
    size_t foundHere;
    while((foundHere = context.find(from, lookHere)) != string::npos)
    {
          context.replace(foundHere, from.size(), to);
          lookHere = foundHere + to.size();
    }
    return context;
}


PlsAPI::PlsAPI()
{
    //Define the verbose and debug flag
    mVerbose=false;
    mDEBUG=false;
    Predictor.SetVerbosity(mVerbose);
    //Assign the default PLS parameters
    //Predictor.setOption("output_the_score","1");
    Predictor.setOption("k","1");
    Predictor.setOption("method","kernel");
    Predictor.setOption("precision","1e-6");
    mSDir = ".";
    mVarNames="";
    mTrained = false;
    mCPUTrainTime = 0;    
    mWallTrainTime = 0;   
    mWallFileProcessTime = 0;
    mOutSymbNum = 0;
}
PlsAPI::~PlsAPI()
{
    Predictor.forget();
}

bool PlsAPI::SetParameter(string ParameterName, string ParameterValue)
{
    struct stat st;

    if (ParameterName.compare("sDir")==0)
    {
        if(stat(ParameterValue.c_str(),&st) == 0)
        {
            mSDir=ParameterValue;
        }
        else
        {
            mSDir=".";
        }
    }
    else if (ParameterName.compare("v")==0)
    {
        if( ParameterValue.compare("0")!=0)
        {
            mVerbose=true;
            Predictor.SetVerbosity(mVerbose);
        }
        else
        {
             mVerbose=false;
             Predictor.SetVerbosity(mVerbose);
        }
    }
    else if (ParameterName.compare("debug")==0)
    {
        if( ParameterValue.compare("0")!=0)
        {
            mDEBUG=true;
        }
        else
        {
             mDEBUG=false;
        }
    }
    else
    {
        try
        {   
            //Assign the desired parameter to the desired value
            Predictor.setOption(ParameterName,ParameterValue);
        }
        catch(...)
        {   
            if (mVerbose) cout<<"ERROR: Parameter '"<<ParameterName<<"' could not be set to '"<<ParameterValue<<"'. Probably it is an invalid parameter!"<<endl;
        return false;
        }
    }

    //Compare the parameter value with the one we wanted to set up
    try
    {
        if (GetParameter(ParameterName).compare(ParameterValue)!=0)
        {
            if (mVerbose) cout<<"Warning: Parameter '"<<ParameterName<<"' could not be set to '"<<ParameterValue<<"'. Its actual value is: "<<GetParameter(ParameterName)<<endl;
            return false;
        }
    }
    catch(...)
    {
        if (mVerbose) cout<<"ERROR: Parameter '"<<ParameterName<<"' could not be set to '"<<ParameterValue<<"'. Probably it is an invalid parameter!"<<endl;
        return false;
    }
    //Parameter has changed, so the model is not trained anymore!
    mTrained=false;
    if (mVerbose) cout<<"PARAMETERS: Parameter '"<<ParameterName<<"' has been set to '"<<ParameterValue<<"'."<<endl;
    return true;
}

string PlsAPI::GetParameter(string ParameterName)
{
    if (ParameterName.compare("sDir")==0)       
            return mSDir;
    if (ParameterName.compare("debug")==0){
        if (mDEBUG) 
            return string("1"); 
        else 
            return string("0");
	}
    if (ParameterName.compare("v")==0){
        if (mVerbose) 
            return string("1"); 
        else 
            return string("0");
	}
    string Value(Predictor.getOption(ParameterName));
    if ((Value[0]=='"') && (Value[Value.length()-1]=='"') )
        return Value.substr(1,Value.length()-2);
    return Value;
}

bool PlsAPI::Train(string TrainFile)
{
    srand ( time(NULL) );
    string ScratchDir(mSDir+string("/scratchdirPLS_API_")+NumToString((int)time(NULL))+"_"+NumToString(rand() % 1000));

    if (mVerbose) cout<<"Actual scratchdir  is:   "<<ScratchDir<<endl;
    string PlearnTrainFile(ScratchDir + string("/PLSTrainData.amat"));

    //Create a temp directory to convert the data into PLearn format 
    system(string(string("mkdir -p ") + ScratchDir).c_str());
    //Convert the data from Orange format to PLearn Format
    if (mVerbose) {time ( &rawtime ); cout<<ctime(&rawtime)<<"Converting train file from Orange to PLearn format (PlsAPI procedure)..."<<flush;}
    mWallFileProcessTime=(double)time(NULL);
    PlearnTrainFile = FileOrngToPLearn(TrainFile,PlearnTrainFile);
    mWallFileProcessTime = (double)time(NULL) - mWallFileProcessTime;
    if (mVerbose) {time ( &rawtime ); cout<<ctime (&rawtime)<<"done."<<endl;}

    if((PlearnTrainFile.substr(0,5)).compare("ERROR")==0)
    {
        //Remove the temp directory
        if(not mDEBUG) system(string(string("/bin/rm -rf ") + ScratchDir).c_str());
        return false;
    }
    mOutSymbNum = 0;
    //Free the memory used by last train if any
    Predictor.forget();
    //Load the PLearn formated file
    if (mVerbose) {time ( &rawtime ); cout<<ctime (&rawtime)<<"Loading train file (PLearn procedure)..."<<flush;}
    mWallFileProcessTime = (double)time(NULL) - mWallFileProcessTime;
    mTrainMatrix = getDataSet(PlearnTrainFile);
    mWallFileProcessTime = (double)time(NULL) - mWallFileProcessTime;
    if (mVerbose) {time ( &rawtime ); cout<<ctime (&rawtime)<<"done."<<endl;}
    if (mVerbose) {time ( &rawtime ); cout<<ctime (&rawtime)<<"Time processing input data: "<<mWallFileProcessTime<<" seconds"<<endl;}
 
    //The Training 
    mCPUTrainTime = clock();  
    mWallTrainTime = time(NULL); 
    try{
        //The build() command makes any changes made to propaate to the parent
        //classes, so it's recomended to call it before any important operation
        Predictor.build();
        //Set the PLS train data
        Predictor.setTrainingSet( mTrainMatrix, false);
        Predictor.build();
        //Call the PLS train methon
        Predictor.train();
    }
    catch (...)
    {
        mCPUTrainTime =((double)clock() - mCPUTrainTime)/(double)CLOCKS_PER_SEC;
        mWallTrainTime =(double)time(NULL) - mWallTrainTime;
        //Remove the temp directory if NOT in debug mode       
        if(not mDEBUG) system(string(string("/bin/rm -rf ") + ScratchDir).c_str());
        mTrained=false;
        return false;
    }
    
    mCPUTrainTime =((double)clock() - mCPUTrainTime)/(double)CLOCKS_PER_SEC;
    mWallTrainTime =(double)time(NULL) - mWallTrainTime;
    if (mVerbose) {time ( &rawtime ); cout<<ctime (&rawtime)<<"Time training PLS: "<<mWallTrainTime<<" seconds"<<endl;}
    //Remove the temp directory if NOT in debug mode       
    if(not mDEBUG) system(string(string("/bin/rm -rf ") + ScratchDir).c_str());
    mTrained=Predictor.IsTrained();
    return mTrained;
}

string PlsAPI::Run(string InputVector, char Token)
{
    if(!mTrained)
        return "ERROR: Cannot run PLS before trained or loaded with a model";
    string* InputValues;                //Array to store the splited input values
    Vec Vin(mTrainMatrix->inputsize()); //Vector input in PLearn format 
    Vec Vout(Predictor.outputsize());   //Vector to store the output
                                        //predicted values
    bool categoricalOut=true;           //Flag indicating the output type
    real res=0;                         //Var to handle singular output value
    string OutputVector="";             //The output vector to return:Size=1

    //The PLS predictor only makes prediction for one output variable
    try{
        categoricalOut=(mTrainMatrix.fieldName(mTrainMatrix->inputsize())).at(0)=='#'?true:false;
        //Convert the input data to PLearn Vector Vi
        int InVecSize=Tokenize(&InputVector, &InputValues , Token);
        if(InVecSize<(mTrainMatrix->inputsize()) )
        {
            delete[] InputValues;
            return string("ERROR: InputVector with too few descriptors(")
                    + NumToString(InVecSize) + string(") the model needs ") 
                    + NumToString(mTrainMatrix->inputsize()) + string(" descriptors");
        }
        for(int i=0;i<(mTrainMatrix->inputsize());i++)
        {
			//Strip and replace all spaces with "_nbsp_"
			strip(InputValues[i]);
            replaceAll(InputValues[i]," ","_nbsp_");
            if((mTrainMatrix.fieldName(i)).at(0)=='#')
            {
                //*InputValues[i] = string("D") + *InputValues[i]; 
                Vin[i] =
                    isnan( mTrainMatrix->getStringVal(i,InputValues[i]))?
                        (real)atof((InputValues[i]).c_str()):
                        mTrainMatrix->getStringVal(i,InputValues[i]);
                if (mDEBUG) cout<<i<<":"<<Vin[i]<<"-"<<InputValues[i].c_str()<<endl;
            }
            else
            {
                Vin[i] = (real)atof((InputValues[i]).c_str());
                if (mDEBUG) cout<<i<<":"<<Vin[i]<<"-"<<InputValues[i].c_str()<<endl;
            }
        }
        //Run the PLS to predict the output
        Predictor.computeOutput(Vin,Vout);
        res=Vout[Vout.size()-1];
        //cout<<"Res: "<<res<<endl;
        for(int idx=0;idx<Vout.size()-1;idx++)
           OutputVector += NumToString(Vout[idx]) + " ";
        if(categoricalOut)
        {
            if (mDEBUG) cout<<"Out is Categ."<<endl;
            //The output is Categorical (Classifier)
            OutputVector = mTrainMatrix->getValString(mTrainMatrix->inputsize(),round(res));
            if (OutputVector=="")
            {
                //cout<<"Res is out!!"<<endl;
                try
                {
                    map<string, real>::iterator itIni;
                    map<string, real>::iterator itEnd;
                    map<string, real>::iterator it;
                    map<string,real> the_map=mTrainMatrix->getStringToRealMapping(mTrainMatrix->inputsize());
                    it = the_map.begin();
                    itIni = the_map.begin();
                    itEnd = the_map.begin();
                    for(it++ ; it!=the_map.end(); it++ )
                    {
                        if(it->second<itIni->second)
                            itIni=it;
                        if(it->second>itEnd->second)
                            itEnd=it;
                    }                    
                     if(mDEBUG) cout<<"Init: "<<itIni->first<<"-"<<itIni->second<<endl;
                     if(mDEBUG) cout<<"End: "<<itEnd->first<<"-"<<itEnd->second<<flush<<endl;
                     if(res < itIni->second)
                        OutputVector =  itIni->first;
                     else if(res > itEnd->second)
                        OutputVector =  itEnd->first;
                     else
                        OutputVector = "ERROR: Unexpected error in iterator";
                    if(mDEBUG) cout<<"Predicted: "<<res<<" -> "<<OutputVector<<endl;
                }
                catch(...)
                {
                    OutputVector = "ERROR: Could not resolve mapping";
                }
            }
            if(OutputVector[0]=='D')
                OutputVector=OutputVector.substr(1,OutputVector.length()-1);
        }
        else
        {
            //The output is Numerical (Regressor)
            if (mDEBUG) cout<<"Numerical"<<endl;
            OutputVector = NumToString(res);
        }
    }
    catch(const PLearnError& e)
    {
        delete[] InputValues;
        return string("ERROR: Error while running the PLS algorithm: ") + e.message();
    }
    catch (...)
    {
        delete[] InputValues;
        return "ERROR: Error while running the PLS algorithm";
    }
    delete[] InputValues;

    if(OutputVector.compare("")==0)
    {
        mTrainMatrix->savePMAT("./ERROR.pmat");   
       return string("ERROR: Returned: ") + NumToString(res) +
                string(" Rounded=") + NumToString(round(res)) + 
                " Metadata saved on ./ERROR.pmat";
    }
    else
		//strip and replace all "_nbsp_" back to spaces
		strip(OutputVector);
		replaceAll(OutputVector," ","\t");
        return replaceAll(OutputVector,"_nbsp_"," ");
}

bool PlsAPI::LoadPLSModel(string PLSModelDirPath)
{
    mCPUTrainTime = clock();
    mWallTrainTime = time(NULL);
    //Check if the specified file exists
    if(!IsFile(PLSModelDirPath + string("/Model.pls")))
        return false;
    //Free any previews used memory
    Predictor.forget();
    mTrained=false;
    mOutSymbNum = 0;
    //Load the model from the file
    PLearn::load(PLSModelDirPath + string("/Model.pls") , Predictor);
    if(!IsFile(PLSModelDirPath + string("/HMatrix.pmat")))
        return false;
    mTrainMatrix = getDataSet(PLSModelDirPath + string("/HMatrix.pmat"));
    if(!IsFile(PLSModelDirPath + string("/HMatrix.pmat.metadata/fieldnames")))
        return false;
    mTrainMatrix->setMetaDataDir(PLSModelDirPath + string("/HMatrix.pmat.metadata"));
    mTrainMatrix->loadFieldInfos();
    if(!mTrainMatrix->hasFieldInfos())
        return false;
    mTrainMatrix->loadAllStringMappings();
    //Check if the PLS is now trained
    mTrained=Predictor.IsTrained();
    mCPUTrainTime =((double)clock() - mCPUTrainTime)/(double)CLOCKS_PER_SEC;
    mWallTrainTime =(double)time(NULL) - mWallTrainTime;
    //Set the number of symbols of output variable if it is descrete

    return mTrained;
}

bool PlsAPI::SavePLSModel(string PLSModelDirPath)
{
    //Only save if the PLS is trained!
    if(!mTrained)
        return false;
    try
    {
        //Create a temporary matrix with only the first row. this is done so
        //that all descriptor names and  mapped strings can be saved
        VMat tmpMatrix(mTrainMatrix.subMatRows(0,1));
        //Create a directory to place all the files needed to specify the model
        system((string("mkdir -p ") + PLSModelDirPath).c_str());
        //Save the PLS model
        PLearn::save(PLSModelDirPath + string("/Model.pls") , Predictor);
        //Check if the file was created
        if(!IsFile(PLSModelDirPath + string("/Model.pls")))
            return false;
        //Save the temporary Header Matrix
        tmpMatrix->savePMAT(PLSModelDirPath + string("/HMatrix.pmat"));
        //Check if the file was created
        if(!IsFile(PLSModelDirPath + string("/HMatrix.pmat")))
            return false;
    }
    catch (...)
    {
        return false;
    }

    return true;
}

bool PlsAPI::IsTrained()
{
    return mTrained;
}

string PlsAPI::NumToString(int Number)
{
    stringstream strmN;
    strmN<<Number;
    return string(strmN.str());
}

string PlsAPI::NumToString(real Number)
{
    stringstream strmN;
    strmN<<Number;
    return string(strmN.str());
}

string PlsAPI::GetClassVarName()
{
    if(!mTrained)
        return string("ERROR: No data loaded");
    string ClassName = mTrainMatrix.fieldName(mTrainMatrix->inputsize());
    
    if(ClassName.at(0)=='#')
        return ClassName.substr(1,ClassName.length()-1);
    else
        return ClassName;
    
}
int PlsAPI::GetVarSymbNum(int col)
{
    if (!IsTrained())
        return 0;
    else if( ((mTrainMatrix.fieldName(col)).at(0)!='#'))
        return 0;
    else
        return mTrainMatrix->getStringToRealMapping(col).size();
    
}

int PlsAPI::GetTrainExamplesNumber()
{
    if (!IsTrained())
        return 0;
    else
        return mTrainMatrix.length();
}

int PlsAPI::GetTrainDescriptorsNumber()
{
    if (!IsTrained())
        return 0;
    else
        return mTrainMatrix->inputsize();
}

int PlsAPI::GetTrainOutputsNumber()
{
    if (!IsTrained())
        return 0;
    else
        return Predictor.outputsize();
}

bool PlsAPI::IsFile(string Path)
{
    FILE* pFile;
    //Check if the specified file exists
    if((pFile = fopen(Path.c_str(),"r")) == NULL)
        return false;
    else
        fclose(pFile);
    return true;
}

string PlsAPI::FileOrngToPLearn(string orngFilePath,string PLrnFilePath)
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
    ofstream PLrnFile(PLrnFilePath.c_str());
	
    if (orngFile.is_open())
    {
        if (PLrnFile.is_open())
        {
            //.amat comment line
			PLrnFile<<"#File automaticaly created from native Orange formated file: "<<orngFilePath<<endl;
            while (! orngFile.eof() )
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
							return "ERROR: File must have at least 2 columns";
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
						}
						//if class column not found, choose the last column as class
						if((classCol<0) || (ignoredRows[classCol]=='1'))
							classCol=ignoredRows.rfind('0');
						//header of data and responses size
						PLrnFile<<"#sizes: "<<dataCols-1<<" "<<"1 0 0"<<endl;
						PLrnFile<<"#: ";
						//header of columns names
                                                mVarNames="";
						for(int i=0;i<inFileCols;i++)
                        {
							if(ignoredRows[i]=='0')
							{
								if( ((firstThreeOrngRows[1][i]).compare("c")==0) || ((firstThreeOrngRows[1][i]).compare("continuous")==0) )
									elementName=(firstThreeOrngRows[0][i]);
								else
									elementName="#"+(firstThreeOrngRows[0][i]);
								//elementName+=(firstThreeOrngRows[0][i]);
								for(unsigned int j=0;j<elementName.length();j++)
									elementName[j] = elementName[j]==' '?'_':elementName[j];                                       
								if(i==classCol)
									className+=elementName;
								else
									PLrnFile<<elementName<<" ";	
							}
                                                mVarNames += elementName + string("\t");
						}
						PLrnFile<<className<<endl;
					}
					if(Rows>=3)//Compute Data Rows
					{
						Cols=0;
						classElement="";
						actPos=-1;
						nxPos=line.find(token,0);
						while((unsigned int)nxPos!=string::npos)//read all row elements but the last
						{
                            element =line.substr(actPos+1,nxPos-actPos-1);
							//Strip and replace all spaces with "_nbsp_"
							strip(element);
                            replaceAll(element," ","_nbsp_");
                            element = element==""?"0":element;
							if(Cols==classCol)
								classElement=element;
							else
								if(ignoredRows[Cols]=='0')
								{
								//if( ((*firstThreeOrngRows[1][Cols]).compare("c")!=0) &&
								//((*firstThreeOrngRows[1][Cols]).compare("continuous")!=0) )
								//   element=string("D")+element;
                                PLrnFile<<element<<" ";
								}
                    	    actPos=nxPos;
                    	    nxPos=line.find(token,actPos+1);
                    	    Cols++;
                    	}
                    	if((unsigned int)actPos<(line.length()-1))//read the last column element
                    	{   
                            element =line.substr(actPos+1,nxPos-actPos-1);
							//strip and replace all spaces with "_nbsp_"
							strip(element);
                            replaceAll(element," ","_nbsp_");
                            element = element==""?"0":element;
                            if(Cols==classCol)
                                classElement = element;
			    			else if(ignoredRows[Cols]=='0')
                            {
								//if(((*firstThreeOrngRows[1][Cols]).compare("c")!=0) &&         
								//((*firstThreeOrngRows[1][Cols]).compare("continuous")!=0) )
								//    element=string("D")+element;
                    	        PLrnFile<<element<<" ";                    
                            }
                    	    Cols++;
                    	}
						if( ((firstThreeOrngRows[1][classCol]).compare("c")!=0) &&         
						         ((firstThreeOrngRows[1][classCol]).compare("continuous")!=0) )
					        classElement=string("D")+classElement;
						PLrnFile<<classElement<<endl;
		    		}   
                    Rows++;
                } 
            }
            PLrnFile.close();
        }
        else 
		{
            delete[] firstThreeOrngRows[0];
            delete[] firstThreeOrngRows[1];
            delete[] firstThreeOrngRows[2];
	    	return "ERROR: Unable to open PLearn file for write";
        }
		orngFile.close();
    }
    else 
    {
        delete[] firstThreeOrngRows[0];
        delete[] firstThreeOrngRows[1];
        delete[] firstThreeOrngRows[2];
    	return "ERROR: Unable to open Orange file for read!";
    }
    delete[] firstThreeOrngRows[0];
    delete[] firstThreeOrngRows[1];
    delete[] firstThreeOrngRows[2];    
    return PLrnFilePath;
}

int PlsAPI::Tokenize(const string* InputStr, string** pSplitedStr ,char Token)
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
		Cols++;				//Count the last element
		
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
		if (mVerbose) cout<<"ERROR:Probably bad formated Orange File!"<<endl;
	}
	return Cols;
}

string NumToString(real Number)
{       
    stringstream strmN;
    strmN<<Number;
    return string(strmN.str());
}

