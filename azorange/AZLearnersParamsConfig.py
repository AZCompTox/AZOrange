"""
This file specifies dicts  of parameters and its possible values for each learner. 
Usage:



       XLearner =                                                               # Name of the Learner Class
                  {str ParameterName :                                          # Name of the parameter recognized by the learner
                              [  str types ParameterType,                       # string indicating the type of the parameter requiered by the Learner. the evaluation of 
                                                                                #this string will be ised to cast the values when assigning to the Learner Ex: "types.IntType"

                                 str valuesRangeType <"interval"|"values">,     # Type of the list of the next parameter. Can be "values" or "interval"
                                                                                #    "Values"   - the valuesRange will be a list of values to be used by the learner
                                                                                #    "interval" - the valuesRange will be a list of 2 values with the limits of the parameter 

                                 str list valuesRange <values|IntervalLimits>,  # String to be evaluated. The evaluations must result in a list! The list must contain the 
                                                                                #values to be used by the optimizer or the limits depending on what choosen on valuesRangeType
                      
                                 list  valuesAlias,                             # a list of strings with the alias for the values specified in valuesRange (can be an empty 
                                                                                #list). This is to be used mainly with the valuesRangeType=values 
                                                                                # To make use of this feature, in order to identify the value for each alias, the values
                                                                                #specified in valuesRange must be separated with ' , ' (i.e., space-comma-space) 

                                 str  DFV,                                      # string with default Fixed Value: Default value when used as fixed parameter

                                 bool DOF,                                      # Default Optimization Flag: Default state of optimization: True: Optimize, False:do Not Optimize (DFV will be used)     

                                 bool ERF,                                      # Editable Range Flag: True: The range is editable, False: The listed values on valuesRange 
                                                                                #are the only accepted for this paramater
                                string RangeToolTip]                             # a string with the explanation to show on a tolltip about the default range defined here in the config file
                    ... (more parameters)
                   }
                   
                     
NOTE: Special KeyWords:
  N_ATTR       Number of attributes in the train set
  N_EX         Number of examples in the train set
     IMPORTANT: The N_EX keyword can and should be used in attributes that have as limit the number of examples, but remember
        that the attributes where this keyword is defined will be set according to the dimensions of the tests done within
        the optimizer, that depends on the user choice. For example, in a learner that one attribute can go from 1 to N_EX,  
        if the dataset used for optimization has 100 examples and the evaluation method choosen is the Classification
        accuracy using 5-fold cross validation (100/5=20  i.e.,  80 ex. to Train and 20 ex. to evaluate the CA), then this parameter 
        can actually go from 1 to 80. The relevance of this parameter is related to the dataset size used in the optimization
        and should be corrected when using the optimized learner with different dataset sizes. Back to the same case, if the 
        resultant optimized parameter of that attribute was 80, that probably means that if we will use the learner with other 
        dataset, this parameter should be set to the number of examples and not to 80!!  
 

To create list of values with decimal or integer intervals use:
   miscUtilities.Range(<FirstValue>,<LastValue>,<Step>)
    Ex:
       miscUtilities.Range(-1,1,0.5)   output: [-1, -0.5, 0.0, 0.5, 1.0] 
       miscUtilities.Range(10,0,-2)    output: [10, 8, 6, 4, 2, 0]

To create list of values of power 2 bases  use:
   miscUtilities.power2Range(<FirstPower>,<LastPower>,<PowerStep>)
    Ex:
       miscUtilities.power2Range(0,10,2)   output: [2^0, 2^2, 2^4, 2^6, 2^8, 2^10] =  [1, 4, 16, 64, 256, 1024]
 
       miscUtilities.power2Range(-3,3,1)   output: [0.125, 0.25, 0.5, 1, 2, 4, 8]
"""

import types
import AZOrangeConfig 

#The version of this file. This version should be change when modifications are made. Also change the requiredParamVer on AZParamOpt,py
version = 12



#CvBoost parameters
"""
AZOrangeConfig.CVBOOSTTYPE = { "DISCRETE":0, "REAL":1, "LOGIT":2, "GENTLE":3 }
AZOrangeConfig.CVBOOSTSPLITCRIT = { "DEFAULT":0, "GINI":1, "MISCLASS":3, "SQERR":4 }
AZOrangeConfig.CVBOOSTDEFAULTDICT = {"boost_type":"DISCRETE","weak_count":100,"split_criteria":"DEFAULT","weight_trim_rate":0.95, "max_depth":1, "use_surrogates":True, "priors":None}
"""
CvBoostLearner = {
             'max_depth':["types.IntType", "values", "[int(round(x)) or 1 for x in miscUtilities.Range(1,20)]",[],AZOrangeConfig.CVBOOSTDEFAULTDICT["max_depth"],True,True,"Integer from 1 to 20.\nMinimum is 1"],\
             'weak_count':["types.IntType", "values", "miscUtilities.Range(1,1000)",[],AZOrangeConfig.CVBOOSTDEFAULTDICT["weak_count"],False,True,"Integer from 1 to 1000"],\
             'weight_trim_rate':["types.FloatType", "interval", "[0 , 1]",[],AZOrangeConfig.CVBOOSTDEFAULTDICT["weight_trim_rate"],False,True,"Continious between"],\
             'boost_type':["types.StringType", "values", "['DISCRETE' , 'REAL' , 'LOGIT' , 'GENTLE']",['DISCRETE','REAL','LOGIT','GENTLE'],AZOrangeConfig.CVBOOSTDEFAULTDICT["boost_type"],False,False,"DISCRETE\nREAL\nLOGIT\nGENTLE"],\
             'split_criteria':["types.StringType", "values", "['DEFAULT' , 'GINI' , 'MISCLASS' , 'SQERR']",['DEFAULT' , 'GINI' , 'MISCLASS' , 'SQERR'],AZOrangeConfig.CVBOOSTDEFAULTDICT["split_criteria"],False,False,"DEFAULT\nGINI\nMISCLASS\nSQERR"],\
             #'priors':["types.StringType", "values", "[None]",[],str(AZOrangeConfig.CVBOOSTDEFAULTDICT["priors"]),False,False,"Do not use for optimization, just to change the default values of priors"],\
             #'useBuiltInMissValHandling':["types.BooleanType", "values", "[False , True]",["No","Yes"],str(AZOrangeConfig.CVBOOSTDEFAULTDICT["useBuiltInMissValHandling"]),False,False,"Do not use for optimization, just for choosing the method to use when there are missing values on datasets."],\
             'use_surrogates':["types.BooleanType", "values", "[False , True]",["No","Yes"],str(AZOrangeConfig.CVBOOSTDEFAULTDICT["use_surrogates"]),False,False,"True or False."]\
            }




# PLS parameters
PLSLearner = {'method' : [ "types.StringType", "values", "['kernel' , 'simpls']",["Kernel","SimPLS"],str(AZOrangeConfig.PLSDEFAULTDICT["method"]),False,False,"Kernel\nsimpls"], \
              'k' : ["types.StringType", "values", "miscUtilities.Range(1,N_ATTR-1)",[],str(AZOrangeConfig.PLSDEFAULTDICT["k"]),True,True,"Integer from 1 to N_ATTR-1"], \
              'precision' : ["types.StringType", "interval", "[1e-7 , 1e-5]",[],str(AZOrangeConfig.PLSDEFAULTDICT["precision"]),False,True,""] \
             }



# SVM parameters
"""
kernel_type:
 0 = "Linear,   x.y"
 1 = "Polynomial,   (g*x.y+c)^d"
 2 = "RBF,   exp(-g*(x-y).(x-y))"
 3 = "Sigmoid,   tanh(g*x.y+c)"

SVM Type:
 0 = C-SVC
 1 = nu-SVC
 2 = one-Class
 3 = epsilon-SVR
 4 =nu-SVR
"""

AZSVMLearner = {'kernel_type' : ["types.IntType", "values", "[0 , 1 , 2 , 3]",['Linear' , 'Polynomial' , 'RBF' , 'Sigmoid'],str(AZOrangeConfig.SVMDEFAULTDICT["kernel_type"]),False,False,"Linear\nPolynomial\nRBF\nSigmoid"],\
                   'svm_type' : ["types.IntType", "values", "[0 , 1 , 3 , 4]",['C-SVC' , 'nu-SVC' , 'epsilon-SVR' , 'nu-SVR'],str(AZOrangeConfig.SVMDEFAULTDICT["svm_type"]),False,False,"C-SVC\nnu-SVC\nepsilon-SVR\nnu-SVR"],\
                   'gamma' : ["types.FloatType", "values", "miscUtilities.power2Range(3,-15,-2)",[],str(AZOrangeConfig.SVMDEFAULTDICT["gamma"]),True,True,"Power of base 2 with power\nfrom -15 to 3 by step of 0.5"],\
                   'C' : ["types.FloatType", "values", "miscUtilities.power2Range(-5,15,2)",[],str(AZOrangeConfig.SVMDEFAULTDICT["C"]),True,True,"Power of base 2 with power\nfrom -5 to 15 by step of 2"],\
                   'eps' : ["types.FloatType", "interval", "[0.001 , 0.01]",[],str(AZOrangeConfig.SVMDEFAULTDICT["epsC"]),False,True,""],\
                   'nu' : ["types.FloatType", "values", "miscUtilities.power2Range(-10,0,0.5)",[],str(AZOrangeConfig.SVMDEFAULTDICT["nu"]),False,True,"Power of base 2 with power\nfrom -10 to 0 by step of 0.5"],\
                   'shrinking' : ["types.IntType", "values", "[0 , 1]",["False","True"],str(AZOrangeConfig.SVMDEFAULTDICT["shrinking"]),False,False,"False and True"],\
                   'probability' : ["types.IntType", "values", "[0 , 1]",["False","True"],str(AZOrangeConfig.SVMDEFAULTDICT["probability"]),False,False,"False and True"],\
                   'scaleData' : ["types.BooleanType", "values", "[0 , 1]",["False","True"],str(AZOrangeConfig.SVMDEFAULTDICT["scaleData"]),False,False,"False and True"],\
                   'p' : ["types.FloatType", "values", "miscUtilities.Range(0,10,0.1)",[],str(AZOrangeConfig.SVMDEFAULTDICT["p"]),False,True,"From 0 to 10 with increments of 0.1"],\
                   'coef0' : ["types.FloatType", "values",  "miscUtilities.Range(0,10,0.5)",[],str(AZOrangeConfig.SVMDEFAULTDICT["coef0"]),False,True,"From 0 to 10 with increments of 0.5"],\
                   'degree' : ["types.FloatType", "values",  "miscUtilities.Range(0,10,0.5)",[],str(AZOrangeConfig.SVMDEFAULTDICT["degree"]),False,True,"From 0 to 10 with increments of 0.5"],\
                  }

#CvBayes parameters
CvBayesLearner = {'scale' : ["types.BooleanType", "values", "[0 , 1]",["False","True"],str(AZOrangeConfig.CVBAYESDEFAULTDICT["scale"]),False,False,"False and True"],\
}

# CvSVM parameters
"""
kernel_type:
 0 = "Linear,   x.y"
 1 = "Polynomial,   (g*x.y+c)^d"
 2 = "RBF,   exp(-g*(x-y).(x-y))"
 3 = "Sigmoid,   tanh(g*x.y+c)"

CvSVM Type:
 100 = C-SVC
 101 = nu-SVC 
 102 = one-Class
 103 = epsilon-SVR 
 104 =nu-SVR
"""
CvSVMLearner = {'kernel_type' : ["types.IntType", "values", "[0 , 1 , 2 , 3]",['Linear' , 'Polynomial' , 'RBF' , 'Sigmoid'],str(AZOrangeConfig.CVSVMDEFAULTDICT["kernel_type"]),False,False,"Linear\nPolynomial\nRBF\nSigmoid"],\
                   'svm_type' : ["types.IntType", "values", "[100 , 101 , 103 , 104]",['C-SVC' , 'nu-SVC' , 'epsilon-SVR' , 'nu-SVR'],str(AZOrangeConfig.CVSVMDEFAULTDICT["svm_type"]),False,False,"C-SVC\nnu-SVC\nepsilon-SVR\nnu-SVR"],\
                   'gamma' : ["types.FloatType", "values", "miscUtilities.power2Range(3,-15,-2)",[],str(AZOrangeConfig.CVSVMDEFAULTDICT["gamma"]),True,True,"Power of base 2 with power\nfrom -15 to 3 by step of 0.5"],\
                   'C' : ["types.FloatType", "values", "miscUtilities.power2Range(-5,15,2)",[],str(AZOrangeConfig.CVSVMDEFAULTDICT["C"]),True,True,"Power of base 2 with power\nfrom -5 to 15 by step of 2"],\
                   'eps' : ["types.FloatType", "interval", "[0.001 , 0.01]",[],str(AZOrangeConfig.CVSVMDEFAULTDICT["epsC"]),False,True,""],\
                   'nu' : ["types.FloatType", "values", "miscUtilities.power2Range(-10,0,0.5)",[],str(AZOrangeConfig.CVSVMDEFAULTDICT["nu"]),False,True,"Power of base 2 with power\nfrom -10 to 0 by step of 0.5"],\
                   'scaleData' : ["types.BooleanType", "values", "[0 , 1]",["False","True"],str(AZOrangeConfig.CVSVMDEFAULTDICT["scaleData"]),False,False,"False and True"],\
                   'scaleClass' : ["types.BooleanType", "values", "[0 , 1]",["False","True"],str(AZOrangeConfig.CVSVMDEFAULTDICT["scaleClass"]),False,False,"False and True"],\
                   'p' : ["types.FloatType", "values", "miscUtilities.Range(0,100,0.1)",[],str(AZOrangeConfig.CVSVMDEFAULTDICT["p"]),False,True,"From 0 to 100 with increments of 0.1"],\
                   'coef0' : ["types.FloatType", "values",  "miscUtilities.Range(0,10,0.5)",[],str(AZOrangeConfig.CVSVMDEFAULTDICT["coef0"]),False,True,"From 0 to 10 with increments of 0.5"],\
                   'degree' : ["types.FloatType", "values",  "miscUtilities.Range(0,10,0.5)",[],str(AZOrangeConfig.CVSVMDEFAULTDICT["degree"]),False,True,"From 0 to 10 with increments of 0.5"],\
                   'priors' : ["types.StringType", "values",  "[None]",[],str(AZOrangeConfig.CVSVMDEFAULTDICT["priors"]),False,False,"Do not use for optimization, just to change the default values of priors"]\
                  }





#RF parameters
"""
maxDepth    Integer from (1/50)N_ATTR to (4/5)N_ATTR with increments of (1/50)N_ATTR
nActVars    Integer from (1/4)*sqrt(N_ATTR) to 1/2*(N_ATTR) with increments of (1/4)*sqrt(N_ATTR)
minSample   Integer from 2 to 50 with increments of 2\nIf N_EX <= 50: Integer from 2 to N_EX/5 with increments of 2
"""
RFLearner = {'maxDepth':["types.StringType", "values", "[int(round(x)) for x in miscUtilities.Range(1.0/50.0*N_ATTR,4.0/5.0*N_ATTR,1.0/50.0*N_ATTR)]",[],AZOrangeConfig.RFDEFAULTDICT["maxDepth"],False,True,"Integer from (1/50)N_ATTR to (4/5)N_ATTR with increments of (1/50)N_ATTR"],\
             'minSample':["types.StringType", "values", "(1.0/5.0)*N_EX<3 and [2] or miscUtilities.Range(2,N_EX >= 50 and 50 or N_EX/5.0,2)",[],AZOrangeConfig.RFDEFAULTDICT["minSample"],False,True,"Integer from 2 to 50 with increments of 2\nIf N_EX <= 50: Integer from 2 to N_EX/5 with increments of 2"],\
             'nActVars':["types.StringType", "values", "[ int(round(x)) for x in miscUtilities.Range(1.0/4.0*sqrt(N_ATTR),(1/2.0)*(N_ATTR),1.0/4.0*sqrt(N_ATTR))]",[],AZOrangeConfig.RFDEFAULTDICT["nActVars"],True,True,"Integer from (1/4)*sqrt(N_ATTR) to 1/2*(N_ATTR) with increments of (1/4)*sqrt(N_ATTR)"],\
             'nTrees':["types.StringType", "values", "miscUtilities.Range(10,1000,10)",[],AZOrangeConfig.RFDEFAULTDICT["nTrees"],False,True,"Integer from 10 to 1000 with increments of 10"],\
             'forestAcc':["types.StringType", "interval", "[0.001 , 0.9]",[],AZOrangeConfig.RFDEFAULTDICT["forestAcc"],False,True,"Interval from 0.001 to 0.9"],\
             'termCrit':["types.StringType", "values", "[0 , 1]",["Number of trees", "Out of bag error"],AZOrangeConfig.RFDEFAULTDICT["termCrit"],False,False,"Number of trees in the forest (0) \nOut of bag error to achive before termination (1)"],\
             'stratify':["types.StringType", "values", "['false' , 'true']",["No","Yes"],AZOrangeConfig.RFDEFAULTDICT["stratify"],False,False,"No and Yes"],\
             'priors':["types.StringType", "values", "[None]",[],str(AZOrangeConfig.RFDEFAULTDICT["priors"]),False,False,"Do not use for optimization, just to change the default values of priors"],\
             'useBuiltInMissValHandling':["types.BooleanType", "values", "[False , True]",["No","Yes"],str(AZOrangeConfig.RFDEFAULTDICT["useBuiltInMissValHandling"]),False,False,"Do not use for optimization, just for choosing the method to use when there are missing values on datasets."],\
             'NumThreads':["types.StringType", "interval", "[0 , 999999]",[],str(AZOrangeConfig.RFDEFAULTDICT["NumThreads"]),False,False,"Do not use for optimization, just for choosing the number of threads to be used by openMP."]\
            }


#RRF parameters
"""
maxDepth    Integer from (1/50)N_ATTR to (4/5)N_ATTR with increments of (1/50)N_ATTR
nActVars    Integer from (1/4)*sqrt(N_ATTR) to 1/2*(N_ATTR) with increments of (1/4)*sqrt(N_ATTR)
minSample   Integer from 2 to 50 with increments of 2\nIf N_EX <= 50: Integer from 2 to N_EX/5 with increments of 2
"""
RRFLearner = {'maxDepth':["types.StringType", "values", "[int(round(x)) for x in miscUtilities.Range(1.0/50.0*N_ATTR,4.0/5.0*N_ATTR,1.0/50.0*N_ATTR)]",[],AZOrangeConfig.RFDEFAULTDICT["maxDepth"],False,True,"Integer from (1/50)N_ATTR to (4/5)N_ATTR with increments of (1/50)N_ATTR"],\
             'minSample':["types.IntType", "values", "(1.0/5.0)*N_EX<3 and [2] or miscUtilities.Range(2,N_EX >= 50 and 50 or N_EX/5.0,2)",[],AZOrangeConfig.RFDEFAULTDICT["minSample"],False,True,"Integer from 2 to 50 with increments of 2\nIf N_EX <= 50: Integer from 2 to N_EX/5 with increments of 2"],\
             'nActVars':["types.IntType", "values", "[ int(round(x)) for x in miscUtilities.Range(1.0/4.0*sqrt(N_ATTR),(1/2.0)*(N_ATTR),1.0/4.0*sqrt(N_ATTR))]",[],AZOrangeConfig.RRFDEFAULTDICT["nActVars"],True,True,"Integer from (1/4)*sqrt(N_ATTR) to 1/2*(N_ATTR) with increments of (1/4)*sqrt(N_ATTR)"],\
             'nTrees':["types.IntType", "values", "miscUtilities.Range(10,1000,10)",[],AZOrangeConfig.RRFDEFAULTDICT["nTrees"],False,True,"Integer from 10 to 1000 with increments of 10"],\
            }

#SignSVM parameters
SignSVMLearner = {
        'startHeight':["types.IntType", "values", "[0]",[],str(AZOrangeConfig.SIGNSVMDEFAULTDICT["startHeight"]),False,True,"The start Height"],\
        'endHeight':["types.IntType", "values", "[3]",[],str(AZOrangeConfig.SIGNSVMDEFAULTDICT["endHeight"]),False,True,"The end Height"],\
        'nrFolds':["types.IntType", "values", "[2, 7, 10]",[],str(AZOrangeConfig.SIGNSVMDEFAULTDICT["nrFolds"]),False,True,"Number of Folds for parameters selection.\n-1 disables optimization"],\
        'nrThreads':["types.IntType", "values", "[2,8]",[],str(AZOrangeConfig.SIGNSVMDEFAULTDICT["nrThreads"]),False,True,"Number of nrThreads to ne\nused in parameters selection"],\
        'priors' : ["types.StringType", "values",  "[None]",[],str(AZOrangeConfig.CVANNDEFAULTDICT["priors"]),False,False,"Do not use for optimization, just to change the default values of priors"],\
        'scaleClass':["types.BooleanType", "values", "[False , True]",[],str(AZOrangeConfig.CVANNDEFAULTDICT["scaleClass"]),False,False,"No and Yes"]
 }

#CvANN parameters
CvANNLearner = {'nHidden':["[types.IntType]", "values", "(1.0/30.0)*N_EX<3 and [2] or [int(round(x)) for x in miscUtilities.Range(2,((1/30.0)*N_EX))]",[],str(AZOrangeConfig.CVANNDEFAULTDICT["nHidden"]),True,True,"2-N_EX/30"],\
              'maxIter':["types.IntType", "values", "[int(x) for x in miscUtilities.Range(100,10000)]",[],str(AZOrangeConfig.CVANNDEFAULTDICT["maxIter"]),False,True,"100-10000"],\
              'optAlg':["types.IntType", "values", "[0 , 1]",["BACKProp","RProp"],str(AZOrangeConfig.CVANNDEFAULTDICT["optAlg"]),False,False,"Incremental\BACKProp\nRProp"],\
              'stopCrit':["types.IntType", "values", "[1 , 2]",["NIter","EPS"],str(AZOrangeConfig.CVANNDEFAULTDICT["stopCrit"]),False,False,"Number of Iterations\nEpsilon"],\
              'eps':["types.FloatType", "interval", "[0.001 , 0.1]",[],str(AZOrangeConfig.CVANNDEFAULTDICT["eps"]),False,True,"0.001-0.1"],\
              'sigmoidAlpha':["types.FloatType", "interval", "[0.0, 1.0]",[],str(AZOrangeConfig.CVANNDEFAULTDICT["sigmoidAlpha"]),False,True,"0.0-1.0"],\
              'sigmoidBeta':["types.FloatType", "interval", "[0.0, 1.0]",[],str(AZOrangeConfig.CVANNDEFAULTDICT["sigmoidBeta"]),False,True,"0.0-1.0"],\
              'rp_dw_plus':["types.FloatType", "interval", "[1.01, 3.0]",[],str(AZOrangeConfig.CVANNDEFAULTDICT["rp_dw_plus"]),False,True,"1.0-3.0"],\
              'rp_dw_minus':["types.FloatType", "interval", "[0.01, 0.99]",[],str(AZOrangeConfig.CVANNDEFAULTDICT["rp_dw_minus"]),False,True,"0.0-1.0"],\
              'rp_dw_max':["types.FloatType", "interval", "[1, 100]",[],str(AZOrangeConfig.CVANNDEFAULTDICT["rp_dw_max"]),False,True,"1-100"],\
              'rp_dw0':["types.FloatType", "interval", "[0.05, 1.0]",[],str(AZOrangeConfig.CVANNDEFAULTDICT["rp_dw0"]),False,True,"0.05-1.0"],\
              'scaleData':["types.BooleanType", "values", "[False , True]",[],str(AZOrangeConfig.CVANNDEFAULTDICT["scaleData"]),False,False,"No and Yes"],\
              'scaleClass':["types.BooleanType", "values", "[False , True]",[],str(AZOrangeConfig.CVANNDEFAULTDICT["scaleClass"]),False,False,"No and Yes"],\
              'priors' : ["types.StringType", "values",  "[None]",[],str(AZOrangeConfig.CVANNDEFAULTDICT["priors"]),False,False,"Do not use for optimization, just to change the default values of priors"]\
             }





# Definition of class for interface with parameters
from copy import deepcopy

class API(object):
    paramFields = {"type":0,            # a STRING with the type required by the Learner ex: "types.IntType"
                                        #    READONLY

                   "rangeType":1,       # a STRING indicating the range type: "values"-(disccrete values)       OR
                                        #    "interval"-(continuous interval between 2 numbers) 
                                        # OBS - interval can only be used with continious attributes

                   "range":2,           # a LIST with the range: "[1 , 2 , 3]" for values                     OR    
                                        #    "[0 , 2.5]" for intervals     

                   "alias":3,           # a LIST of alias for the respective values in 'range' Ex: ["No","Yes"]

                   "default":4,         # the VALUE to use  when NOT optimising. it must be a valid value for the learner

                   "optimize":5,        # the FLAG indicating either the parameter is to be optimized or not i.e.  True or Flase 

                   "editableRange":6,   # a FLAG indicating if the values specified in the 'range' can be edited. 
                                        #    If False, this parameter can only be one of the values defined in the 'range'
                                        #    READ ONLY

                   "description":7}     # a STRING with a description of the 'range'
 
    def __init__(self, learner, **kwds):
        self.verbose = 0                        # Verbose Flag
        self.__dict__.update(kwds)              # Update the input keywords
        self.parametersDict = None              # The user changed parameters in the correct configuration format
        self.originalParametersDict = None      # The original changed parameters in the correct configuration format
        self.parameters = None                  # A list with the name of the available parameters for this learner
                               #paramFields     # A list with the available field names 
        if (not learner) or (type(learner) not in types.StringTypes):
            return
        self.learner = learner
        try:
            self.originalParametersDict = eval(self.learner)
            self.parametersDict = deepcopy(self.originalParametersDict)
        except:
            if self.verbose > 0: print "Parameter '"+self.learner+"' not defined in the optimizer configuration file"
            return
        self.parameters = [param for param in self.parametersDict]
    def setOptimizeAllParameters(self,state):
        """setOptimizeAllParameters(state)   -  Set all parameters 'optimize' flag to 'state'
           state must be True or False
        """
        if state: 
            state = True
        else:
            state = False
        for par in self.parameters:
            self.setParameter(par,"optimize",state)


    def getParameterNames(self):
        """Returns a list with the parameter names of the current learner"""
        return deepcopy(self.parameters)

    def getParameterFields(self):
        """Returns a list with the available parameter fields """
        return [field for field in self.paramFields]

    def getParametersDict(self):
        """Returns the dictionary of the user changed parameters according to the correct configuration format"""
        return deepcopy(self.parametersDict)

    def getOriginalParametersDict(self):
        """Returns the dictionary of the original parameters according to the correct configuration format"""
        return deepcopy(self.originalParametersDict)

    def setParameter(self,name,field,value):
        """setParameter(name,field,value)  ->  Set the 'field' of parameter 'name' with 'value'
           Returns True if the assigned or False if not"""

        if name not in self.parameters:
            if self.verbose > 0: print "Parameter '"+name+"' is not valid for learner  '"+self.learner
            return False
        if field not in self.paramFields:
            if self.verbose > 0: print "Field '"+field+"' is not valid"
            return False
        if field == "default":
            realValue = value 
            try:
                if type(realValue) in types.StringTypes:
                    if realValue in self.originalParametersDict[name][self.paramFields["alias"]]:
                        idx = self.originalParametersDict[name][self.paramFields["alias"]].index(realValue)
                        realValue = eval(self.originalParametersDict[name][self.paramFields["range"]])[idx]
            except:
                if self.verbose > 0: print "Field '"+field+"' was not able to be set up via alias!"
                realValue = value
            try:
                if type(eval(self.parametersDict[name][self.paramFields["type"]])) == types.ListType:
                    self.parametersDict[name][self.paramFields[field]] = eval(self.parametersDict[name][self.paramFields["type"][0]])(realValue)
                else:
                    self.parametersDict[name][self.paramFields[field]] = eval(self.parametersDict[name][self.paramFields["type"]])(realValue)
                # special cases where the parameter is defined as string, it cannot have spaces
                try:
                    if eval(self.parametersDict[name][self.paramFields["type"]]) == types.StringType:
                        self.parametersDict[name][self.paramFields[field]] = self.parametersDict[name][self.paramFields[field]].replace(" ","")
                except:
                    pass
            except:
                if self.verbose > 0: print "Field '"+field+"' was not able to be set up"
                return False
        elif (field == "editableRange"):
            if self.verbose > 0: print "The field '"+field+"' is read only"
            return False
        elif (field == "optimize"):
            if type(value) == types.BooleanType:
                self.parametersDict[name][self.paramFields[field]] = value
            elif type(value) == types.IntType and (value == 0 or value == 1):
                self.parametersDict[name][self.paramFields[field]] = types.BooleanType(value)
            else:
                if self.verbose > 0: print "Field '"+field+"' was not able to be set up. Wrong type"
                return False
        elif (field == "description"):
            self.parametersDict[name][self.paramFields[field]] = types.StringType(value)     
        elif (field == "alias"):
            if type(value) == types.ListType:
                self.parametersDict[name][self.paramFields[field]] = [str(alias) for alias in values]
            else:
                if self.verbose > 0: print "Field '"+field+"' was not able to be set up. Wrong type"
                return False
        elif (field == "type"):
            if self.verbose > 0: print "The field '"+field+"' is read only"
            return False 
        elif (field == "rangeType"):
            if type(value) in types.StringTypes:
                if value.lower() == "values" or value.lower() == "interval":
                    self.parametersDict[name][self.paramFields[field]] = value.lower()
                else:
                    if self.verbose > 0: print "Field '"+field+"' was not able to be set up. Wrong keyword. Use 'values' or 'interval'"
                    return False
            else:
                if self.verbose > 0: print "Field '"+field+"' was not able to be set up. Wrong type"
                return False
        else: # field == "range"
            if type(value) == types.ListType:
                self.parametersDict[name][self.paramFields[field]] = str(value).replace(", "," , ")
            elif type(value) in types.StringTypes:
                self.parametersDict[name][self.paramFields[field]] = value.strip()
            else:
                if self.verbose > 0: print "Field '"+field+"' was not able to be set up. Wrong type"
                return False
        return True


    def getParameter(self,name,field):
        """getParameter(name,field)  ->  Get the value of the 'field' of the parameter 'name'
           Returned the value of the parameter or None if there is no such 'name' or 'field'"""
        if name not in self.parameters:
            if self.verbose > 0: print "Parameter '"+name+"' is not valid for learner  '"+self.learner
            return None
        if field not in self.paramFields:
            if self.verbose > 0: print "Field '"+field+"' is not valid"
            return None
        return self.parametersDict[name][self.paramFields[field]]


    def getOriginalParameter(self,name,field):
        """getOriginalParameter(name,field)  ->  Get the value of the 'field' of the parameter 'name'
           Returned the value of the parameter or None if there is no such 'name' or 'field'"""
        if name not in self.originalParametersDict:
            if self.verbose > 0: print "Parameter '"+name+"' is not valid for learner  '"+self.learner
            return None
        if field not in self.paramFields:
            if self.verbose > 0: print "Field '"+field+"' is not valid"
            return None
        return self.originalParametersDict[name][self.paramFields[field]]


