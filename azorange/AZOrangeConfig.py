"""
This file should contain global configuration variables for AZOrange.
Some of these variables may have to be changed for a new installation. 
"""
import os

# Data Load Options
#   Unsupported symbols that need conversion 
# NOTE: For testing in python the symbol code, use for ex: unichr(181)
#       Then in CVT_SYM use the chr method.
# 
# unichr(181)- micro sign
#   Sometimes the \xc2 is preceding the synmbol. We just remove it!
CVT_SYM = {"\xc2":"", chr(181):"u"} 

# General environement
AZORANGEHOME = os.environ["AZORANGEHOME"]
SCRATCHDIR = "/tmp"
NFS_SCRATCHDIR = os.path.join(os.environ["HOME"],"AZO_NFS_scratchDir")

if os.path.isfile(os.path.join(AZORANGEHOME,"azorange","AZOrangeExtraConfig.py")):
    from AZOrangeExtraConfig import *

# These strings will automatically be identified as the smiles attribute of the data.
SMILESNAMES = ["SMILES", "Molecule SMILES", "SMILES_", "Compound Structure", "glf_smiles", "MolSmiles", "Smiles", "Molecule Structure", "Structure", "SMILEStoPred", "SMILES_1", "smiles", "SMILEStoPred"]
AZIDNAMES = ["Compound Name", "MolName"]



OPTALGDICT = {"FANN_TRAIN_INCREMENTAL":0, "FANN_TRAIN_BATCH":1, "FANN_TRAIN_RPROP":2, \
              "FANN_TRAIN_QUICKPROP":3}
ANNDEFAULTDICT = {"nHidden":5, "randomWeights":True, "scale":True ,"nEpochs":1000, "optAlg":"FANN_TRAIN_QUICKPROP", "MSE":0.001}

#CvBayes parameters
CVBAYESDEFAULTDICT = {"scale":False}

#CvBoost parameters
CVBOOSTTYPE = { "DISCRETE":0, "REAL":1, "LOGIT":2, "GENTLE":3 }
CVBOOSTSPLITCRIT = { "DEFAULT":0, "GINI":1, "MISCLASS":3, "SQERR":4 }
CVBOOSTDEFAULTDICT = {"boost_type":"GENTLE","weak_count":100,"split_criteria":"SQERR","weight_trim_rate":0.95, "max_depth":1, "use_surrogates":False, "priors":None}

# CvANN specific parameters
# Backprop does not seem to work very well, in general no use trying
CVANNOPTALGDICT = {"FANN_TRAIN_BACKPROP":0, "FANN_TRAIN_RPROP":1} # RPROP = resilient-propagation
CVANNSTOPCRITDICT = {"ITER":1, "EPS":2}
#from opencv for the activationFunction:
#ml.CvANN_MLP.GAUSSIAN = 2    # Not completly supported by the moment
#ml.CvANN_MLP.SIGMOID_SYM = 1
#ml.CvANN_MLP.IDENTITY = 0    # this was tested and 'nan' was always returmed
CVANNDEFAULTDICT = {
        "activationFunction"    : 1    ,
        "sigmoidAlpha"          : 0    , # Alpha for the sigmoid activation function
        "sigmoidBeta"           : 0    , # Beta for the sigmoid activation function
        "nHidden"               : 5    ,     
        "scaleData"             : True , 
        "scaleClass"            : True ,
        "optAlg"                : 1    , 
        "bp_dw_scale"           : 0.1  ,  
        "bp_moment_scale"       : 0.1  ,
        "rp_dw0"                : 0.1  ,  # Has a great impact on accuracy, default seems ok, add to optimizer?
        "rp_dw_plus"            : 1.2  ,  # Has a great impact on accuracy, default seems ok, add to optimizer?
        "rp_dw_minus"           : 0.5  ,  # Has a great impact on accuracy, default seems ok, add to optimizer?
        #params.rp_dw_min = ##default is the minimum float value 
        "rp_dw_max"             : 50   ,
        'stopCrit'              : 1    ,# Exist in opencv only
                                        #     cv.CV_TERMCRIT_EPS = 2
                                        #     cv.CV_TERMCRIT_ITER = 1
        'priors'                : None ,
        'maxIter'               : 1000 ,  
        'eps'                   : 0.001,
        "nDiffIniWeights"       : 0,    #default 10
        "stopUPs"               : 5}    #default 5. If > 0, use early stopping with 5 epochs inbetween evaluations




# Random forest specific parameters
# Claire's recommendations below. Using shallower tree to reduce the default training time. 20 is the OpenCV default. 
# Using also OpenCV default of minSample. 
#TermCrti:  0          :  cv.CV_TERMCRIT_ITER = 1
#           other value:  cv.CV_TERMCRIT_EPS = 2   
#RFDEFAULTDICT = {"maxDepth":"1000", "minSample":"2", "useSurrogates":"false", "getVarVariance":"false", "nActVars":"0",
#                 "nTrees":"50", "forestAcc":"0.1", "termCrit":"0"}
RFDEFAULTDICT = {"maxDepth":"20", "minSample":"5", "useSurrogates":"false", "getVarVariance":"false", "nActVars":"0",
                 "nTrees":"100", "forestAcc":"0.1", "termCrit":"0", "stratify":"false", "priors":None, "useBuiltInMissValHandling":False, "NumThreads":"1"}

##scPA
#R-RF default parameters
RRFDEFAULTDICT = {"nTrees":500, "minSample":1, "maxDepth":"NULL", "nActVars":0}
#nActVars = 0 => nActVars = sqrt(dataAttrs) 
#nActVars = -1 => nActVars (mTry) will be optimized by it's internal algorithm 

# PLS defaul parameters
OPTMETHODDICT = {"kernel":0, "simpls":1, "pls1":2}
PLSDEFAULTDICT = {"precision":"1e-6", "k":"10", "method":"kernel"}



# SVM defaul parameters
KERNELTYPEDICT = {"Linear":0, "Polynomial":1, "RBF":2, "Sigmoid":3}
SVMTYPEDICT = {"C-SVC":0, "nu-SVC":1, "One-Class":2, "epsilon-SVR":3,"nu-SVR":4}
SVMDEFAULTDICT = {
	'kernel_type'  : 2     ,
	'svm_type'     : 0     ,
	'gamma'        : 1     ,
	'C'            : 1     ,
	'p'            : 0.1   ,
	'epsC'         : 0.001 ,
	'epsR'         : 0.001 ,
	'probability'  : 0     ,
	'shrinking'    : 1     ,
	'nu'           : 0.5   ,
	'scaleData'    : 1     ,    

        'useNu'        : 0     ,
	'coef0'        : 0     , 
	'degree'       : 5     }


#CV-SVM defaul parameters
#From opnncv: 
#                      SVM_TYPE
#CvSVM.C_SVC = 100
#CvSVM.NU_SVC = 101
#CvSVM.ONE_CLASS = 102
#ml.CvSVM.EPS_SVR = 103
#ml.CvSVM.NU_SVR = 104
#                      KERNEL_TYPE
#ml.CvSVM.LINEAR = 0
#ml.CvSVM.POLY = 1
#ml.CvSVM.RBF = 2
#ml.CvSVM.SIGMOID = 3
CVSVMKERNELTYPEDICT = {"Linear":0, "Polynomial":1, "RBF":2, "Sigmoid":3}
CVSVMTYPEDICT = {"C-SVC":100, "nu-SVC":101, "One-Class":102, "epsilon-SVR":103,"nu-SVR":104}
CVSVMSTOPCRITDICT = {"ITER":1, "EPS":2}
CVSVMDEFAULTDICT = {
        'kernel_type'  : 2     ,
        'svm_type'     : 100   ,
        'gamma'        : 0.05     ,
        'C'            : 1     ,
        'p'            : 0.1     ,

        'stopCrit'     : 2     , # Exist in opencv only
                                 #     cv.CV_TERMCRIT_EPS = 2
                                 #     cv.CV_TERMCRIT_ITER = 1
        'maxIter'      : 50    , # Exist in opencv only
        'epsC'         : 0.001 ,
        'epsR'         : 0.1  ,

        'nu'           : 0.5   ,
        'scaleData'    : 1     ,
        'scaleClass'   : 1     ,

        'priors'       : None  ,

        'coef0'        : 5     ,
        'degree'       : 5     }


