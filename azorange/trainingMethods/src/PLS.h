// -*- C++ -*-
// PLS.h

/*! \file PLS.h */

#ifndef PLS_INC
#define PLS_INC

#include <plearn_learners/generic/PLearner.h>

namespace PLearn {
using namespace std;


class PLS: public PLearner
{

public:

    typedef PLearner inherited;
  
protected:

    // *********************
    // * protected options *
    // *********************

    Mat B;
    int m;
    Vec mean_input;
    Vec mean_target;
    int p;
    Vec stddev_input;
    Vec stddev_target;
    Mat W;

    /// Estimate of the residual variance for each output variable.  Saved as a
    /// learned option to allow outputting confidence intervals when model is
    /// reloaded and used in test mode.
    Vec resid_variance;
    
    //variable used as aux. for time stamps
    time_t rawtime;
  
public:

    // ************************
    // * public build options *
    // ************************

    int k;
    string method;
    real precision;
    bool output_the_score;
    bool output_the_target;
    bool compute_confidence;

    // ****************
    // * Constructors *
    // ****************

    // Default constructor, make sure the implementation in the .cc
    // initializes all fields to reasonable default values.
    PLS();


    // ********************
    // * PLearner methods *
    // ********************

private: 

    //! This does the actual building. 
    void build_();
    bool mVerbose;

protected: 
  
    //! Declares this class' options.
    static void declareOptions(OptionList& ol);

    //! Compute the variance of residuals on the specified dataset
    void computeResidVariance(VMat dataset, Vec& resid_variance);
    
public:

    // ************************
    // **** Object methods ****
    // ************************
	//Added by Pedro
	//Get teh train status of PLS
	bool IsTrained() { if( (W.width()==0) || (stage==0)) return false; else return true; };
    //! Simply calls inherited::build() then build_().
    virtual void build();

    //! Transforms a shallow copy into a deep copy.
    virtual void makeDeepCopyFromShallowCopy(CopiesMap& copies);

    // Declares other standard object methods.
    // If your class is not instantiatable (it has pure virtual methods)
    // you should replace this by PLEARN_DECLARE_ABSTRACT_OBJECT_METHODS.
    PLEARN_DECLARE_OBJECT(PLS);


    // **************************
    // **** PLearner methods ****
    // **************************

    //! Returns the size of this learner's output, (which typically
    //! may depend on its inputsize(), targetsize() and set options).
    virtual int outputsize() const;

    //! (Re-)initializes the PLearner in its fresh state (that state may depend on the 'seed' option)
    //! And sets 'stage' back to 0 (this is the stage of a fresh learner!).
    virtual void forget();
    
    //! The role of the train method is to bring the learner up to stage==nstages,
    //! updating the train_stats collector with training costs measured on-line in the process.
    virtual void train();

////scPA
    // SET VERBOSITY FLAG
    void SetVerbosity(bool verbose=true){mVerbose=verbose; this->verbosity=verbose;}
////ecPA

    //! Computes the output from the input.
    virtual void computeOutput(const Vec& input, Vec& output) const;

    //! Computes the costs from already computed output. 
    virtual void computeCostsFromOutputs(const Vec& input, const Vec& output, 
                                         const Vec& target, Vec& costs) const;

    //! Compute confidence intervals from already-computed outputs.
    virtual bool computeConfidenceFromOutput(const Vec&, const Vec& output, real probability,
                                             TVec< pair<real,real> >& intervals) const;
    
    //! Returns the names of the costs computed by computeCostsFromOutpus (and thus the test method).
    virtual TVec<string> getTestCostNames() const;

    //! Returns the names of the objective costs that the train method computes and 
    //! for which it updates the VecStatsCollector train_stats.
    virtual TVec<string> getTrainCostNames() const;

    // *** SUBCLASS WRITING: ***
    // While in general not necessary, in case of particular needs 
    // (efficiency concerns for ex) you may also want to overload
    // some of the following methods:
    // virtual void computeOutputAndCosts(const Vec& input, const Vec& target, Vec& output, Vec& costs) const;
    // virtual void computeCostsOnly(const Vec& input, const Vec& target, Vec& costs) const;
    // virtual void test(VMat testset, PP<VecStatsCollector> test_stats, VMat testoutputs=0, VMat testcosts=0) const;
    // virtual int nTestCosts() const;
    // virtual int nTrainCosts() const;


    //***** STATIC METHODS *****

    //! Compute the largest eigenvector of m with the NIPALS algorithm:
    //!   (1) v <- random initialization (but normalized)
    //!   (2) v = m.v, normalize v
    //!   (3) if there is a v[i] that has changed by more than 'preicision',
    //!       go to (2), otherwise return v.
    static void NIPALSEigenvector(const Mat& m, Vec& v, real precision);

};

// Declares a few other classes and functions related to this class.
DECLARE_OBJECT_PTR(PLS);
  
} // end of namespace PLearn

#endif


/*
  Local Variables:
  mode:c++
  c-basic-offset:4
  c-file-style:"stroustrup"
  c-file-offsets:((innamespace . 0)(inline-open . 0))
  indent-tabs-mode:nil
  fill-column:79
  End:
*/
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:encoding=utf-8:textwidth=79 :
