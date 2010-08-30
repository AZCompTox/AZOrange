// -*- C++ -*-
/*! \file PLS.cc */

#define PL_LOG_MODULE_NAME "PLS"

// From PLearn
#include "PLS.h"
#include <plearn/io/pl_log.h>
#include <plearn/math/TMat_maths.h>    //!< For dist.
#include <plearn/math/pl_erf.h>
#include <plearn/math/plapack.h>
#include <plearn/vmat/ShiftAndRescaleVMatrix.h>
#include <plearn/vmat/SubVMatrix.h>
#include <plearn/vmat/VMat_linalg.h>
//#include <plearn/db/getDataSet.h>
//#include <time.h>
//#include <plearn/io/load_and_save.h>
//#include <iostream>
//#include <fstream>

namespace PLearn {
	using namespace std;
	PLS::PLS() 
    : m(-1),
	p(-1),
	k(1),
	method("kernel"),
	precision(1e-6),
	output_the_score(false),
	output_the_target(true),
	compute_confidence(false),
        mVerbose(false)
{}
	
	PLEARN_IMPLEMENT_OBJECT(
							PLS,
							"Partial Least Squares Regression (PLSR).",
							"You can use this learner to perform regression, and / or dimensionality\n"
							"reduction.\n"
							"PLS regression assumes the target Y and the data X are linked through:\n"
							" Y = T.Q' + E\n"
							" X = T.P' + F\n"
							"The underlying coefficients T (the 'scores') and the loading matrices\n"
							"Q and P are seeked. It is then possible to compute the prediction y for\n"
							"a new input x, as well as its score vector t (its representation in\n"
							"lower-dimensional coordinates).\n"
							"The available algorithms to perform PLS (chosen by the 'method' option) are:\n"
							"\n"
							" ====  PLS1  ====\n"
							"The classical PLS algorithm, suitable only for a 1-dimensional target. The\n"
							"following algorithm is taken from 'Factor Analysis in Chemistry', with an\n"
							"additional loop that (I believe) was missing:\n"
							" (1) Let X (n x p) = the centered and normalized input data\n"
							"     Let y (n x 1) = the centered and normalized target data\n"
							"     Let k be the number of components extracted\n"
							" (2) s = y\n"
							" (3) lx' = s' X, s = X lx (normalized)\n"
							" (4) If s has changed by more than 'precision', loop to (3)\n"
							" (5) ly = s' y\n"
							" (6) lx' = s' X\n"
							" (7) Store s, lx and ly in the columns of respectively T, P and Q\n"
							" (8) X = X - s lx', y = y - s ly, loop to (2) k times\n"
							" (9) Set W = (T P')^(+) T, where the ^(+) is the right pseudoinverse\n"
							"\n"
							" ==== Kernel ====\n"
							"The code implements a NIPALS-PLS-like algorithm, which is a so-called\n"
							"'kernel' algorithm (faster than more classical implementations).\n"
							"The algorithm, inspired from 'Factor Analysis in Chemistry' and above all\n"
							"www.statsoftinc.com/textbook/stpls.html, is the following:\n"
							" (1) Let X (n x p) = the centered and normalized input data\n"
							"     Let Y (n x m) = the centered and normalized target data\n"
							"     Let k be the number of components extracted\n"
							" (2) Initialize A_0 = X'Y, M_0 = X'X, C_0 = Identity(p), and h = 0\n"
							" (3) q_h = largest eigenvector of B_h = A_h' A_h, found by the NIPALS method:\n"
							"       (3.a) q_h = a (normalized) randomn column of B_h\n"
							"       (3.b) q_h = B_h q_h\n"
							"       (3.c) normalize q_h\n"
							"       (3.d) if q_h has changed by more than 'precision', go to (b)\n"
							" (4) w_h = C_h A_h q_h, normalize w_h and store it in a column of W (p x k)\n"
							" (5) p_h = M_h w_h, c_h = w_h' p_h, p_h = p_h / c_h and store it in a column\n"
							"     of P (p x k)\n"
							" (6) q_h = A_h' w_h / c_h, and store it in a column of Q (m x k)\n"
							" (7) A_h+1 = A_h - c_h p_h q_h'\n"
							"     M_h+1 = M_h - c_h p_h p_h',\n"
							"     C_h+1 = C_h - w_h p_h\n"
							" (8) h = h+1, and if h < k, go to (3)\n"
							"\n"
							"The result is then given by:\n"
							" - Y = X B, with B (p x m) = W Q'\n"
							" - T = X W, where T is the score (reduced coordinates)\n"
							"\n"
							"You can choose to have the score (T) and / or the target (Y) in the output\n"
							"of the learner (default is target only, i.e. regression)."
							);
	
	////////////////////
	// declareOptions //
	////////////////////
	void PLS::declareOptions(OptionList& ol)
	{
		// Build options.
		
		declareOption(ol, "k", &PLS::k, OptionBase::buildoption,
					  "The number of components (factors) computed.");
		
		declareOption(ol, "method", &PLS::method, OptionBase::buildoption,
					  "The PLS algorithm used ('pls1' or 'kernel', see help for more details).\n");
		
		declareOption(ol, "output_the_score", &PLS::output_the_score, OptionBase::buildoption,
					  "If set to 1, then the score (the low-dimensional representation of the input)\n"
					  "will be included in the output (before the target).");
		
		declareOption(ol, "output_the_target", &PLS::output_the_target, OptionBase::buildoption,
					  "If set to 1, then (the prediction of) the target will be included in the\n"
					  "output (after the score).");
		
		declareOption(ol, "compute_confidence", &PLS::compute_confidence,
					  OptionBase::buildoption,
					  "If set to 1, the variance of the residuals on the training set is\n"
					  "computed after training in order to allow the computation of confidence\n"
					  "intervals.  In the current implementation, this entails performing another\n"
					  "traversal of the training set.");
		
		// Learnt options.
		
		declareOption(ol, "B", &PLS::B, OptionBase::learntoption,
					  "The regression matrix in Y = X.B + E.");
		
		declareOption(ol, "m", &PLS::m, OptionBase::learntoption,
					  "Used to store the target size.");
		
		declareOption(ol, "mean_input", &PLS::mean_input, OptionBase::learntoption,
					  "The mean of the input data X.");
		
		declareOption(ol, "mean_target", &PLS::mean_target, OptionBase::learntoption,
					  "The mean of the target data Y.");
		
		declareOption(ol, "p", &PLS::p, OptionBase::learntoption,
					  "Used to store the input size.");
		
		declareOption(ol, "precision", &PLS::precision, OptionBase::buildoption,
					  "The precision to which we compute the eigenvectors.");
		
		declareOption(ol, "stddev_input", &PLS::stddev_input, OptionBase::learntoption,
					  "The standard deviation of the input data X.");
		
		declareOption(ol, "stddev_target", &PLS::stddev_target, OptionBase::learntoption,
					  "The standard deviation of the target data Y.");
		
		declareOption(ol, "W", &PLS::W, OptionBase::learntoption,
					  "The regression matrix in T = X.W.");
		
		declareOption(ol, "resid_variance", &PLS::resid_variance, OptionBase::learntoption,
					  "Estimate of the residual variance for each output variable.  Saved as a\n"
					  "learned option to allow outputting confidence intervals when model is\n"
					  "reloaded and used in test mode.  These are saved only if the option\n"
					  "'compute_confidence' is true at train-time.");
		
		// Now call the parent class' declareOptions
		inherited::declareOptions(ol);
	}
	
	///////////
	// build //
	///////////
	void PLS::build()
	{
		inherited::build();
		build_();
	}
	
	////////////
	// build_ //
	////////////
	void PLS::build_()
	{
		if (train_set) {
			this->m = train_set->targetsize();
			this->p = train_set->inputsize();
			mean_input.resize(p);
			stddev_input.resize(p);
			mean_target.resize(m);
			stddev_target.resize(m);
			if (train_set->weightsize() > 0) {
				PLWARNING("In PLS::build_ - The train set has weights, but the optimization algorithm won't use them");
			}
			// Check method consistency.
			if (method == "pls1") {
				// Make sure the target is 1-dimensional.
				if (m != 1) {
					PLERROR("In PLS::build_ - With the 'pls1' method, target should be 1-dimensional");
				}
			} else if (method == "kernel") {
				// Everything should be ok.
			} else if (method == "simpls"){
                                //New Algorithm added by Pedro
                        } else{
				PLERROR("In PLS::build_ - Unknown value for option 'method'");
			}
		}
		if (!output_the_score && !output_the_target) {
			// Weird, we don't want any output ??
			PLWARNING("In PLS::build_ - There will be no output");
		}
	}
	
	/////////////////////////////
	// computeCostsFromOutputs //
	/////////////////////////////
	void PLS::computeCostsFromOutputs(const Vec& input, const Vec& output, 
									  const Vec& target, Vec& costs) const
	{
		// No cost computed.
	}
	
	///////////////////
	// computeOutput //
	///////////////////
	void PLS::computeOutput(const Vec& input, Vec& output) const
	{
		static Vec input_copy;
		if (W.width()==0)
			PLERROR("PLS::computeOutput but model was not trained!");
		// Compute the output from the input
		int nout = outputsize();
		output.resize(nout);
		// First normalize the input.
		input_copy.resize(this->p);
		input_copy << input;
		input_copy -= mean_input;
		input_copy /= stddev_input;
		int target_start = 0;
		if (output_the_score) {
			transposeProduct(output.subVec(0, this->k), W, input_copy);
			target_start = this->k;
		}
		if (output_the_target) {
			if (this->m > 0) {
				Vec target = output.subVec(target_start, this->m);
				transposeProduct(target, B, input_copy);
				target *= stddev_target;
				target += mean_target;
			} else {
				// This is just a safety check, since it should never happen.
				PLWARNING("In PLS::computeOutput - You ask to output the target but the target size is <= 0");
			}
		}
	}
	
	
	/////////////////////////////////
	// computeConfidenceFromOutput //
	/////////////////////////////////
	
	bool PLS::computeConfidenceFromOutput(const Vec&, const Vec& output, real probability,
										  TVec< pair<real,real> >& intervals) const
	{
		// Must figure out where the real output starts within the output vector
		if (! output_the_target)
			PLERROR("PLS::computeConfidenceFromOutput: the option 'output_the_target' "
					"must be enabled in order to compute confidence intervals");
		int ostart = (output_the_score? k : 0);
		Vec regr_output = output.subVec(ostart, m);
		
		if (m != resid_variance.size())
			PLERROR("PLS::computeConfidenceFromOutput: residual variance not yet computed "
					"or its size (= %d) does not match the output size (= %d)",
					resid_variance.size(), m);
		
		// two-tailed
		const real multiplier = gauss_01_quantile((1+probability)/2);
		intervals.resize(m);
		for (int i=0; i<m; ++i) {
			real half_width = multiplier * sqrt(resid_variance[i]);
			intervals[i] = std::make_pair(output[i] - half_width,
										  output[i] + half_width);
		}
		return true;
	}
	
	
	////////////
	// forget //
	////////////
	void PLS::forget()
	{
		stage = 0;
		// Free memory.
		B = Mat();
		W = Mat();
	}
	
	//////////////////////
	// getTestCostNames //
	//////////////////////
	TVec<string> PLS::getTestCostNames() const
	{
		// No cost computed.
		TVec<string> t;
		return t;
	}
	
	///////////////////////
	// getTrainCostNames //
	///////////////////////
	TVec<string> PLS::getTrainCostNames() const
	{
		// No cost computed.
		TVec<string> t;
		return t;
	}
	
	/////////////////////////////////
	// makeDeepCopyFromShallowCopy //
	/////////////////////////////////
	void PLS::makeDeepCopyFromShallowCopy(CopiesMap& copies)
	{
		inherited::makeDeepCopyFromShallowCopy(copies);
		
		// ### Call deepCopyField on all "pointer-like" fields 
		// ### that you wish to be deepCopied rather than 
		// ### shallow-copied.
		// ### ex:
		deepCopyField(B, copies);
		deepCopyField(mean_input, copies);
		deepCopyField(mean_target, copies);
		deepCopyField(stddev_input, copies);
		deepCopyField(stddev_target, copies);
		deepCopyField(W, copies);
		deepCopyField(resid_variance, copies);
	}
	
	///////////////////////
	// NIPALSEigenvector //
	///////////////////////
	void PLS::NIPALSEigenvector(const Mat& m, Vec& v, real precision) {
		int n = v.length();
		Vec w(n);
		v << m.column(0);
		normalize(v, 2.0);
		bool ok = false;
		while (!ok) {
			w << v;
			product(v, m, w);
			normalize(v, 2.0);
			ok = true;
			for (int i = 0; i < n && ok; i++) {
				if (fabs(v[i] - w[i]) > precision) {
					ok = false;
				}
			}
		}
	}
	
	////////////////
	// outputsize //
	////////////////
	int PLS::outputsize() const
	{
		int os = 0;
		if (output_the_score) {
			os += this->k;
		}
		if (output_the_target && m >= 0) {
			// If m < 0, this means we don't know yet the target size, thus we
			// shouldn't report it here.
			os += this->m;
		}
		return os;
	}
	
	///////////
	// train //
	///////////
	void PLS::train()
	{
                report_progress=mVerbose;
		if (stage == 1) {
			// Already trained.
			MODULE_LOG << "Skipping PLS training" << endl;
			return;
		}
		MODULE_LOG << "PLS training started" << endl;
		
		// Construct the centered and normalized training set, for the input
		// as well as the target part.
		DBG_MODULE_LOG << "Normalizing of the data" << endl;
		VMat input_part = new SubVMatrix(train_set,
										 0, 0,
										 train_set->length(),
										 train_set->inputsize());
                VMat target_part = new SubVMatrix( train_set,
										   0, train_set->inputsize(),
										   train_set->length(),
										   train_set->targetsize());
                PP<ShiftAndRescaleVMatrix> X_vmat =
			new ShiftAndRescaleVMatrix(input_part, true);
                //Force the algorithm to accept the matrix as it is making NO
                //  scale (change also ShiftAndRescaleVMatrix(input_part, true)
                //  to FALSE)
                //X_vmat->setOption("no_scale","1");
                //X_vmat->build();
                X_vmat->verbosity = this->verbosity;
		mean_input << X_vmat->shift;
		stddev_input << X_vmat->scale;
		negateElements(mean_input);
		invertElements(stddev_input);
		PP<ShiftAndRescaleVMatrix> Y_vmat =
			new ShiftAndRescaleVMatrix(target_part, target_part->width(), true);
                Y_vmat->verbosity = this->verbosity;
		mean_target << Y_vmat->shift;
		stddev_target << Y_vmat->scale;
		negateElements(mean_target);
		invertElements(stddev_target);
		// Some common initialization.
		W.resize(p, k);
		Mat P(p, k);
		Mat Q(m, k);
		int n = X_vmat->length();
		VMat X_vmatrix = static_cast<ShiftAndRescaleVMatrix*>(X_vmat);
		VMat Y_vmatrix = static_cast<ShiftAndRescaleVMatrix*>(Y_vmat);
		if (method == "kernel") {
			// Initialize the various coefficients.
			DBG_MODULE_LOG << "Initialization of the coefficients" << endl;
			Vec ph(p);
			Vec qh(m);
			Vec wh(p);
			Vec tmp(p);
			real ch;
                        if (mVerbose) {time ( &rawtime ); cout<<ctime(&rawtime)<<"Computing A0=X'Y ... "<<flush;}
			Mat Ah = transposeProduct(X_vmatrix, Y_vmatrix);
                        if (mVerbose) {time ( &rawtime ); cout<<ctime(&rawtime)<<"done"<<endl<<"Computing M0=X'X ... "<<flush;}
			Mat Mh = transposeProduct(X_vmatrix, X_vmatrix);
                        if (mVerbose) {time ( &rawtime ); cout<<ctime(&rawtime)<<"done"<<endl;}
			Mat Ch(p,p);    // Initialized to Identity(p).
			Mat Ah_t_Ah;
			Mat update_Ah(p,m);
                        Mat update_Mh(p,p);
			Mat update_Ch(p,p);
			for (int i = 0; i < p; i++) {
				for (int j = i+1; j < p; j++) {
					Ch(i,j) = Ch(j,i) = 0;
				}
				Ch(i,i) = 1;
			}
			// Iterate k times to find the k first factors.
			PP<ProgressBar> pb(
							   report_progress? new ProgressBar("Computing the PLS-kernel components", k)
											  : 0);
			
			for (int h = 0; h < this->k; h++) {
                                if (mVerbose) {time ( &rawtime ); cout<<ctime(&rawtime)<<" Computing component "<<h<<"..."<<flush;}
				Ah_t_Ah = transposeProduct(Ah,Ah);
				if (m == 1) {
					// No need to compute the eigenvector.
					qh[0] = 1;
				} else {
					NIPALSEigenvector(Ah_t_Ah, qh, precision);
				}
				product(tmp, Ah, qh);
				product(wh, Ch, tmp);
				normalize(wh, 2.0);
				W.column(h) << wh;
				product(ph, Mh, wh);
				ch = dot(wh, ph);
				ph /= ch;
				P.column(h) << ph;
				transposeProduct(qh, Ah, wh);
				qh /= ch;
				Q.column(h) << qh;
				Mat ph_mat(p, 1, ph);
				Mat qh_mat(m, 1, qh);
				Mat wh_mat(p, 1, wh);
				update_Ah = productTranspose(ph_mat, qh_mat);
				update_Ah *= ch;
				Ah -= update_Ah;
				update_Mh = productTranspose(ph_mat, ph_mat);
				update_Mh *= ch;
				Mh -= update_Mh;
				update_Ch = productTranspose(wh_mat, ph_mat);
				Ch -= update_Ch;
				if (pb)
					pb->update(h + 1);
                                if (mVerbose) {time ( &rawtime ); cout<<ctime(&rawtime)<<"done"<<endl;}
			}
		} else if (method == "pls1") {
			Vec s(n);
			Vec old_s(n);
			Vec y(n);
			Vec lx(p);
			Vec ly(1);
			Mat T(n,k);
			Mat X = X_vmatrix->toMat();
			y << Y_vmatrix->toMat();
			
			PP<ProgressBar> pb(
							   report_progress? new ProgressBar("Computing the PLS1 components", k)
											  : 0);
			
			for (int h = 0; h < k; h++) {
				if (pb)
					pb->update(h);
				s << y;
				normalize(s, 2.0);
				bool finished = false;
				while (!finished) {
					old_s << s;
					transposeProduct(lx, X, s);
					product(s, X, lx);
					normalize(s, 2.0);
					if (dist(old_s, s, 2) < precision) {
						finished = true;
					}
				}
				ly[0] = dot(s, y);
				transposeProduct(lx, X, s);
				T.column(h) << s;
				P.column(h) << lx;
				Q.column(h) << ly;
				// X = X - s lx'
				// y = y - s ly
				for (int i = 0; i < n; i++) {
					for (int j = 0; j < p; j++) {
						X(i,j) -= s[i] * lx[j];
					}
					y[i] -= s[i] * ly[0];
				}
			}
			DBG_MODULE_LOG << " Computation of the corresponding coefficients" << endl;
			Mat tmp(n, p);
			productTranspose(tmp, T, P);
			Mat U, Vt;
			Vec D;
			real safeguard = 1.1; // Because the SVD may crash otherwise.
			SVD(tmp, U, D, Vt, 'A', safeguard);
			for (int i = 0; i < D.length(); i++) {
				if (abs(D[i]) < precision) {
					D[i] = 0;
				} else {
					D[i] = 1.0 / D[i];
				}
			}
			Mat tmp2(n,p);
			tmp2.fill(0);
			for (int i = 0; i < D.length(); i++) {
				if (!fast_exact_is_equal(D[i], 0)) {
					tmp2(i) << D[i] * Vt(i);
				}
			}
			product(tmp, U, tmp2);
			transposeProduct(W, tmp, T);
		}else if (method == "simpls") {
			// Initialize the various coefficients.
			DBG_MODULE_LOG << "Initialization of the coefficients" << endl;
			Vec ph(p);
			Vec qh(m);
                        Vec vh(p);
			Vec wh(p);
			Vec tmp(p);
			real ch;
                        if (mVerbose) {time ( &rawtime ); cout<<ctime(&rawtime)<<"Computing A0=X'Y ..."<<flush;}
                        Mat Ah = transposeProduct(X_vmatrix, Y_vmatrix);
                        if (mVerbose)  {time ( &rawtime ); cout<<ctime(&rawtime)<<"done"<<endl<<"Computing M0=X'X ..."<<flush;}
                        Mat Mh = transposeProduct(X_vmatrix, X_vmatrix);
                        if (mVerbose) {time ( &rawtime ); cout<<ctime(&rawtime)<<"done"<<endl;}
			Mat Ch(p,p);    // Initialized to Identity(p).
			Mat Ah_t_Ah;
			Mat update_Ah(p,m);
			Mat update_Mh(p,p);
			Mat update_Ch(p,p);
			for (int i = 0; i < p; i++) {
				for (int j = i+1; j < p; j++) {
					Ch(i,j) = Ch(j,i) = 0;
				}
				Ch(i,i) = 1;
			}
			
			// Iterate k times to find the k first factors.
			PP<ProgressBar> pb(
							   report_progress? new ProgressBar("Computing the PLS-simpls components", k)
											  : 0);
			for (int h = 0; h < this->k; h++) {
                                //1
				Ah_t_Ah = transposeProduct(Ah,Ah);
				if (m == 1) {
					// No need to compute the eigenvector.
					qh[0] = 1;
				} else {
					NIPALSEigenvector(Ah_t_Ah, qh, precision);
				}
                                //2
				product(wh, Ah, qh);
				product(tmp, Mh, wh);
				ch = dot(wh, tmp);
				wh /= sqrt(ch);
				W.column(h) << wh;
                                //3
                                product(ph, Mh, wh);
				P.column(h) << ph;
                                //4
                                transposeProduct(qh, Ah, wh);
				Q.column(h) << qh;
                                //5
                                product(vh, Ch, ph);
                                normalize(vh , 2.0);
                                //6 
                                Mat ph_mat(p, 1, ph);
                                Mat vh_mat(p, 1, vh);

                                update_Ch = productTranspose(vh_mat, vh_mat);
                                Ch -= update_Ch;
                                update_Mh = productTranspose(ph_mat, ph_mat);
                                Mh -= update_Mh;
                                //7
				Ah = product(Ch, Ah);

				if (pb)
					pb->update(h + 1);
			}
		}
		B.resize(p,m);
                if (mVerbose) {time ( &rawtime ); cout<<ctime(&rawtime)<<"Computing B=WQ' ... "<<flush;}
		productTranspose(B, W, Q);
                if (mVerbose) {time ( &rawtime ); cout<<ctime(&rawtime)<<"done"<<endl;}
		// If we requested confidence intervals, compute the variance of the
		// residuals on the training set
		if (compute_confidence)
			computeResidVariance(train_set, resid_variance);
		else
			resid_variance.resize(0);
		
		MODULE_LOG << "PLS training ended" << endl;
		stage = 1;
	}
	
	
	//#####  computeResidVariance  ################################################
	
	void PLS::computeResidVariance(VMat dataset, Vec& resid_variance)
	{
		PLASSERT( dataset.isNotNull() && m >= 0 );
		bool old_output_score = output_the_score;
		bool old_output_target= output_the_target;
		output_the_score  = false;
		output_the_target = true;
		
		resid_variance.resize(m);
		resid_variance.fill(0.0);
		Vec input, target, output(m);
		real weight;
		for (int i=0, n=dataset.length() ; i<n ; ++i) {
			dataset->getExample(i, input, target, weight);
			computeOutput(input, output);
			target -= output;
			target *= target;                    // Square of residual
			resid_variance += target;
		}
		resid_variance /= (dataset.length() - inputsize());
		
		output_the_score  = old_output_score;
		output_the_target = old_output_target;
	}
	
	
} // end of namespace PLearn


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
