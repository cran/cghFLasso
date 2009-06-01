///#include "FL.h"

//-----------------------------------------------------
// Some constants and utility functions
//-----------------------------------------------------

#define FLV_X 0
#define FLV_DY 1
#define FL_SEGSZ 2

#define FL_RESID_COEF 0.5
#define FL_ENDPT_KNOT 1
#define FL_ENDPT_KNOT_FUDGE 1e-4

//////////////////// added by Pei 
#ifndef FUSED_LASSO_C_SRC_FL_H__

#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <R_ext/Print.h>

SEXP L2L1Vit(SEXP obsSeq, SEXP obsWts, SEXP lambda2, SEXP lambda1,
			 SEXP retPath, SEXP maxSegs, SEXP nSegs, SEXP backPtrs,
			 SEXP primBds);

SEXP L2L1VitPath(SEXP obsSeq, SEXP lambda2,  SEXP lambda1, SEXP retPath, SEXP maxSegs,
				 SEXP segmentVec, SEXP primBds);

SEXP L2L1ExpandFit(SEXP betaPath, SEXP betaSegs, SEXP lam1, SEXP lam2Inds, SEXP retFit);

#endif
/////////////////// End of Pei's edit

int GetInt(SEXP p, int default_val, int* err_code){
	if(p == R_NilValue){
		if(err_code) *err_code = 1;
		return default_val;
	}else if(IS_INTEGER(p)){
		return INTEGER(p)[0];
	}else if(IS_LOGICAL(p)){
		if(LOGICAL(p)[0]) return 1;
		else return 0;
	}else if(IS_NUMERIC(p)){
		return (int)(REAL(p)[0]);
	}else{
		if(err_code) *err_code = 2;
		return default_val;
	}
}

double GetNumeric(SEXP p, double default_val, int* err_code){
	if(p == R_NilValue){
		if(err_code) *err_code = 1;
		return default_val;
	}else if(IS_INTEGER(p)){
		return INTEGER(p)[0];
	}else if(IS_LOGICAL(p)){
		if(LOGICAL(p)[0]) return 1.0;
		else return 0.0;
	}else if(IS_NUMERIC(p)){
		return REAL(p)[0];
	}else{
		if(err_code) *err_code = 2;
		return default_val;
	}
}

double * AllocProtectReal(int n){

  SEXP x = R_NilValue;
  PROTECT(x = NEW_NUMERIC(n));
  return( REAL(x) );
}

int * AllocProtectInt(int n){

  SEXP x = R_NilValue;
  PROTECT(x = NEW_INTEGER(n));
  return( INTEGER(x) );
}

void AddQuadPcwsQuad(double * inp_segs, int n_inp_segs,
					 double qd_root, double qd_scale)
{
  double a = 2.0*qd_scale;
  double b = -2.0 * qd_root * qd_scale;

  for(int i = 0; i < n_inp_segs; i++, inp_segs += FL_SEGSZ){
    double d1 = a * inp_segs[FLV_X] + b;
	inp_segs[FLV_DY] += d1;
  }
}

int L2L1VitArgmax(double * inp_segs, int n_inp_segs,
				 double * ret_segs, int * n_ret_segs,
				 double lambda2, double * mid_seg)
{
  int n_seg = n_inp_segs;
  int n_rs = 0;
  double * rs = ret_segs;

  double neg_lam2 = -lambda2;

  //some temporary variables used in the computations below
  double d1, seg_width, a1, x_new;

  int left_added = -1;

  mid_seg[0] = inp_segs[0];
  mid_seg[1] = inp_segs[ FL_SEGSZ*(n_seg-1) ];

  double * is_lst = inp_segs + ( FL_SEGSZ * (n_seg-1) );

  int j = 0;
  while(j < n_seg){
    double * is = inp_segs + (FL_SEGSZ*j);

	double dy = is[FLV_DY];
	double x = is[FLV_X];

    if(left_added < 0){

      if(j == 0 && dy <= lambda2){
	    left_added = j;
	    // The derivative of lambda2 is never acheived
		//copy the left knot entirely
		for(int k = 0; k < FL_SEGSZ; k++) rs[k] = is[k];
		//increment the pointers
		rs += FL_SEGSZ;
		n_rs++;

		j--;

	  }else if(j+1 < n_seg && lambda2 >= is[FL_SEGSZ + FLV_DY] ){
        //If the left-derivative the next knot to the right is
	    //less than lambda2, then the derivative lambda2 will
	    //be acheived somewhere between

        left_added = j;
        d1 = is[FL_SEGSZ + FLV_DY] - dy;
        seg_width = is[FL_SEGSZ+FLV_X] - x;

        a1 = seg_width * (lambda2 - dy) / d1;
        mid_seg[0] = x_new = x + a1;

#ifdef FL_ENDPT_KNOT
		//insert a knot at the left-most position
        rs[FLV_X] = inp_segs[FLV_X];
        rs[FLV_DY] = lambda2;

        n_rs++;
		rs += FL_SEGSZ;
#endif

        rs[0] = x_new;
        rs[FLV_DY] = lambda2;

        n_rs++;
		rs += FL_SEGSZ;

        j--;
      }
    }else{

      if(j+1 < n_seg && is[FL_SEGSZ + FLV_DY] <= neg_lam2 ){
        //If the left-derivative the next knot to the right is
	    //less than -lambda2, then the derivative -lambda2 will
	    //be achieved somewhere between

        d1 = is[FL_SEGSZ + FLV_DY] - dy;
        seg_width = is[FL_SEGSZ+FLV_X] - x;

        a1 = seg_width * (neg_lam2 - dy) / d1;
        x_new = x + a1;

		if(left_added != j){
			//copy the left knot entirely
			for(int k = 0; k < FL_SEGSZ; k++) rs[k] = is[k];
			//increment the pointers
			rs += FL_SEGSZ;
			n_rs++;
		}

        mid_seg[1] = rs[FLV_X] = x_new;
        rs[FLV_DY] = neg_lam2;

		rs += FL_SEGSZ;
        n_rs++;

#ifdef FL_ENDPT_KNOT
        rs[FLV_X] = is_lst[FLV_X];
        rs[FLV_DY] = neg_lam2;

        n_rs++;
		rs += FL_SEGSZ;
#endif

        break;
      }else if(left_added != j){
        // This is not the last segment.  We add this knot unchanged.

        for(int k = 0; k < FL_SEGSZ; k++) rs[k] = is[k];

		rs += FL_SEGSZ;
        n_rs++;
      }
    }
    j++;
  }
  if(left_added < 0){
    return -1; //check concavity
  }

  *n_ret_segs = n_rs;
  return 1;
}

void L2L1VitMsgMax(double * inp_segs, int n_inp_segs,
				   double * x_opt)
{
  double * is = inp_segs;
  for(int j = 0; j < n_inp_segs; j++, is += FL_SEGSZ){

	double dy = is[FLV_DY];
	double x = is[FLV_X];

    if( dy == 0.0){
	  if(x_opt) *x_opt = x;
      return;

    }else if(j+1 < n_inp_segs && dy >= 0.0 && 0.0 >= is[FL_SEGSZ + FLV_DY] ){

      double d1 = is[FL_SEGSZ + FLV_DY] - dy;
      double seg_width = is[FL_SEGSZ+FLV_X] - x;

      double x_new = x + seg_width * (0.0 - dy) / d1;

	  if(x_opt) *x_opt = x_new;
      return;
    }
  }
  return;
}

int L2L1VitFwd(double lam2, double * o, double * wts,
			   double ** msg_buf, int * buf_len, int max_segs,
			   double * back_ptrs, int * nsegs, int n_o, int vit_msg_len,
			   double obs_min, double obs_max, double * last_beta)
{
  double * mbuf = *msg_buf;
  int mbuf_len = *buf_len;
  double * vit_msg1 = mbuf;

  vit_msg1[0] = obs_min - FL_ENDPT_KNOT_FUDGE;
  vit_msg1[1] = 0.0;
  vit_msg1[2] = obs_max + FL_ENDPT_KNOT_FUDGE;
  vit_msg1[3] = 0.0;

  int vit1_len = 2, vit2_len = -1;

  if(nsegs) nsegs[0] = vit1_len;

  if(R_FINITE(o[0])){
	  AddQuadPcwsQuad(vit_msg1, vit1_len, o[0],
		              FL_RESID_COEF * ((wts) ? -wts[0] : -1.0) );
  }

  for(int i = 1; i < n_o; i++){
	  double * bp = back_ptrs + (2*i);

	  double * vm2 = (vit_msg1 == mbuf) ? (mbuf + (vit1_len+3)*FL_SEGSZ) : mbuf;

	  int r1 = L2L1VitArgmax(vit_msg1, vit1_len,
						     vm2, &vit2_len,
						     lam2, bp);

	  if(r1 != 1){
		  return r1;
	  }

	  // Check to see if the buffer is almost full.  If it is,
	  // allocate a larger one and switch pointers.

	  if((vit1_len + 3)*2*FL_SEGSZ > mbuf_len && vm2 == mbuf){

		  int new_len = (vit1_len + 20)*2*FL_SEGSZ;

		  if(vit1_len + 20 > max_segs){
			  free(mbuf);
			  // We will not allocate more than this
			  error("Viterbi message at index %d has %d segments (max is %d).  exiting\n",
				    i, vit1_len, vit1_len);
		  }

		  double * tmp_ptr = malloc(new_len * sizeof(double) );
		  if(tmp_ptr == NULL){
			  return -100;
		  }else{

			memcpy(tmp_ptr, mbuf, mbuf_len*sizeof(double) );
			free(mbuf);

			mbuf_len = *buf_len = new_len;
			vm2 = mbuf = *msg_buf = tmp_ptr;
		  }
	  }

	  if(R_FINITE(o[i])){
		  AddQuadPcwsQuad(vm2, vit2_len, o[i],
			              FL_RESID_COEF * ((wts) ? -wts[i] : -1.0) );
	  }

	  vit_msg1 = vm2;
	  vit1_len = vit2_len;

	  if(nsegs) nsegs[i] = vit1_len;
  }

  if(last_beta){
	  L2L1VitMsgMax(vit_msg1, vit1_len, last_beta);
  }

  return 1;
}

double soft_thresh(double x, double amt){
	if(amt == 0.0) return x;
	return ( fabs(x) < amt ) ? 0.0 :
	                           ( (x > 0.0) ? (x-amt) : (x+amt) );

}

SEXP L2L1ExpandFit(SEXP betaPath, SEXP betaSegs, SEXP lam1, SEXP lam2Inds, SEXP retFit)
{
	double lm1 = GetNumeric(lam1, 0, 0);

	double * ret_fit = REAL(retFit);

	int n_lam2i = LENGTH(lam2Inds);
	int * lam2i = INTEGER(lam2Inds);

	for(int i = 0; i < n_lam2i; i++){
		SEXP bp_sxp = VECTOR_ELT(betaPath, lam2i[i]);
		SEXP bsg_sxp = VECTOR_ELT(betaSegs, lam2i[i]);

		int n_seg = LENGTH(bp_sxp);
		double * bv = REAL(bp_sxp);

		int * sgv = INTEGER(bsg_sxp);

		double cur_b = soft_thresh(bv[0], lm1);

		int n_fill = 1 + sgv[1] - sgv[0];
		for(int j = 0; j < n_fill; j++, ret_fit++){
			*ret_fit = cur_b;
		}

		for(int k = 1, k2 = 2; k < n_seg; k++, k2 += 2){
			cur_b = soft_thresh(bv[k], lm1);
			n_fill = 1 + sgv[k2 + 1] - sgv[k2];

			for(int j = 0; j < n_fill; j++, ret_fit++){
				*ret_fit = cur_b;
			}
		}
	}

	return R_NilValue;
}

SEXP L2L1GetPrimBds(SEXP betaPath, SEXP betaSegs, SEXP lam1, SEXP lam2Inds, SEXP retBds, SEXP retFit)
{
	int n_lam1 = LENGTH(lam1);
	double * lam1v = REAL(lam1);

	double * bds = REAL(retBds);

	int n_lam2i = LENGTH(lam2Inds);
	int * lam2i = INTEGER(lam2Inds);

	for(int i = 0; i < n_lam2i; i++){
		SEXP bp_sxp = VECTOR_ELT(betaPath, lam2i[i]);
		SEXP bsg_sxp = VECTOR_ELT(betaSegs, lam2i[i]);

		int n_seg = LENGTH(bp_sxp);
		double * bv = REAL(bp_sxp);

		int * sgv = INTEGER(bsg_sxp);

		for(int j = 0; j < n_lam1; j++){

			double lm1 = lam1v[j];
			double bd1 = 0.0, bd2 = 0.0;

			double cur_b = soft_thresh(bv[0], lm1);

			bd1 += fabs(cur_b) * (double)(1 + sgv[1] - sgv[0]);

			for(int k = 1, k2 = 2; k < n_seg; k++, k2 += 2){
				double next_b = soft_thresh(bv[k], lm1);

				bd1 += fabs(next_b) * (double)(1 + sgv[k2 + 1] - sgv[k2]);
				bd2 += fabs(next_b - cur_b);

				cur_b = next_b;
			}
			bds[0] = bd1;
			bds[1] = bd2;

			bds += 2;
		}
	}

	return R_NilValue;
}

int L2L1GetNFused(double beta_hat, int n_o, double * back_ptrs)
{
	double btht = beta_hat;
	int nfsd2 = 0;

	if(n_o == 1){
	  nfsd2 = 1;
	}else{
	  double * bp = back_ptrs + (2*(n_o-1));

	  for(int i = n_o-2; i >= 0; i--, bp -= 2){

		  if(btht > bp[1]){
			  btht = bp[1];
			  nfsd2++;
		  }else if(btht < bp[0]){
			  btht = bp[0];
			  nfsd2++;
		  }
		  if(i == 0){
			  nfsd2++;
		  }
	  }
	}
	return nfsd2;
}

void L2L1BackTrace(double beta_last, double lam1, double * out_beta, int n_obs,
				   double * back_ptrs, double * bd1, double * bd2)
{
  double prev_b = beta_last, cur_b = -1.0, cur_b_shr = -1.0, prev_b_shr = -1;

  double * bp = back_ptrs + (n_obs-1)*2;

  prev_b_shr = out_beta[n_obs-1] = soft_thresh(prev_b, lam1);

  double rbd1 = fabs(prev_b_shr), rbd2 = 0.0;

  for(int i = n_obs-2; i >= 0; i--, bp -= 2){

	  if(prev_b > bp[1]){
		  cur_b = bp[1];
		  cur_b_shr = soft_thresh(cur_b, lam1);
		  rbd2 += fabs(cur_b_shr - prev_b_shr);
	  }else if(prev_b < bp[0]){
		  cur_b = bp[0];
		  cur_b_shr = soft_thresh(cur_b, lam1);
		  rbd2 += fabs(cur_b_shr - prev_b_shr);
	  }else{
		  cur_b = prev_b;
		  cur_b_shr = prev_b_shr;
	  }
	  out_beta[i] = cur_b_shr;

	  rbd1 += fabs(cur_b_shr);

	  prev_b_shr = cur_b_shr;
	  prev_b = cur_b;
  }

  if(bd1) *bd1 = rbd1;
  if(bd2) *bd2 = rbd2;
}

SEXP L2L1Vit(SEXP obsSeq, SEXP obsWts, SEXP lambda2, SEXP lambda1,
			 SEXP retPath, SEXP maxSegs, SEXP nSegs, SEXP backPtrs,
			 SEXP primBds)
{
  int max_segs = GetInt(maxSegs, 0, 0);

  double * o   = REAL(obsSeq);
  double * wts = REAL(obsWts);
  double lam2 = GetNumeric(lambda2, 0, 0);
  double lam1 = GetNumeric(lambda1, 0, 0);

  int n_obs = LENGTH(obsSeq);
  int n_protect = 0;

  double * back_ptrs = REAL(backPtrs);

  int msg_buf_len = FL_SEGSZ*2*30;
  double * msg_buf = malloc( msg_buf_len*sizeof(double) );

  int * n_segs = INTEGER(nSegs);

  double obs_min = R_PosInf, obs_max = R_NegInf;
  for(int i = 0; i < n_obs; i++){
	  if(R_FINITE(o[i])){
		  if(o[i] < obs_min) obs_min = o[i];
		  else if(o[i] > obs_max) obs_max = o[i];
	  }
  }

  SEXP ret_sxp;
  PROTECT(ret_sxp = NEW_INTEGER(1)); n_protect++;

  double * rp = REAL(retPath);

  int r1 = L2L1VitFwd(lam2, o, wts,
	                  &msg_buf, &msg_buf_len, max_segs,
	                  back_ptrs, n_segs, n_obs, max_segs, obs_min, obs_max,
					  (rp + (n_obs-1)) );

  if(r1 != 1){
	  INTEGER(ret_sxp)[0] = r1;
	  UNPROTECT(n_protect);
	  return ret_sxp;
  }

  double * bd1 = NULL, *bd2 = NULL;
  if(primBds != R_NilValue){
	  bd1 = REAL(primBds);
	  bd2 = bd1 + 1;
  }

  L2L1BackTrace(rp[n_obs-1], lam1, rp, n_obs, back_ptrs, bd1, bd2);

  free(msg_buf);

  INTEGER(ret_sxp)[0]  = 1;
  if(n_protect > 0) UNPROTECT(n_protect);
  return ret_sxp;
}


SEXP L2L1VitPath(SEXP obsSeq, SEXP lambda2, SEXP lambda1, SEXP retPath, SEXP maxSegs,
				 SEXP segmentVec, SEXP primBds)
{
  int segmented_ret = (segmentVec != R_NilValue) ? 1 : 0;

  int max_segs = GetInt(maxSegs, 0, 0);

  double * all_obs   = REAL(obsSeq);

  int n_obs = LENGTH(obsSeq);
  int n_protect = 0;

  double * back_ptrs     = AllocProtectReal(2*n_obs);  n_protect++;
  int * fused_segs1   = NULL;
  int * fused_segs2   = NULL;

  double *o2 = NULL, *wts2 = NULL, *o3 = NULL, *wts3 = NULL;

  int msg_buf_len = FL_SEGSZ*2*30;
  double * msg_buf = malloc( msg_buf_len*sizeof(double) );

  SEXP ret_sxp;
  PROTECT(ret_sxp = NEW_INTEGER(1)); n_protect++;

  double obs_min = R_PosInf, obs_max = R_NegInf;
  for(int i = 0; i < n_obs; i++){

	  if(R_FINITE(all_obs[i])){
		  if(all_obs[i] < obs_min) obs_min = all_obs[i];
		  if(all_obs[i] > obs_max) obs_max = all_obs[i];
	  }
  }

  double lam1 = GetNumeric(lambda1, 0, 0);

  int n_lam2 = LENGTH(lambda2);
  int n_o = n_obs;
  double * o = all_obs;

  double * wts = NULL;

  double * prim_bds = (primBds == R_NilValue) ? NULL : REAL(primBds);

  for(int lam2i = 0; lam2i < n_lam2; lam2i++){

	  double lam2 = REAL(lambda2)[lam2i];

	  double beta_hat = 0.0;

	  int r1 = L2L1VitFwd(lam2, o, wts,
	                      &msg_buf, &msg_buf_len, max_segs,
						  back_ptrs, NULL, n_o, max_segs, obs_min, obs_max,
						  &beta_hat);

	  if(r1 != 1){
		  INTEGER(ret_sxp)[0] = r1;
		  UNPROTECT(n_protect);
		  return ret_sxp;
	  }

	  int * fs = fused_segs1;

	  int nfsd2 = 0;
	  if(o2 == NULL || segmented_ret){
		  //We haven't allocated the buffers for the
		  //fused observations yet
		  nfsd2 = L2L1GetNFused(beta_hat, n_o, back_ptrs);

		  o2 = AllocProtectReal(nfsd2);  n_protect++;
		  wts2 = AllocProtectReal(nfsd2);  n_protect++;

		  fused_segs1 = AllocProtectInt(2*(nfsd2+1));  n_protect++;
		  fused_segs2 = AllocProtectInt(2*(nfsd2+1));  n_protect++;
	  }

	  double * fit_v = NULL;

	  if(segmented_ret){
		  SEXP tmp_sxp;
		  PROTECT(tmp_sxp = NEW_NUMERIC(nfsd2));
		  SET_VECTOR_ELT(retPath, lam2i, tmp_sxp);
		  UNPROTECT(1);

		  fit_v = REAL(VECTOR_ELT(retPath, lam2i));
	  }else{
		  fit_v = REAL(retPath) + n_obs * lam2i;
	  }

	  int seg_R = (fs) ? fs[0] : (n_obs-1);
	  int seg_L = (fs) ? fs[1] : (n_obs-1);

	  int n_fused2 = 0;
	  fused_segs2[0] = seg_R;

	  if(fs) fs += 2;

	  double bd1 = 0.0, bd2 = 0.0;
	  double beta_hat_shr = beta_hat;

	  if(segmented_ret){
		  fit_v[(nfsd2-1) - n_fused2] = beta_hat_shr;
		  bd1 += fabs(beta_hat_shr);
	  }else{
		  for(int k = seg_L; k <= seg_R; k++){
			  fit_v[k] = beta_hat_shr;
		  }
		  bd1 += fabs(beta_hat_shr) * (double)(1+seg_R - seg_L);
	  }
	  if( !R_FINITE(o[n_o-1]) ){
		  o2[n_fused2] = wts2[n_fused2] = 0;
	  }else if(wts){
		  o2[n_fused2] = o[n_o-1]*wts[n_o-1];
		  wts2[n_fused2] = wts[n_o-1];
	  }else{
		  o2[n_fused2] = o[n_o-1];
		  wts2[n_fused2] = 1.0;
	  }

	  if(n_o == 1){
		  n_fused2 = 1;
		  fused_segs2[0] = n_obs - 1;
		  fused_segs2[1] = 0;

		  o2[0] = o[0] * wts[0];
		  wts2[0] = wts[0];
	  }

	  for(int i = n_o-2; i >= 0; i--){
		  seg_R = (fs) ? fs[0] : i;
		  seg_L = (fs) ? fs[1] : i;

		  double * bp = back_ptrs + (2*(i+1));

		  if(beta_hat > bp[1]){
			  bd2 += fabs(beta_hat - bp[1]);

			  beta_hat = bp[1];
			  beta_hat_shr = beta_hat;

			  fused_segs2[2*n_fused2 + 1] = seg_R+1;
			  n_fused2++;

			  o2[n_fused2] = wts2[n_fused2] = 0.0;
			  fused_segs2[2*n_fused2] = seg_R;


		  }else if(beta_hat < bp[0]){
			  bd2 += fabs(beta_hat - bp[0]);

			  beta_hat = bp[0];
			  beta_hat_shr = beta_hat;

			  fused_segs2[2*n_fused2 + 1] = seg_R+1;
			  n_fused2++;

			  o2[n_fused2] = wts2[n_fused2] = 0.0;
			  fused_segs2[2*n_fused2] = seg_R;
		  }

		  if(R_FINITE(o[i])){
			  if(wts){
				  o2[n_fused2] += o[i]*wts[i];
				  wts2[n_fused2] += wts[i];
			  }else{
				  o2[n_fused2] += o[i];
				  wts2[n_fused2] += 1.0;
			  }
		  }

		  if(segmented_ret){
			  fit_v[(nfsd2-1) - n_fused2] = beta_hat_shr;
		  }else{
			  for(int k = seg_L; k <= seg_R; k++){
				  fit_v[k] = beta_hat_shr;
			  }
		  }
		  bd1 += fabs(beta_hat_shr) * (double)(1+seg_R - seg_L);

		  if(i == 0){
			  fused_segs2[2*n_fused2 + 1] = seg_L;
			  n_fused2++;
		  }

		  if(fs) fs += 2;
	  }
	  if(prim_bds){
		  double * bdv = prim_bds + 2*lam2i;
		  bdv[0] = bd1;
		  bdv[1] = bd2;
	  }
	  // We have stored the fitted parameters.  Now we collapse
	  // observations and fit on the new sequence at the next
	  // iteration

	  obs_min = R_PosInf;
	  obs_max = R_NegInf;

	  if(o3 == NULL){
		  o3 = AllocProtectReal(n_fused2);  n_protect++;
		  wts3 = AllocProtectReal(n_fused2);  n_protect++;
	  }

	  for(int i = 0; i < n_fused2; i++){
		  if( wts2[n_fused2-1-i] > 0.0 ){

			  double z = o2[n_fused2-1-i] / wts2[n_fused2-1-i];
			  if(z < obs_min) obs_min = z;
			  if(z > obs_max) obs_max = z;

			  o3[i] = z;
		  }else{
			  o3[i] = NA_REAL;
		  }
		  wts3[i] = wts2[n_fused2-1-i];
	  }

	  if(n_o == 1){
		  obs_max = obs_min + FL_ENDPT_KNOT_FUDGE;
		  obs_min -= FL_ENDPT_KNOT_FUDGE;
	  }

	  if(segmented_ret){
		  SEXP tmp_sxp, seg_dim;
		  PROTECT(tmp_sxp = NEW_INTEGER(2*nfsd2));

		  PROTECT(seg_dim=NEW_INTEGER(2));
		  INTEGER(seg_dim)[0] = 2;
		  INTEGER(seg_dim)[1] = nfsd2;

		  SET_DIM(tmp_sxp,seg_dim);

		  SET_VECTOR_ELT(segmentVec, lam2i, tmp_sxp);
		  UNPROTECT(2);

		  int * seg_v = INTEGER(VECTOR_ELT(segmentVec, lam2i));
		  for(int k = 0; k < nfsd2; k++){
			  seg_v[1+2*k] = fused_segs2[(nfsd2-1-k)*2]+1;
			  seg_v[2*k] = fused_segs2[1+(nfsd2-1-k)*2]+1;
		  }
	  }

	  o = o3;
	  wts = wts3;

	  fs = fused_segs2;
	  fused_segs2 = fused_segs1;
	  fused_segs1 = fs;

	  n_o = n_fused2;
  }

  free(msg_buf);

  if(segmented_ret){
	  for(int lam2i = 0; lam2i < n_lam2; lam2i++){
		   double * bv = REAL(VECTOR_ELT(retPath, lam2i));
		   int m = LENGTH(VECTOR_ELT(retPath, lam2i));

		   for(int i = 0; i < m; i++){
			   bv[i] = soft_thresh(bv[i], lam1);
		   }
	  }
  }else{
	   double * bv = REAL(retPath);
	   int m = LENGTH(retPath);

	   for(int i = 0; i < m; i++){
		   bv[i] = soft_thresh(bv[i], lam1);
	   }
  }



  INTEGER(ret_sxp)[0]  = 1;
  UNPROTECT(n_protect);
  return ret_sxp;
}


