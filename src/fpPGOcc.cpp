#define USE_FC_LEN_T
#include <string>
#include "util.h"
#include "rpg.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#define R_NO_REMAP
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#ifndef FCONE
# define FCONE
#endif

extern "C" {
  SEXP fpPGOcc(SEXP y_r, SEXP X_r, SEXP Xmu_r, SEXP Xomega_r, SEXP XRE_r, 
               SEXP XmuRE_r, SEXP XomegaRE_r, SEXP d_r, SEXP nVal_r, 
               SEXP consts_r, SEXP nOccRELong_r, SEXP nTPDetRELong_r, 
               SEXP nFPDetRELong_r, SEXP betaStarting_r, SEXP gammaStarting_r, 
               SEXP deltaStarting_r, SEXP sigmaSqPsiStarting_r, SEXP sigmaSqMuStarting_r, 
               SEXP sigmaSqOmegaStarting_r, SEXP betaStarStarting_r, SEXP gammaStarStarting_r, 
               SEXP deltaStarStarting_r, SEXP zStarting_r, SEXP zLongIndx_r, 
               SEXP betaStarIndx_r, SEXP betaLevelIndx_r, SEXP gammaStarIndx_r, 
               SEXP gammaLevelIndx_r, SEXP deltaStarIndx_r, SEXP deltaLevelIndx_r, 
               SEXP muBeta_r, SEXP muGamma_r, SEXP muDelta_r, 
               SEXP SigmaBeta_r, SEXP SigmaGamma_r, SEXP SigmaDelta_r, 
               SEXP sigmaSqPsiA_r, SEXP sigmaSqPsiB_r, SEXP sigmaSqMuA_r, 
               SEXP sigmaSqMuB_r, SEXP sigmaSqOmegaA_r, SEXP sigmaSqOmegaB_r, 
               SEXP tuning_r, SEXP nBatch_r, SEXP batchLength_r, SEXP acceptRate_r, 
               SEXP nThreads_r, SEXP verbose_r, SEXP nReport_r, SEXP samplesInfo_r,
               SEXP chainInfo_r){
   
    /**********************************************************************
     * Initial constants
     * *******************************************************************/
    int i, j, s, r, l, ll, info, nProtect=0, q, k;
    const int inc = 1;
    const double one = 1.0;
    const double zero = 0.0;
    char const *lower = "L";
    char const *ntran = "N";
    char const *ytran = "T";
    
    /**********************************************************************
     * Get Inputs
     * *******************************************************************/
    double *y = REAL(y_r);
    double *X = REAL(X_r);
    double *d = REAL(d_r);
    double *nVal = REAL(nVal_r);
    double *Xmu = REAL(Xmu_r);
    double *Xomega = REAL(Xomega_r);
    int *XRE = INTEGER(XRE_r); 
    int *XmuRE = INTEGER(XmuRE_r); 
    int *XomegaRE = INTEGER(XomegaRE_r); 
    // Load constants
    int J = INTEGER(consts_r)[0];
    int nObs = INTEGER(consts_r)[1]; 
    int pOcc = INTEGER(consts_r)[2];
    int pOccRE = INTEGER(consts_r)[3];
    int nOccRE = INTEGER(consts_r)[4];
    int pTPDet = INTEGER(consts_r)[5];
    int pTPDetRE = INTEGER(consts_r)[6];
    int nTPDetRE = INTEGER(consts_r)[7];
    int pFPDet = INTEGER(consts_r)[8];
    int pFPDetRE = INTEGER(consts_r)[9];
    int nFPDetRE = INTEGER(consts_r)[10];
    int tpDetInt = INTEGER(consts_r)[11];
    int fpDetInt = INTEGER(consts_r)[12];
    int ppTPDet = pTPDet * pTPDet;
    int ppFPDet = pFPDet * pFPDet;
    int ppOcc = pOcc * pOcc; 
    double *muBeta = (double *) R_alloc(pOcc, sizeof(double));   
    F77_NAME(dcopy)(&pOcc, REAL(muBeta_r), &inc, muBeta, &inc);
    double *muGamma = (double *) R_alloc(pTPDet, sizeof(double));   
    F77_NAME(dcopy)(&pTPDet, REAL(muGamma_r), &inc, muGamma, &inc);
    double *muDelta = (double *) R_alloc(pFPDet, sizeof(double));   
    F77_NAME(dcopy)(&pFPDet, REAL(muDelta_r), &inc, muDelta, &inc);
    double *SigmaBetaInv = (double *) R_alloc(ppOcc, sizeof(double));   
    F77_NAME(dcopy)(&ppOcc, REAL(SigmaBeta_r), &inc, SigmaBetaInv, &inc);
    double *SigmaGamma = (double *) R_alloc(ppTPDet, sizeof(double));   
    F77_NAME(dcopy)(&ppTPDet, REAL(SigmaGamma_r), &inc, SigmaGamma, &inc);
    double *SigmaDelta = (double *) R_alloc(ppFPDet, sizeof(double));   
    F77_NAME(dcopy)(&ppFPDet, REAL(SigmaDelta_r), &inc, SigmaDelta, &inc);
    double *sigmaSqPsiA = REAL(sigmaSqPsiA_r); 
    double *sigmaSqPsiB = REAL(sigmaSqPsiB_r); 
    double *sigmaSqMuA = REAL(sigmaSqMuA_r); 
    double *sigmaSqMuB = REAL(sigmaSqMuB_r); 
    double *sigmaSqOmegaA = REAL(sigmaSqOmegaA_r); 
    double *sigmaSqOmegaB = REAL(sigmaSqOmegaB_r); 
    int *nOccRELong = INTEGER(nOccRELong_r); 
    int *nTPDetRELong = INTEGER(nTPDetRELong_r); 
    int *nFPDetRELong = INTEGER(nFPDetRELong_r); 
    int *zLongIndx = INTEGER(zLongIndx_r); 
    int *gammaStarIndx = INTEGER(gammaStarIndx_r); 
    int *gammaLevelIndx = INTEGER(gammaLevelIndx_r);
    int *deltaStarIndx = INTEGER(deltaStarIndx_r); 
    int *deltaLevelIndx = INTEGER(deltaLevelIndx_r);
    int *betaStarIndx = INTEGER(betaStarIndx_r); 
    int *betaLevelIndx = INTEGER(betaLevelIndx_r);
    double *tuning = REAL(tuning_r); 
    int nBatch = INTEGER(nBatch_r)[0]; 
    int batchLength = INTEGER(batchLength_r)[0]; 
    int nSamples = nBatch * batchLength; 
    int nBurn = INTEGER(samplesInfo_r)[0]; 
    int nThin = INTEGER(samplesInfo_r)[1];
    int nPost = INTEGER(samplesInfo_r)[2]; 
    int currChain = INTEGER(chainInfo_r)[0];
    int nChain = INTEGER(chainInfo_r)[1];
    double acceptRate = REAL(acceptRate_r)[0];
    int nThreads = INTEGER(nThreads_r)[0];
    int verbose = INTEGER(verbose_r)[0];
    int nReport = INTEGER(nReport_r)[0];

    int status = 0; 
    int thinIndx = 0;
    int sPost = 0;  

#ifdef _OPENMP
    omp_set_num_threads(nThreads);
#else
    if(nThreads > 1){
      Rf_warning("n.omp.threads > 1, but source not compiled with OpenMP support.");
      nThreads = 1;
    }
#endif
    
    /**********************************************************************
     * Print Information 
     * *******************************************************************/
    if(verbose){
      if (currChain == 1) {
        Rprintf("----------------------------------------\n");
        Rprintf("\tModel description\n");
        Rprintf("----------------------------------------\n");
        Rprintf("False positive occupancy model with Polya-Gamma latent\nvariable fit with %i sites.\n\n", J);
        Rprintf("Samples per chain: %i (%i batches of length %i)\n", nSamples, nBatch, batchLength);
        Rprintf("Burn-in: %i \n", nBurn); 
        Rprintf("Thinning Rate: %i \n", nThin); 
        Rprintf("Number of Chains: %i \n", nChain);
        Rprintf("Total Posterior Samples: %i \n\n", nPost * nChain); 
#ifdef _OPENMP
        Rprintf("Source compiled with OpenMP support and model fit using %i thread(s).\n\n", nThreads);
#else
        Rprintf("Source not compiled with OpenMP support.\n\n");
#endif
      }
      Rprintf("----------------------------------------\n");
      Rprintf("\tChain %i\n", currChain);
      Rprintf("----------------------------------------\n");
      Rprintf("Sampling ... \n");
      #ifdef Win32
        R_FlushConsole();
      #endif
    }

    /**********************************************************************
     * Parameters
     * *******************************************************************/
    // Occupancy covariates
    double *beta = (double *) R_alloc(pOcc, sizeof(double));   
    F77_NAME(dcopy)(&pOcc, REAL(betaStarting_r), &inc, beta, &inc);
    // Occupancy random effect variances
    double *sigmaSqPsi = (double *) R_alloc(pOccRE, sizeof(double)); 
    F77_NAME(dcopy)(&pOccRE, REAL(sigmaSqPsiStarting_r), &inc, sigmaSqPsi, &inc); 
    // Latent occupancy random effects
    double *betaStar = (double *) R_alloc(nOccRE, sizeof(double)); 
    F77_NAME(dcopy)(&nOccRE, REAL(betaStarStarting_r), &inc, betaStar, &inc); 

    // True positive detection covariates
    double *gamma = (double *) R_alloc(pTPDet, sizeof(double));   
    F77_NAME(dcopy)(&pTPDet, REAL(gammaStarting_r), &inc, gamma, &inc);
    // True positive detection random effect variances
    double *sigmaSqMu = (double *) R_alloc(pTPDetRE, sizeof(double)); 
    F77_NAME(dcopy)(&pTPDetRE, REAL(sigmaSqMuStarting_r), &inc, sigmaSqMu, &inc); 
    // Latent true positive detection random effects
    double *gammaStar = (double *) R_alloc(nTPDetRE, sizeof(double)); 
    F77_NAME(dcopy)(&nTPDetRE, REAL(gammaStarStarting_r), &inc, gammaStar, &inc); 

    // False positive detection covariates
    double *delta = (double *) R_alloc(pFPDet, sizeof(double));   
    F77_NAME(dcopy)(&pFPDet, REAL(deltaStarting_r), &inc, delta, &inc);
    // False positive detection random effect variances
    double *sigmaSqOmega = (double *) R_alloc(pFPDetRE, sizeof(double)); 
    F77_NAME(dcopy)(&pFPDetRE, REAL(sigmaSqOmegaStarting_r), &inc, sigmaSqOmega, &inc); 
    // Latent false positive detection random effects
    double *deltaStar = (double *) R_alloc(nFPDetRE, sizeof(double)); 
    F77_NAME(dcopy)(&nFPDetRE, REAL(deltaStarStarting_r), &inc, deltaStar, &inc); 
    
    // Latent Occurrence
    double *z = (double *) R_alloc(J, sizeof(double));   
    F77_NAME(dcopy)(&J, REAL(zStarting_r), &inc, z, &inc);

    // Auxiliary variables
    double *omegaPGOcc = (double *) R_alloc(J, sizeof(double)); zeros(omegaPGOcc, J);
    double *kappaOcc = (double *) R_alloc(J, sizeof(double)); zeros(kappaOcc, J);

    /**********************************************************************
     * Return Stuff
     * *******************************************************************/
    SEXP betaSamples_r;
    PROTECT(betaSamples_r = Rf_allocMatrix(REALSXP, pOcc, nPost)); nProtect++;
    zeros(REAL(betaSamples_r), pOcc * nPost);
    SEXP gammaSamples_r; 
    PROTECT(gammaSamples_r = Rf_allocMatrix(REALSXP, pTPDet, nPost)); nProtect++;
    zeros(REAL(gammaSamples_r), pTPDet * nPost);
    SEXP deltaSamples_r; 
    PROTECT(deltaSamples_r = Rf_allocMatrix(REALSXP, pFPDet, nPost)); nProtect++;
    zeros(REAL(deltaSamples_r), pFPDet * nPost);
    SEXP zSamples_r; 
    PROTECT(zSamples_r = Rf_allocMatrix(REALSXP, J, nPost)); nProtect++; 
    zeros(REAL(zSamples_r), J * nPost);
    SEXP psiSamples_r; 
    PROTECT(psiSamples_r = Rf_allocMatrix(REALSXP, J, nPost)); nProtect++; 
    zeros(REAL(psiSamples_r), J * nPost);
    // True positive detection random effects
    SEXP sigmaSqMuSamples_r; 
    SEXP gammaStarSamples_r; 
    if (pTPDetRE > 0) {
      PROTECT(sigmaSqMuSamples_r = Rf_allocMatrix(REALSXP, pTPDetRE, nPost)); nProtect++;
      zeros(REAL(sigmaSqMuSamples_r), pTPDetRE * nPost);
      PROTECT(gammaStarSamples_r = Rf_allocMatrix(REALSXP, nTPDetRE, nPost)); nProtect++;
      zeros(REAL(gammaStarSamples_r), nTPDetRE * nPost);
    }
    // False positive detection random effects
    SEXP sigmaSqOmegaSamples_r; 
    SEXP deltaStarSamples_r;
    if (pFPDetRE > 0) {
      PROTECT(sigmaSqOmegaSamples_r = Rf_allocMatrix(REALSXP, pFPDetRE, nPost)); nProtect++;
      zeros(REAL(sigmaSqOmegaSamples_r), pFPDetRE * nPost);
      PROTECT(deltaStarSamples_r = Rf_allocMatrix(REALSXP, nFPDetRE, nPost)); nProtect++;
      zeros(REAL(deltaStarSamples_r), nFPDetRE * nPost);
    }
    // Occurrence random effects
    SEXP sigmaSqPsiSamples_r; 
    SEXP betaStarSamples_r; 
    if (pOccRE > 0) {
      PROTECT(sigmaSqPsiSamples_r = Rf_allocMatrix(REALSXP, pOccRE, nPost)); nProtect++;
      zeros(REAL(sigmaSqPsiSamples_r), pOccRE * nPost);
      PROTECT(betaStarSamples_r = Rf_allocMatrix(REALSXP, nOccRE, nPost)); nProtect++;
      zeros(REAL(betaStarSamples_r), nOccRE * nPost);
    }
    // Likelihood samples for WAIC. 
    SEXP likeSamples_r;
    PROTECT(likeSamples_r = Rf_allocMatrix(REALSXP, J, nPost)); nProtect++;
    zeros(REAL(likeSamples_r), J * nPost);
    
    /********************************************************************
      Some constants and temporary variables to be used later
    ********************************************************************/
    int JpOcc = J * pOcc; 
    int JpOccRE = J * pOccRE; 
    int nObspTPDetRE = nObs * pTPDetRE;
    int nObspFPDetRE = nObs * pFPDetRE;
    double tmp_0, tmp_02; 
    double *tmp_one = (double *) R_alloc(inc, sizeof(double)); 
    double *tmp_ppOcc = (double *) R_alloc(ppOcc, sizeof(double)); 
    double *tmp_pOcc = (double *) R_alloc(pOcc, sizeof(double));
    double *tmp_pOcc2 = (double *) R_alloc(pOcc, sizeof(double));
    int *tmp_J = (int *) R_alloc(J, sizeof(int));
    for (j = 0; j < J; j++) {
      tmp_J[j] = 0; 
    }
    double *tmp_nObs = (double *) R_alloc(nObs, sizeof(double)); 
    double *tmp_nObs2 = (double *) R_alloc(nObs, sizeof(double)); 
    double *tmp_JpOcc = (double *) R_alloc(JpOcc, sizeof(double));
    double *tmp_J1 = (double *) R_alloc(J, sizeof(double));
   
    // For latent occupancy and WAIC
    double logitPsi, logLikeMax;
    double *mu = (double *) R_alloc(nObs, sizeof(double)); zeros(mu, nObs);
    double *omega = (double *) R_alloc(nObs, sizeof(double)); zeros(omega, nObs);
    double *yWAIC = (double *) R_alloc(J, sizeof(double)); zeros(yWAIC, J);
    double *psi = (double *) R_alloc(J, sizeof(double));
    zeros(psi, J);
    // Site-level log-likelihood of y and d when z = 1 and when z = 0
    double *logLikeOcc = (double *) R_alloc(J, sizeof(double)); zeros(logLikeOcc, J);
    double *logLikeUnocc = (double *) R_alloc(J, sizeof(double)); zeros(logLikeUnocc, J);
    // Site-level number of validated true positives
    double *dSum = (double *) R_alloc(J, sizeof(double)); zeros(dSum, J);

    // For normal priors
    // Occupancy regression coefficient priors. 
    F77_NAME(dpotrf)(lower, &pOcc, SigmaBetaInv, &pOcc, &info FCONE); 
    if(info != 0){Rf_error("c++ error: dpotrf SigmaBetaInv failed\n");}
    F77_NAME(dpotri)(lower, &pOcc, SigmaBetaInv, &pOcc, &info FCONE); 
    if(info != 0){Rf_error("c++ error: dpotri SigmaBetaInv failed\n");}
    double *SigmaBetaInvMuBeta = (double *) R_alloc(pOcc, sizeof(double)); 
    F77_NAME(dsymv)(lower, &pOcc, &one, SigmaBetaInv, &pOcc, muBeta, &inc, &zero, 
        	    SigmaBetaInvMuBeta, &inc FCONE);

    /**********************************************************************
     * Prep for random effects
     * *******************************************************************/
    // Site-level sums of the occurrence random effects
    double *betaStarSites = (double *) R_alloc(J, sizeof(double)); 
    zeros(betaStarSites, J); 
    int *betaStarLongIndx = (int *) R_alloc(JpOccRE, sizeof(int));
    // Initial sums
    for (j = 0; j < J; j++) {
      for (l = 0; l < pOccRE; l++) {
        betaStarLongIndx[l * J + j] = which(XRE[l * J + j], betaLevelIndx, nOccRE);
        betaStarSites[j] += betaStar[betaStarLongIndx[l * J + j]];
      }
    }
    // Observation-level sums of the true positive detection random effects
    double *gammaStarObs = (double *) R_alloc(nObs, sizeof(double)); 
    zeros(gammaStarObs, nObs); 
    int *gammaStarLongIndx = (int *) R_alloc(nObspTPDetRE, sizeof(int));
    // Get sums of the current REs for each site/visit combo
    for (i = 0; i < nObs; i++) {
      for (l = 0; l < pTPDetRE; l++) {
        gammaStarLongIndx[l * nObs + i] = which(XmuRE[l * nObs + i], gammaLevelIndx, nTPDetRE);
        gammaStarObs[i] += gammaStar[gammaStarLongIndx[l * nObs + i]];
      }
    }
    // Observation-level sums of the false positive detection random effects
    double *deltaStarObs = (double *) R_alloc(nObs, sizeof(double)); 
    zeros(deltaStarObs, nObs); 
    int *deltaStarLongIndx = (int *) R_alloc(nObspFPDetRE, sizeof(int));
    // Get sums of the current REs for each site/visit combo
    for (i = 0; i < nObs; i++) {
      for (l = 0; l < pFPDetRE; l++) {
        deltaStarLongIndx[l * nObs + i] = which(XomegaRE[l * nObs + i], deltaLevelIndx, nFPDetRE);
        deltaStarObs[i] += deltaStar[deltaStarLongIndx[l * nObs + i]];
      }
    }
    // Starting index for occurrence random effects
    int *betaStarStart = (int *) R_alloc(pOccRE, sizeof(int)); 
    for (l = 0; l < pOccRE; l++) {
      betaStarStart[l] = which(l, betaStarIndx, nOccRE); 
    }
    // Starting index for true positive detection random effects
    int *gammaStarStart = (int *) R_alloc(pTPDetRE, sizeof(int)); 
    for (l = 0; l < pTPDetRE; l++) {
      gammaStarStart[l] = which(l, gammaStarIndx, nTPDetRE); 
    }
    // Starting index for false positive detection random effects
    int *deltaStarStart = (int *) R_alloc(pFPDetRE, sizeof(int)); 
    for (l = 0; l < pFPDetRE; l++) {
      deltaStarStart[l] = which(l, deltaStarIndx, nFPDetRE); 
    }

    /**********************************************************************
     * Set up MH stuff
     *********************************************************************/
    double logPostGammaCurr = 0.0, logPostGammaCand = 0.0;
    double logPostDeltaCurr = 0.0, logPostDeltaCand = 0.0;
    double *logPostGammaStarCand = (double *) R_alloc(nTPDetRE, sizeof(double));
    double *logPostGammaStarCurr = (double *) R_alloc(nTPDetRE, sizeof(double));
    for (j = 0; j < nTPDetRE; j++) {
      logPostGammaStarCurr[j] = R_NegInf;
      logPostGammaStarCand[j] = logPostGammaStarCurr[j];
    }
    double *logPostDeltaStarCand = (double *) R_alloc(nFPDetRE, sizeof(double));
    double *logPostDeltaStarCurr = (double *) R_alloc(nFPDetRE, sizeof(double));
    for (j = 0; j < nFPDetRE; j++) {
      logPostDeltaStarCurr[j] = R_NegInf;
      logPostDeltaStarCand[j] = logPostDeltaStarCurr[j];
    }
    int nAMCMC = pTPDet + pFPDet;
    if (pTPDetRE > 0) {
      nAMCMC += nTPDetRE;
    }
    if (pFPDetRE > 0) {
      nAMCMC += nFPDetRE;
    }
    int gammaAMCMCIndx = 0;
    int deltaAMCMCIndx = gammaAMCMCIndx + pTPDet;
    int gammaStarAMCMCIndx = deltaAMCMCIndx + pFPDet;
    int deltastarAMCMCIndx = gammaStarAMCMCIndx + nTPDetRE;
    double *accept = (double *) R_alloc(nAMCMC, sizeof(double)); zeros(accept, nAMCMC);
    // Set the initial candidate values for everything to the inital values.
    double *gammaCand = (double *) R_alloc(pTPDet, sizeof(double));
    for (j = 0; j < pTPDet; j++) {
      gammaCand[j] = gamma[j];
    }
    double *deltaCand = (double *) R_alloc(pFPDet, sizeof(double));
    for (j = 0; j < pFPDet; j++) {
      deltaCand[j] = delta[j];
    }
    double *gammaStarCand = (double *) R_alloc(nTPDetRE, sizeof(double));
    for (j = 0; j < nTPDetRE; j++) {
      gammaStarCand[j] = gammaStar[j];
    }
    double *deltaStarCand = (double *) R_alloc(nFPDetRE, sizeof(double));
    for (j = 0; j < nFPDetRE; j++) {
      deltaStarCand[j] = deltaStar[j];
    }
    SEXP acceptSamples_r;
    PROTECT(acceptSamples_r = Rf_allocMatrix(REALSXP, nAMCMC, nBatch)); nProtect++;
    zeros(REAL(acceptSamples_r), nAMCMC * nBatch);
    SEXP tuningSamples_r;
    PROTECT(tuningSamples_r = Rf_allocMatrix(REALSXP, nAMCMC, nBatch)); nProtect++;
    zeros(REAL(tuningSamples_r), nAMCMC * nBatch);

    // Initialize omega for use in gamma update.
    for (i = 0; i < nObs; i++) {
      omega[i] = exp(F77_NAME(ddot)(&pFPDet, &Xomega[i], &nObs, delta, &inc) + deltaStarObs[i]);
    }

    GetRNGstate(); 

    for (s = 0, q = 0; s < nBatch; s++) {
      for (r = 0; r < batchLength; r++, q++) {
        /********************************************************************
         *Update Occupancy Auxiliary Variables 
         *******************************************************************/
        for (j = 0; j < J; j++) {
          omegaPGOcc[j] = rpg(1.0, F77_NAME(ddot)(&pOcc, &X[j], &J, beta, &inc) + betaStarSites[j]);
        } // j

        /********************************************************************
         *Update Occupancy Regression Coefficients
         *******************************************************************/
        for (j = 0; j < J; j++) {
          kappaOcc[j] = z[j] - 1.0 / 2.0; 
          tmp_J1[j] = kappaOcc[j] - omegaPGOcc[j] * betaStarSites[j]; 
        } // j

        /********************************
         * Compute b.beta
         *******************************/
        F77_NAME(dgemv)(ytran, &J, &pOcc, &one, X, &J, tmp_J1, &inc, &zero, tmp_pOcc, &inc FCONE); 	 
        for (j = 0; j < pOcc; j++) {
          tmp_pOcc[j] += SigmaBetaInvMuBeta[j]; 
        } // j 

        /********************************
         * Compute A.beta
         * *****************************/
        // tmp_JpOcc is X %*% omegaPGOcc. 
        for(j = 0; j < J; j++){
          for(i = 0; i < pOcc; i++){
            tmp_JpOcc[i*J+j] = X[i*J+j]*omegaPGOcc[j];
          }
        }

        // This finishes off A.beta
        // 1 * X * tmp_JpOcc + 0 * tmp_ppOcc = tmp_ppOcc
        F77_NAME(dgemm)(ytran, ntran, &pOcc, &pOcc, &J, &one, X, &J, tmp_JpOcc, &J, &zero, tmp_ppOcc, &pOcc FCONE FCONE);
        for (j = 0; j < ppOcc; j++) {
          tmp_ppOcc[j] += SigmaBetaInv[j]; 
        } // j

        F77_NAME(dpotrf)(lower, &pOcc, tmp_ppOcc, &pOcc, &info FCONE); 
        if(info != 0){Rf_error("c++ error: dpotrf here failed\n");}
        F77_NAME(dpotri)(lower, &pOcc, tmp_ppOcc, &pOcc, &info FCONE); 
        if(info != 0){Rf_error("c++ error: dpotri here failed\n");}
        F77_NAME(dsymv)(lower, &pOcc, &one, tmp_ppOcc, &pOcc, tmp_pOcc, &inc, &zero, tmp_pOcc2, &inc FCONE);
        F77_NAME(dpotrf)(lower, &pOcc, tmp_ppOcc, &pOcc, &info FCONE); 
        if(info != 0){Rf_error("c++ error: dpotrf here failed\n");}
        mvrnorm(beta, tmp_pOcc2, tmp_ppOcc, pOcc);
        
        /********************************************************************
         *Update True Positive Detection Regression Coefficients 
         *******************************************************************/
        // Proposal
        for (k = 0; k < pTPDet; k++) {
          logPostGammaCand = 0.0;
          logPostGammaCurr = 0.0;
          gammaCand[k] = rnorm(gamma[k], exp(tuning[gammaAMCMCIndx + k]));
          logPostGammaCand += dnorm(gammaCand[k], muGamma[k], sqrt(SigmaGamma[k * pTPDet + k]), 1);
          logPostGammaCurr += dnorm(gamma[k], muGamma[k], sqrt(SigmaGamma[k * pTPDet + k]), 1);
          for (i = 0; i < nObs; i++) {
            // TP covariates only inform when the site is actually occupied. 
            if (z[zLongIndx[i]] == 1) {
              // Just mu
              tmp_nObs[i] = exp(F77_NAME(ddot)(&pTPDet, &Xmu[i], &nObs, gammaCand, &inc) + 
                                gammaStarObs[i]);
              // mu + gamma
              tmp_nObs2[i] = tmp_nObs[i] + omega[i]; 
              logPostGammaCand += poisson_logpost(y[i], tmp_nObs2[i], 1.0);
              logPostGammaCand += Rf_dbinom(d[i], nVal[i], tmp_nObs[i] / tmp_nObs2[i], 1);


              tmp_nObs[i] = exp(F77_NAME(ddot)(&pTPDet, &Xmu[i], &nObs, gamma, &inc) + 
                                gammaStarObs[i]);
              // mu + gamma
              tmp_nObs2[i] = tmp_nObs[i] + omega[i]; 
              logPostGammaCurr += poisson_logpost(y[i], tmp_nObs2[i], 1.0);
              logPostGammaCurr += Rf_dbinom(d[i], nVal[i], tmp_nObs[i] / tmp_nObs2[i], 1);
            }
          }
          if (runif(0.0, 1.0) <= exp(logPostGammaCand - logPostGammaCurr)) {
            gamma[k] = gammaCand[k];
            accept[gammaAMCMCIndx + k]++;
          } else {
            gammaCand[k] = gamma[k];
          }
        }
        // Put updated values in mu for use in other updates.
        for (i =0; i < nObs; i++) {
          mu[i] = exp(F77_NAME(ddot)(&pTPDet, &Xmu[i], &nObs, gamma, &inc) + 
                            gammaStarObs[i]); 
        }
        

        /********************************************************************
         *Update False Positive Detection Regression Coefficients 
         *******************************************************************/
        // Proposal
        for (k = 0; k < pFPDet; k++) {
          logPostDeltaCand = 0.0;
          logPostDeltaCurr = 0.0;
          deltaCand[k] = rnorm(delta[k], exp(tuning[deltaAMCMCIndx + k]));
          logPostDeltaCand += dnorm(deltaCand[k], muDelta[k], sqrt(SigmaDelta[k * pFPDet + k]), 1);
          logPostDeltaCurr += dnorm(delta[k], muDelta[k], sqrt(SigmaDelta[k * pFPDet + k]), 1);
          for (i = 0; i < nObs; i++) {
            // Notice all observations are used here since FPs can occur 
            // everywhere, not just at occupied sites. 
            tmp_nObs[i] = exp(F77_NAME(ddot)(&pFPDet, &Xomega[i], &nObs, deltaCand, &inc) + deltaStarObs[i]) + 
                          mu[i] * z[zLongIndx[i]];
            logPostDeltaCand += poisson_logpost(y[i], tmp_nObs[i], 1.0);
            logPostDeltaCand += Rf_dbinom(d[i], nVal[i], mu[i] * z[zLongIndx[i]] / tmp_nObs[i], 1);
            
            tmp_nObs[i] = exp(F77_NAME(ddot)(&pFPDet, &Xomega[i], &nObs, delta, &inc) + deltaStarObs[i]) + 
                          mu[i] * z[zLongIndx[i]];
            logPostDeltaCurr += poisson_logpost(y[i], tmp_nObs[i], 1.0);
            logPostDeltaCurr += Rf_dbinom(d[i], nVal[i], mu[i] * z[zLongIndx[i]] / tmp_nObs[i], 1);

          }
          if (runif(0.0, 1.0) <= exp(logPostDeltaCand - logPostDeltaCurr)) {
            delta[k] = deltaCand[k];
            accept[deltaAMCMCIndx + k]++;
          } else {
            deltaCand[k] = delta[k];
          }
        }
        // Put updated values in omega for use in other updates.
        for (i = 0; i < nObs; i++) {
          omega[i] = exp(F77_NAME(ddot)(&pFPDet, &Xomega[i], &nObs, delta, &inc) + deltaStarObs[i]);
        }

        /********************************************************************
         *Update Occupancy random effects variance
         *******************************************************************/
        for (l = 0; l < pOccRE; l++) {
          tmp_0 = F77_NAME(ddot)(&nOccRELong[l], &betaStar[betaStarStart[l]], &inc, &betaStar[betaStarStart[l]], &inc); 
          tmp_0 *= 0.5; 
          sigmaSqPsi[l] = rigamma(sigmaSqPsiA[l] + nOccRELong[l] / 2.0, sigmaSqPsiB[l] + tmp_0); 
        }

        /********************************************************************
         *Update True Positive Detection Random Effects Variance
         *******************************************************************/
        for (l = 0; l < pTPDetRE; l++) {
          tmp_0 = F77_NAME(ddot)(&nTPDetRELong[l], &gammaStar[gammaStarStart[l]], &inc, &gammaStar[gammaStarStart[l]], &inc); 
          tmp_0 *= 0.5; 
          sigmaSqMu[l] = rigamma(sigmaSqMuA[l] + nTPDetRELong[l] / 2.0, sigmaSqMuB[l] + tmp_0); 
        }
        
        /********************************************************************
         *Update False Positive Detection Random Effects Variance
         *******************************************************************/
        for (l = 0; l < pFPDetRE; l++) {
          tmp_0 = F77_NAME(ddot)(&nFPDetRELong[l], &deltaStar[deltaStarStart[l]], &inc, &deltaStar[deltaStarStart[l]], &inc); 
          tmp_0 *= 0.5; 
          sigmaSqOmega[l] = rigamma(sigmaSqOmegaA[l] + nFPDetRELong[l] / 2.0, sigmaSqOmegaB[l] + tmp_0); 
        }

        /********************************************************************
         *Update Occupancy random effects
         *******************************************************************/
        if (pOccRE > 0) {
          // Update each individual random effect one by one. 
          for (l = 0; l < nOccRE; l++) {
            /********************************
             * Compute b.beta.star
             *******************************/
            zeros(tmp_one, inc);
            tmp_0 = 0.0;	      
            // Only allow information to come from when XRE == betaLevelIndx[l]. 
            // aka information only comes from the sites with any given level 
            // of a random effect. 
            for (j = 0; j < J; j++) {
              if (XRE[betaStarIndx[l] * J + j] == betaLevelIndx[l]) {
                tmp_02 = 0.0;
                for (ll = 0; ll < pOccRE; ll++) {
                  tmp_02 += betaStar[betaStarLongIndx[ll * J + j]];
	        } 
                tmp_one[0] += kappaOcc[j] - (F77_NAME(ddot)(&pOcc, &X[j], &J, beta, &inc) + 
          		    tmp_02 - betaStar[l]) * omegaPGOcc[j];
                tmp_0 += omegaPGOcc[j];
              }
            }
            /********************************
             * Compute A.beta.star
             *******************************/
            tmp_0 += 1.0 / sigmaSqPsi[betaStarIndx[l]]; 
            tmp_0 = 1.0 / tmp_0; 
            betaStar[l] = rnorm(tmp_0 * tmp_one[0], sqrt(tmp_0)); 
          }
          // Update the RE sums for the current species
          zeros(betaStarSites, J);
          for (j = 0; j < J; j++) {
            for (l = 0; l < pOccRE; l++) {
              betaStarSites[j] += betaStar[betaStarLongIndx[l * J + j]];
            }
          }
        }

        /********************************************************************
         *Update True Positive Detection random effects
         *******************************************************************/
        if (pTPDetRE > 0) {
          // Update each random effect level one at a time.
          for (l = 0; l < nTPDetRE; l++) {
            logPostGammaStarCand[l] = 0.0;
            logPostGammaStarCurr[l] = 0.0;
            gammaStarCand[l] = rnorm(gammaStar[l], exp(tuning[gammaStarAMCMCIndx + l]));
            logPostGammaStarCand[l] += dnorm(gammaStarCand[l], 0.0,
                                             sqrt(sigmaSqMu[gammaStarIndx[l]]), 1);
            logPostGammaStarCurr[l] += dnorm(gammaStar[l], 0.0,
                                             sqrt(sigmaSqMu[gammaStarIndx[l]]), 1);
            // Only observations at this level of the random effect, and only
            // at occupied sites, carry information about gammaStar[l].
            for (i = 0; i < nObs; i++) {
              if ((XmuRE[gammaStarIndx[l] * nObs + i] == gammaLevelIndx[l]) &&
                  (z[zLongIndx[i]] == 1)) {
                // mu[i] holds the current TP rate. Swap gammaStar[l] for the
                // candidate on the log scale.
                tmp_nObs[i] = mu[i] * exp(gammaStarCand[l] - gammaStar[l]);
                tmp_nObs2[i] = tmp_nObs[i] + omega[i];
                logPostGammaStarCand[l] += poisson_logpost(y[i], tmp_nObs2[i], 1.0);
                logPostGammaStarCand[l] += Rf_dbinom(d[i], nVal[i], tmp_nObs[i] / tmp_nObs2[i], 1);

                tmp_nObs2[i] = mu[i] + omega[i];
                logPostGammaStarCurr[l] += poisson_logpost(y[i], tmp_nObs2[i], 1.0);
                logPostGammaStarCurr[l] += Rf_dbinom(d[i], nVal[i], mu[i] / tmp_nObs2[i], 1);
              }
            }
            if (runif(0.0, 1.0) <= exp(logPostGammaStarCand[l] - logPostGammaStarCurr[l])) {
              // Keep the observation-level RE sums and mu in sync so later
              // levels (possibly from other RE factors) see the new value.
              tmp_0 = gammaStarCand[l] - gammaStar[l];
              for (i = 0; i < nObs; i++) {
                if (XmuRE[gammaStarIndx[l] * nObs + i] == gammaLevelIndx[l]) {
                  gammaStarObs[i] += tmp_0;
                  mu[i] *= exp(tmp_0);
                }
              }
              gammaStar[l] = gammaStarCand[l];
              accept[gammaStarAMCMCIndx + l]++;
            } else {
              gammaStarCand[l] = gammaStar[l];
            }
          }
        }

        /********************************************************************
         *Joint shift of TP intercept and TP random effects
         *******************************************************************/
        // gamma[0] + c and gammaStar - c (for all levels of one RE factor)
        // leaves every mu[i] unchanged, so the likelihood cancels and c is
        // drawn exactly from its Gaussian conditional given the priors
        // (generalized Gibbs move; Liu and Sabatti 2000).
        if (pTPDetRE > 0 && tpDetInt) {
          for (ll = 0; ll < pTPDetRE; ll++) {
            tmp_0 = 0.0;
            for (l = gammaStarStart[ll]; l < gammaStarStart[ll] + nTPDetRELong[ll]; l++) {
              tmp_0 += gammaStar[l];
            }
            tmp_02 = 1.0 / (1.0 / SigmaGamma[0] + nTPDetRELong[ll] / sigmaSqMu[ll]);
            tmp_02 = rnorm(tmp_02 * ((muGamma[0] - gamma[0]) / SigmaGamma[0] +
                                     tmp_0 / sigmaSqMu[ll]), sqrt(tmp_02));
            gamma[0] += tmp_02;
            gammaCand[0] = gamma[0];
            for (l = gammaStarStart[ll]; l < gammaStarStart[ll] + nTPDetRELong[ll]; l++) {
              gammaStar[l] -= tmp_02;
              gammaStarCand[l] = gammaStar[l];
            }
            for (i = 0; i < nObs; i++) {
              gammaStarObs[i] -= tmp_02;
            }
          }
        }

        /********************************************************************
         *Update False Positive Detection random effects
         *******************************************************************/
        if (pFPDetRE > 0) {
          // Update each random effect level one at a time.
          for (l = 0; l < nFPDetRE; l++) {
            logPostDeltaStarCand[l] = 0.0;
            logPostDeltaStarCurr[l] = 0.0;
            deltaStarCand[l] = rnorm(deltaStar[l], exp(tuning[deltastarAMCMCIndx + l]));
            logPostDeltaStarCand[l] += dnorm(deltaStarCand[l], 0.0,
                                             sqrt(sigmaSqOmega[deltaStarIndx[l]]), 1);
            logPostDeltaStarCurr[l] += dnorm(deltaStar[l], 0.0,
                                             sqrt(sigmaSqOmega[deltaStarIndx[l]]), 1);
            // All observations at this level of the random effect carry
            // information, since FPs can occur at both occupied and unoccupied
            // sites.
            for (i = 0; i < nObs; i++) {
              if (XomegaRE[deltaStarIndx[l] * nObs + i] == deltaLevelIndx[l]) {
                // omega[i] holds the current FP rate. Swap deltaStar[l] for the
                // candidate on the log scale.
                tmp_nObs[i] = omega[i] * exp(deltaStarCand[l] - deltaStar[l]) +
                              mu[i] * z[zLongIndx[i]];
                logPostDeltaStarCand[l] += poisson_logpost(y[i], tmp_nObs[i], 1.0);
                logPostDeltaStarCand[l] += Rf_dbinom(d[i], nVal[i], mu[i] * z[zLongIndx[i]] / tmp_nObs[i], 1);

                tmp_nObs[i] = omega[i] + mu[i] * z[zLongIndx[i]];
                logPostDeltaStarCurr[l] += poisson_logpost(y[i], tmp_nObs[i], 1.0);
                logPostDeltaStarCurr[l] += Rf_dbinom(d[i], nVal[i], mu[i] * z[zLongIndx[i]] / tmp_nObs[i], 1);
              }
            }
            if (runif(0.0, 1.0) <= exp(logPostDeltaStarCand[l] - logPostDeltaStarCurr[l])) {
              // Keep the observation-level RE sums and omega in sync so later
              // levels (possibly from other RE factors) see the new value.
              tmp_0 = deltaStarCand[l] - deltaStar[l];
              for (i = 0; i < nObs; i++) {
                if (XomegaRE[deltaStarIndx[l] * nObs + i] == deltaLevelIndx[l]) {
                  deltaStarObs[i] += tmp_0;
                  omega[i] *= exp(tmp_0);
                }
              }
              deltaStar[l] = deltaStarCand[l];
              accept[deltastarAMCMCIndx + l]++;
            } else {
              deltaStarCand[l] = deltaStar[l];
            }
          }
        }

        /********************************************************************
         *Joint shift of FP intercept and FP random effects
         *******************************************************************/
        // Same move as for the TP random effects: omega[i] is unchanged.
        if (pFPDetRE > 0 && fpDetInt) {
          for (ll = 0; ll < pFPDetRE; ll++) {
            tmp_0 = 0.0;
            for (l = deltaStarStart[ll]; l < deltaStarStart[ll] + nFPDetRELong[ll]; l++) {
              tmp_0 += deltaStar[l];
            }
            tmp_02 = 1.0 / (1.0 / SigmaDelta[0] + nFPDetRELong[ll] / sigmaSqOmega[ll]);
            tmp_02 = rnorm(tmp_02 * ((muDelta[0] - delta[0]) / SigmaDelta[0] +
                                     tmp_0 / sigmaSqOmega[ll]), sqrt(tmp_02));
            delta[0] += tmp_02;
            deltaCand[0] = delta[0];
            for (l = deltaStarStart[ll]; l < deltaStarStart[ll] + nFPDetRELong[ll]; l++) {
              deltaStar[l] -= tmp_02;
              deltaStarCand[l] = deltaStar[l];
            }
            for (i = 0; i < nObs; i++) {
              deltaStarObs[i] -= tmp_02;
            }
          }
        }

        /********************************************************************
         *Update Latent Occupancy
         *******************************************************************/
        // Accumulate the site-level log-likelihood of (y, d) under each value
        // of z. When z = 1: y ~ Poisson(mu + omega), d ~ Bin(n, mu / (mu + omega)).
        // When z = 0: y ~ Poisson(omega), and d must be 0 (Bin(d; n, 0)).
        // Visits with y = 0 are included, as they are evidence for absence.
        for (j = 0; j < J; j++) {
          logLikeOcc[j] = 0.0;
          logLikeUnocc[j] = 0.0;
          dSum[j] = 0.0;
        } // j
        for (i = 0; i < nObs; i++) {
          tmp_0 = mu[i] + omega[i];
          logLikeOcc[zLongIndx[i]] += poisson_logpost(y[i], tmp_0, 1.0) +
                                      Rf_dbinom(d[i], nVal[i], mu[i] / tmp_0, 1);
          logLikeUnocc[zLongIndx[i]] += poisson_logpost(y[i], omega[i], 1.0);
          dSum[zLongIndx[i]] += d[i];
        } // i

        // Sample z and compute the z-marginalized likelihood for WAIC
        for (j = 0; j < J; j++) {
          logitPsi = F77_NAME(ddot)(&pOcc, &X[j], &J, beta, &inc) + betaStarSites[j];
          psi[j] = logitInv(logitPsi, zero, one);
          if (dSum[j] > 0.0) {
            // A validated true positive means the site is occupied.
            z[j] = one;
            yWAIC[j] = psi[j] * exp(logLikeOcc[j]);
          } else {
            // P(z = 1 | .) = psi L1 / (psi L1 + (1 - psi) L0), which on the
            // log scale is logitInv(logit(psi) + log L1 - log L0). Working
            // with the difference avoids underflow of L1 and L0.
            // L1 is the likelihood if site j is occupied. L0 is the likelihood
            // if site j is unoccupied. 
            z[j] = rbinom(one, logitInv(logitPsi + logLikeOcc[j] - logLikeUnocc[j],
                                        zero, one));
            // log(psi L1 + (1 - psi) L0) via log-sum-exp
            tmp_0 = log(psi[j]) + logLikeOcc[j];
            tmp_02 = log(1.0 - psi[j]) + logLikeUnocc[j];
            logLikeMax = std::max(tmp_0, tmp_02);
            yWAIC[j] = exp(logLikeMax + log(exp(tmp_0 - logLikeMax) +
                                            exp(tmp_02 - logLikeMax)));
          }
        } // j

        /********************************************************************
         *Save samples
         *******************************************************************/
        if (q >= nBurn) {
          thinIndx++; 
          if (thinIndx == nThin) {
            F77_NAME(dcopy)(&pOcc, beta, &inc, &REAL(betaSamples_r)[sPost*pOcc], &inc);
            F77_NAME(dcopy)(&pTPDet, gamma, &inc, &REAL(gammaSamples_r)[sPost*pTPDet], &inc);
            F77_NAME(dcopy)(&pFPDet, delta, &inc, &REAL(deltaSamples_r)[sPost*pFPDet], &inc);
            F77_NAME(dcopy)(&J, psi, &inc, &REAL(psiSamples_r)[sPost*J], &inc); 
            F77_NAME(dcopy)(&J, z, &inc, &REAL(zSamples_r)[sPost*J], &inc); 
            if (pOccRE > 0) {
              F77_NAME(dcopy)(&pOccRE, sigmaSqPsi, &inc, 
          		    &REAL(sigmaSqPsiSamples_r)[sPost*pOccRE], &inc);
              F77_NAME(dcopy)(&nOccRE, betaStar, &inc, 
          		    &REAL(betaStarSamples_r)[sPost*nOccRE], &inc);
            }
            if (pTPDetRE > 0) {
              F77_NAME(dcopy)(&pTPDetRE, sigmaSqMu, &inc, 
          		    &REAL(sigmaSqMuSamples_r)[sPost*pTPDetRE], &inc);
              F77_NAME(dcopy)(&nTPDetRE, gammaStar, &inc, 
          		    &REAL(gammaStarSamples_r)[sPost*nTPDetRE], &inc);
            }
            if (pFPDetRE > 0) {
              F77_NAME(dcopy)(&pFPDetRE, sigmaSqOmega, &inc, 
          		    &REAL(sigmaSqOmegaSamples_r)[sPost*pFPDetRE], &inc);
              F77_NAME(dcopy)(&nFPDetRE, deltaStar, &inc, 
          		    &REAL(deltaStarSamples_r)[sPost*nFPDetRE], &inc);
            }
            F77_NAME(dcopy)(&J, yWAIC, &inc, 
          		  &REAL(likeSamples_r)[sPost*J], &inc);
            sPost++; 
            thinIndx = 0; 
          }
        }
        R_CheckUserInterrupt();
      } // r (end batch)

      /********************************************************************
       *Adjust tuning
       *******************************************************************/
      for (j = 0; j < nAMCMC; j++) {
        REAL(acceptSamples_r)[s * nAMCMC + j] = accept[j]/batchLength;
        REAL(tuningSamples_r)[s * nAMCMC + j] = tuning[j];
        if (accept[j] / batchLength > acceptRate) {
          tuning[j] += std::min(0.01, 1.0/sqrt(static_cast<double>(s)));
        } else{
            tuning[j] -= std::min(0.01, 1.0/sqrt(static_cast<double>(s)));
          }
        accept[j] = 0;
      }

      /********************************************************************
       *Report
       *******************************************************************/
      if (verbose) {
        if (status == nReport) {
          Rprintf("Batch: %i of %i, %3.2f%%\n", s, nBatch, 100.0*s/nBatch);
          Rprintf("\tParameter\tAcceptance\tTuning\n");
          for (j = 0; j < pTPDet; j++) {
            Rprintf("\tgamma[%i]\t\t%3.1f\t\t%1.5f\n", j + 1, 100.0*REAL(acceptSamples_r)[s * nAMCMC + gammaAMCMCIndx + j], exp(tuning[gammaAMCMCIndx + j]));
          }
          for (j = 0; j < pFPDet; j++) {
            Rprintf("\tdelta[%i]\t\t%3.1f\t\t%1.5f\n", j + 1, 100.0*REAL(acceptSamples_r)[s * nAMCMC + deltaAMCMCIndx + j], exp(tuning[deltaAMCMCIndx + j]));
          }
          Rprintf("-------------------------------------------------\n");
          #ifdef Win32
          R_FlushConsole();
          #endif
          status = 0;
        }
      }
      status++;

    } // s all batches

    if (verbose) {
      Rprintf("Batch: %i of %i, %3.2f%%\n", s, nBatch, 100.0*s/nBatch);
    }
    PutRNGstate();

    SEXP result_r, resultName_r;
    int nResultListObjs = 12;

    PROTECT(result_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;

    // Setting the components of the output list.
    SET_VECTOR_ELT(result_r, 0, betaSamples_r);
    SET_VECTOR_ELT(result_r, 1, gammaSamples_r);
    SET_VECTOR_ELT(result_r, 2, deltaSamples_r);
    SET_VECTOR_ELT(result_r, 3, zSamples_r); 
    SET_VECTOR_ELT(result_r, 4, psiSamples_r);
    SET_VECTOR_ELT(result_r, 5, likeSamples_r);
    if (pTPDetRE > 0) {
      SET_VECTOR_ELT(result_r, 6, sigmaSqMuSamples_r);
      SET_VECTOR_ELT(result_r, 7, gammaStarSamples_r);
    }
    if (pOccRE > 0) {
      SET_VECTOR_ELT(result_r, 8, sigmaSqPsiSamples_r);
      SET_VECTOR_ELT(result_r, 9, betaStarSamples_r);
    }
    if (pFPDetRE > 0) {
      SET_VECTOR_ELT(result_r, 10, sigmaSqOmegaSamples_r);
      SET_VECTOR_ELT(result_r, 11, deltaStarSamples_r);
    }

    SET_VECTOR_ELT(resultName_r, 0, Rf_mkChar("beta.samples")); 
    SET_VECTOR_ELT(resultName_r, 1, Rf_mkChar("gamma.samples")); 
    SET_VECTOR_ELT(resultName_r, 2, Rf_mkChar("delta.samples")); 
    SET_VECTOR_ELT(resultName_r, 3, Rf_mkChar("z.samples")); 
    SET_VECTOR_ELT(resultName_r, 4, Rf_mkChar("psi.samples"));
    SET_VECTOR_ELT(resultName_r, 5, Rf_mkChar("like.samples"));
    if (pTPDetRE > 0) {
      SET_VECTOR_ELT(resultName_r, 6, Rf_mkChar("sigma.sq.mu.samples")); 
      SET_VECTOR_ELT(resultName_r, 7, Rf_mkChar("gamma.star.samples")); 
    }
    if (pOccRE > 0) {
      SET_VECTOR_ELT(resultName_r, 8, Rf_mkChar("sigma.sq.psi.samples")); 
      SET_VECTOR_ELT(resultName_r, 9, Rf_mkChar("beta.star.samples")); 
    }
    if (pFPDetRE > 0) {
      SET_VECTOR_ELT(resultName_r, 10, Rf_mkChar("sigma.sq.omega.samples")); 
      SET_VECTOR_ELT(resultName_r, 11, Rf_mkChar("delta.star.samples")); 
    }
   
    Rf_namesgets(result_r, resultName_r);
    
    UNPROTECT(nProtect);
    
    return(result_r);
  }
}



