#include <stdlib.h>
#include <math.h>
#include "mex.h"
#include "matrix.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double *InData, *tic, *init_double;
    int NRow, NCol, Ntic, Ndim;
    double lambda1, lambda2, lambda3, mu1, mu2, mu3, INIT;

    int *Freq0, *Freq1, *Freq2, *fundM, *init;
    double *FreqOut1, *FreqOut2, *FreqOut0;
    // double *FValOut;
    double *Energy, *FVal;
    double minval, val;

    double multi1, multi2, multi3;
    double highfund, lowfund;
    int lowIdx = 0, highIdx = 0, bw2, bw3, bw1;
    int i,j,k,l;

    /* Portal to matlab */
    InData = mxGetPr(prhs[0]);
    NRow = mxGetM(prhs[0]);
    NCol = mxGetN(prhs[0]);

    tic = mxGetPr(prhs[1]);
    Ntic = mxGetM(prhs[1]);

    lambda1 = mxGetScalar(prhs[2]);
    lambda2 = mxGetScalar(prhs[3]);
    lambda3 = mxGetScalar(prhs[4]);

    mu1 = mxGetScalar(prhs[5]);
    mu2 = mxGetScalar(prhs[6]);

    multi1 = mxGetScalar(prhs[7]);
    multi2 = mxGetScalar(prhs[8]);

    lowfund = mxGetScalar(prhs[9]);
    highfund = mxGetScalar(prhs[10]);
    for (k = Ntic-1; k >= 1; --k) {
        // search for the true indices for multiple frequency in tfrtic
        if (tic[k]>=lowfund && tic[k-1]<lowfund)
            lowIdx = k;
        if (tic[k]>=highfund && tic[k-1]<highfund)
            highIdx = k;
        if ((lowIdx != 0) && (highIdx != 0))
            break;
    }
    bw1 = highIdx - lowIdx + 1;
    bw2 = mxGetScalar(prhs[11]);
    bw3 = mxGetScalar(prhs[12]);

    INIT = mxGetScalar(prhs[13]);

    init_double = mxGetPr(prhs[14]);
    Ndim = mxGetM(prhs[14]);

    double FUND;    //1.0 or 0.0: Consider the fundamental's energy or not
    FUND = mxGetScalar(prhs[15]);

    init = (int *)mxMalloc(sizeof(int)*(Ndim));
    for (i = 0; i < Ndim; ++i) {
        init[i] = (int) init_double[i];
        init[i] -= 1;   //Matlab starts from 1
    }
    printf("\n");
    
    plhs[0] = mxCreateNumericMatrix(NRow, 1, mxDOUBLE_CLASS, mxREAL);
    plhs[1] = mxCreateNumericMatrix(NRow, 1, mxDOUBLE_CLASS, mxREAL);
    plhs[2] = mxCreateNumericMatrix(NRow, 1, mxDOUBLE_CLASS, mxREAL);

    FreqOut1 = mxGetPr(plhs[1]);
    FreqOut2 = mxGetPr(plhs[2]);
    FreqOut0 = mxGetPr(plhs[0]);
    printf("NRow, NCol = %d, %d\n", NRow, NCol);

    /* Main operations start here */
    const double eps = 1e-8;
    double sum = 0, tmp;

    printf("BW1, BW2, BW3, INIT, low, high = %d, %d, %d, %d, %d, %d\n", bw1, bw2, bw3, INIT, lowIdx, highIdx);
    printf("Multiples are %.1f, %.1f, %.1f\n", multi1, multi2, multi3);

    Energy = (double *)mxMalloc(sizeof(double)*NRow*NCol);
    FVal = (double *)mxMalloc(sizeof(double)*NRow*(bw1)*(2*bw2)*(2*bw3));
    Freq0 = (int *)mxMalloc(sizeof(int)*NRow);
    Freq1 = (int *)mxMalloc(sizeof(int)*NRow);
    Freq2 = (int *)mxMalloc(sizeof(int)*NRow);

    // BaseLine array for the multiples: fundM
    fundM = (int *)mxMalloc(sizeof(int)*(2*bw1));
    double freqVal1, freqVal2;
    for (int j1 = 0; j1 < bw1; ++j1) {
        freqVal1 = tic[lowIdx+j1]*multi1;
        freqVal2 = tic[lowIdx+j1]*multi2;
        fundM[j1] = Ntic-1; fundM[bw1+j1] = Ntic-1;
        for (k = Ntic-1; k >= 1; --k) {
            // search for the true indices for multiple frequency in tfrtic
            if (tic[k]>=freqVal1 && tic[k-1]<freqVal1)
                fundM[j1] = k;
            if (tic[k]>=freqVal2 && tic[k-1]<freqVal2)
                fundM[bw1+j1] = k;
        }
    }

    // TFR matrix
    for (i=0; i<NRow; ++i) {
        for (j=0; j<NCol; ++j) {
            Energy[i*NCol+j] = InData[i+j*NRow];
            sum += Energy[i*NCol+j];
        }
    }
    for (i=0;i<NRow;++i) {
        for (j=0;j<NCol;++j)
            Energy[i*NCol+j] = -log(Energy[i*NCol+j]/sum+eps);
    }

    // Pre-allocation
    // for (i=0; i<NRow; ++i)
    
    for (int j1 = 0; j1 < bw1; ++j1) {
        for (int j2 = 0; j2 < 2*bw2; ++j2) {
            for (int j3 = 0; j3 < 2*bw3; ++j3) {
                FVal[j1*(2*bw2*2*bw3) + j2*(2*bw3) + j3] = FUND*Energy[j1+lowIdx]\
                 + Energy[j2 + fundM[j1] - bw2]\
                 + Energy[j3 + fundM[bw1+j1] - bw3]\
                 + INIT * (lambda1 * (j1+lowIdx-init[0]) * (j1+lowIdx-init[0])\
                 + lambda2 * ((j2+fundM[j1]-bw2)-init[1]) * ((j2+fundM[j1]-bw2)-init[1])\
                 + lambda3 * ((j3+fundM[bw1+j1]-bw3)-init[2]) * ((j3+fundM[bw1+j1]-bw3)-init[2])\
                 + mu1 * fabs(multi1*tic[j1+lowIdx]-tic[j2+fundM[j1]-bw2]) * fabs(multi1*tic[j1+lowIdx]-tic[j2+fundM[j1]-bw2])\
                 + mu2 * fabs(multi2*tic[j1+lowIdx]-tic[j3+fundM[bw1+j1]-bw3]) * fabs(multi2*tic[j1+lowIdx]-tic[j3+fundM[bw1+j1]-bw3]));
            }
        }
    }

    int j2val, j3val, k2val, k3val;
    double k1freq, k2freq;
    for (i = 1; i < floor(NRow); ++i) {
        for (int j1 = 0; j1 < bw1; ++j1) {
            for (int j2 = 0; j2 < 2*bw2; ++j2) {
                j2val = j2 + fundM[j1];
                for (int j3 = 0; j3 < 2*bw3; ++j3) {
                    minval = 1e16; //FVal[i*(bw1*2*bw2*2*bw3)+j1*(2*bw2*2*bw3)+j2*2*bw3+j3];
                    j3val = j3 + fundM[bw1+j1];
                    for (int k1 = 0; k1 < bw1; ++k1) {
                        k1freq = tic[k1+lowIdx];
                        k2val = fundM[k1];
                        k3val = fundM[bw1+k1];
                        for (int k2 = 0; k2 < 2*bw2; ++k2) {
                            k2freq = tic[k2+k2val-bw2];
                            for (int k3 = 0; k3 < 2*bw3; ++k3) {
                                tmp = FVal[(i-1)*(bw1*2*bw2*2*bw3)+k1*(2*bw2*2*bw3)+k2*(2*bw3)+k3]\
                                 + lambda1 * (k1-j1) * (k1-j1)\
                                 + lambda2 * ((k2+k2val)-j2val) * ((k2+k2val)-j2val)\
                                 + lambda3 * ((k3+k3val)-j3val) * ((k3+k3val)-j3val)\
                                 + mu1 * fabs(multi1*k1freq-k2freq) * fabs(multi1*k1freq-k2freq)\
                                 + mu2 * fabs(multi2*k1freq-tic[k3+k3val-bw3]) * fabs(multi2*k1freq-tic[k3+k3val-bw3]);
                                if(tmp < minval)
                                    minval = tmp;
                            }
                        }
                    }
                    FVal[i*(bw1*2*bw2*2*bw3)+j1*(2*bw2*2*bw3)+j2*(2*bw3)+j3] = minval\
                     + FUND*Energy[i*NCol+(j1+lowIdx)]\
                     + Energy[i*NCol+(j2+fundM[j1]-bw2)]\
                     + Energy[i*NCol+(j3+fundM[bw1+j1]-bw3)];
                }
            }
        }
    }

    minval = FVal[(NRow-1)*(bw1*2*bw2*2*bw3)];
    Freq0[NRow-1] = 0;
    Freq1[NRow-1] = 0 + fundM[0] - bw2;
    Freq2[NRow-1] = 0 + fundM[bw1] - bw3;
    // Freq3[NRow-1] = 0 + fundM[3*NRow-1] - bw3;
    for (int j1 = 0; j1 < bw1; ++j1) {
        for (int j2 = 0; j2 < 2*bw2; ++j2) {
            for (int j3 = 0; j3 < 2*bw3; ++j3) {
                tmp = FVal[(NRow-1)*(bw1*2*bw2*2*bw3) + j1*(2*bw2*2*bw3) + j2*(2*bw3) + j3];
                if (tmp < minval) {
                    minval = tmp;
                    // Freq1[NRow-1] = j1 + fundM[NRow-1] - bw1;
                    Freq0[NRow-1] = j1;
                    Freq1[NRow-1] = j2 + fundM[j1] - bw2;
                    Freq2[NRow-1] = j3 + fundM[bw1+j1] - bw3;
                }
            }
        }
    }

    for (i = NRow-2; i >= 0; --i) {
        val = FVal[(i+1)*(bw1*2*bw2*2*bw3) + Freq0[i+1]*(2*bw2*2*bw3)\
         + (Freq1[i+1]-fundM[Freq0[i+1]]+bw2)*(2*bw3) + (Freq2[i+1]-fundM[bw1+Freq0[i+1]]+bw3)]\
         - FUND*Energy[(i+1)*NCol+Freq0[i+1]+lowIdx] - Energy[(i+1)*NCol+Freq1[i+1]] - Energy[(i+1)*NCol+Freq2[i+1]];

        for (int j1 = 0; j1 < bw1; ++j1) {
            for (int j2 = 0; j2 < 2*bw2; ++j2) {
                for (int j3 = 0; j3 < 2*bw3; ++j3) {
                    tmp = FVal[i*(bw1*2*bw2*2*bw3) + j1*(2*bw2*2*bw3) + j2*(2*bw3) + j3]\
                     + lambda1 * (j1-Freq0[i+1]) * (j1-Freq0[i+1])\
                     + lambda2 * ((j2+fundM[j1]-bw2) - Freq1[i+1]) * ((j2+fundM[j1]-bw2) - Freq1[i+1])\
                     + lambda3 * ((j3+fundM[bw1+j1]-bw3) - Freq2[i+1]) * ((j3+fundM[bw1+j1]-bw3) - Freq2[i+1])\
                     + mu1 * fabs(multi1*tic[j1+lowIdx]-tic[j2+fundM[j1]-bw2]) * fabs(multi1*tic[j1+lowIdx]-tic[j2+fundM[j1]-bw2])\
                     + mu2 * fabs(multi2*tic[j1+lowIdx]-tic[j3+fundM[bw1+j1]-bw3]) * fabs(multi2*tic[j1+lowIdx]-tic[j3+fundM[bw1+j1]-bw3]);
                    if(fabs(val-tmp) < eps) {
                        Freq0[i] = j1;
                        Freq1[i] = j2 + fundM[j1] - bw2;
                        Freq2[i] = j3 + fundM[bw1+j1] - bw3;
                        j1 = bw1; j2 = 2*bw2; j3 = 2*bw3; // get out of the loops
                    }
                }
            }
        }
    }

    for (i=0; i<NRow; ++i) {
        FreqOut0[i] = (double)Freq0[i]+lowIdx+1.;
        FreqOut1[i] = (double)Freq1[i]+1.;
        FreqOut2[i] = (double)Freq2[i]+1.;
    }

    mxFree(Freq1); mxFree(Freq2); mxFree(Freq0);
    mxFree(fundM);
    mxFree(Energy);
    mxFree(FVal); 
    return;
}