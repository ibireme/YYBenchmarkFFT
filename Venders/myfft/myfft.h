//
//  myfft.h
//  TestFFT
//
//  Created by ibireme on 14/11/26.
//  Copyright (c) 2014 ibireme. All rights reserved.
//

/*
 My dft/fft.
 Simple code, low performance.
 */

#ifndef __TestFFT__myfft__
#define __TestFFT__myfft__

#ifdef __cplusplus
extern "C" {
#endif
    
    
#include <stdint.h>
#include <stdbool.h>
    
    
    
    /**
     fft
     @param sign   forward:-1, inverse:1
     */
    void mydft(float *in_r, float *in_i, float *out_r, float *out_i, int length, int sign);
    void mydftd(double *in_r, double *in_i, double *out_r, double *out_i, int length, int sign);
    
    
    /**
     @param cpx    [real0,imag0,real1,imag1...realN,imagN]
     @param sign         forward:-1, inverse:1
     */
    void myfft(float *cpx, int length, int sign);
    void myfftd(double *cpx, int length, int sign);
    
    
#ifdef __cplusplus
}
#endif

#endif /* defined(__TestFFT__yyfft__) */
