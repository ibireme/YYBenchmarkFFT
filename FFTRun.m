//
//  FFTProfile.m
//  TestFFT
//
//  Created by ibireme on 14/12/2.
//  Copyright (c) 2014 ibireme. All rights reserved.
//

#import "FFTRun.h"
#import <Accelerate/Accelerate.h>

#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>


#include "fftw3.h"
#include "kiss_fft.h"
#include "NE10_fft.h"
#include "ckfft.h"
#include "pffft.h"
#include "fftn.h"
#include "nsfft.h"
#include "myfft.h"


static inline void ProfileTime(void (^block)(void), void (^complete)(double ms)) {
    struct timeval t0, t1;
    gettimeofday(&t0, NULL);
    block();
    gettimeofday(&t1, NULL);
    double ms = (double)(t1.tv_sec - t0.tv_sec) * 1e3 + (double)(t1.tv_usec - t0.tv_usec) * 1e-3;
    complete(ms);
}


bool eq_(float a, float b, int ulps) {
    if (fabsf(a) < 0.001 && fabsf(b) < 0.001) return true;
    
    assert(ulps > 0 && ulps < 4 * 1024 * 1024);
    int ai = *(int*)&a;
    if (ai < 0) ai = 0x80000000 - ai;
    int bi = *(int*)&b;
    if (bi < 0) bi = 0x80000000 - bi;
    int diff = abs(ai - bi);
    if (diff <= ulps) return true;
    //printf("diff:%d ai:%d  bi:%d\n",diff, ai, bi);
    return false;
}

bool eq(float a, float b) {
    bool isEqual = eq_(a, b, 1 << 20);
    if (!isEqual) {
        printf("%f, %f\n", a,b);
    }
    return isEqual;
}


@implementation FFTRun


- (void)run {
    
    
#ifdef __ARM_NEON__
#if ULONG_MAX > 0xFFFFFFFF
    char *env = "arm64";
#else
    char *env = "armv7";
#endif
#else
    #if ULONG_MAX > 0xFFFFFFFF
    char *env = "x86_64";
#else
    char *env = "i386";
#endif
#endif
    
    printf("run fft in %s\n",env);
    printf("----------------\n");
    
    for (int i = 4; i < 16; i++) {
        int length = pow(2, i);
        [self run:length];
    }
}


- (void)run:(int)length {
    printf("%d-------\n",length);
    long repeat = 16777216 / length;
    
    
    double *data_real = (double *)calloc(length, sizeof(double));
    double *data_imag = (double *)calloc(length, sizeof(double));
    {   /// generate audio data
        double freq, rate, amp, phase, dc;
        freq = 100; rate = 44100; amp = 2; phase = M_PI / 4; dc = -0.5;
        for (int i = 0; i < length; i++) {
            data_real[i] += cos(i * (2 * M_PI) / (rate / freq) + phase) * amp + dc;
        }
        freq = 500; rate = 44100; amp = 1; phase = M_PI / 3; dc = 1;
        for (int i = 0; i < length; i++) {
            data_real[i] += cos(i * (2 * M_PI) / (rate / freq) + phase) * amp + dc;
        }
        freq = 10000; rate = 44100; amp = 0.5; phase = M_PI / 2; dc = 0;
        for (int i = 0; i < length; i++) {
            data_real[i] += cos(i * (2 * M_PI) / (rate / freq) + phase) * amp + dc;
        }
        freq = 20000; rate = 44100; amp = 0.3; phase = M_PI / 9; dc = 0;
        for (int i = 0; i < length; i++) {
            data_real[i] += cos(i * (2 * M_PI) / (rate / freq) + phase) * amp + dc;
        }
    }
    
    double *idft_real = (double *)calloc(length, sizeof(double));
    double *idft_imag = (double *)calloc(length, sizeof(double));
    mydftd(data_real, data_imag, idft_real, idft_imag, length, -1);
    
    
    { /// myfft
        double *in = (double *)calloc(length * 2, sizeof(double));
        
        { // validate
            for (int i = 0; i < length; i++) {
                in[i*2] = data_real[i];
                in[i*2+1] = data_imag[i];
            }
            myfftd(in, length, -1);
            for (int i = 0; i < length; i++) {
                if (!eq(in[i*2], idft_real[i])) {printf("data err!\n"); break;}
                if (!eq(in[i*2+1], idft_imag[i])) {printf("data err!\n"); break;}
            }
            myfftd(in, length, 1);
            for (int i = 0; i < length; i++) {
                if (!eq(in[i*2]/length, data_real[i])) {printf("data err!\n"); break;}
                if (!eq(in[i*2+1]/length, data_imag[i])) {printf("data err!\n"); break;}
            }
        }
        
        { // profile
            ProfileTime(^{
                for (int r = 0 ; r < repeat; r++) {
                    for (int i = 0; i < length; i++) {
                        in[i*2] = data_real[i];
                        in[i*2+1] = data_imag[i];
                    }
                    myfftd(in, length, -1);
                    myfftd(in, length, 1);
                }
            }, ^(double ms) {
                printf("myfft:%10d ms  %10lldM/s\n",(int)ms, (uint64_t)(1000/ms/1024.0/1024.0 * length * repeat));
            });
        }
        
        free(in);
    }
    
    
    { /// fftn
        float *in_real = (float *)calloc(length, sizeof(float));
        float *in_imag = (float *)calloc(length, sizeof(float));
        
        {// validate
            for (int i = 0; i < length; i++) {
                in_real[i] = data_real[i];
                in_imag[i] = data_imag[i];
            }
            
            int dim[] = {length};
            fftnf(1, dim, in_real, in_imag, -1, 1);
            for (int i = 0; i < length; i++) {
                if (!eq(in_real[i], idft_real[i])) {printf("data err!\n"); break;}
                if (!eq(in_imag[i], idft_imag[i])) {printf("data err!\n"); break;}
            }
            fftnf(1, dim, in_real, in_imag, 1, 1);
            for (int i = 0; i < length; i++) {
                if (!eq(in_real[i] / length, data_real[i])) {printf("data err!\n"); break;}
                if (!eq(in_imag[i] / length, data_imag[i])) {printf("data err!\n"); break;}
            }
        }
        
        { //profile
            ProfileTime(^{
                for (int r = 0 ; r < repeat; r++) {
                    for (int i = 0; i < length; i++) {
                        in_real[i] = data_real[i];
                        in_imag[i] = data_imag[i];
                    }
                    int dim[] = {length};
                    fftnf(1, dim, in_real, in_imag, -1, 1);
                    fftnf(1, dim, in_real, in_imag, 1, 1);
                }
            }, ^(double ms) {
                printf("fftn: %10d ms  %10lldM/s\n",(int)ms, (uint64_t)(1000/ms/1024.0/1024.0 * length * repeat));
            });
        }
        
        free(in_real);
        free(in_imag);
    }
    
    
    { /// kissfft
        kiss_fft_cpx *in = (kiss_fft_cpx *)calloc(length, sizeof(kiss_fft_cpx));
        kiss_fft_cpx *out = (kiss_fft_cpx *)calloc(length, sizeof(kiss_fft_cpx));
        
        kiss_fft_cfg cfg = kiss_fft_alloc(length, 0, NULL, NULL);
        kiss_fft_cfg icfg = kiss_fft_alloc(length, 1, NULL, NULL);
        
        { // validate
            for (int i = 0; i < length; i++) {
                in[i].r = data_real[i];
                in[i].i = data_imag[i];
            }
            
            kiss_fft(cfg, in, out);
            for (int i = 0; i < length; i++) {
                if (!eq(out[i].r, idft_real[i])) {printf("data err!\n"); break;}
                if (!eq(out[i].i, idft_imag[i])) {printf("data err!\n"); break;}
            }
            kiss_fft(icfg, out, in);
            for (int i = 0; i < length; i++) {
                if (!eq(in[i].r / length, data_real[i])) {printf("data err!\n"); break;}
                if (!eq(in[i].i / length, data_imag[i])) {printf("data err!\n"); break;}
            }
        }
        
        { // profile
            ProfileTime(^{
                for (int r = 0 ; r < repeat; r++) {
                    for (int i = 0; i < length; i++) {
                        in[i].r = data_real[i];
                        in[i].i = data_imag[i];
                    }
                    kiss_fft(cfg, in, out);
                    kiss_fft(icfg, out, in);
                }
            }, ^(double ms) {
                printf("kiss: %10d ms  %10lldM/s\n",(int)ms, (uint64_t)(1000/ms/1024.0/1024.0 * length * repeat));
            });
        }
        
        free(cfg);
        free(icfg);
        
        free(in);
        free(out);
    }
    
    { /// nsfft
        int mode = SIMDBase_chooseBestMode(SIMDBase_TYPE_FLOAT);
        //int veclen = SIMDBase_getModeParamInt(SIMDBase_PARAMID_VECTOR_LEN, mode); //4
        int sizeOfVect = SIMDBase_getModeParamInt(SIMDBase_PARAMID_SIZE_OF_VECT, mode); // 16
        
        DFT *cfg = DFT_init(mode, length, 0);
        float *in = (float *)SIMDBase_alignedMalloc(sizeOfVect * length * 2);
        memset(in, 0, sizeOfVect * length * 2);
        
        { // validate
            for (int i = 0; i < length; i++) {
                in[(i * 2) * sizeOfVect / sizeof(float)] = data_real[i];
                in[(i * 2 + 1) * sizeOfVect / sizeof(float)] = data_imag[i];
            }
            
            DFT_execute(cfg, mode, in, -1);
            
            for (int i = 0; i < length; i++) {
                if (!eq(in[(i * 2) * sizeOfVect / sizeof(float)], idft_real[i])) {printf("data err!\n"); break;}
                if (!eq(in[(i * 2 + 1) * sizeOfVect / sizeof(float)], idft_imag[i])) {printf("data err!\n"); break;}
            }
            
            DFT_execute(cfg, mode, in, 1);
            for (int i = 0; i < length; i++) {
                if (!eq(in[(i * 2) * sizeOfVect / sizeof(float)] / length, data_real[i])) {printf("data err!\n"); break;}
                if (!eq(in[(i * 2 + 1) * sizeOfVect / sizeof(float)] / length, data_imag[i])) {printf("data err!\n"); break;}
            }
        }
        
        { // profile
            ProfileTime(^{
                for (int k = 0 ; k < repeat; k++) {
                    for (int i = 0; i < length; i++) {
                        in[(i * 2) * sizeOfVect / sizeof(float)] = data_real[i];
                        in[(i * 2 + 1) * sizeOfVect / sizeof(float)] = data_imag[i];
                    }
                    DFT_execute(cfg, mode, in, -1);
                    DFT_execute(cfg, mode, in, 1);
                }
            }, ^(double ms) {
                printf("nsfft:%10d ms  %10lldM/s\n",(int)ms, (uint64_t)(1000/ms/1024.0/1024.0 * length * repeat));
            });
        }
        
        free(cfg);
        free(in);
    }
    
    { /// pffft
        float *in = (float *)calloc(length * 2, sizeof(float));
        float *out = (float *)calloc(length * 2, sizeof(float));
        PFFFT_Setup *setup = pffft_new_setup(length, PFFFT_COMPLEX); // valid length >= 16
        
        { // validate
            for (int i = 0; i < length; i++) {
                in[i*2] = data_real[i];
                in[i*2+1] = data_imag[i];
            }
            pffft_transform_ordered(setup, in, out, NULL, PFFFT_FORWARD);
            for (int i = 0; i < length; i++) {
                if (!eq(out[i*2], idft_real[i])) {printf("data err!\n"); break;}
                if (!eq(out[i*2+1], idft_imag[i])) {printf("data err!\n"); break;}
            }
            pffft_transform_ordered(setup, out, in, NULL, PFFFT_BACKWARD);
            for (int i = 0; i < length; i++) {
                if (!eq(in[i*2]/length, data_real[i])) {printf("data err!\n"); break;}
                if (!eq(in[i*2+1]/length, data_imag[i])) {printf("data err!\n"); break;}
            }
        }
        
        { // profile
            ProfileTime(^{
                for (int p = 0 ; p < repeat; p++) {
                    for (int i = 0; i < length; i++) {
                        in[i*2] = data_real[i];
                        in[i*2+1] = data_imag[i];
                    }
                    pffft_transform_ordered(setup, in, out, NULL, PFFFT_FORWARD);
                    pffft_transform_ordered(setup, out, in, NULL, PFFFT_BACKWARD);
                }
            }, ^(double ms) {
                printf("pffft:%10d ms  %10lldM/s\n",(int)ms, (uint64_t)(1000/ms/1024.0/1024.0 * length * repeat));
            });
        }
        
        pffft_destroy_setup(setup);
        free(in);
        free(out);
    }
    
    
    
    { /// ckfft
        
        CkFftComplex* in = (CkFftComplex*)calloc(length, sizeof(CkFftComplex));
        CkFftComplex* out = (CkFftComplex*)calloc(length, sizeof(CkFftComplex));
        CkFftComplex* buf = (CkFftComplex*)calloc(length, sizeof(CkFftComplex));
        CkFftContext* context = CkFftInit(length, kCkFftDirection_Both, NULL, NULL);
        
        { // validate
            for (int i = 0; i < length; i++) {
                in[i].real = data_real[i];
                in[i].imag = data_imag[i];
            }
            CkFftComplexForward(context, length, in, out);
            for (int i = 0; i < length; i++) {
                if (!eq(out[i].real, idft_real[i])) {printf("data err!\n"); break;}
                if (!eq(out[i].imag, idft_imag[i])) {printf("data err!\n"); break;}
            }
            CkFftComplexInverse(context, length, out, in);
            for (int i = 0; i < length; i++) {
                if (!eq(in[i].real / length, data_real[i])) {printf("data err!\n"); break;}
                if (!eq(in[i].imag / length, data_imag[i])) {printf("data err!\n"); break;}
            }
        }
        
        { // profile
            ProfileTime(^{
                for (int p = 0 ; p < repeat; p++) {
                    for (int i = 0; i < length; i++) {
                        in[i].real = data_real[i];
                        in[i].imag = data_imag[i];
                    }
                    CkFftComplexForward(context, length, in, out);
                    CkFftComplexInverse(context, length, out, in);
                }
            }, ^(double ms) {
                printf("ckfft:%10d ms  %10lldM/s\n",(int)ms, (uint64_t)(1000/ms/1024.0/1024.0 * length * repeat));
            });
        }
        
        CkFftShutdown(context);

        free(in);
        free(out);
        free(buf);
    }
    
    
    { /// fftw
        fftwf_complex *in = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * length);
        fftwf_complex *out = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * length);
        fftwf_plan plan = fftwf_plan_dft_1d(length, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
        fftwf_plan iplan = fftwf_plan_dft_1d(length, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);
        
        { // validate
            for (int i = 0; i < length; i++) {
                in[i][0] = data_real[i];
                in[i][1] = data_imag[i];
            }
            fftwf_execute(plan);
            for (int i = 0; i < length; i++) {
                if (!eq(out[i][0], idft_real[i])) {printf("data err!\n"); break;}
                if (!eq(out[i][1], idft_imag[i])) {printf("data err!\n"); break;}
            }
            fftwf_execute(iplan);
            for (int i = 0; i < length; i++) {
                if (!eq(in[i][0]/length, data_real[i])) {printf("data err!\n"); break;}
                if (!eq(in[i][1]/length, data_imag[i])) {printf("data err!\n"); break;}
            }
        }
        
        { // profile
            ProfileTime(^{
                for (int r = 0 ; r < repeat; r++) {
                    for (int i = 0; i < length; i++) {
                        in[i][0] = data_real[i];
                        in[i][1] = data_imag[i];
                    }
                    fftwf_execute(plan);
                    fftwf_execute(iplan);
                }
            }, ^(double ms) {
                printf("fftw: %10d ms  %10lldM/s\n",(int)ms, (uint64_t)(1000/ms/1024.0/1024.0 * length * repeat));
            });
        }
        
        fftwf_destroy_plan(plan);
        fftwf_destroy_plan(iplan);
        
        fftwf_free(in);
        fftwf_free(out);
    }
    
    
    { /// ne10
        ne10_fft_cpx_float32_t *in = (ne10_fft_cpx_float32_t *)calloc(length, sizeof(ne10_fft_cpx_float32_t));
        ne10_fft_cpx_float32_t *out = (ne10_fft_cpx_float32_t *)calloc(length, sizeof(ne10_fft_cpx_float32_t));
        ne10_fft_cfg_float32_t cfg = ne10_fft_alloc_c2c_float32(length);
        
        { // validate
            for (int i = 0; i < length; i++) {
                in[i].r = data_real[i];
                in[i].i = data_imag[i];
            }
            ne10_fft_c2c_1d_float32(out, in, cfg, 0);
            for (int i = 0; i < length; i++) {
                if (!eq(out[i].r, idft_real[i])) {printf("data err!\n"); break;}
                if (!eq(out[i].i, idft_imag[i])) {printf("data err!\n"); break;}
            }
            ne10_fft_c2c_1d_float32(in, out, cfg, 1);
            for (int i = 0; i < length; i++) {
                if (!eq(in[i].r, data_real[i])) {printf("data err!\n"); break;}
                if (!eq(in[i].i, data_imag[i])) {printf("data err!\n"); break;}
            }
        }
        
        { // profile
            ProfileTime(^{
                for (int r = 0 ; r < repeat; r++) {
                    for (int i = 0; i < length; i++) {
                        in[i].r = data_real[i];
                        in[i].i = data_imag[i];
                    }
                    ne10_fft_c2c_1d_float32(out, in, cfg, 0);
                    ne10_fft_c2c_1d_float32(in, out, cfg, 1);
                }
            }, ^(double ms) {
                printf("ne10: %10d ms  %10lldM/s\n",(int)ms, (uint64_t)(1000/ms/1024.0/1024.0 * length * repeat));
            });
        }
        
        free(cfg);
        free(in);
        free(out);
    }
    
    
    { /// vdsp out-place cpx
        COMPLEX *in = (COMPLEX *)calloc(length, sizeof(COMPLEX));
        for (int i = 0; i < length; i++) {
            in[i].real = data_real[i];
            in[i].imag = data_imag[i];
        }
        
        uint32_t log2n = log2f((float)length);
        FFTSetup setup = vDSP_create_fftsetup(log2n, FFT_RADIX2);
        
        COMPLEX_SPLIT cpx;
        cpx.realp = (float *)malloc(sizeof(float) * length);
        cpx.imagp = (float *)malloc(sizeof(float) * length);
        
        COMPLEX_SPLIT icpx;
        icpx.realp = (float *)malloc(sizeof(float) * length);
        icpx.imagp = (float *)malloc(sizeof(float) * length);
        
        { // validate
            vDSP_ctoz(in, 2, &cpx, 1, length);
            vDSP_fft_zop(setup, &cpx, 1, &icpx, 1, log2n, FFT_FORWARD);
            for (int i = 0 ; i < length; i++) {
                if (!eq(icpx.realp[i], idft_real[i])) {printf("data err!\n"); break;}
                if (!eq(icpx.imagp[i], idft_imag[i])) {printf("data err!\n"); break;}
            }
            vDSP_fft_zop(setup, &icpx, 1, &cpx, 1, log2n, FFT_INVERSE);
            for (int i = 0 ; i < length; i++) {
                if (!eq(cpx.realp[i] / length, data_real[i])) {printf("data err!\n"); break;}
                if (!eq(cpx.imagp[i] / length, data_imag[i])) {printf("data err!\n"); break;}
            }
        }
        
        ProfileTime ( ^{
            for (int p = 0; p < repeat; p++) {
                vDSP_ctoz(in, 2, &cpx, 1, length);
                vDSP_fft_zop(setup, &cpx, 1, &icpx, 1, log2n, FFT_FORWARD);
                vDSP_fft_zop(setup, &icpx, 1, &cpx, 1, log2n, FFT_INVERSE);
            }
        }, ^(double ms) {
            printf("vdsp: %10d ms  %10lldM/s  (out-place cpx)\n", (int)ms, (uint64_t)(1000 / ms / 1024.0 /1024.0 * length * repeat));
        });
        
        free(in);
        free(cpx.realp);
        free(cpx.imagp);
        free(icpx.realp);
        free(icpx.imagp);
        vDSP_destroy_fftsetup(setup);
    }
    
    
    
    { /// vdsp in-place cpx with tmp buf
        COMPLEX *in = (COMPLEX *)calloc(length, sizeof(COMPLEX));
        for (int i = 0; i < length; i++) {
            in[i].real = data_real[i];
            in[i].imag = data_imag[i];
        }
        
        uint32_t log2n = log2f((float)length);
        FFTSetup setup = vDSP_create_fftsetup(log2n, FFT_RADIX2);
        
        COMPLEX_SPLIT cpx;
        cpx.realp = (float *)malloc(sizeof(float) * length);
        cpx.imagp = (float *)malloc(sizeof(float) * length);
        
        COMPLEX_SPLIT tmp;
        tmp.realp = (float *)malloc(sizeof(float) * length);
        tmp.imagp = (float *)malloc(sizeof(float) * length);
        
        { // validate
            vDSP_ctoz(in, 2, &cpx, 1, length);
            vDSP_fft_zipt(setup, &cpx, 1, &tmp, log2n, FFT_FORWARD);
            for (int i = 0 ; i < length; i++) {
                if (!eq(cpx.realp[i], idft_real[i])) {printf("data err!\n"); break;}
                if (!eq(cpx.imagp[i], idft_imag[i])) {printf("data err!\n"); break;}
            }
            vDSP_fft_zipt(setup, &cpx, 1, &tmp, log2n, FFT_INVERSE);
            for (int i = 0 ; i < length; i++) {
                if (!eq(cpx.realp[i] / length, data_real[i])) {printf("data err!\n"); break;}
                if (!eq(cpx.imagp[i] / length, data_imag[i])) {printf("data err!\n"); break;}
            }
        }
        
        ProfileTime ( ^{
            for (int p = 0; p < repeat; p++) {
                vDSP_ctoz(in, 2, &cpx, 1, length);
                vDSP_fft_zipt(setup, &cpx, 1, &tmp, log2n, FFT_FORWARD);
                vDSP_fft_zipt(setup, &cpx, 1, &tmp, log2n, FFT_INVERSE);
            }
        }, ^(double ms) {
            printf("vdsp: %10d ms  %10lldM/s  (in-place cpx with tmp buf)\n", (int)ms, (uint64_t)(1000 / ms / 1024.0/1024.0 * length * repeat));
        });
        
        free(in);
        free(cpx.realp);
        free(cpx.imagp);
        free(tmp.realp);
        free(tmp.imagp);
        vDSP_destroy_fftsetup(setup);
    }
    
    { /// vdsp in-place cpx
        COMPLEX *in = (COMPLEX *)calloc(length, sizeof(COMPLEX));
        for (int i = 0; i < length; i++) {
            in[i].real = data_real[i];
            in[i].imag = data_imag[i];
        }
        
        uint32_t log2n = log2f((float)length);
        FFTSetup setup = vDSP_create_fftsetup(log2n, FFT_RADIX2);
        
        COMPLEX_SPLIT A;
        A.realp = (float *)malloc(sizeof(float) * length);
        A.imagp = (float *)malloc(sizeof(float) * length);
        
        { // validate
            vDSP_ctoz(in, 2, &A, 1, length);
            vDSP_fft_zip(setup, &A, 1, log2n, FFT_FORWARD);
            for (int i = 0 ; i < length; i++) {
                if (!eq(A.realp[i], idft_real[i])) {printf("data err!\n"); break;}
                if (!eq(A.imagp[i], idft_imag[i])) {printf("data err!\n"); break;}
            }
            vDSP_fft_zip(setup, &A, 1, log2n, FFT_INVERSE);
            for (int i = 0 ; i < length; i++) {
                if (!eq(A.realp[i] / length, data_real[i])) {printf("data err!\n"); break;}
                if (!eq(A.imagp[i] / length, data_imag[i])) {printf("data err!\n"); break;}
            }
        }
        
        ProfileTime ( ^{
            for (int p = 0; p < repeat; p++) {
                vDSP_ctoz(in, 2, &A, 1, length);
                vDSP_fft_zip(setup, &A, 1, log2n, FFT_FORWARD);
                vDSP_fft_zip(setup, &A, 1, log2n, FFT_INVERSE);
            }
        }, ^(double ms) {
            printf("vdsp: %10d ms  %10lldM/s  (in-place cpx)\n", (int)ms, (uint64_t)(1000 / ms / 1024.0/1024.0 * length * repeat));
        });
        
        free(in);
        free(A.realp);
        free(A.imagp);
        vDSP_destroy_fftsetup(setup);
    }
    
    
    
    
    { /// vdsp in-place real
        float *in = (float *)calloc(length, sizeof(float));
        for (int i = 0; i < length; i++) {
            in[i] = data_real[i];
        }
        
        float *buf = (float *)calloc(length, sizeof(float));
        COMPLEX_SPLIT scpx = {buf, buf + length / 2};
        
        uint32_t log2n = log2f((float)length);
        FFTSetup setup = vDSP_create_fftsetup(log2n, FFT_RADIX2);
        
        { // validate
            vDSP_ctoz((COMPLEX *)in, 2, &scpx, 1, length / 2);
            vDSP_fft_zrip(setup, &scpx, 1, log2n, FFT_FORWARD);
            float scale = 0.5;
            vDSP_vsmul(buf, 1, &scale, buf, 1, length);
            if (!eq(scpx.realp[0], idft_real[0]) ) {printf("data err!\n");}
            if (!eq(scpx.imagp[0], idft_real[length/2]) ) {printf("data err!\n");}
            for (int i = 1; i < length / 2; i++) {
                if (!eq(scpx.realp[i], idft_real[i])) {printf("data err!\n"); break;}
                if (!eq(scpx.imagp[i], idft_imag[i])) {printf("data err!\n"); break;}
            }
            
            vDSP_fft_zrip(setup, &scpx, 1, log2n, FFT_INVERSE);
            vDSP_ztoc(&scpx, 1, (COMPLEX *)in, 2, length / 2);
            
            for (int i = 0; i < length; i++) {
                if (!eq(in[i] / length, data_real[i])) {printf("data err!\n"); break;}
            }
        }
        
        { // profile
            ProfileTime(^{
                for (int p = 0 ; p < repeat; p++) {
                    vDSP_ctoz((COMPLEX *)in, 2, &scpx, 1, length / 2); // set value
                    vDSP_fft_zrip(setup, &scpx, 1, log2n, FFT_FORWARD); // fft
                    float scale = 0.5;
                    vDSP_vsmul(buf, 1, &scale, buf, 1, length); // fix scale
                    vDSP_fft_zrip(setup, &scpx, 1, log2n, FFT_INVERSE); // ifft
                    vDSP_ztoc(&scpx, 1, (COMPLEX *)in, 2, length / 2); // copy out
                }
            }, ^(double ms) {
                printf("vdsp: %10d ms  %10lldM/s  (in-place real)\n",(int)ms, (uint64_t)(1000/ms/1024.0/1024.0 * length * repeat));
            });
        }
        
        free(in);
        free(buf);
        vDSP_destroy_fftsetup(setup);
    }

    
    
    free(data_real);
    free(data_imag);
    free(idft_real);
    free(idft_imag);
    
}


@end
