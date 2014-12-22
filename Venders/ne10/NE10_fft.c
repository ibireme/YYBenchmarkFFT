//
//  NE10_fft.c
//  TestFFT
//
//  Created by ibireme on 14/11/25.
//  Copyright (c) 2014 ibireme. All rights reserved.
//

#include "NE10_fft.h"

extern void ne10_fft_c2c_1d_float32_c(ne10_fft_cpx_float32_t *fout,
                                      ne10_fft_cpx_float32_t *fin,
                                      ne10_fft_cfg_float32_t  cfg,
                                      ne10_int32_t            inverse_fft);

extern void ne10_fft_r2c_1d_float32_c(ne10_fft_cpx_float32_t *   fout,
                                      ne10_float32_t *           fin,
                                      ne10_fft_r2c_cfg_float32_t cfg);

extern void ne10_fft_c2r_1d_float32_c(ne10_float32_t *           fout,
                                      ne10_fft_cpx_float32_t *   fin,
                                      ne10_fft_r2c_cfg_float32_t cfg);


#ifdef __ARM_NEON__
extern void ne10_fft_c2c_1d_float32_neon(ne10_fft_cpx_float32_t *fout,
                                         ne10_fft_cpx_float32_t *fin,
                                         ne10_fft_cfg_float32_t  cfg,
                                         ne10_int32_t            inverse_fft);

extern void ne10_fft_r2c_1d_float32_neon(ne10_fft_cpx_float32_t *   fout,
                                         ne10_float32_t *           fin,
                                         ne10_fft_r2c_cfg_float32_t cfg);

extern void ne10_fft_c2r_1d_float32_neon(ne10_float32_t *           fout,
                                         ne10_fft_cpx_float32_t *   fin,
                                         ne10_fft_r2c_cfg_float32_t cfg);

#endif



void ne10_fft_c2c_1d_float32(ne10_fft_cpx_float32_t *fout,
                             ne10_fft_cpx_float32_t *fin,
                             ne10_fft_cfg_float32_t  cfg,
                             ne10_int32_t            inverse_fft) {
#ifdef __ARM_NEON__
    ne10_fft_c2c_1d_float32_neon(fout, fin, cfg, inverse_fft);
#else
    ne10_fft_c2c_1d_float32_c(fout, fin, cfg, inverse_fft);
#endif
}

void ne10_fft_r2c_1d_float32(ne10_fft_cpx_float32_t *   fout,
                             ne10_float32_t *           fin,
                             ne10_fft_r2c_cfg_float32_t cfg) {
#ifdef __ARM_NEON__
    ne10_fft_r2c_1d_float32_neon(fout, fin, cfg);
#else
    ne10_fft_r2c_1d_float32_c(fout, fin, cfg);
#endif
}

void ne10_fft_c2r_1d_float32(ne10_float32_t *           fout,
                             ne10_fft_cpx_float32_t *   fin,
                             ne10_fft_r2c_cfg_float32_t cfg) {
#ifdef __ARM_NEON__
    ne10_fft_c2r_1d_float32_neon(fout, fin, cfg);
#else
    ne10_fft_c2r_1d_float32_c(fout, fin, cfg);
#endif
}
