/*
 *  Copyright 2013-14 ARM Limited
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *    * Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimer.
 *    * Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions and the following disclaimer in the
 *      documentation and/or other materials provided with the distribution.
 *    * Neither the name of ARM Limited nor the
 *      names of its contributors may be used to endorse or promote products
 *      derived from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY ARM LIMITED AND CONTRIBUTORS "AS IS" AND
 *  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 *  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 *  DISCLAIMED. IN NO EVENT SHALL ARM LIMITED BE LIABLE FOR ANY
 *  DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 *  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 *  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 *  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 *  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/*
 * NE10 Library : dsp/NE10_fft.h
 */


#include "NE10_types.h"

#ifndef NE10_FFT_H
#define NE10_FFT_H

#ifdef __cplusplus
extern "C" {
#endif
    
#define NE10_FFT_BYTE_ALIGNMENT 8
    

    extern ne10_fft_cfg_float32_t ne10_fft_alloc_c2c_float32(ne10_int32_t nfft);
    
    extern ne10_fft_r2c_cfg_float32_t ne10_fft_alloc_r2c_float32(ne10_int32_t nfft);
    
    
    extern void ne10_fft_c2c_1d_float32(ne10_fft_cpx_float32_t *fout,
                                        ne10_fft_cpx_float32_t *fin,
                                        ne10_fft_cfg_float32_t  cfg,
                                        ne10_int32_t            inverse_fft);
    
    extern void ne10_fft_r2c_1d_float32(ne10_fft_cpx_float32_t *   fout,
                                        ne10_float32_t *           fin,
                                        ne10_fft_r2c_cfg_float32_t cfg);
    
    extern void ne10_fft_c2r_1d_float32(ne10_float32_t *           fout,
                                        ne10_fft_cpx_float32_t *   fin,
                                        ne10_fft_r2c_cfg_float32_t cfg);
    
#ifdef __cplusplus
}
#endif

#endif
