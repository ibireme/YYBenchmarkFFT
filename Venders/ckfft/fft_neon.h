#pragma once
#include "ckfft.h"


namespace ckfft
{

void fft_neon(
        CkFftContext* context, 
        const CkFftComplex* input, 
        CkFftComplex* output, 
        int count, 
        bool inverse,
        int stride, 
        const CkFftComplex* expTable,
        int expTableStride);

}


