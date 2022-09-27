/**
 * @file BatchMandelCalculator.cc
 * @author René Rešetár <xreset00@stud.fit.vutbr.cz>
 * @brief Implementation of Mandelbrot calculator that uses SIMD paralelization over small batches
 * @date 14.11.2021
 */

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include <stdlib.h>
#include <stdexcept>
#include <random>

#include "BatchMandelCalculator.h"

BatchMandelCalculator::BatchMandelCalculator (unsigned matrixBaseSize, unsigned limit) :
	BaseMandelCalculator(matrixBaseSize, limit, "BatchMandelCalculator")
{
    data = (int *)(malloc(height * width * sizeof(int)));
    zReal = (float *)(_mm_malloc(width * sizeof(float), 256));
    zImag = (float *)(_mm_malloc(width * sizeof(float), 256));
    real2 = (float *)(_mm_malloc(width * sizeof(float), 256));
    imag2 = (float *)(_mm_malloc(width * sizeof(float), 256));
}

BatchMandelCalculator::~BatchMandelCalculator() {
    free(data);
    data = NULL;
    _mm_free(zReal);
    zReal = NULL;
    _mm_free(zImag);
    zImag = NULL;
    _mm_free(real2);
    zReal = NULL;
    _mm_free(imag2);
    zImag = NULL;
}

int * BatchMandelCalculator::calculateMandelbrot () {
    float *zR = zReal;
    float *zI = zImag;
    float *r2 = real2;
    float *i2 = imag2;
    const int tile = 64;
    int control = 0;

    for (int tj = 0; tj < height / tile; tj++) {
        for (int ti = 0; ti < width / tile; ti++) {
            for (int lj = 0; lj < tile; lj++) {
                int j = tj * tile + lj;
                // aby sme nemuseli počítať index vždy posunieme pdata na začiatok práve prechádzaného úseku
                int *pdata = data;
                pdata = pdata + j * width;
                float y = y_start + j * dy; // current imaginary value
                for (int l = 0; control < tile && l < limit; ++l) {
//                for (int l = 0; l < limit; ++l) {
#pragma omp simd reduction(+:control) aligned(zI, zR, i2, r2, pdata:64) simdlen(64)
//#pragma omp simd reduction(+:l) aligned(pdata:64) simdlen(64)
                    for (int li = 0; li < tile; li++) {
                        int i = ti * tile + li;
                        if (l == 0) {
                            pdata[i] = 200;
//                            pdata[j * width + i] = 0;
                            zI[i] = y;
                            zR[i] = x_start + i * dx;
                        }
                        r2[i] = zR[i] * zR[i];
                        i2[i] = zI[i] * zI[i];
                        float x = x_start + i * dx; // current real value

//                        if (r2[i] + i2[i] <= 4.0f) {
                        if (r2[i] + i2[i] > 4.0f) {
                            if(pdata[i] > limit) {
//                          pdata[j * width + i] += 1;
                                pdata[i] = l;
                                control++;
                            }
                        }
                        zI[i] = 2.0f * zR[i] * zI[i] + y;
                        zR[i] = r2[i] - i2[i] + x;
                    }
                }
                control = 0;
            }
        }
    }
    return data;
}
