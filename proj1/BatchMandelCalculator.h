/**
 * @file BatchMandelCalculator.h
 * @author René Rešetár <xreset00@stud.fit.vutbr.cz>
 * @brief Implementation of Mandelbrot calculator that uses SIMD paralelization over small batches
 * @date 14.11.2021
 */
#ifndef BATCHMANDELCALCULATOR_H
#define BATCHMANDELCALCULATOR_H

#include <BaseMandelCalculator.h>

class BatchMandelCalculator : public BaseMandelCalculator
{
public:
    BatchMandelCalculator(unsigned matrixBaseSize, unsigned limit);
    ~BatchMandelCalculator();
    int * calculateMandelbrot();

private:
    // @TODO add all internal parameters
    int *data;
    float *zReal;
    float *zImag;
    float *real2;
    float *imag2;
};

#endif