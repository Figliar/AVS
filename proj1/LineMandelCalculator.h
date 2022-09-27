/**
 * @file LineMandelCalculator.h
 * @author René Rešetár <xreset00@stud.fit.vutbr.cz>
 * @brief Implementation of Mandelbrot calculator that uses SIMD paralelization over lines
 * @date 14.11.2021
 */

#include <BaseMandelCalculator.h>

class LineMandelCalculator : public BaseMandelCalculator
{
public:
    LineMandelCalculator(unsigned matrixBaseSize, unsigned limit);
    ~LineMandelCalculator();
    int *calculateMandelbrot();

private:
    // @TODO add all internal parameters
    int *data;
    float *zReal;
    float *zImag;
    float *real2;
    float *imag2;
};