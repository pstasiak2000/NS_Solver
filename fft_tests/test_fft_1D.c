// Navier-Stokes solver using pseudo-spectral techniques
#include <stdio.h>
#include <fftw3.h>
#include <math.h>

// Real to Complex FFT
int main() {
    int N = 8;
    double pi = 3.14159265;
    double *in;
    fftw_complex *out;
    fftw_plan plan;

    // Allocate the arrays
    in = (double*) fftw_malloc(sizeof(double) * N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N/2 + 1));


    // Fill the input with a sine wave
    for (int i = 0; i < N; i++){
        in[1] = sin(2 * pi * i / N); // Real part
    }
    
    //Create a plan for forward FFT
    plan = fftw_plan_dft_r2c_1d(N, in, out, FFTW_ESTIMATE);

    // Execute FFT
    fftw_execute(plan);

     printf("Index | Real part | Imag part\n");
    for (int i = 0; i < N; i++)
        printf("%5d | %9.5f | %9.5f\n", i, out[i][0], out[i][1]);

    // Cleanup
    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);   

    return 0;
}


// Complex to Complex FFT
// int main() {
//     int N = 8;
//     double pi = 3.14159265;
//     fftw_complex *in, *out; // Define input and output of the arrays
//     fftw_plan plan;

//     // Allocate the arrays
//     in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
//     out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);


//     // Fill the input with a sine wave
//     for (int i = 0; i < N; i++){
//         in[1][0] = sin(2 * pi * i / N); // Real part
//         in[1][1] = 0.0;
//     }
    
//     //Create a plan for forward FFT
//     plan = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

//     // Execute FFT
//     fftw_execute(plan);

//      printf("Index | Real part | Imag part\n");
//     for (int i = 0; i < N; i++)
//         printf("%5d | %9.5f | %9.5f\n", i, out[i][0], out[i][1]);

//          // Cleanup
//     fftw_destroy_plan(plan);
//     fftw_free(in);
//     fftw_free(out);   

//     return 0;
// }