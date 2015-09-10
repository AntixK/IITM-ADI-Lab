/**************************************************************************

               C O O L E Y - T U K E Y   A L G O R I T H M

                                 F O R

                F A S T   F O U R I E R   T R A N S F O R M

        By S Anand Krishnamoorthy, Project Associate, ADI DSP Lab
                            IIT MADRAS, India

----------------------------------------------------------------------------

Reference(s):
[1]   Cooley, James W.; Tukey, John W. (1965).
      "An algorithm for the machine calculation of complex Fourier series".
	  Math. Comput. 19: 297–301. doi:10.2307/2003354.

External Links:
[1]   Fast Fourier transform http://www.librow.com/articles/article-10

[2]   Understanding the FFT Algorithm
      https://jakevdp.github.io/blog/2013/08/28/understanding-the-fft/

Interface and Usage:

=> Set the N value for the N-point FFT (must be a power of 2)

=> Run the calc_twiddle() function to pre-compute the phase factors

=> To compute the FFT of a sequence through Decimation in Time, call the DIT_FFT(); function
   Parameters: Input array of Complex data type (In place computation)

=> To compute the FFT of a sequence through Decimation in Frequency,
   call the DIF_FFT(); function
   Parameters: Input array of Complex data type (In place computation)

=> To compute the IFFT of a sequence through Decimation in Time, call the DIT_IFFT(); function
   Parameters: Input array of Complex data type (In place computation)

=> To compute the IFFT of a sequence through Decimation in Frequency,
   call the DIF_IFFT(); function
   Parameters: Input array of Complex data type (In place computation)


Start Date: 27 June 2015
End Date:   27 June 2015
Last Modified: 28 June 2015

Copyright (c) 2015 IIT-MADRAS
All rights reserved.
****************************************************************************/

#ifndef FFT_H_
#define FFT_H_

#include <math.h>

#define N 16         // Must be a power of 2
#define PI 3.1415925535897

const int Stages = log2(N);
using namespace std;

typedef struct
{
    double Re;
    double Im;
}Complex;

static Complex twiddle_factors[N/2];

/** Declaration of functions**/

static void calc_twiddle(bool flag);
static void DIT_FFT(Complex input[]);
static void DIF_FFT(Complex input[]);
static void DIT_IFFT(Complex input[]);
static void DIF_IFFT(Complex input[]);
static void bit_reverse(Complex input[]);
static void fft_main(Complex input[],bool flag);


/*-------------------------------------------------------------------------*/

/** Function Definitions **/

static void calc_twiddle(bool flag)
{

    for(uint8_t n=0; n<(N/2); ++n)
    {
        twiddle_factors[n].Re = cos(2*PI*n/N);
        twiddle_factors[n].Im = sin(2*PI*n/N);

         if(flag ==0)
            twiddle_factors[n].Im = (-1)*twiddle_factors[n].Im;
    }
}

static void bit_reverse(Complex input[])
{
    uint8_t rev_indx, temp, bit_mask = 0x00000001;

    for(uint8_t i=0; i<N; ++i)
    {
        temp = i;rev_indx = 0;
        for(uint8_t j=0; j<Stages; ++j)
        {
            rev_indx <<= 1;
            rev_indx += temp & bit_mask;
            temp >>= 1;
        }

        if( rev_indx > i)
        {
            input[rev_indx].Re += input[i].Re;
            input[rev_indx].Im += input[i].Im;

            input[i].Re = input[rev_indx].Re - input[i].Re;
            input[i].Im = input[rev_indx].Im - input[i].Im;

            input[rev_indx].Re -= input[i].Re;
            input[rev_indx].Im -= input[i].Im;

        }
    }

}

static void fft_main(Complex input[], bool flag)
{
    uint8_t Blocks, Butterfly;

    if(flag ==0)
    {
        Blocks = (N/2);
        Butterfly = 1;
    }
    else
    {
        Blocks = 1;
        Butterfly = (N/2);
    }

    uint8_t indx1 = 0, indx2 =0;
    double temp_real = 0.0, temp_comp = 0.0;

    for(uint8_t i=0; i<Stages; ++i)
    {
        for(uint8_t j=0; j<Blocks; ++j)
        {
            for(uint8_t k=0; k<Butterfly; ++k)
            {
                indx2 = indx1 + Butterfly;

                if(flag ==0)
                {

                    /** (a + ib) * (c + id) = (ac - bd)+i(ad + bc) **/
                    temp_real = (twiddle_factors[k*Blocks].Re * input[indx2].Re)-(twiddle_factors[k*Blocks].Im*input[indx2].Im);
                    temp_comp = (twiddle_factors[k*Blocks].Re * input[indx2].Im)+(twiddle_factors[k*Blocks].Im*input[indx2].Re);

                    input[indx2].Re = input[indx1].Re - temp_real;
                    input[indx2].Im = input[indx1].Im - temp_comp;

                    input[indx1].Re += temp_real;
                    input[indx1].Im += temp_comp;
                }
                else
                {
                    double t_real = input[indx2].Re;
                    double t_comp = input[indx2].Im;

                    input[indx2].Re = input[indx1].Re - input[indx2].Re;
                    input[indx2].Im = input[indx1].Im - input[indx2].Im;

                    input[indx1].Re += t_real;
                    input[indx1].Im += t_comp;

                    t_real = input[indx2].Re;

                    input[indx2].Re = (twiddle_factors[k*Blocks].Re * input[indx2].Re)-(twiddle_factors[k*Blocks].Im*input[indx2].Im);
                    input[indx2].Im = (twiddle_factors[k*Blocks].Re * input[indx2].Im)+(twiddle_factors[k*Blocks].Im*t_real);
                }


                ++indx1;
            }
            indx1 = indx2+1;
        }

         if(flag ==0)
         {
            Blocks/=2;
            Butterfly*=2;
         }
         else
         {
            Blocks*=2;
            Butterfly/=2;
         }

        indx1 =indx2=0;
    }

}

static void DIT_FFT(Complex input[])
{
    calc_twiddle(0);
    bit_reverse(input);
    fft_main(input,0);
}

static void DIF_FFT(Complex input[])
{
   calc_twiddle(0);
   fft_main(input,1);
   bit_reverse(input);
}

static void DIT_IFFT(Complex input[])
{
    calc_twiddle(1);
    bit_reverse(input);
    fft_main(input,0);

    for(uint8_t i=0; i<N; ++i)
    {
        input[i].Re/=(double)N;
        input[i].Im/=(double)N;

    }
}

static void DIF_IFFT(Complex input[])
{
   calc_twiddle(1);
   fft_main(input,1);
   bit_reverse(input);

     for(uint8_t i=0; i<N; ++i)
    {
        input[i].Re/=(double)N;
        input[i].Im/=(double)N;

    }
}

#endif

/** End of Header File -----------------------------------------------**/
