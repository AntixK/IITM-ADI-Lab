/**************************************************************************

               G O E R T Z E L 'S   A L G O R I T H M

                                 F O R

        D I S C R E T E   F O U R I E R   T R A N S F O R M

        By S Anand Krishnamoorthy, Project Associate, ADI DSP Lab
                            IIT MADRAS, India

----------------------------------------------------------------------------

Reference(s):
[1]  Proakis, John G., and Dimitris G. Manolakis.
     "Digital Signal Processing: Principles, Algorithms, and Applications".
     3rd Edition. Upper Saddle River, NJ: Prentice Hall, 1996, pp. 480–481.

External Links:
[1]  "Goertzel's Algorithm". Cnx.org. 2006-09-12. Retrieved 2014-02-03.

[2]   C Code implementation of iterative Goertzel algorithm
      http://netwerkt.wordpress.com/2011/08/25/goertzel-filter/

Interface and Usage:

=> Set the N value for the N-point DFT

=> Run the phase_factors() function to pre-compute the phase factors

=> To compute the DFT of a sequence, call the Goertzel_FT();
   Parameters: Input and Output arrays of Complex data type

=> To compute the DFT at a particular index, call the function Goertzel();
   Parameters: Input array of Complex data type and
   index (k) of unsigned 16 bit integer field

 => To compute the magnitude of DFT at a test frequency, call the function
    Goertzel_get_mag();
    Parameter: Input array of complex data type,
    sampling frequency and test frequency of unsigned 16 bit integer field.

Start Date: 17 June 2015
End Date:   18 June 2015
Last Modified: 4 July 2015

Copyright (c) 2015 IIT-MADRAS
All rights reserved.
****************************************************************************/

#ifndef Goertzel_H_
#define Goertzel_H_

#include <math.h>

#define PI 3.1415926535
#define N 205            //Length of the array

using namespace std;

typedef struct
{
    double Re, Im;

} Complex;

//Declation of Global Variables
static Complex W[N];

/** Declaration of functions**/
static void phase_factors();
static Complex Goertzel(double input[], uint16_t k);
static void Goertzel_FT(Complex output[],double input[]);
static double Goertzel_get_mag(double input[], uint16_t smpl_freq, uint16_t test_freq);

/*-------------------------------------------------------------------------*/

/** Function Definitions **/
static void phase_factors()
{
    for(uint16_t k =0; k < N; ++k)
    {
        W[k].Re = cos(2*PI*k/N);
        W[k].Im = sin(2*PI*k/N);
    }

}

static Complex Goertzel(double input[], uint16_t k)
{
    double sequence[N+1];
    Complex output;
    uint16_t n;
    /* Generating the Sequence */
    sequence[0] = input[0];
    sequence[1] = 2*cos(2*PI*k/N)*sequence[0] + input[1];

    for(n =2; n < N; ++n)
    {
        sequence[n] = 2*cos(2*PI*k/N)*sequence[n-1] - sequence[n-2] + input[n];

    }
     sequence[n] = 2*cos(2*PI*k/N)*sequence[n-1] - sequence[n-2];

    output.Re = sequence[N] - W[k].Re*sequence[N-1];
    output.Im = W[k].Im*sequence[N-1];

    return output;

}

static void Goertzel_FT(Complex output[], double input[])
{
    for (uint16_t k = 0; k < N; ++k)
    {
        output[k] = Goertzel(input, k);
    }
}

static double Goertzel_get_mag(double input[], uint16_t smpl_freq, uint16_t test_freq )
{
    uint16_t k = (int) (0.5 + ((double)(N*test_freq)/(double)smpl_freq));
    Complex out;
    out = Goertzel(input, k);
    double mag = sqrt((out.Re*out.Re)+(out.Im*out.Im));
    return mag;
}

#endif

/** End of Header File -----------------------------------------------**/
