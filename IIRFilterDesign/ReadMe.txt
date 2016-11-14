========================================================================
    The IIR namespace files:  IIR.h and IIR.cpp
========================================================================

The IIR namespace includes a set of classes to design and test Infinite Impulse Response
Digital Filters.

________________________________________________________________________
IIR::BiquadSection 

This class represents a biquadratic digital section
H(z) =  (b0 + b1 * z^-1 + b2 * z ^-2) / (1 + a1 * z^-1 + a2 * z ^-2)

________________________________________________________________________
IIR::PolarBiquadSection 

This class represents a biquadratic digital section with the poles and zeros
in polar form

________________________________________________________________________
IIR::BiquadsCascade

This class represents an array of IIR::BiquadSections implemented in Direct Form 2
- to compute the output of the filter: use the ComputeOutput function 

________________________________________________________________________
IIR::PolarBiquadsCascade

This class represents an array of IIR::PolarBiquadSections
- to compute the Magnitude of the Frequency Response: use the GetMagnitude function
- to compute the Group Delay: use the GetGroupDelay function 

________________________________________________________________________
IIR::BiquadFilter

This class represents an IIR filter with a set of biquad sections in cascade

________________________________________________________________________
IIR::ButterworthFilter

This class represents a Butterworth IIR filter

________________________________________________________________________
IIR::ChebyshevFilter

This class represents a Chebyshev IIR filter

________________________________________________________________________
IIR::EllipticFilter

This class represents a Elliptic IIR filter
________________________________________________________________________
 EXAMPLE of a Chebyshev lowpass filter

          Cut frequency:  1 rad
          Stop frequency:  2.0 rad
          Passband ripples:  1.0 dB
          Stopband gain:  -30 dB

// IIRFilterDesign.cpp : Defines the entry point for the console application.
//
#include "stdafx.h"
#include <iostream>
#include "IIR.h"

int main()
{
	IIR::ChebyshevFilter filter;
	if (filter.CreateLowPass(1.0, 1.0, 2.0, -30.0) == false)
	{
		std::cout << "Unable to create low pass filter" << std::endl;
		return 0;
	}
	IIR::BiquadsCascade biquadsCascade = filter.biquadsCascade;
	const int size = biquadsCascade.GetSize();
	for(int i = 0; i < size; i++)
	{
		std::cout << "b0 = " << biquadsCascade[i].b0 << std::endl;
		std::cout << "b1 = " << biquadsCascade[i].b1 << std::endl;
		std::cout << "b2 = " << biquadsCascade[i].b2 << std::endl;
		std::cout << "a1 = " << biquadsCascade[i].a1 << std::endl;
		std::cout << "a2 = " << biquadsCascade[i].a2 << std::endl;
		std::cout << "____________________"  << std::endl;
	}
	getchar();
    return 0;
}