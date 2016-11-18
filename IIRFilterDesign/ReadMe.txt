========================================================================
    The IIR Library:  IIR.h and IIR.cpp
    Version 1.0
========================================================================

The IIR namespace includes a set of classes to design and test Infinite 
Impulse Response Digital Filters.

________________________________________________________________________
LICENSE

The IIR library (files IIR.h and IIR.cpp) is licensed under the GNU 
General Pulbic License (GPL), version 3

________________________________________________________________________
IIR::BiquadSection 

This class represents a biquadratic digital section
H(z) =  (b0 + b1 * z^-1 + b2 * z ^-2) / (1 + a1 * z^-1 + a2 * z ^-2)

________________________________________________________________________
IIR::PolarBiquadSection 

This class represents a biquadratic digital section with the poles 
and zeros in polar form

________________________________________________________________________
IIR::BiquadsCascade

This class represents an array of IIR::BiquadSections implemented in 
Direct Form 2

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
 EXAMPLE 1:
	Elliptic lowpass filter
	Cut frequency:  2 rad
	Stop frequency:  2.5 rad
	Passband ripples:  1.0 dB
	Stopband gain:  -50 dB

// IIRFilterDesign.cpp : Defines the entry point for the console application.
//
#include "stdafx.h"
#include <iostream>
#include "IIR.h"

int main()
{
	IIR::EllipticFilter filter;
	if (filter.CreateLowPass(2.0, 1.0, 2.5, -50.0) == false)
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

________________________________________________________________________
 EXAMPLE 2:
	Butterworth lowpass filter
	Cut frequency:  1 rad
	Stop frequency:  2.0 rad
	Passband ripples:  1.0 dB
	Stopband gain:  -30 dB

// IIRFilterDesign.cpp : Defines the entry point for the console application.
//
#include "stdafx.h"
#include <iostream>
#include "IIR.h"
#define FREQ_COUNT 16

int main()
{
	IIR::ButterworthFilter filter;
	double frequency;
	if (filter.CreateLowPass(1.0, 1.0, 2.0, -30.0) == false)
	{
		std::cout << "Unable to create low pass filter" << std::endl;
		return 0;
	}
	for (int i = 0; i < FREQ_COUNT; i++)
	{
		frequency = i*M_PI/(FREQ_COUNT-1.0);
		std::cout << frequency << " rads -> ";
		std::cout << filter.biquadsCascade.GetMagnitude(frequency);
		std::cout << std::endl;
	}
	getchar();
	return 0;
}

________________________________________________________________________
 EXAMPLE 3: Real time Digital Signal Processing
	Elliptic highpass filter
	Cut frequency:  2.5 rad
	Stop frequency:  2.0 rad
	Passband ripples:  1.0 dB
	Stopband gain:  -50 dB

// IIRFilterDesign.cpp : Defines the entry point for the console application.
//
#include "stdafx.h"
#include <iostream>
#include "IIR.h"
#define BUFFER_SIZE 1024

int main()
{
	IIR::EllipticFilter filter;
	if (filter.CreateHighPass(2.5, 1.0, 2.0, -50.0) == false)
	{
		std::cout << "Unable to create high pass filter" << std::endl;
		return 0;
	}
	IIR::BiquadsCascade biquadsCascade = filter.biquadsCascade;
	double buffer[BUFFER_SIZE];
	while(...)
	{
		//______________________________________________ 1. Fill the buffer
		for (int i = 0; i < BUFFER_SIZE; i++) buffer[i] = rand()/(RAND_MAX-1.0);
		//______________________________________________ 2. DSP
		biquadsCascade.ComputeOutput(buffer, BUFFER_SIZE);
		//______________________________________________ 3. Do something with buffer
		// ...
	}
	return 0;
}