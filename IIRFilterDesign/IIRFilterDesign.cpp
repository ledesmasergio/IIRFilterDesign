// IIRFilterDesign.cpp : Defines the entry point for the console application.
//
#include "stdafx.h"
#include <iostream>
#include "IIR.h"

//_________________________________________ Design an Elliptic Lowpass filter wc=2.0, ws=2.5
//int main()
//{
//	IIR::EllipticFilter filter;
//	if (filter.CreateLowPass(2.0, 1.0, 2.5, -50.0) == false)
//	{
//		std::cout << "Unable to create low pass filter" << std::endl;
//		return 0;
//	}
//	IIR::BiquadsCascade biquadsCascade = filter.biquadsCascade;
//	const int size = biquadsCascade.GetSize();
//	for(int i = 0; i < size; i++)
//	{
//		std::cout << "b0 = " << biquadsCascade[i].b0 << std::endl;
//		std::cout << "b1 = " << biquadsCascade[i].b1 << std::endl;
//		std::cout << "b2 = " << biquadsCascade[i].b2 << std::endl;
//		std::cout << "a1 = " << biquadsCascade[i].a1 << std::endl;
//		std::cout << "a2 = " << biquadsCascade[i].a2 << std::endl;
//		std::cout << "____________________"  << std::endl;
//	}
//	getchar();
//    return 0;
//}

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
	int i;
	//_____________________ Fill the buffer with noise
	for (i = 0; i < BUFFER_SIZE; i++) buffer[i] = rand()/(RAND_MAX-1.0);
	biquadsCascade.ComputeOutput(buffer, BUFFER_SIZE);
	//_____________________ Do something with buffer
	// ...
	return 0;
}





