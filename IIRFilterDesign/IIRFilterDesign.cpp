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
		cout << "Unable to create low pass filter" << endl;
		return 0;
	}
	IIR::BiquadsCascade biquadsCascade = filter.biquadsCascade;
	const int size = biquadsCascade.GetSize();
	for(int i = 0; i < size; i++)
	{
		cout << "b0 = " << biquadsCascade[i].b0 << endl;
		cout << "b1 = " << biquadsCascade[i].b1 << endl;
		cout << "b2 = " << biquadsCascade[i].b2 << endl;
		cout << "a1 = " << biquadsCascade[i].a1 << endl;
		cout << "a2 = " << biquadsCascade[i].a2 << endl;
		cout << "____________________"  << endl;
	}
	getchar();
    return 0;
}

