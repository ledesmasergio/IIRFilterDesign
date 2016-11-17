#pragma once
//____________________________________________________________________
// IIR Library Version 1.0 (Files: IIR.h and IIR.cpp)
// The IIR library has a set of classes to design, analyze and test IIR digital filters
// © Copyright (C) 2016.  Sergio Ledesma (selo@ugto.mx)
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//____________________________________________________________________ IIR.h
// The IIR namespace includes a set of classes to design and test Infinite Impulse Response Digital Filters
// Date: November 16, 2016

#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define FILTER_TYPE_LOWPASS 0
#define FILTER_TYPE_HIGHPASS 1
namespace IIR //________________________________________ namespace IIR::Ini
{

class PolarBiquadSection;
//_____________________________________________________________________ IIR::BiquadSection
// H(z) =  (b0 + b1 * z^-1 + b2 * z ^-2) / (1 + a1 * z^-1 + a2 * z ^-2)
class BiquadSection
{
public:
	BiquadSection();
	BiquadSection(const IIR::BiquadSection& init);
	BiquadSection(const IIR::PolarBiquadSection& init);
	~BiquadSection();
	void Clear();
	void ResetTaps();
	double b0;
	double b1;
	double b2;
	double a1;
	double a2;
	IIR::BiquadSection& operator=(const IIR::BiquadSection& init);
	IIR::BiquadSection& operator=(const IIR::PolarBiquadSection& init);
	//
	// It returns y[n] given a value of x[n]
	// w[n] = x[n] - a1*w[n-1] - a2*w[n-2]
	// y[n] = b0*w[n] + b1*w[n-1] + b2*w[n-2]
	inline double ComputeOutput(double x);
	//
	// It modifies the samples in buffer, the operation is in place (that is the output is written in the same buffer)
	void ComputeOutput(double* buffer, const int numSamples);
	//
	// It modifies the samples in buffer, the operation is in place (that is the output is written in the same buffer)
	// buffer has two channels, i.e., ch1, ch2, ch1, ch2, ch1, ch2 (numSamples is 3)
	void ComputeOutput_2Channels(double* buffer, const int numSamples);
private:
	void Copy(const IIR::BiquadSection& init);
	void Copy(const IIR::PolarBiquadSection& init);
	//
	double w0;
	double w11;
	double w12;
	double w21;
	double w22;
};

class PolarBiquadsCascade;
//_____________________________________________________________________ IIR::BiquadsCascade
class BiquadsCascade
{
public:
	BiquadsCascade();
	BiquadsCascade(const IIR::BiquadsCascade& init);
	BiquadsCascade(const IIR::PolarBiquadsCascade& init);
	~BiquadsCascade();
	bool SetSize(int numberBiquadSections);
	int GetSize();
	void Delete();
	//
	// It returns y[n] given a value of x[n]
	double ComputeOutput(double x);
	//
	// It modifies the samples in buffer, the operation is in place (that is the output is written in the same buffer)
	void ComputeOutput(double* buffer, const int numSamples);
	//
	// It modifies the samples in buffer, the operation is in place (that is the output is written in the same buffer)
	// buffer has two channels, i.e., ch1, ch2, ch1, ch2, ch1, ch2 (numSamples is 3)
	void ComputeOutput_2Channels(double* buffer, const int numSamples);
	//
	IIR::BiquadsCascade& operator=(const IIR::BiquadsCascade& init);
	IIR::BiquadsCascade& operator=(const IIR::PolarBiquadsCascade& init);
	IIR::BiquadSection& operator[](long index);
	const IIR::BiquadSection& operator[](long index) const;
private:
	IIR::BiquadSection* biquadSection;
	int size;
	bool Copy(const IIR::BiquadsCascade& init);
	bool Copy(const IIR::PolarBiquadsCascade& init);
	friend class PolarBiquadsCascade;
};

//_____________________________________________________________________ IIR::PolarBiquadSection
// Andrew G. Deczky 1972 Eq. 10
//k (z - zeroR * exp(j * zeroAngle)) (z - zeroR * exp(-j * zeroAngle) ) / (z - poleR *exp (j * poleAngle) ) (z - poleR *exp (-j * poleAngle) )
class PolarBiquadSection
{
public:
	PolarBiquadSection();
	PolarBiquadSection(const IIR::PolarBiquadSection& init);
	PolarBiquadSection(const IIR::BiquadSection& init);
	~PolarBiquadSection();
	void Clear();
	double zeroR; // alpha
	double zeroAngle; // alpha
	double poleR; // beta
	double poleAngle; // beta
	double k;
	IIR::PolarBiquadSection& operator=(const IIR::PolarBiquadSection& init);
	IIR::PolarBiquadSection& operator=(const IIR::BiquadSection& init);
	//
	void GetPoles(double& out_reP1, double& out_imP1, double& out_reP2, double& out_imP2);
	void GetZeros(double& out_reZ1, double& out_imZ1, double& out_reZ2, double& out_imZ2);
private:
	void Copy(const IIR::PolarBiquadSection& init);
	void Copy(const IIR::BiquadSection& init);
};

//_____________________________________________________________________ IIR::PolarBiquadsCascade
class PolarBiquadsCascade
{
public:
	PolarBiquadsCascade();
	PolarBiquadsCascade(const IIR::PolarBiquadsCascade& init);
	PolarBiquadsCascade(const IIR::BiquadsCascade& init);
	~PolarBiquadsCascade();
	bool SetSize(int numberBiquadSections);
	int GetSize() const;
	void Delete();
	//
	IIR::PolarBiquadsCascade& operator=(const IIR::PolarBiquadsCascade& init);
	IIR::PolarBiquadsCascade& operator=(const IIR::BiquadsCascade& init);
	IIR::PolarBiquadSection& operator[](long index);
	const IIR::PolarBiquadSection& operator[](long index) const;
	//
	double GetMagnitude(double freq_rads);
	double GetGroupDelay(double freq_rads);
private:
	IIR::PolarBiquadSection* biquadSection;
	int size;
	bool Copy(const IIR::PolarBiquadsCascade& init);
	bool Copy(const IIR::BiquadsCascade& init);
	friend class BiquadsCascade;
};

//_____________________________________________________________________ BiquadFilter
// It provides support to implement an infinite impulse response filter
class BiquadFilter
{
public:
	BiquadFilter(void);
	~BiquadFilter(void);
	bool CreateLowPass(double cutFreq_rads, double passBand_dB, double stopFreq_rads, double stopBand_dB);
	bool CreateHighPass(double cutFreq_rads, double passBand_dB, double stopFreq_rads, double stopBand_dB);
	IIR::PolarBiquadsCascade biquadsCascade;
protected:
	//
	// input:  A*s*s + B*s + C
	// output:  a  + b * z^-1 + c * z^-2
	// filter_type: FILTER_TYPE_LOWPASS s = 2*(1 - z^-1)/(1 + z^-1), FILTER_TYPE_HIGHPASS s = 2*(1 + z^-1)/(1 - z^-1)
	static void BilinearTransformation(int filter_type, double A, double B, double C, double& out_a, double& out_b, double& out_c);
	//
	virtual bool Create(int filter_type, double cutFreq_rads, double passBand_dB, double stopFreq_rads, double stopBand_dB);
};

//_____________________________________________________________________ ButterworthFilter
class ButterworthFilter: public IIR::BiquadFilter
{
public:
	ButterworthFilter();
	~ButterworthFilter();
private:
	bool Create(int filter_type, double cutFreq_rads, double passBand_dB, double stopFreq_rads, double stopBand_dB);
};

//_____________________________________________________________________ ChebyshevFilter
class ChebyshevFilter: public IIR::BiquadFilter
{
public:
	ChebyshevFilter();
	~ChebyshevFilter();
private:
	bool Create(int filter_type, double cutFreq_rads, double passBand_dB, double stopFreq_rads, double stopBand_dB);
};

//_____________________________________________________________________ EllipticFilter
class EllipticFilter: public IIR::BiquadFilter
{
public:
	EllipticFilter();
	~EllipticFilter();
private:
	bool Create(int filter_type, double cutFreq_rads, double passBand_dB, double stopFreq_rads, double stopBand_dB);
};

}
