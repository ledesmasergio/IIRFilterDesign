#include "stdafx.h"  //_____________________________________________ IIR.cpp
#include "IIR.h" 

namespace IIR //________________________________________ namespace IIR::Ini
{

//_____________________________________________________________________ IIR::BiquadSection
BiquadSection::BiquadSection()
{
	Clear();
}

BiquadSection::~BiquadSection()
{
}

void BiquadSection::Clear()
{
	b0 = 1.0;
	b1 = 0.0;
	b2 = 0.0;
	a1 = 0.0;
	a2 = 0.0;
	//
	ResetTaps();
}

void BiquadSection::ResetTaps()
{
	w11 = 0.0;
	w12 = 0.0;
	w21 = 0.0;
	w22 = 0.0;
}

BiquadSection::BiquadSection(const IIR::BiquadSection& init)
{
	Clear();
	Copy(init);
}

BiquadSection::BiquadSection(const IIR::PolarBiquadSection& init)
{
	Clear();
	Copy(init);
}

IIR::BiquadSection& BiquadSection::operator=(const IIR::BiquadSection& init)
{
	Copy(init);
	return *this;
}

IIR::BiquadSection& BiquadSection::operator=(const IIR::PolarBiquadSection& init)
{
	Copy(init);
	return *this;
}

void BiquadSection::Copy(const IIR::BiquadSection& init)
{
	b0 = init.b0;
	b1 = init.b1;
	b2 = init.b2;
	a1 = init.a1;
	a2 = init.a2;
	//
	w11 = init.w11;
	w12 = init.w12;
	w21 = init.w21;
	w22 = init.w22;
}

void BiquadSection::Copy(const IIR::PolarBiquadSection& init)
{
	ResetTaps();
	b0 = init.k;
	//___________________________________________________________ zeros
	double re = init.zeroR * cos(init.zeroAngle);
	b1 = -2.0*re*b0; // b1/b0 = - zero - (zero)* = -2 Re (zero)
	b2 = init.zeroR*init.zeroR*b0; // b2/b0 = (zero)(zero)* =  |z|^2
	//___________________________________________________________ poles
	re = init.poleR * cos(init.poleAngle);
	a1 = -2.0*re; // a1 = - pole - (pole)* = -2 Re (pole)
	a2 = init.poleR*init.poleR; // a2 = (pole)(pole)* =  |p|^2
}

// w[n] = x[n] - a1*w[n-1] - a2*w[n-2]
// y[n] = b0*w[n] + b1*w[n-1] + b2*w[n-2]
double BiquadSection::ComputeOutput(double x)
{
	w0 = x - a1*w11 - a2*w12;
	x = b0*w0 + b1*w11 + b2*w12;
	w12 = w11;
	w11 = w0;
	return x;
}

void BiquadSection::ComputeOutput(double* buffer, const int numSamples)
{
	for (int i = 0; i < numSamples; i++, buffer++)
	{
		w0 = *buffer - a1*w11 - a2*w12;
		*buffer = b0*w0 + b1*w11 + b2*w12;
		w12 = w11;
		w11 = w0;
	}
}

void BiquadSection::ComputeOutput_2Channels(double* buffer, const int numSamples)
{
	for (int i = 0; i < numSamples; i++)
	{
		w0 = *buffer - a1*w11 - a2*w12;
		*buffer = b0*w0 + b1*w11 + b2*w12;
		w12 = w11;
		w11 = w0;
		buffer++;
		//
		w0 = *buffer - a1*w21 - a2*w22;
		*buffer = b0*w0 + b1*w21 + b2*w22;
		w22 = w21;
		w21 = w0;
		buffer++;
	}
}

//_____________________________________________________________________ IIR::BiquadsCascade
BiquadsCascade::BiquadsCascade()
{
	biquadSection = NULL;
	size = 0;
}

BiquadsCascade::BiquadsCascade(const IIR::BiquadsCascade& init)
{
	biquadSection = NULL;
	size = 0;
	Copy(init);
}

BiquadsCascade::BiquadsCascade(const IIR::PolarBiquadsCascade& init)
{
	biquadSection = NULL;
	size = 0;
	Copy(init);
}

BiquadsCascade::~BiquadsCascade()
{
	if (biquadSection != NULL) delete [] biquadSection;
}

bool BiquadsCascade::SetSize(int numberBiquadSections)
{
	if (numberBiquadSections <= 0) return false;
	if (size == numberBiquadSections) return true;
	Delete();
	biquadSection = new IIR::BiquadSection[numberBiquadSections];
	if (biquadSection == NULL) return false;
	size = numberBiquadSections;
	return true;
}

int BiquadsCascade::GetSize()
{
	return size;
}

void BiquadsCascade::Delete()
{
	if (biquadSection != NULL) delete[] biquadSection;
	biquadSection = NULL;
	size = 0;
}

double BiquadsCascade::ComputeOutput(double x)
{
	int count = 0;
	for (IIR::BiquadSection *p = biquadSection; count < size; count++, p++)
	{
		x = p->ComputeOutput(x);
	}
	return x;
}

void BiquadsCascade::ComputeOutput(double* buffer, const int numSamples)
{
	int count = 0;
	for (IIR::BiquadSection *p = biquadSection; count < size; count++, p++)
	{
		p->ComputeOutput(buffer, numSamples);
	}
}

void BiquadsCascade::ComputeOutput_2Channels(double* buffer, const int numSamples)
{
	int count = 0;
	for (IIR::BiquadSection *p = biquadSection; count < size; count++, p++)
	{
		p->ComputeOutput_2Channels(buffer, numSamples);
	}
}

IIR::BiquadsCascade& BiquadsCascade::operator=(const IIR::BiquadsCascade& init)
{
	Delete();
	Copy(init);
	return *this;
}

IIR::BiquadsCascade& BiquadsCascade::operator=(const IIR::PolarBiquadsCascade& init)
{
	Delete();
	Copy(init);
	return *this;
}

IIR::BiquadSection& BiquadsCascade::operator[](long index)
{
#ifdef _DEBUG
	if (index < 0) throw L"BiquadsCascade: The index must be positive";
	if (index >= size) throw L"BiquadsCascade: The index was too large";
#endif
	return biquadSection[index];
}

const IIR::BiquadSection& BiquadsCascade::operator[](long index) const
{
	#ifdef _DEBUG
	if (index < 0) throw L"BiquadsCascade: The index must be positive";
	if (index >= size) throw L"BiquadsCascade: The index was too large";
	#endif
	return biquadSection[index];
}

bool BiquadsCascade::Copy(const IIR::BiquadsCascade& init)
{
	//________________________________________________ Empty Cascade
	if (init.size == 0)
	{
		Delete();
		return true;
	}
	if (init.biquadSection == NULL)
	{
		Delete();
		return true;
	}
	//________________________________________________ Cascade with sections
	if (SetSize(init.size) == false) return false;
	for (int i = 0; i < size; i++)
	{
		biquadSection[i] = init.biquadSection[i];
	}
	return true;
}

bool BiquadsCascade::Copy(const IIR::PolarBiquadsCascade& init)
{
	//________________________________________________ Empty Cascade
	if (init.size == 0)
	{
		Delete();
		return true;
	}
	if (init.biquadSection == NULL)
	{
		Delete();
		return true;
	}
	//________________________________________________ Cascade with sections
	if (SetSize(init.size) == false) return false;
	for (int i = 0; i < size; i++)
	{
		biquadSection[i] = init.biquadSection[i];
	}
	return true;
}

//_____________________________________________________________________ IIR::PolarBiquadSection
PolarBiquadSection::PolarBiquadSection()
{
	Clear();
}

PolarBiquadSection::PolarBiquadSection(const IIR::PolarBiquadSection& init)
{
	Clear();
	Copy(init);
}

PolarBiquadSection::PolarBiquadSection(const IIR::BiquadSection& init)
{
	Clear();
	Copy(init);
}

PolarBiquadSection::~PolarBiquadSection()
{
}

void PolarBiquadSection::Clear()
{
	k = 1.0;
	zeroR = 0.0; // alpha
	zeroAngle = 0.0; // alpha
	poleR = 0.0; // beta
	poleAngle = 0.0; // beta
}

IIR::PolarBiquadSection& PolarBiquadSection::operator=(const IIR::PolarBiquadSection& init)
{
	Copy(init);
	return *this;
}

IIR::PolarBiquadSection& PolarBiquadSection::operator=(const IIR::BiquadSection& init)
{
	Copy(init);
	return *this;
}

void PolarBiquadSection::GetPoles(double& out_reP1, double& out_imP1, double& out_reP2, double& out_imP2)
{
	out_reP1 = poleR * cos(poleAngle);
	out_imP1 = poleR * sin(poleAngle);
	//_____________________________________________ Complex Conjugate
	out_reP2 = out_reP1;
	out_imP2 = -out_imP1;
}

void PolarBiquadSection::GetZeros(double& out_reZ1, double& out_imZ1, double& out_reZ2, double& out_imZ2)
{
	out_reZ1 = zeroR * cos(zeroAngle);
	out_imZ1 = zeroR * sin(zeroAngle);
	//_____________________________________________ Complex Conjugate
	out_reZ2 = out_reZ1;
	out_imZ2 = -out_imZ1;
}

void PolarBiquadSection::Copy(const IIR::PolarBiquadSection& init)
{
	k = init.k;
	zeroR = init.zeroR;
	zeroAngle = init.zeroAngle;
	poleR = init.poleR;
	poleAngle = init.poleAngle;
}

void PolarBiquadSection::Copy(const IIR::BiquadSection& init)
{
	k = init.b0;
	double re, im;
	//___________________________________________________________ zeros
	re = -init.b1/(2.0*init.b0); // b1/b0 = - zero - (zero)* = -2 Re (zero) ;      Re(zero) = - (b1/b0) / 2
	im = sqrt(init.b2/init.b0 - re*re); // b2/b0 = (zero)(zero)* =  |z|^2;        Im(zero) = sqrt(b2/b0 - re*re);
	//
	zeroR = sqrt(re*re + im*im);
	zeroAngle = atan2(im, re); // -pi <= angle < pi
	//___________________________________________________________ poles
	re = -init.a1/2.0; // c = - pole - (pole)* = -2 Re (pole) ;      Re(pole) = - c / 2
	im = sqrt(init.a2 - re*re); // d = (pole)(pole)* =  |p|^2;        Im(pole) = sqrt(d - re*re);
	//
	poleR = sqrt(re*re + im*im);
	poleAngle = atan2(im, re); //  -pi <= angle < pi
}

//_____________________________________________________________________ IIR::PolarBiquadsCascade
PolarBiquadsCascade::PolarBiquadsCascade()
{
	biquadSection = NULL;
	size = 0;
}

PolarBiquadsCascade::PolarBiquadsCascade(const IIR::PolarBiquadsCascade& init)
{
	biquadSection = NULL;
	size = 0;
	Copy(init);
}

PolarBiquadsCascade::PolarBiquadsCascade(const IIR::BiquadsCascade& init)
{
	biquadSection = NULL;
	size = 0;
	Copy(init);
}

PolarBiquadsCascade::~PolarBiquadsCascade()
{
	if (biquadSection != NULL) delete[] biquadSection;
}

bool PolarBiquadsCascade::SetSize(int numberBiquadSections)
{
	if (numberBiquadSections <= 0) return false;
	if (size == numberBiquadSections) return true;
	Delete();
	biquadSection = new IIR::PolarBiquadSection[numberBiquadSections];
	if (biquadSection == NULL) return false;
	size = numberBiquadSections;
	return true;
}

int PolarBiquadsCascade::GetSize() const
{
	return size;
}

void PolarBiquadsCascade::Delete()
{
	if (biquadSection != NULL) delete[] biquadSection;
	biquadSection = NULL;
	size = 0;
}

IIR::PolarBiquadsCascade& PolarBiquadsCascade::operator=(const IIR::PolarBiquadsCascade& init)
{
	Delete();
	Copy(init);
	return *this;
}

IIR::PolarBiquadsCascade& PolarBiquadsCascade::operator=(const IIR::BiquadsCascade& init)
{
	Delete();
	Copy(init);
	return *this;
}

IIR::PolarBiquadSection& PolarBiquadsCascade::operator[](long index)
{
#ifdef _DEBUG
	if (index < 0) throw L"PolarBiquadsCascade: The index must be positive";
	if (index >= size) throw L"PolarBiquadsCascade: The index was too large";
#endif
	return biquadSection[index];
}

const IIR::PolarBiquadSection& PolarBiquadsCascade::operator[](long index) const
{
#ifdef _DEBUG
	if (index < 0) throw L"PolarBiquadsCascade: The index must be positive";
	if (index >= size) throw L"PolarBiquadsCascade: The index was too large";
#endif
	return biquadSection[index];
}

// Andrew G. Deczky 1972 Eq. 11
double PolarBiquadsCascade::GetMagnitude(double freq_rads)
{
	const int N = size;
	double alpha = 1.0;
	double tmp;
	double fap, fan;
	double fbp, fbn;
	int i;
	IIR::PolarBiquadSection* pbs = biquadSection;
	for (i = 0; i < N; i++, pbs++)
	{
		tmp = 1.0 + pbs->zeroR * pbs->zeroR;
		fap = tmp - 2.0*pbs->zeroR*cos(freq_rads + pbs->zeroAngle);
		fan = tmp - 2.0*pbs->zeroR*cos(freq_rads - pbs->zeroAngle);
		//
		tmp = 1.0 + pbs->poleR * pbs->poleR;
		fbp = tmp - 2.0*pbs->poleR*cos(freq_rads + pbs->poleAngle);
		fbn = tmp - 2.0*pbs->poleR*cos(freq_rads - pbs->poleAngle);
		//
		alpha *= ((fap*fan) / (fbp*fbn));
		alpha *= (pbs->k*pbs->k);
	}
	return sqrt(alpha);
}

// Andrew G. Deczky 1972 Eq. 12
double PolarBiquadsCascade::GetGroupDelay(double freq_rads)
{
	if (freq_rads == M_PI) freq_rads = 0.9999*M_PI;
	if (freq_rads == -M_PI) freq_rads = -0.9999*M_PI;
	const int N = size;
	double tao = 0.0;
	double tmp;
	IIR::PolarBiquadSection* pbs = biquadSection;
	for (int i = 0; i < N; i++, pbs++)
	{
		tmp = 1.0 - 2*pbs->poleR*cos(freq_rads - pbs->poleAngle) + pbs->poleR*pbs->poleR;
		if (tmp == 0.0) tmp = 1.0e-10;
		tao += ((1.0 - pbs->poleR*cos(freq_rads - pbs->poleAngle)) / tmp);
		//
		tmp = 1.0 - 2*pbs->poleR*cos(freq_rads + pbs->poleAngle) + pbs->poleR*pbs->poleR;
		if (tmp == 0.0) tmp = 1.0e-10;
		tao += ((1.0 - pbs->poleR*cos(freq_rads + pbs->poleAngle)) / tmp);
		//
		tmp = 1.0 - 2*pbs->zeroR*cos(freq_rads - pbs->zeroAngle) + pbs->zeroR*pbs->zeroR;
		if (tmp == 0.0) tmp = 1.0e-10;
		tao += ((1.0 - pbs->zeroR*cos(freq_rads - pbs->zeroAngle)) / tmp);
		//
		tmp = 1.0 - 2*pbs->zeroR*cos(freq_rads + pbs->zeroAngle) + pbs->zeroR*pbs->zeroR;
		if (tmp == 0.0) tmp = 1.0e-10;
		tao += ((1.0 - pbs->zeroR*cos(freq_rads + pbs->zeroAngle)) / tmp);
	}
	return tao;
}

bool PolarBiquadsCascade::Copy(const IIR::PolarBiquadsCascade& init)
{
	//________________________________________________ Empty Cascade
	if (init.size == 0)
	{
		Delete();
		return true;
	}
	if (init.biquadSection == NULL)
	{
		Delete();
		return true;
	}
	//________________________________________________ Cascade with sections
	if (SetSize(init.size) == false) return false;
	for (int i = 0; i < size; i++)
	{
		biquadSection[i] = init.biquadSection[i];
	}
	return true;
}

bool PolarBiquadsCascade::Copy(const IIR::BiquadsCascade& init)
{
	//________________________________________________ Empty Cascade
	if (init.size == 0)
	{
		Delete();
		return true;
	}
	if (init.biquadSection == NULL)
	{
		Delete();
		return true;
	}
	//________________________________________________ Cascade with sections
	if (SetSize(init.size) == false) return false;
	for (int i = 0; i < size; i++)
	{
		biquadSection[i] = init.biquadSection[i];
	}
	return true;
}


//_____________________________________________________________________ IIR::BiquadFilter
BiquadFilter::BiquadFilter(void)
{
}

BiquadFilter::~BiquadFilter(void)
{
}

void BiquadFilter::BilinearTransformation(int filter_type, double A, double B, double C, double& out_a, double& out_b, double& out_c)
{
	out_a = 4.0*A + 2.0*B + C;
	out_c = 4.0*A - 2.0*B + C;
	//
	if (filter_type == FILTER_TYPE_LOWPASS) // (4*A + 2*B + C) + (-8*A+2*C) z^-1 + (4*A - 2*B + C) z^-2
	{
		out_b = -8.0*A + 2.0*C;
	}
	else if (filter_type == FILTER_TYPE_HIGHPASS) // (4*A + 2*B + C) + (8*A-2*C) z^-1 + (4*A - 2*B + C) z^-2
	{
		out_b = 8.0*A - 2.0*C;
	}
}

bool BiquadFilter::CreateLowPass(double cutFreq_rads, double passBand_dB, double stopFreq_rads, double stopBand_dB)
{
	if (cutFreq_rads >= stopFreq_rads) return false;
	return Create(FILTER_TYPE_LOWPASS, cutFreq_rads, passBand_dB, stopFreq_rads, stopBand_dB);
}

bool BiquadFilter::CreateHighPass(double cutFreq_rads, double passBand_dB, double stopFreq_rads, double stopBand_dB)
{
	if (cutFreq_rads <= stopFreq_rads) return false;
	cutFreq_rads = M_PI - cutFreq_rads;
	stopFreq_rads = M_PI - stopFreq_rads;
	return Create(FILTER_TYPE_HIGHPASS, cutFreq_rads, passBand_dB, stopFreq_rads, stopBand_dB);
}

bool BiquadFilter::Create(int filter_type, double cutFreq_rads, double passBand_dB, double stopFreq_rads, double stopBand_dB)
{
	return false;
}

//_____________________________________________________________________ ButterworthFilter
ButterworthFilter::ButterworthFilter()
{
}

ButterworthFilter::~ButterworthFilter()
{
}

// Oppenheim Example 7.4 (page 420)
bool ButterworthFilter::Create(int filter_type, double cutFreq_rads, double passBand_dB, double stopFreq_rads, double stopBand_dB)
{
	passBand_dB = - fabs(passBand_dB);
	const double gainPassBand = pow(10.0, passBand_dB/20.0);
	const double gainStopBand = pow(10.0, stopBand_dB/20.0);
	const double wCut = 2.0*tan(cutFreq_rads/2.0);
	const double wStop = 2.0*tan(stopFreq_rads/2.0);
	const double pass = 1.0/(gainPassBand*gainPassBand) - 1.0;
	const double stop = 1.0/(gainStopBand*gainStopBand) - 1.0;
	const double n = log10(stop/pass)/ (2.0*log10(wStop/wCut));
	if (n <= 0) return false;
	const int in = (int)ceil(n);
	const int N = (in %2 == 0) ? in : in +1; // Force even
	if (N < 2) return false;
	if (N > 150) return false;
	const double omegaCut = wStop / pow(stop, 1.0/(2.0*N));
	const double sq_omegaCut = omegaCut*omegaCut;
	const int half = N/2;
	IIR::BiquadsCascade biquads;
	if (biquads.SetSize(half) == false) return false;
	double A, B, C, a, b, c;
	int i;
	const double delta = M_PI/N;
	double poleAngle;
	for (i = 0; i < half; i++)
	{
		//______________ Compute pole
		poleAngle = M_PI/2.0 + i*delta + delta/2.0;
		//______________________________________________________________ Bilinear Transformation
		//______________ Numerator: A*s*s + B*s + C which is  only a constant value
		A = 0.0;
		B = 0.0;
		C = sq_omegaCut;
		BilinearTransformation(filter_type, A, B, C, a, b, c);
		biquads[i].b0 = a;
		biquads[i].b1 = b;
		biquads[i].b2 = c;
		//______________ Denominator: A*s*s + B*s + C wich is s*s - 2*s*omegaCut*cos(poleAngle)+omegaCut*omegaCut
		A = 1.0;
		B = -2.0*omegaCut*cos(poleAngle);
		C = sq_omegaCut;
		BilinearTransformation(filter_type, A, B, C, a, b, c);
		biquads[i].a1 = b/a;
		biquads[i].a2 = c/a;
		//
		biquads[i].b0 /= a;
		biquads[i].b1 /= a;
		biquads[i].b2 /= a;
	}
	this->biquadsCascade = biquads;
	return true;
}

//_____________________________________________________________________ ChebyshevFilter
ChebyshevFilter::ChebyshevFilter()
{
}

ChebyshevFilter::~ChebyshevFilter()
{
}

// Oppenheim  B2 Chebyshev Filters, page 847
bool ChebyshevFilter::Create(int filter_type, double cutFreq_rads, double passBand_dB, double stopFreq_rads, double stopBand_dB)
{
	if (fabs(passBand_dB) > 3.0) passBand_dB = 3.0;
	const double epsilon = sqrt(pow(10.0, fabs(passBand_dB)/10.0) - 1.0);
	const double gainStopBand = pow(10.0, stopBand_dB/20.0);
	const double wCut = 2.0*tan(cutFreq_rads/2.0);
	const double wStop = 2.0*tan(stopFreq_rads/2.0);
	const double stop = 1.0/(gainStopBand*gainStopBand) - 1.0;
	const double n = acosh(sqrt(stop/(epsilon*epsilon)))/acosh(wStop/wCut);
	if (n <= 0) return false;
	const int in = (int)ceil(n);
	const int N = (in %2 == 0) ? in : in +1; // Force even
	if (N > 150) return false;
	if (N < 2) return false;
	const double omegaCut = wCut*cosh(acosh(1.0/epsilon)/N);
	const double alpha = 1.0/epsilon + sqrt(1.0+1.0/(epsilon*epsilon));
	const double beta = (pow(alpha, 1.0/N) - pow(alpha, -1.0/N))/2.0;
	const double gamma = (pow(alpha, 1.0/N) + pow(alpha, -1.0/N))/2.0;
	const double zero = 0.0;
	const int half = N/2;
	IIR::BiquadsCascade biquads;
	if (biquads.SetSize(half) == false) return false;
	double A, B, C, a, b, c;
	int i;
	const double delta = M_PI/N;
	double angle;
	double re, im;
	for (i = 0; i < half; i++)
	{
		//______________ Compute pole
		angle = M_PI/2.0 + i*delta + delta/2.0;
		re = beta*omegaCut*cos(angle);
		im = gamma*omegaCut*sin(angle);
		//______________________________________________________________ Bilinear Transformation
		//______________ Numerator: A*s*s + B*s + C which is  only a constant value
		A = 0.0;
		B = 0.0;
		C = re*re + im*im;
		BilinearTransformation(filter_type, A, B, C, a, b, c);
		biquads[i].b0 = a;
		biquads[i].b1 = b;
		biquads[i].b2 = c;
		//______________ Denominator: A*s*s + B*s + C
		A = 1.0;
		B = -2.0*re;
		C = re*re + im*im;
		BilinearTransformation(filter_type, A, B, C, a, b, c);
		biquads[i].a1 = b/a;
		biquads[i].a2 = c/a;
		//
		biquads[i].b0 /= a;
		biquads[i].b1 /= a;
		biquads[i].b2 /= a;
	}
	this->biquadsCascade = biquads;
	return true;
}

//_____________________________________________________________________ EllipticFilter
EllipticFilter::EllipticFilter()
{
}

EllipticFilter::~EllipticFilter()
{
}

//Digital Filter Designer'S Handbook: Featuring C Routines C. Britton Rorabaugh
bool EllipticFilter::Create(int filter_type, double cutFreq_rads, double passBand_dB, double stopFreq_rads, double stopBand_dB)
{
	const double Ap = fabs(passBand_dB);
	const double As = fabs(stopBand_dB);
	const double wp = 2.0*tan(cutFreq_rads/2.0);
	const double ws = 2.0*tan(stopFreq_rads/2.0);
	const double k = wp/ws;
	const double kk = sqrt(sqrt(1.0-k*k));
	const double u = (1.0 - kk)/(2.0*(1.0 + kk));  // EQ. (5.2)
	const double u4 = u*u*u*u;
	const double u5 = u4*u;
	const double u9 = u4*u5;
	const double u13 = u4*u9;
	const double q = u + 2.0*u5 + 15.0*u9 + 150.0*u13; // EQ. (5.1)
	const double D = (pow(10.0, As/10.0) -1.0)/(pow(10.0, Ap/10.0) -1.0); // EQ. (5.3)
	const double dn = (log10(16.0*D))/(log10(1.0/q));
	const int in = (int)ceil(dn); // EQ. (5.4)
	const int n = (in %2 == 0) ? in : in +1; // Force even
	const double actual_As = 10.0*log10(1.0 + (pow(10.0, Ap/10.0)-1.0)/(16.0*pow(q, (double)n))); // EQ. (5.5)
	const double alpha = sqrt(wp*ws); // EQ. (5.8)
	const double vv = pow(10.0, Ap/20.0);
	const double V = (1.0/(2.0*n))*log((vv + 1.0)/(vv - 1.0)); // EQ. (5.12)
	//______________________________________________________________ EQ. (5.13)
	double pp;
	int m;
	//____________________ Numerator
	double sum_numerator = 0.0;
	for (m = 0; m < 100; m++)
	{
		pp = pow(q, m*(m+1))*sinh((2*m+1)*V);
		if (m%2 != 0) pp = -pp;
		sum_numerator += pp;
		if (fabs(pp)/fabs(sum_numerator) < 1.0e-10*fabs(sum_numerator)) break; // done
	}
	//____________________ Denominator
	double sum_denominator = 0.5;
	for (m = 1; m < 100; m++)
	{
		pp = pow(q, m*m)*cosh(2*m*V);
		if (m%2 != 0) pp = -pp;
		sum_denominator += pp;
		if (fabs(pp)/fabs(sum_denominator) < 1.0e-10*fabs(sum_denominator)) break; // done
	}
	const double p0 = fabs(sqrt(sqrt(q))*sum_numerator/sum_denominator);
	const double W = sqrt((1.0 + (p0*p0)/k)  *  (1.0 + k*p0*p0)); //EQ. (5.14)
	const bool is_n_even = (n%2 == 0);
	const int r = is_n_even ? n/2: (n-1)/2;
	//______________________________________________________________ EQ. (5.15)
	double * X = new double[r];
	int i;
	double mu;
	for (i = 1; i <= r; i++)
	{
		mu = is_n_even ? (i-0.5) : (double)i;
		//____________________ Numerator
		sum_numerator = 0.0;
		for (m = 0; m < 100; m++)
		{
			pp = pow(q, m*(m+1))*sin((2*m+1)*mu*M_PI/n);
			if (m%2 != 0) pp = -pp;
			sum_numerator += pp;
			if (fabs(pp)/fabs(sum_numerator) < 1.0e-10*fabs(sum_numerator)) break; // done
		}
		//____________________ Denominator
		sum_denominator = 0.0;
		for (m = 1; m < 100; m++)
		{
			pp = pow(q, m*m)*cos(2.0*m*mu*M_PI/n);
			if (m%2 != 0) pp = -pp;
			sum_denominator += pp;
			if (fabs(pp)/fabs(sum_denominator) < 1.0e-10*fabs(sum_denominator)) break; // done
		}
		X[i-1] = (2.0*sqrt(sqrt(q))*sum_numerator) / (1.0 + 2.0*sum_denominator);
	}
	//______________________________________________________________ EQ. (5.16)
	double * Y= new double[r];
	for (i = 1; i <= r; i++)
	{
		Y[i-1] = sqrt((1.0 - (X[i-1]*X[i-1]/k))  *  (1.0 - k*X[i-1]*X[i-1]));
	}
	//______________________________________________________________ EQ. (5.17)
	double * a_= new double[r];
	for (i = 1; i <= r; i++)
	{
		a_[i-1] = 1.0/(X[i-1]*X[i-1]);
	}
	//______________________________________________________________ EQ. (5.18)
	double * b_= new double[r];
	for (i = 1; i <= r; i++)
	{
		b_[i-1] = 2.0*p0*Y[i-1]/(1.0 + p0*p0*X[i-1]*X[i-1]);
	}
	//______________________________________________________________ EQ. (5.19)
	double * c_= new double[r];
	double tmp;
	for (i = 1; i <= r; i++)
	{
		pp = p0*Y[i-1];
		tmp = X[i-1]*W;
		c_[i-1] = pp*pp + tmp*tmp;
		//
		tmp = 1.0 + p0*p0*X[i-1]*X[i-1];
		c_[i-1] /= (tmp*tmp);
	}
	//______________________________________________________________ EQ. (5.20)
	tmp = 1.0;
	for (i = 1; i <= r; i++) tmp *= (c_[i-1]/a_[i-1]);
	const double H0 = (is_n_even == true) ? (pow(10.0, -Ap/20.0)*tmp) : (p0*tmp);
	//______________________________________________________________ EQ. (5.22)
	const double k0 = H0; //(is_n_even == true) ? H0 : ( (H0*alpha)/(s + alpha*p0));
	for (i = 0; i < r; i++)
	{
		a_[i] = alpha*alpha*a_[i];
		b_[i] = alpha*b_[i];
		c_[i] = alpha*alpha*c_[i];
	}
	//______________________________________________________________ Bilinear Transformation
	const int N = n;
	const int half = N/2;
	IIR::BiquadsCascade biquads;
	if (biquads.SetSize(half))
	{
		double A, B, C, a, b, c;
		const double delta = M_PI/N;
		for (i = 0; i < half; i++)
		{
			//______________ Numerator: A*s*s + B*s + C
			A = 1.0;
			B = 0.0;
			C = a_[i];
			BilinearTransformation(filter_type, A, B, C, a, b, c);
			biquads[i].b0 = a;
			biquads[i].b1 = b;
			biquads[i].b2 = c;
			//______________ Denominator: A*s*s + B*s + C
			A = 1.0;
			B = b_[i];
			C = c_[i];
			BilinearTransformation(filter_type, A, B, C, a, b, c);
			biquads[i].a1 = b/a;
			biquads[i].a2 = c/a;
			//
			biquads[i].b0 /= a;
			biquads[i].b1 /= a;
			biquads[i].b2 /= a;
		}
		//_______________________________ Incorporate k0 in first section
		biquads[0].b0 *= k0;
		biquads[0].b1 *= k0;
		biquads[0].b2 *= k0;
		//________________
		this->biquadsCascade = biquads;
	}
	if (X != NULL) delete [] X;
	if (Y != NULL) delete[] Y;
	if (a_ != NULL) delete[] a_;
	if (b_ != NULL) delete[] b_;
	if (c_ != NULL) delete[] c_;
	return true;
}

}
