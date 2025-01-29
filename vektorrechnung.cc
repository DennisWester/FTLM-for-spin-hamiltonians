#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <complex>

using namespace std;

inline long double normierung(complex<double> *psi, const int index, const unsigned long hilbertraumDim) 
{
	long double normierung = 0;
	for (unsigned long i = 0; i < hilbertraumDim; i++)
	{
		normierung += norm(psi[i + (index*hilbertraumDim)]);	
		
	}
	normierung = sqrt(normierung);
	return normierung;
}

//Komplexe Vektor-Vektor-Multiplikation (Skalarprodukt)
inline complex<double> complexPsiProduct(const complex<double> *psi, const unsigned long hilbertraumDim,
										 const int index1, const int index2)
{
	complex<double> sum = 0;
	for (unsigned long i = 0; i < hilbertraumDim; i++)
	{
		sum += conj(psi[i + (index1 * hilbertraumDim)]) * (psi[i + (index2 * hilbertraumDim)]);
	}
	return sum;
}

//Gram-Schmidtsche Orthogonalisierung
// inline void orthogonalisiere(complex<float> *psi, const complex<float> alpha, const complex<float> beta,
// 						const int index_berechnung, const unsigned long hilbertraumDim)
// {
// 	unsigned long index_tochter = hilbertraumDim * ((index_berechnung + 2) % 3);
// 	unsigned long index_enkel = hilbertraumDim * ((index_berechnung + 1) % 3);
// 	for (unsigned long i = 0; i < hilbertraumDim; i++)
// 	{
// 		psi[i + (index_berechnung * hilbertraumDim)] = psi[i + (index_berechnung * hilbertraumDim)] - (psi[index_tochter + i] * alpha + psi[index_enkel + i] * beta);
// 	}
// }

// void evnormtest(const double *evmat, int size)
// {
	
// 	for (int i = 0; i < size; ++i)
// 	{
// 		double sum = 0;
// 		for (int j = 0; j < size; ++j)
// 		{
// 			sum += evmat[i*size + j] * evmat[i*size + j];
// 		}
// 		cout << setprecision(3);
// 		cout << "Norm(EV[" << i << "]: " << sum << endl;
// 	}
// }

// void normtest(const complex<double> *psi, const unsigned long index, const unsigned long hilbertraumDim)
// {
// 	complex<double> sum = 0;
// 	complex<double> temp;
// 	for (unsigned long i = 0; i < hilbertraumDim; ++i)
// 	{
// 		temp = psi[(index * hilbertraumDim) + i];
// 		sum += conj(temp) * psi[(index * hilbertraumDim) + i];
// 	}
// 	cout << "||psi>| = " << sum << "\n";
// }

// void orthogonaltest(const complex<double> *psi, const unsigned long index1, const unsigned long index2, const unsigned long hilbertraumDim)
// {
// 	complex<double> sum = 0;
// 	complex<double> temp;
// 	for (unsigned long i = 0; i < hilbertraumDim; i++)
// 	{
// 		temp = psi[(index1 * hilbertraumDim) + i];
// 		sum += conj(temp) * psi[(index2 * hilbertraumDim) + i];
// 	}
// 	if (sum.real() < 0.0000001 && sum.imag() < 0.0000001)
// 	{
// 		cout << "Die Vektoren sind (ziemlich) orthogonal\n";
// 	}
// 	else
// 	{
// 		cout << sum << " Die Vektoren sind nicht orthogonal\n";
// 	}
// }

// int testEVorder(const double *trimatrix,const double *evmatrix,const double *ew, const int size)
// {
// 	double matrix[size*size];
// 	for (int i = 0; i < size * size; i++)
// 	{
// 		matrix[i] = trimatrix[i];
// 	}
// 	//Aufbauen der gesamten Matrix zur Matrix - Vektor Multiplikation
// 	int matrixEntry = 0, tempEntry = 0;
// 	for (int i = 0; i < size; i++)
// 	{
// 		matrix[matrixEntry] = matrix[matrixEntry];
// 		matrixEntry++;
// 		matrix[matrixEntry] = matrix[matrixEntry];
// 	 	tempEntry = matrixEntry;
// 		matrixEntry++;
// 		for (int p = 0; p < size-1; p++)
// 		{
// 			matrix[matrixEntry] = 0;
// 			matrixEntry++;
// 		}
// 		matrix[matrixEntry] = matrix[tempEntry];
// 	}
// 	matrix[matrixEntry] = matrix[matrixEntry];

// 	//Eigentliche Testroutine
// 	double ergev[size];
// 	int result = 0;
// 	for (int k = 0; k < size; k++)
// 	{
// 		for (int i = 0; i < size; i++)
// 		{
// 			ergev[i] = 0;
// 			for (int j = 0; j < size; j++)
// 			{
// 				ergev[i] += matrix[j + i * size] * evmatrix[j * size + k];
// 			}
// 			ergev[i] /= ew[i];
// 			if (ergev[i] == evmatrix[i * size + k])
// 			{
// 				result++;
// 			}
// 		}
// 	}
// 	if (result == size * size)
// 	{
// 		cout << "juhu ";
// 		return 1;
// 	} else
// 	{
// 		cout << "EV-System in der falschen Reihenfolge\n";
// 		return 0;
// 	}

// }