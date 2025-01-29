#include <iostream>
#include <fstream>
#include <math.h>
#include <complex>

using namespace std;

/*
/Erstellt das struct Kopplungsmatrix, welches die Spins s_i, s_j miteinander Verbindet.
/Die Variablen spini und spinj stehen für die jeweiligen Spins,
/für die die Kopplungsmatrix gilt.
/Die 3x3-Matrizen werden als 9-elementiges Array gespeichert.
/Die Kopplungsmatrizen werden durch den Konstruktor direkt als
/J_ij(+,-,z) Matrix abgespeichert.
*/

struct Kopplungsmatrix 
{
	short spini,spinj;
	complex<double> jMatrix[9];
	explicit Kopplungsmatrix(double *readInFileJ);
	explicit Kopplungsmatrix();
};

Kopplungsmatrix::Kopplungsmatrix(double *readInFileJ) 
{
	spini = (short) (readInFileJ[0] - 1);
	spinj = (short) (readInFileJ[1] - 1);
	//Berechnen der Kombinationen: ++,+-,-+,--
	//++
	jMatrix[0].real(0.25*(readInFileJ[2]-readInFileJ[6]));
	jMatrix[0].imag(-0.25*(readInFileJ[3]+readInFileJ[5]));
	//+-
	jMatrix[1].real(0.25*(readInFileJ[2]+readInFileJ[6]));
	jMatrix[1].imag(0.25*(readInFileJ[3]-readInFileJ[5]));
	//-+
	jMatrix[3].real(0.25*(readInFileJ[2]+readInFileJ[6]));
	jMatrix[3].imag(-0.25*(readInFileJ[3]-readInFileJ[5]));
	//--
	jMatrix[4].real(0.25*(readInFileJ[2]-readInFileJ[6]));
	jMatrix[4].imag(0.25*(readInFileJ[3]+readInFileJ[5]));
	//Berechnen der Kombinationen: +z,-z,z+,z-
	//+z
	jMatrix[2].real(0.5*(readInFileJ[4]));
	jMatrix[2].imag(-0.5*(readInFileJ[7]));
	//-z
	jMatrix[5].real(0.5*(readInFileJ[4]));
	jMatrix[5].imag(0.5*(readInFileJ[7]));
	//z+
	jMatrix[6].real(0.5*(readInFileJ[8]));
	jMatrix[6].imag(-0.5*(readInFileJ[9]));
	//z-
	jMatrix[7].real(0.5*(readInFileJ[8]));
	jMatrix[7].imag(0.5*(readInFileJ[9]));
	//Berechnen von zz
	jMatrix[8].real(readInFileJ[10]);
}

Kopplungsmatrix::Kopplungsmatrix()
{
	spini=0;
	spinj=0;
	for (int i = 0; i < 9; ++i)
	{
		*(jMatrix+i) = 0;
		
	}
}

