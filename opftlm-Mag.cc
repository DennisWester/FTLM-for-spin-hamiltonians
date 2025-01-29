//ANNAHME: Die Einzelspin-Anisotropie ist symmetrisch (Also die jeweiligen Tensoren)
//Die Inhalte des Hauptprogramms werden am Ende Kommentiert.
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <complex>
#include <time.h>
#include <random>
#include <iomanip>
#include <cstring>
#include <omp.h>
#include <boost/scoped_array.hpp>
#include <boost/shared_array.hpp>
#include "nummer-basis.cc"
#include "vektorrechnung.cc"
#include "kopplungsmatrizenJ.cc"
#include "holeInfosAusDatei.cc"

using namespace std;

#ifdef __cplusplus
extern "C" {
#endif
/*#ifdef USE_NAG
extern void f06qff_(char* MATRIX, int* M, int* N, double* A,
int *LDA, double *B, int* LDB);
extern void f02faf_(char* JOB, char* UPLO, int* N, double* A,
int* LDA, double* W, double* WORK, int* LWORK, int* IFAIL);
#else*/
extern void dcopy_(int* N, double* DX, int* INCX, double* DY, int* INCY);
extern void dsyev_(char* JOBZ, char* UPLO, int* N, double* A, int* LDA,
double* W, double* WORK, int* LWORK, int* INFO);
//#endif
#ifdef __cplusplus
}
#endif
const double muB_durch_kB = 0.6717141002;
const double k_Bdurchmu_B = 1.48872872505340; 
const double cm_in_kelvin = 1.;
const double PI = 3.14159265358979323846;
const int threadzahl = 10;



/*
/Speichert die in der Datei befindlichen Matrizen J_ij(x,y,z)
/in das Array kopplungen von Kopplungsmatrizen [J_ij(+,-,z)].
/Dabei wird die struct Kopplugsmatrix aus dem Programmteil spinwirkungJ.cc benutzt.
*/
void speicherJMatrizen(Kopplungsmatrix *kopplungen, double *xyzkopp, int jMatrizenAnzahl, ifstream &ifs)
{
 	ifs.seekg(ios::beg);
	char ctemp;
	int itemp;
	while(!ifs.eof())
	{
		ifs >> ctemp;
		if (ctemp == 'j' || ctemp == 'J')
		 {
		 	break;
		 }else if (ifs.bad())
		 {
		 	exit(0);
		 }
	}
	ifs >> itemp;

	double hilfsArray[11];
	for (int i = 0; i < jMatrizenAnzahl; i++)//Schleife zum durchlaufen aller J-Matrizen
	{
		for (int j = 0; j < 11; j++)// Schleife zum durchlaufen der Zeilen einer J-Matrix
		{
			ifs >> hilfsArray[j];
		}
		xyzkopp[i*11] = hilfsArray[0];
		xyzkopp[i*11 + 1] = hilfsArray[1];
		for (int j = 2; j < 11; ++j)
		{
			hilfsArray[j] *= cm_in_kelvin;
			xyzkopp[j + i*11] = hilfsArray[j];
		}
		kopplungen[i] = Kopplungsmatrix(hilfsArray);
	}
}

void infoausgabe(const unsigned short *spins, const double *gFaktoren, const int random_vectors, const int lanczosSteps, const unsigned long hilbertraumDimension, const int spinAnzahl,
				 const int jMatrizenAnzahl, const double *xyzkopp, const vector<double> &vecbfeld, const vector<double> &temperatur, const string text, const string magtext,
				 const string inputname, const long samen, const double exec_time, const double cpu_time)
{

	string infodat = "FTLM-Informationen_";
	infodat += magtext;
	ofstream outinfo(infodat);
	if (!outinfo.is_open() )
	{
		cerr << "Keine Outputdatei geöffnet" << endl;
	}
	outinfo << "Die Inputdatei heißt: " << inputname << endl;
	outinfo << "Die Outputdateien heißen: " << text << "  " << magtext << endl;
	outinfo << "Das Programm hat: " << exec_time << " Minuten gebraucht und die CPUs " << cpu_time << " min gebraten" << endl;
	outinfo << endl;
	outinfo << "Die Pseudozufallszahlen wurden mit dem Seed: " << samen << " erstellt." << endl;
	outinfo << "Es wurden R=" << random_vectors << " Zufallsvektoren verwendet und N_L=" << lanczosSteps << " lanczosSteps. Die Hilbertraumdimension beträgt: " << hilbertraumDimension << endl;
	outinfo <<  "Der Algorithmus wurde für die Felder ";
	for (double feldstaerke : vecbfeld)
	{
		outinfo << "T, B= " << feldstaerke; 
	}
	outinfo << "T durchgeführt." << endl;
	outinfo << "Die Spins:" << setprecision(3) << endl;
	for (int i = 0; i < spinAnzahl; ++i)
	{
		double s = 0.5 * spins[i];
		outinfo << s << " ";
	}
	outinfo << endl;
	outinfo << "Die g-Faktoren (skalar) lauten:" << endl;
	for (int i = 0; i < spinAnzahl; ++i)
	{
		outinfo << gFaktoren[i] << " ";
	}
	outinfo << endl;
	outinfo << "Die Kopplungsmatrizen in xyz-Darstellung: ";
	if (cm_in_kelvin > 1.1)
	{
		outinfo << "(in cm^-1 angegeben)" << endl;
	}
	else{
		outinfo << "(in K angegeben)" << endl;
	}
	for (int i = 0; i < jMatrizenAnzahl; ++i)
	{
		int s = xyzkopp[i*11];
		outinfo << s << " ";
		s = xyzkopp[i*11 + 1];
		outinfo << s << endl;
		for (int j = 0; j < 3; ++j)
		{
			for (int k = 0; k < 3; ++k)
			{
				outinfo << xyzkopp[k + 2 + j*3 + i*11] << " ";
			}
			outinfo << endl;
		}
	}
	outinfo << "Folgende Temperaturen wurden dabei eingestellt: (in K)" << endl;
	for (double temp : temperatur)
	{
		outinfo << temp << endl;
	}
	outinfo.close();
}

//new_handler zum behandeln von Ausnahmen der Speicherallokierung von new
void outOfMemory()
{
	cerr << "Speicherplatzanfrage konnte nicht erfüllt werden\n";
	abort();
}

inline double s_plus(const unsigned short *spins, const unsigned short *a, const int index)
{
	double spin = spins[index] * 0.5;
	double mi = spin - a[index];
	double plus = sqrt(spin * (spin + 1.) - mi * (mi + 1.));
	return plus;
}

inline double s_minus(const unsigned short *spins, const unsigned short *a, const int index)
{
	double spin = spins[index] * 0.5;
	double mi = spin - a[index];
	double plus = sqrt(spin * (spin + 1.) - mi * (mi - 1.));
	return plus;
}

inline double s_z(const unsigned short *spins, const unsigned short *a, const int index)
{
	double spin = spins[index] * 0.5;
	double mi = spin - a[index];
	return mi;
}

long double zufallspsi(complex<double> *psi, const unsigned long hilbertraumDimension, const int r, const long samen, const bool zeitumkehrvec)
{
	double speicher_random;
	//Der Lanczos Algorithmus
	//Einspeichern der Random-Zahlen in |psi_1>
	if (zeitumkehrvec == 0)
	{
		auto zufall = bind(normal_distribution<double>(0,1.0), mt19937(samen));
		for (unsigned long counter = 0; counter < 2 * hilbertraumDimension * (r/2); ++counter)
		{ //Das Zufallszahlen-Array wieder auf die richtige Stelle setzen
			speicher_random = zufall();
		}
		for (unsigned long counter = 0; counter < hilbertraumDimension; counter++)
		{
			psi[counter].real(zufall());
			psi[counter].imag(zufall());
		}
	} else
	{
		auto zufall = bind(normal_distribution<double>(0,1.0), mt19937(samen));
		for (unsigned long counter = 0; counter < 2 * hilbertraumDimension * (r/2); ++counter)
		{ //Das Zufallszahlen-Array wieder auf die richtige Stelle setzen
			speicher_random = zufall();
		}
		for (unsigned long counter = 1; counter <= hilbertraumDimension; counter++)
		{
			psi[hilbertraumDimension - counter].real(zufall());
			psi[hilbertraumDimension - counter].imag(-zufall());
		}
	}			
	//Normierung des Startvektors
	long double normierungZufall = normierung(psi, 0, hilbertraumDimension);
	normierungZufall = 1. / normierungZufall;
	for (unsigned long counter = 0; counter < hilbertraumDimension; counter++)
	{
		psi[counter] *= normierungZufall;
	}
	return normierungZufall;
}

/*
/Implementierung der Funktion die den Hamiltonian in Form der J_ij
/auf einen Eintrag |a_1,...,a_N> wirken lässt.
/Die J_ij befinden sich schon in der gewünschten Form, d.h.
/mit den Richtungen +/-/z anstatt x/y/z.
/Dabei muss zwischen allen Möglichkeiten (9) unterschieden werden,
/dass s_i oder s_j entweder Maximal oder Minimal sind und somit Terme wegfallen.
*/
void HamiltonMalEintrag(unsigned short *a, const Kopplungsmatrix *kopplungen, complex<double> *psi, const unsigned short *spins, 
						const double *gFaktoren, const int spinAnzahl, const double *b_Feld, const unsigned long ergebnisshift, 
						const unsigned long tochtershift, const unsigned long h, const int jMatrizenAnzahl, const unsigned long *spindimtable,
						complex<double> *hdimserw, const complex<double> *rvec, const unsigned long hilbertraumDimension)
{
	unsigned long hilfsNummer;//Hilfsnummer für VecToNum also Nummer aus Psi
	double sj, si;//Zum einspeichern der Faktoren durch s_plus, s_minus, s_z
	complex<double> zeemanTerm;//Hilfsvariable zur Berechnung der Wirkung des Zeeman-Terms 
	for (int i = 0; i < spinAnzahl; i++)
	{

		complex<double> serwhelpp = 0, serwhelpm = 0;
		// Berechnung des Zeeman-Term Anteils der Wirkung
		// Berechnung s_+ Term des Zeeman Terms der durch (B_1/2)-(B_2*i/2) gegeben ist
		if (a[i] != spins[i])
		{
			a[i] ++;
			si = s_plus(spins, a, i);
			hilfsNummer = vecToNum(a, spins, spinAnzahl, spindimtable);
			//Wirkung s_+ auf Zustand
			//Multiplikation mit (B_1/2)-(i*B_2/2)
			zeemanTerm.real(si * b_Feld[i*3 + 0] * 0.5 * gFaktoren[i]);
			zeemanTerm.imag(-si * b_Feld[i*3 + 1] * 0.5 * gFaktoren[i]);
			psi[ergebnisshift+h] += (complex<double>) zeemanTerm * muB_durch_kB * psi[tochtershift + hilfsNummer];
			serwhelpp = (complex<double>) conj(psi[tochtershift + hilfsNummer]) * si * gFaktoren[i] * rvec[h];
			a[i] --;
		}
		//Berechnung s_- Term des Zeeman Terms der durch (B_1/2)+(i*B_2/2) gegeben ist
		if (a[i] != 0)
		{
			a[i] --;
			si = s_minus(spins, a, i);
			hilfsNummer = vecToNum(a, spins, spinAnzahl, spindimtable);
			//Wirkung s_+ auf Zustand
			//Multiplikation mit (B_1/2)-(B_2*i/2)
			zeemanTerm.real(si * b_Feld[i*3 + 0] * 0.5 * gFaktoren[i]);
			zeemanTerm.imag(si * b_Feld[i*3 + 1] * 0.5 * gFaktoren[i]);
			psi[ergebnisshift+h] += (complex<double>) zeemanTerm * muB_durch_kB * psi[tochtershift + hilfsNummer];
			serwhelpm = (complex<double>) conj(psi[tochtershift + hilfsNummer]) * si * gFaktoren[i] * rvec[h];
			a[i] ++;
		}
		//Berechnung s_z Term des Zeeman Terms der durch B_3 gegeben ist
		si = s_z(spins, a, i);
		sj = (double) si * b_Feld[i*3 + 2] * gFaktoren[i];
		psi[ergebnisshift + h] += (complex<double>) sj * muB_durch_kB * psi[tochtershift + h];

		hdimserw[h + 2*hilbertraumDimension] += (complex<double>) conj(psi[tochtershift + h]) * si * gFaktoren[i] * rvec[h];
		hdimserw[h + 1*hilbertraumDimension] += (complex<double>) serwhelpm;
		hdimserw[h] += (complex<double>) serwhelpp;

	}
			
	for (int kopplungsIndex = 0; kopplungsIndex < jMatrizenAnzahl; kopplungsIndex++)
	{
		int i = kopplungen[kopplungsIndex].spini;
		int j = kopplungen[kopplungsIndex].spinj;
		if (i == j)//Hier beginnen die Auswahlregeln für die Einzelspin-Anisotropie
		{ // Es wird im Folgenden angenommen, dass der Anisotropietensor symmetrisch ist !!!
			//s_+ s_+ Term wenn m_i<(max-1)
			if (a[i] < spins[i] - 1)
			{
				a[i] ++;
				si = s_plus(spins, a, i);
				a[i] ++;
				si *= s_plus(spins, a, i);
				hilfsNummer = vecToNum(a, spins, spinAnzahl, spindimtable);
				psi[ergebnisshift + h] += (complex<double>) kopplungen[kopplungsIndex].jMatrix[0] * si * psi[tochtershift + hilfsNummer];
				a[i] -= 2;
			}
			//s_- s_- Term wenn m_i > (min+1)
			if (a[i] > 1)
			{
				a[i] --;
				si = s_minus(spins, a, i);
				a[i] --;
				si *= s_minus(spins, a, i);
				hilfsNummer = vecToNum(a, spins, spinAnzahl, spindimtable);
				psi[ergebnisshift + h] += (complex<double>) kopplungen[kopplungsIndex].jMatrix[4] * si * psi[tochtershift + hilfsNummer];
				a[i] += 2;
			}
			//s_+ s_z und s_z s_+ Term wenn m_i<max
			if (a[i] != spins[i])
			{
				si = s_z(spins, a, i);//s_+ s_z
				a[i] ++;
				si *= s_plus(spins, a, i); //s_+ s_z
				sj = s_plus(spins, a, i); // s_z s_+
				sj *= s_z(spins, a, i); //s_z s_+
				hilfsNummer = vecToNum(a, spins, spinAnzahl, spindimtable);
				psi[ergebnisshift+h] += (complex<double>) kopplungen[kopplungsIndex].jMatrix[2] * (si + sj) * psi[tochtershift + hilfsNummer];
				a[i] --;
			}else { //s_- s_+ Term wenn s_+ nicht anwendbar -> kommutieren zu s_+ s_- (zusätzlich noch normale Anwendung von s_+ s_-)
				si = spins[i];
				psi[ergebnisshift+h] += (complex<double>) kopplungen[kopplungsIndex].jMatrix[1] * si * psi[tochtershift + h];
			}
			//s_- s_z Term wenn m_i>min
			if (a[i] != 0)
			{
				si = s_z(spins, a, i);
				a[i] --;
				si *= s_minus(spins, a, i);
				sj = s_minus(spins, a, i);
				sj *= s_z(spins, a, i);
				hilfsNummer = vecToNum(a, spins, spinAnzahl, spindimtable);
				psi[ergebnisshift+h] += (complex<double>) kopplungen[kopplungsIndex].jMatrix[5] * (si + sj) * psi[tochtershift + hilfsNummer];
				a[i] ++;
			} else { //s_+ s_- Term wenn s_- nicht anwendbar -> kommutieren zu s_- s_+ (zusätzlich noch normale Anwendung von s_- s_+)
				si = spins[i];
				psi[ergebnisshift+h] += (complex<double>) kopplungen[kopplungsIndex].jMatrix[3] * si * psi[tochtershift + h];
			}
			if (a[i] != 0 && a[i] != spins[i])
			{
				//s_+ s_- und s_- s_+ Terme 
				a[i] --;
				si = s_minus(spins, a, i);
				a[i] ++;
				si *= s_plus(spins, a, i);
				a[i] ++;
				sj = s_plus(spins, a, i);
				a[i] --;
				sj *= s_minus(spins, a, i);
				hilfsNummer = vecToNum(a, spins, spinAnzahl, spindimtable);
				psi[ergebnisshift+h] += (complex<double>) kopplungen[kopplungsIndex].jMatrix[1] * (si + sj) * psi[tochtershift + hilfsNummer];
			}
		} 
		else // Hier beginnen die Kopplungsterme
		{
			//s_+ s_+ Term wenn m_i<max und m_j<max
			if (a[i] != spins[i] && a[j] != spins[j])
			{
				a[i] ++;
				a[j] ++;
				si = s_plus(spins, a, i);
				sj = s_plus(spins, a, j);
				hilfsNummer = vecToNum(a, spins, spinAnzahl, spindimtable);
				psi[ergebnisshift+h] += (complex<double>) kopplungen[kopplungsIndex].jMatrix[0] * (si * sj) * psi[tochtershift + hilfsNummer];
				a[i] --;
				a[j] --;
			}
			//s_+ s_- Term wenn m_i<max und m_j>min
			if (a[i] != spins[i] && a[j] != 0)
			{
				a[i] ++;
				a[j] --;
				si = s_plus(spins, a, i);
				sj = s_minus(spins, a, j);
				hilfsNummer = vecToNum(a, spins, spinAnzahl, spindimtable);
				psi[ergebnisshift+h] += (complex<double>) kopplungen[kopplungsIndex].jMatrix[1] * (si * sj) * psi[tochtershift + hilfsNummer];
				a[i] --;
				a[j] ++;
			}
			//s_+ s_z Term wenn m_i<max
			if (a[i] != spins[i])
			{
				a[i] ++;
				si = s_plus(spins, a, i);
				sj = s_z(spins, a, j);
				hilfsNummer = vecToNum(a, spins, spinAnzahl, spindimtable);
				psi[ergebnisshift+h] += (complex<double>) kopplungen[kopplungsIndex].jMatrix[2] * (si * sj) * psi[tochtershift + hilfsNummer];
				a[i] --;
			}
			//s_- s_+ Term wenn m_i>min und m_j<max
			if (a[i] != 0 && a[j] != spins[j])
			{
				a[i] --;
				a[j] ++;
				si = s_minus(spins, a, i);
				sj = s_plus(spins, a, j);
				hilfsNummer = vecToNum(a, spins, spinAnzahl, spindimtable);
				psi[ergebnisshift+h] += (complex<double>) kopplungen[kopplungsIndex].jMatrix[3] * (si * sj) * psi[tochtershift + hilfsNummer];
				a[i] ++;
				a[j] --;
			}
			//s_- s_- Term wenn m_i>min und m_j>min
			if (a[i] != 0 && a[j] != 0)
			{
				a[i] --;
				a[j] --;
				si = s_minus(spins, a, i);
				sj = s_minus(spins, a, j);
				hilfsNummer = vecToNum(a, spins, spinAnzahl, spindimtable);
				psi[ergebnisshift+h] += (complex<double>) kopplungen[kopplungsIndex].jMatrix[4] * (si * sj) * psi[tochtershift + hilfsNummer];
				a[i] ++;
				a[j] ++;
			}
			//s_- s_z Term wenn m_i>min
			if (a[i] != 0)
			{
				a[i] --;
				si = s_minus(spins, a, i);
				sj = s_z(spins, a, j);
				hilfsNummer = vecToNum(a, spins, spinAnzahl, spindimtable);
				psi[ergebnisshift+h] += (complex<double>) kopplungen[kopplungsIndex].jMatrix[5] * (si * sj) * psi[tochtershift + hilfsNummer];
				a[i] ++;
			}
			//s_z s_+ Term wenn m_j<max
			if (a[j] != spins[j])
			{
				a[j] ++;
				si = s_z(spins, a, i);
				sj = s_plus(spins, a, j);
				hilfsNummer = vecToNum(a, spins, spinAnzahl, spindimtable);
				psi[ergebnisshift+h] += (complex<double>) kopplungen[kopplungsIndex].jMatrix[6] * (si * sj) * psi[tochtershift + hilfsNummer];
				a[j] --;
			}
			//s_z s_- Term wenn m_j>min
			if (a[j] != 0)
			{
				a[j] --;
				si = s_z(spins, a, i);
				sj = s_minus(spins, a, j);
				hilfsNummer = vecToNum(a, spins, spinAnzahl, spindimtable);
				psi[ergebnisshift+h] += (complex<double>) kopplungen[kopplungsIndex].jMatrix[7] * (si * sj) * psi[tochtershift + hilfsNummer];
				a[j] ++;
			}
		}
		//s_z s_z Term der in jedem Fall vorkommt
		si = s_z(spins, a, i);
		sj = s_z(spins, a, j);
		psi[ergebnisshift+h] += (complex<double>) kopplungen[kopplungsIndex].jMatrix[8] * (si * sj) * psi[tochtershift + h];
	}
	
}



/*
/Implementierung der Funktion zur Berechnung von H|psi>
/Diese Funktion steuert die anderen Subroutinen zum Nummer-Basiswechsel
/und zur Wirkung der einzelnen J_ij auf |a_1,...,a_N>
*/
void hamiltonMalPsi(complex<double> *psi, const Kopplungsmatrix *kopplungen, const unsigned short *spins, 
					const double *gFaktoren, const int spinAnzahl, const unsigned long hilbertraumDimension, 
					const double *b_Feld, const unsigned long k, const int jMatrizenAnzahl, const unsigned long *spindimtable,
					complex<double> *hdimserw, const complex<double> *rvec)
{
	unsigned long ergebnisshift = hilbertraumDimension * ((k+1) % 3);
	unsigned long tochtershift = hilbertraumDimension * (k % 3);
	for (unsigned long h = 0; h < hilbertraumDimension; h++)
	{
		psi[ergebnisshift+h] = 0;
	}
	#pragma omp parallel
	{
		#pragma omp for 
		for (unsigned long h = 0; h < hilbertraumDimension; h++)
		{
			unsigned short a[spinAnzahl];//Initialisieren des |a_i> für einen Eintrag von |Psi>
			numToVec(a, spins, spinAnzahl, h, spindimtable);
			HamiltonMalEintrag(a, kopplungen, psi, spins, gFaktoren, spinAnzahl, b_Feld, ergebnisshift, tochtershift, h, jMatrizenAnzahl, spindimtable, hdimserw, rvec, hilbertraumDimension);
		}
	}
}

// void erwartungswertSpin(const complex<double> *psi, const complex<double> *rvec, complex<double> *serw, const unsigned short *spins,
// 						const int spinAnzahl, const unsigned long hilbertraumDimension, const int k,
// 						const unsigned long *spindimtable, complex<double> *hdimserw)
// {
// 	const unsigned long tochtershift = (unsigned long) (k % 3) * hilbertraumDimension;// n ist sowohl für den die Auswahl des richtigen Psi als auch für die Speicherung der Werte für verschiedene l zuständig
// 	#pragma omp parallel
// 	{
// 		#pragma omp for 
// 		for (int i = 0; i < spinAnzahl; ++i)
// 		{
// 			for (unsigned long m = 0; m < hilbertraumDimension; ++m)
// 			{
// 				unsigned short a[spinAnzahl];
// 				double si;
// 				unsigned long m_strich;
// 				numToVec(a, spins, spinAnzahl, m, spindimtable);
// 				if (a[i] > 0)
// 				{
// 					si = s_plus(spins, a, i);
// 					a[i]--;
// 					m_strich = vecToNum(a, spins, spinAnzahl, spindimtable);
// 					a[i]++;
// 					serw[i*3 + k*3*spinAnzahl] += (complex<double>) conj(psi[tochtershift + m_strich]) * si * rvec[m];
// 				}
// 				if (a[i] < spins[i])
// 				{
// 					si = s_minus(spins, a, i);
// 					a[i]++;
// 					m_strich = vecToNum(a, spins, spinAnzahl, spindimtable);
// 					a[i]--;
// 					serw[1 + i*3 + k*3*spinAnzahl] += (complex<double>) conj(psi[tochtershift + m_strich]) * si * rvec[m];
// 				}
// 				si = s_z(spins, a, i);
// 				serw[2 + i*3 + k*3*spinAnzahl] += (complex<double>) conj(psi[tochtershift + m]) * si * rvec[m];
// 			}
// 		}
// 	}
// }



void lanczos_method(const Kopplungsmatrix *kopplungen, double *allEW, double *allsquarednr,
					const unsigned short *spins, const double *gFaktoren, const int spinAnzahl,
					const unsigned long hilbertraumDimension, const double *b_Feld,
					int lanczosSteps, complex<double> *psi, const int jMatrizenAnzahl, 
					complex<double> *SerwBrichtung, const int richtungen, const int r, const int b,
					const int random_vectors, const long samen, const bool zeitumkehrvec,
					const unsigned long *spindimtable)
{
	//Die Koeffizienten <m|r> werden gespeichert:
	boost::shared_array<complex<double>> rvec{new complex<double>[hilbertraumDimension]};
	memcpy(rvec.get(), psi, sizeof(complex<double>) * hilbertraumDimension);
	//Erstellen der Koeffizienten-Arrays alpha_k, beta_k
	complex<double> alpha[lanczosSteps];
	complex<double> beta[lanczosSteps - 1];
	//Ab hier beginnt der eigentliche Lanzcos-Algorithmus
	//Der Algorithmus ist nach Dissertation-Hanebaum S.17 1.3.2 Algorithmus 1.3.1 implementiert
	complex<double> serw[3 * lanczosSteps];//Für die inneren Summen der Berechnung
	memset(serw, 0, 3 * sizeof(complex<double>) * lanczosSteps);
	boost::shared_array<complex<double>> hdimserw{new complex<double>[3 * hilbertraumDimension]};
	for (int k = 0; k < lanczosSteps; k++)
	{	
		memset(hdimserw.get(), 0, 3 * sizeof(complex<double>) * hilbertraumDimension);
		//erwartungswertSpin(psi, rvec.get(), serw, spins, spinAnzahl, hilbertraumDimension, k, spindimtable);
		// index_enkel ist der Index des Vektors 2 vorher
		// index_tochter ist der Index des vorherigen Vektors
		// index_berechnung ist der Index des zu berechnenden Vektors
		unsigned long index_tochter = k % 3;
		unsigned long index_berechnung = (k + 1) % 3;
		unsigned long index_enkel = (k + 2) % 3;
		//Mit psi wird ein Zeiger auf (das erste Element von) |psi> 
		//und mit index_tochter der Vektor der als letztes berechnet wurde.
		hamiltonMalPsi(psi, kopplungen, spins, gFaktoren, spinAnzahl, hilbertraumDimension, b_Feld, k, jMatrizenAnzahl, spindimtable, hdimserw.get(), rvec.get());
		for (int j = 0; j < 3; ++j)
		{
			for (unsigned long i = 0; i < hilbertraumDimension; ++i)
			{
				serw[k + j*lanczosSteps] += hdimserw[i + j*hilbertraumDimension];
			}
		}
		alpha[k] = complexPsiProduct(psi, hilbertraumDimension, index_tochter, index_berechnung);
		//Beta wird für den 0-ten Durchlauf noch nicht benötigt (zum Orthonormalisieren)
		if (k != 0)
		{
			beta[k-1] = complexPsiProduct(psi, hilbertraumDimension, index_enkel, index_berechnung);
			//orthogonalisiere(psi, alpha[k], beta[k - 1], index_berechnung, hilbertraumDimension);
			#pragma omp parallel
			{
				#pragma omp for
				for (unsigned long i = 0; i < hilbertraumDimension; i++)
				{
					psi[i + (index_berechnung * hilbertraumDimension)] = psi[i + (index_berechnung * hilbertraumDimension)] - (psi[index_tochter*hilbertraumDimension + i] * alpha[k] + psi[index_enkel*hilbertraumDimension + i] * beta[k-1]);
				}
			}
		}
		else
		{ 
			#pragma omp parallel
			{
				#pragma omp for
				for (unsigned long i = 0; i < hilbertraumDimension; i++)
				{	
					psi[(index_berechnung * hilbertraumDimension) + i] -= psi[i] * alpha[k];
				}
			}
		}
		long double eta = normierung(psi, index_berechnung, hilbertraumDimension);
		//Abfragen der Toleranzgrenze
		eta = 1. / eta;
		#pragma omp parallel
		{
			#pragma omp for
			for (unsigned long i = 0; i < hilbertraumDimension; i++)//Normierung auf 1
			{
				psi[(index_berechnung * hilbertraumDimension) + i] *= eta;
			}
		}
	}
	boost::scoped_array<double> trimatrix{new double[lanczosSteps*lanczosSteps]};
	int triMaEntry = 0;
	//Einspeichern der Elemente in die Tridiagonalmatrix
	//Dabei wir nur die obere Nebendiagoale mit angegeben
	for (int i = 0; i < lanczosSteps - 1; i++)
	{
		trimatrix[triMaEntry] = alpha[i].real();
		triMaEntry++;
		trimatrix[triMaEntry] = beta[i].real();
		triMaEntry++;
		for (int p = 0; p < lanczosSteps - 2; p++)
		{
			trimatrix[triMaEntry] = 0;
			triMaEntry++;
		}
		trimatrix[triMaEntry] = beta[i].real();
		triMaEntry++;
	}
	trimatrix[triMaEntry] = alpha[lanczosSteps-1].real();
	int lwork = -1;
	double tempwork[1];
	int bugreport;
	char eigenvektorenIncluded = 'V';
	char obereDreiecksmatrix = 'L';	
	boost::scoped_array<double> eigenwerte{new double[lanczosSteps]};
	memset(eigenwerte.get(), 0, sizeof(double) * lanczosSteps);
	dsyev_(&eigenvektorenIncluded, &obereDreiecksmatrix, &lanczosSteps, trimatrix.get(), &lanczosSteps, eigenwerte.get(), tempwork, &lwork, &bugreport);
	lwork = tempwork[0];
	boost::scoped_array<double> work{new double[lwork]};
	memset(work.get(), 0, sizeof(double) * lwork);
	//Die LAPACK-Routine zum Berechnen von Eigenwerten (und Eigenvektoren) einer Tridiagonalmatrix
	dsyev_(&eigenvektorenIncluded, &obereDreiecksmatrix, &lanczosSteps, trimatrix.get(), &lanczosSteps, eigenwerte.get(), work.get(), &lwork, &bugreport);		
	//Speichert die Eigenwerte und Eigenvektoren in die Arrays ew und ev, wobei die Eigenvektoren
	//in column-major Order eingespeichert werden.
	for (int l = 0; l < lanczosSteps; l++)
	{
		allEW[l + r*lanczosSteps + b*random_vectors*lanczosSteps] = eigenwerte[l];
		allsquarednr[l + r*lanczosSteps + b*random_vectors*lanczosSteps] = trimatrix[l*lanczosSteps] * trimatrix[l*lanczosSteps];
	}
	//Hier wird die Berechnung von <n|S|r> vervollständigt...
	for (int n = 0; n < lanczosSteps; ++n)
	{
		for (int l = 0; l < lanczosSteps; ++l)
		{
			for (int j = 0; j < 3; ++j)
			{
				SerwBrichtung[n + r*lanczosSteps + b*random_vectors*lanczosSteps + j*random_vectors*lanczosSteps*richtungen] += serw[l + j*lanczosSteps] * trimatrix[n*lanczosSteps + l] * trimatrix[n*lanczosSteps];
			}
		}
	}
	eigenwerte.reset();
	work.reset();
	trimatrix.reset();
	rvec.reset();
}


long double zustandssumme(const double *allEW, const double *allsquarednr, const int lanczosSteps, const int random_vectors, const int b, const double temperatur) 
{	
	long double zsum = 0;
	for (int i = lanczosSteps-1; i >= 0; --i)
	{
		for (int r = 0; r < random_vectors; ++r)
		{
			zsum += (long double) allsquarednr[i + r*lanczosSteps + b*random_vectors*lanczosSteps] * exp((long double)(- allEW[i + r*lanczosSteps + b*random_vectors*lanczosSteps] / temperatur));
		}
	}
	return zsum;
}

long double magnetisierung(const double *allEW, const double *allsquarednr, const complex<double> *SerwBrichtung, const int lanczosSteps, const int random_vectors, const int richtungen,
						   const unsigned long hilbertraumDimension, const double gFaktor, const double temperatur, double *orientierungsmittel)
{
	complex<long double> Mag[3 * richtungen];
	double magproj[richtungen];
	memset(Mag, 0, 3 * sizeof(complex<long double>) * richtungen);
	memset(magproj, 0, sizeof(double) * richtungen);
	for (int b = 0; b < richtungen; ++b)
	{
		long double zsum = zustandssumme(allEW, allsquarednr, lanczosSteps, random_vectors, b, temperatur);
		zsum *= (long double) hilbertraumDimension/random_vectors;
		for (int i = lanczosSteps-1; i >= 0; --i)
		{
			for (int r = 0; r < random_vectors; ++r)
			{
				long double expo = - allEW[i + r*lanczosSteps + b*random_vectors*lanczosSteps] / temperatur;
				for (int j = 0; j < 3; ++j)
				{
					Mag[b + j*richtungen] += ((complex<long double>) SerwBrichtung[i + r*lanczosSteps + b*random_vectors*lanczosSteps + j*random_vectors*lanczosSteps*richtungen]) * (complex<long double>) exp(expo);
				}				
			}
		}
		for (int j = 0; j < 3; ++j)
		{
			Mag[b + j*richtungen] *= (long double) hilbertraumDimension / (random_vectors * richtungen);
			Mag[b + j*richtungen] *= (long double) -1.0 / zsum;
		}
		// cout << "z: " << Mag[b + 2*richtungen] << " +: " << Mag[b] << " -: " << Mag[b + richtungen] << endl;
		magproj[b] += Mag[b + 2*richtungen].real() * orientierungsmittel[b*3 + 2];
		magproj[b] += (0.5 * (Mag[b].real() + Mag[b + richtungen].real() )) * orientierungsmittel[b*3 + 0];
		magproj[b] += (0.5 * (Mag[b +richtungen].imag() - Mag[b].imag() )) * orientierungsmittel[b*3 + 1];

	}
	for (int b = 1; b < richtungen; ++b)
	{
		magproj[0] += magproj[b];
	}
	return magproj[0];
}






//Die main-Methode braucht zum Aufruf eine Textdatei mit den Spins etc. (Aufbau siehe unten)
//und die Anzahl der Lanczos Iterationen die maximal gemacht werden sollen
int main(int argc, char *argv[]) 
{
	omp_set_num_threads(threadzahl);
	if (argc!=3)//Zum Abfangen von Programmaufrufen mit einer unzulässigen Anzahl an Argumenten
	{
		cerr << "Der Programmaufruf benötigt nur die Anzahl der Lanzcos Steps und die Anzahl der zufälligen Vektoren,\n";
		 return 0;
	}

	// set_new_handler(outOfMemory);
	//Anzahl der Lanzcos Steps die gemacht werden sollen
	int lanczosSteps = atoi(argv[1]);
	//Anzahl der zufälligen Vektoren über die summiert werden soll
	int random_vectors = atoi(argv[2]);



	char inputname[60];
	sprintf(inputname, "input-example.txt");
	ifstream indat(inputname);
	if (!indat.is_open() )
	{
		cerr << "Keine Inputdatei geöffnet" << endl;
		return 0;
	}
	//Erstellen der ersten Variablen
	unsigned int spinJanzahlen[3];
	int einleseAbfang;
	einleseAbfang = dateizahlen(indat, spinJanzahlen);
	if(einleseAbfang == -1)
	{
		cerr << "Da ist was beim Einlesen schiefgelaufen!\n";
		return 0;
	}
	const unsigned int spinAnzahl = spinJanzahlen[0];
	const int jMatrizenAnzahl = spinJanzahlen[2];
	//Spin-Array erstellen
	unsigned short spins[spinAnzahl];
	speicherSpins(spins, spinAnzahl, indat);
	//gFaktorArray erstellen
	double gFaktoren[spinAnzahl];
	speicherGFaktoren(gFaktoren, spinAnzahl, indat);
	// J-Matrizen speichern
	double xyzkopp[jMatrizenAnzahl*11];
	Kopplungsmatrix kopplungen[jMatrizenAnzahl + 1];
	speicherJMatrizen(kopplungen, xyzkopp, jMatrizenAnzahl, indat);
	indat.close();
	
	unsigned long spindimtable[spinAnzahl];
	for (int i = 0; i < spinAnzahl; ++i)
	{
		spindimtable[i] = spindimension(i, spins, spinAnzahl);
	}

	//Abfragen der Hilbertraumdimension
	const unsigned long hilbertraumDimension = hilbertDimension(spins, spinAnzahl);


	const int richtungen = 10;

	//Zufallszahlen-Seed setzen
	long samen;
	time(&samen);
	double exec_time = 0, cpu_time;
 	// Start measuring time
    struct timespec begin, end, beg_cpu, end_cpu; 
    clock_gettime(CLOCK_REALTIME, &begin);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &beg_cpu);
	

	//Initialisiren der Vektoren |psi_1>, |psi_2>, |psi_3> in einem Array  
	boost::shared_array<complex<double>> psi{new complex<double>[3 * hilbertraumDimension]};
	//memset(psi.get(), 0, sizeof(complex<double>) * 3 * hilbertraumDimension);

	boost::shared_array<complex<double>> SerwBrichtung{new complex<double>[3 * lanczosSteps* richtungen * random_vectors]};
	boost::shared_array<double> allEW{new double[lanczosSteps* richtungen * random_vectors]};
	boost::shared_array<double> allsquarednr{new double[lanczosSteps* richtungen * random_vectors]};

	vector<double> temperatur;
	temperatur.push_back(2.0);
	temperatur.push_back(3.0);
	temperatur.push_back(5.0);

	//speichern des B-Feldes in kartesischen Koordinaten
	//https:www.etsy.com/de/market/dodecahedron_3d
	double orientierungsmittel[30] = {0.577350269189625764, 0.577350269189625764, 0.577350269189625764,
									 -0.577350269189625764, 0.577350269189625764, 0.577350269189625764,
									 0.577350269189625764, -0.577350269189625764, 0.577350269189625764,
									 0.577350269189625764, 0.577350269189625764, -0.577350269189625764,
									 0.0, 0.356822089773089932, 0.934172358962715698,
									 0.0, -0.356822089773089932, 0.934172358962715698,
									 0.934172358962715698, 0.0, 0.356822089773089932,
									 -0.934172358962715698, 0.0, 0.356822089773089932,
									 0.356822089773089932, 0.934172358962715698, 0.0,
									 -0.356822089773089932, 0.934172358962715698, 0.0};
	
	
	//Durchlaufen der B-Felder > 0 mit Mittelung über Kugeloberfläche nach Lebedev-Laikov (Dissertation Hanebaum)
	vector<double> vecbfeld;
	vecbfeld.push_back(0.1);
	vecbfeld.push_back(0.5);
	vecbfeld.push_back(0.8);
	vecbfeld.push_back(1.3);
	vecbfeld.push_back(1.9);
	vecbfeld.push_back(3.0);
	vecbfeld.push_back(5.0);
	vecbfeld.push_back(7.0);
	string magstring;
	char magtext[65];
	sprintf(magtext, "mag-");
	magstring = magtext;
	magstring += inputname;
	ofstream outmag(magstring);
	if (!outmag.is_open() )
	{
		cerr << "Keine Outputdatei geöffnet" << endl;
	} 
	outmag << "Mag-pulv Temperatur B-Feld" << endl;
	string ftlmstring;
	for (auto feldcounter : vecbfeld)
	{
		
		double feldstaerke = feldcounter;
		char text[65];
		sprintf(text, "eigenvals_b%2.2f-", feldstaerke);
		ftlmstring = text;
		ftlmstring += inputname;
		ofstream outdat(ftlmstring);
		if (!outdat.is_open() )
		{
			cerr << "Keine Outputdatei geöffnet" << endl;
		}
		//Kopfzeile in die Outputdatei schreiben
		outdat << "a " << lanczosSteps << " " << random_vectors << " " << hilbertraumDimension << endl;
		outdat << "f " << feldstaerke << endl;
		long double opsum[lanczosSteps];
		memset(SerwBrichtung.get(), 0, 3 * sizeof(complex<double>) * lanczosSteps * richtungen * random_vectors);
		//allEw und allsquarednr werden sowieso gesetzt und benötigen daher kein memset

		for (int b = 0; b < richtungen; ++b)//Diese Schleife über die Richtungsmittelung wird nur für Pulverproben benötigt
		{		
			//Es wird ein Feld benötigt, das für jeden Spin eine Unterschiedliche Richtung haben kann
			double b_Feld[3 * spinAnzahl];
			memset(b_Feld, 0, sizeof(double) * 3 * spinAnzahl);
			for (int i = 0; i < spinAnzahl; ++i)//Der Fall des homegenen Magnetfeldes
			{
				b_Feld[i*3] = feldstaerke * orientierungsmittel[b*3];
				b_Feld[i*3 + 1] = feldstaerke * orientierungsmittel[b*3 + 1];
				b_Feld[i*3 + 2] = feldstaerke * orientierungsmittel[b*3 + 2];
			}
			bool zeitumkehrvec = 0;
			double speicher_random = 0;
			
			for (int r = 0; r < random_vectors; ++r)
			{
				long double normierungZufall = zufallspsi(psi.get(), hilbertraumDimension, r, samen, zeitumkehrvec);
				lanczos_method(kopplungen, allEW.get(), allsquarednr.get(), spins, gFaktoren, spinAnzahl, hilbertraumDimension, b_Feld, lanczosSteps,
							   psi.get(), jMatrizenAnzahl, SerwBrichtung.get(), richtungen, r, b, random_vectors, samen, zeitumkehrvec, spindimtable);
				if (zeitumkehrvec == 1)
				{
					zeitumkehrvec = 0;
				}else {
					zeitumkehrvec = 1;
				}
				outdat << fixed << setprecision(3);
				outdat << "bidx" << endl;
				outdat << feldstaerke << "  " << b_Feld[0] << " " << b_Feld[1] << " " << b_Feld[2];
				outdat << fixed << setprecision(13);
				outdat << " e" << endl; 
				for (int j = 0; j < lanczosSteps; j++)
				{
					outdat << allEW[j + r*lanczosSteps + b*random_vectors*lanczosSteps] << " " << allsquarednr[j + r*lanczosSteps + b*random_vectors*lanczosSteps];
					outdat << " " << SerwBrichtung[j + r*lanczosSteps + b*random_vectors*lanczosSteps].real() << " " << SerwBrichtung[j + r*lanczosSteps + b*random_vectors*lanczosSteps].imag();
					outdat << " " << SerwBrichtung[j + r*lanczosSteps + b*random_vectors*lanczosSteps + 1*random_vectors*lanczosSteps*richtungen].real() << " " << SerwBrichtung[j + r*lanczosSteps + b*random_vectors*lanczosSteps + 1*random_vectors*lanczosSteps*richtungen].imag();
					outdat << " " << SerwBrichtung[j + r*lanczosSteps + b*random_vectors*lanczosSteps + 2*random_vectors*lanczosSteps*richtungen].real() << " " << SerwBrichtung[j + r*lanczosSteps + b*random_vectors*lanczosSteps + 2*random_vectors*lanczosSteps*richtungen].imag();
					outdat << endl;
				}
			}
			
		}
		for (auto t : temperatur)
		{	
			long double mag = magnetisierung(allEW.get(), allsquarednr.get(), SerwBrichtung.get(), lanczosSteps, random_vectors, richtungen, hilbertraumDimension, gFaktoren[0], t, orientierungsmittel);
			outmag << setprecision(11);
			outmag << mag << "  " << t << "  " << feldstaerke << endl;
		}
		outdat.close();
	}
	clock_gettime(CLOCK_REALTIME, &end);
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end_cpu);

	long seconds = end.tv_sec - begin.tv_sec;
	exec_time = (double) seconds / 60.0;	
	seconds = end_cpu.tv_sec - beg_cpu.tv_sec;
	cpu_time = (double) seconds / 60.0;	

	outmag << "z " << 'z' << endl;
	outmag.close();
	infoausgabe(spins, gFaktoren, random_vectors, lanczosSteps, hilbertraumDimension, spinAnzahl, jMatrizenAnzahl, xyzkopp, vecbfeld, temperatur, ftlmstring, magstring, inputname, samen, exec_time, cpu_time);
	psi.reset();
	allEW.reset();
	allsquarednr.reset();
	SerwBrichtung.reset();
	return 0;
}


/*Speichern der Dreiecksmatrix in einer Datei zum Überprüfen mit Mathematica

	ofstream trimaout("dreiecksdatei.txt");
	if (!trimaout.is_open())
	{
		cerr << "Keine Ausgabedatei geöffnet.!.!." << endl;
	}
	trimaout << setprecision(12);
	trimaout << "{";
	for (int i = 0; i < lanczosSteps; i++)
	{
		trimaout << "{";
		for (int j = 0; j < lanczosSteps - 1; j++)
		{
			trimaout << trimatrix[j + i * lanczosSteps] << ",";
		}
		trimaout << trimatrix[lanczosSteps - 1 + i * lanczosSteps];
		trimaout << "}," << endl;
	}
	trimaout << "}";
	trimaout.close();
*/



/*
Dies ist die main-Routine des FTLM-Algorithmus.
In diesem Programmteil wird der Hauptteil des Programms, 
die Wirkung von H auf |Psi>, gesteuert.
Außerdem werden die Werte aus der Inputdatei geholt 
und dann beginnt der eigentliche FTLM-Algorithmus:
Es wird für eine Anzahl zufälliger Vektoren (random_vectors)
der Laczos Algorithmus durchgeführt...
*/
/*
Der Hauptteil des Lanczos-Algorithmus berechnet H|psi> LanczosSteps mal
*/
/*
Die Textdatei in der die Daten gespeichert sind sollte folgendermaßen aufgebaut sein:
/l1:Vorangestelltes s und dann die Spins mit Lehrzeichen getrennt  (unsigned short)
/l2:Vorangestelltes g und dann die g-Faktoren mit Lehrzeichen getrennt  (double)
/l3:Vorangestelltes j und dann die Anzahl der noch folgenden J-Matrizen  (unsigned int)
Nun folgen die Blöcke für die J-Matrizen [1 Block = 4lines]
/l4: spini und spinj mit Lehrzeichen getrennt  (short)
/l5: Einträge xx,xy,xz mit Lehrzeichen getrennt  (double)
/l6: Einträge yx,yy,yz mit Lehrzeichen getrennt  (double)
/l7: Einträge zx,zy,zz mit Lehrzeichen getrennt  (double)
Die letzten 4 Lines wiederholen sich nun für alle i,j
*/

/* Comment1:
Dieser Teil (in der if-Bedingung) berechnet die Eigenwerte der Tridiagonalmatrix,
die auf der Hauptdiagonalen die alpha_k und auf den Nebendiagonalen die beta_k hat.
Damit ist diese Matrix reell symmetrisch und wird deswegen mit der dsyev Routine aus
Lapack berechnet: http://www.netlib.org/lapack/explore-html/d2/d8a/group__double_s_yeigen_ga442c43fca5493590f8f26cf42fed4044.html
Dazu wird zuerst immer die Größe des work-Arrays bestimmt um dann damit die EW zu berechnen.
Alle relevanten Arrays werden als scoped_arrays gespeichert um Speicherlecks zu vermeiden
*/

/*
Der Algorithmus im Überblick:
l 43-76: Das Speichern der J_ij und der new_handler
l 78-99: Hilfsfunktionen zum berechnen der Wirkung von s_+,-,z auf die Zustände
l 122-331: Die Wirkung des Hamiltonians auf die einzelnen Zustände.
l 230-346: Die Wirkung des Hamiltonians auf einen der Vektoren im Hilbertraum
l 348-490: Der Lanczos-Algorithmus
l 496-Ende: Der FTLM-Algorithmus in der main Methode
*/

//kompilieren mit:
//g++ -g -lm -lblas -llapack

//Aufruf über:
// ./<name> <Anzahl Lanczos Steps> <Anzahl zufälliger Vektoren> <Wert des größten B-Feldes> <Anzahl Schritte in denen das b_Feld durchlaufen wird>
// natürlich ohne die Dreiecksklammern.