#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;


inline unsigned long spindimension(const int j, const unsigned short *spins, const int spinAnzahl) {
	unsigned long actualdim = 1;

	for (int i = spinAnzahl - 1 ; i > j; i--) // Dimension bis zum i-ten Spin berechnen
	{
		actualdim *= (unsigned long) (spins[i] + 1); 
	}
	return actualdim;
}


//Funktion, die die Nummer des Basisvektors ausgibt
inline unsigned long vecToNum (const unsigned short *a, const unsigned short *spins,
							  const unsigned int spinAnzahl, const unsigned long *spindimtable)
{
	unsigned long vecnumber = a[spinAnzahl - 1];
	for (int j = spinAnzahl - 2; j >= 0; j--)
	{
		vecnumber += a[j] * spindimtable[j];
	}
	return vecnumber;
}

//Funktion die eine Nummer in einen Basisvektor |a> umwandelt
inline void numToVec(unsigned short *a, const unsigned short *spins, const unsigned int spinAnzahl,
					 const unsigned long inputNum, const unsigned long *spindimtable)
{
	unsigned long num = inputNum;
	for (int j = 0; j < spinAnzahl - 1; j++)
	{
		a[j] = (unsigned long) num / spindimtable[j];
		num -= a[j] * spindimtable[j];
	}
	a[spinAnzahl - 1] = num;
}
