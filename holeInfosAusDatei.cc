#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <locale>
#include <math.h>

using namespace std;
/*
Dies ist das Unterprogramm zum einlesen der Daten aus einer Datei.
In dieser wird die Ursprungsdatei gelesen und es werden die Spins, g-Faktoren,
Anzahl der J-Matrizen aber nicht die eigentlichen J-Matrizen (da diese das struct benötigen)
aus der Datei gelesen und in seperaten Arrays gespeichert.
*/
/*
Dabei sollte die Textdatei folgendermaßen aufgebaut sein:
/l1:Vorangestelltes s und dann die Spins mit Lehrzeichen getrennt  (unsigned short)
/l2:Vorangestelltes g und dann die g-Faktoren mit Lehrzeichen getrennt  (double)
/l3:Vorangestelltes j und dann die Anzahl der noch folgenden J-Matrizen  (unsigned int)
Nun folgen die Blöcke für die J-Matrizen [1 Block= 4lines]
/l4: spini und spinj mit Lehrzeichen getrennt  (short)
/l5: Einträge xx,xy,xz mit Lehrzeichen getrennt  (double)
/l6: Einträge yx,yy,yz mit Lehrzeichen getrennt  (double)
/l7: Einträge zx,zy,zz mit Lehrzeichen getrennt  (double)
Die letzten 4 Lines wiederholen sich nun für alle i,j
*/

/*
/Zählt die Spins und die g-Faktoren in der Inputdatei und vergleicht die Werte
/Außerdem wird die Anzahl der J-Matrizen gespeichert 
/und bei falschen Textdateien ein Fehler ausgegeben
*/
int dateizahlen(ifstream &ifs, unsigned int *spinsJausgabe)
{
	
	int sperre = 5000; //Die Anzahl der maximalen Einträge in der Inputdatei
	int counter[3] = {0,0,0};
	int index=0;
	bool jbool = false;
	char ctemp;
	double dtemp;
	ifs.seekg(ios::beg);
	while(!ifs.eof())
	{
		ifs >> ctemp;
		sperre--;
		if(ifs.bad() || sperre <=0)
		{
			cerr << "Im Input ist was schief gelaufen\n";
			return -1;
		}
		if (isdigit(ctemp))
		{
			ifs.unget();
			ifs >> dtemp;
			counter[index]++;

		}else if (ctemp == 's' || ctemp == 'S'){
			index = 0;
		}else if (ctemp == 'g' || ctemp == 'G')
		{
			index = 1;
		} else if (ctemp == 'j' || ctemp == 'J')
		{
			jbool = true;
			break;
		} else
		{
			break;
		}
	}
	ifs.clear();
	if (jbool)
	{
		ifs >> dtemp;
		counter[2] = (int) dtemp;
	} else
	{
		cerr << "In der Datei fehlt ein j für die Anzahl der Kopplungen\n";
		return -1;
	}
	if (counter[0] != counter[1])
	{
		cerr << "Die Anzahl der Spins stimmt nicht mit der Anzahl der g-Faktoren überein\n";
		return -1;
	}
	for (int i = 0; i < 3; i++)
	{
		spinsJausgabe[i] = counter[i];
	}
	return 0;
}


/*
/Speichert die Spins aus der Datei "spindatei"
/in dem Array spins der Länge N.
/Dabei werden die Spins schon ganzzahlig gemacht um short zu erhalten
*/
void speicherSpins(unsigned short *spins,const int spinAnzahl, ifstream &spindatei)
{
	spindatei.seekg(ios::beg);	
	char ctemp;
	spindatei >> ctemp;
	for (int i = 0; i < spinAnzahl; i++)
	{
		float hilfsVar;
		spindatei >> hilfsVar;
		hilfsVar *= 2.001;
		*(spins+i) = (unsigned short) hilfsVar;
		
	}

}
	
/*
/Speichert die g-Faktoren in einem Array
*/
void speicherGFaktoren(double *gFaktoren, const int spinAnzahl, ifstream &ifs)
{
	ifs.seekg(ios::beg);
	char ctemp;
	while(!ifs.eof())
	{
		ifs >> ctemp;
		if (ctemp == 'g' || ctemp == 'G')
		 {
		 	break;
		 }else if (ifs.bad())
		 {
		 	exit(0);
		 }
	}
	for (int i = 0; i < spinAnzahl; i++)
	{
		ifs >> *(gFaktoren+i);
	}
}

unsigned int hilbertDimension(const unsigned short *spinArray, const unsigned int spinAnzahl) 
{
	unsigned int actualdim = 1;

	for (unsigned int j = 0; j < spinAnzahl; j++) //Dimension des Hilbertraums bestimmen
	{
		actualdim *= *(spinArray+j) + 1; 
	}

	return actualdim;
}