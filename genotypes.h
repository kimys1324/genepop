/***************************************************************************
© F. Rousset 2005-2006

rousset@isem.univ-montp2.fr

This file is part of Genepop'007
This software is a computer program whose purpose is to perform statistical analyses.

This software is governed by the CeCILL license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.

 ***************************************************************************/
#ifndef H_GENOTYPES
#define H_GENOTYPES
#include <string>
#include <iterator> //for iterators in genetic programming (not much effect on problem?
#include <map>



// classe de création et de comptage de génotypes
template <typename Type>
class TGenotypes {
	private:
		map<Type, int>genotypes; // compteur de tous les génotypes existants
		int genotypeExists(Type genotype); // renvoie vrai si le génotype <genotype> existe
//something wrong with Type in this context; if this is really needed some, replace iterators (cf Primer+ p.894)		
		map<Type, int>::const_iterator iter; // itérateur sur les clés
		int _sum;
	public:
		int declareGenotype(Type genotype); // crée ou incrémente le génotype <genotype> ; renvoie le nb de génotypes
		void clear();
		int getEffective(Type genotype); // renvoie l'effectif du génotype <genotype>
		void fillGenotypes(int iLoc, CPopulation *pop); // remplit le truc avec les machins fournis à la construction
		void resetIterator(); // reset l'itérateur "iter" de "genotypes" à begin()
		int getNext(); // renvoie la clé suivant, -1 si end() atteint
		int getSum(); // renvoie le nombre total de génotypes
		int getNumber(); // nombre de génotypes différents
		void printKeys(int rank, string *output, int keyWidth, int foldLines); // imprime les génotypes, le rang indique si on doit imprimer l'allèle fort ou faible
		void printValues(TGenotypes<Type> *model, string *output, int effectiveWidth, int foldLines); // imprime les effectifs des génotypes
};

// ------------------------------------------------------------------------------------------
// TGenotypes -------------------------------------------------------------------------------

// renvoie vrai si le génotype <genotype> existe
template <typename Type>
int TGenotypes<Type>::genotypeExists(Type genotype) {
	map<Type, int>::iterator p = genotypes.find(genotype);
	return(p != genotypes.end() ? -1 : 0);
}

// renvoie l'effectif du génotype <genotype>
template <typename Type>
int TGenotypes<Type>::getEffective(Type genotype) {
	if (genotypeExists(genotype)) { return genotypes[genotype]; } else { return 0; }
}



// crée ou incrémente le génotype <genotype> ; renvoie le nb de génotypes
template <typename Type>
int TGenotypes<Type>::declareGenotype(Type genotype) {
	if (genotypeExists(genotype)) {
		genotypes[genotype]++;
	} else {
		genotypes[genotype] = 1;
	}
	_sum ++;
	return(genotypes.size());
}

// remplit les génotypes avec les informations de typages issues du fichier genepop
template <typename Type>
void TGenotypes<Type>::fillGenotypes(int geno, CPopulation *pop) {
	for (int iInd=0; iInd < pop->inds.size(); iInd ++) {
		if (pop->inds[iInd]->isValid(geno)) {
			declareGenotype(pop->inds[iInd]->getMaxAllele(geno) * 100 + pop->inds[iInd]->getMinAllele(geno));
		}
	}
}

// purge
template <typename Type>
void TGenotypes<Type>::clear() {
	genotypes.clear();
	resetIterator();
	_sum = 0;
}

// itérateur sur les clés
template <typename Type>
void TGenotypes<Type>::resetIterator() { iter = genotypes.begin(); }

template <typename Type>
int TGenotypes<Type>::getNext() {
	Type resu;
	if (iter == genotypes.end()) {
		resu = "";
	} else {
		resu = iter->first;
		iter ++;
	}
	return resu;
}

// renvoie le nombre total de génotypes (doit être = au nb d'individus si pas de trous 0000)
template <typename Type>
int TGenotypes<Type>::getSum() {return _sum;}

// renvoie le nombre de génotypes différents
template <typename Type>
int TGenotypes<Type>::getNumber() {return genotypes.size();}

// sortie d'une ligne (éventuellement repliée) au format "genotip2" des génotypes,
//     le rang définit lequel des deux allèles (faible ou fort) on écrit.
template <typename Type>
void TGenotypes<Type>::printKeys(int rank, string *output, int keyWidth, int foldLines) {
	Type genotype; // génotype courant
	int colCount; // comptage des colonnes en sortie
	char buf[10];
	char fmt[10];

	resetIterator();
	sprintf(fmt, "%%-%dd", keyWidth);
	colCount = output->length();
	while (!(strcmp(genotype = getNext(),""))) {
		if (foldLines && (colCount > 252)) {
			rtrim(output);
			*output += "\n";
			colCount = 0;
		} else {
			*output += " "; 
			colCount ++;
		}
		if (rank == 0) {
			sprintf(buf, fmt, genotype % 100);
		} else {
			sprintf(buf, fmt, int(genotype / 100));
		}
		*output += buf;
		colCount += 3;
	}
	rtrim(output);
}
// sortie d'une ligne (éventuellement repliée) des effectifs des génotypes, le modèle donne la liste des génotypes
//     dont on doit donner l'effectif
template <typename Type>
void TGenotypes<Type>::printValues(TGenotypes<Type> *model, string *output, int effectiveWidth, int foldLines) {
	Type genotype; // génotype courant
	int colCount; // comptage des colonnes en sortie
	char buf[10];
	
	colCount = output->length();
	model->resetIterator();
	while (!(strcmp(genotype = getNext(),""))) {   //probleme de comparaison selon le type...
		if ((colCount > 240) && foldLines) {
			rtrim(output);
			*output += "\n"; 
			colCount = 0;
		} else {
			*output += " ";
			colCount ++;
		}
		if (effectiveWidth == 2) {sprintf(buf, "%-2d", getEffective(genotype));
		} else {sprintf(buf, "%-3d", getEffective(genotype));}
		*output += buf;
		colCount += effectiveWidth;
	}
	rtrim(output);
}




#endif
