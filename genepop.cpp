/***************************************************************************
© R. Leblois 2001-2004
© F. Rousset 2005-2006

rousset@isem.univ-montp2.fr

This file is part of Genepop'007
This software is a computer program whose purpose is to perform statistical analyses.

This software is governed by the CeCILL license under French law and
abiding by the rules of distribution of free software.  You can  use,
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt".

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
#include <cmath> // rajout pour Dev-C++ ...
#include <iostream>
#include <iomanip>
#include <string>
#include <ctime>
#include <cstdlib>
#include <list>
#include <map>
#include <fstream>
#include <sstream>
#include <cctype>
#include <cstdio>
#include <algorithm> // for replace in string
#include "genepop.h"

//PAS d'EXTERN dans ce fichier !!

bool cinGetOnError=true; //to pause on cerr messages in cinGetOnError mode; overridden by explicit call to Batch mode
bool inputCheckBool=false; //debugging aid; there's a keyword too but coming a bit late in execution...
CFichier_genepop *fichier_genepop; // fichier genepop courant


/*
Le bordel:
analyses par individus utilise les effectifs determin['e]s par parsefile.
Les alleles pas dans des g['e]nos diploides di-typ['e]s sont perdus
tests HW utilisent la meme info => reconstitution des effectifs g['e]nos par fillGenotypes utilise
le nombre d'all[`e]les compt['e]s par parsefile;

*/
bool debug=false;

using namespace std;
string EOLtype="";
const string fichierIn="fichier.in";
const string strEMPTY = "";


//---------------------------Rafraichir l'['e]cran---------------------------------
//------------------------------------------------------------------------------
//using namespace NS_GP;
/*   #if (defined(linux) || defined(__CYGWIN__))
    system("clear");
   #elif (defined(WIN32))
    system("cls");
   #endif */



#ifdef WIN32
#include <windows.h>
void _gotoxy(int x,int y) {
    COORD mescoord;
    mescoord.X=x;
    mescoord.Y=y;
    SetConsoleCursorPosition(GetStdHandle(STD_OUTPUT_HANDLE),mescoord);
}/**/
void effacer_ecran() {
    system("cls");
}
// wherex/y() added 2015/02/04
int wherex()  {
  CONSOLE_SCREEN_BUFFER_INFO csbi;
  if (!GetConsoleScreenBufferInfo(
         GetStdHandle( STD_OUTPUT_HANDLE ),
         &csbi
         ))
    return -1;
  return csbi.dwCursorPosition.X;
}

int wherey()  {
  CONSOLE_SCREEN_BUFFER_INFO csbi;
  if (!GetConsoleScreenBufferInfo(
         GetStdHandle( STD_OUTPUT_HANDLE ),
         &csbi
         ))
    return -1;
  return csbi.dwCursorPosition.Y;
}
#else //vt100 terminals
void _gotoxy(int x,int y) {
//    printf("\033[%d;%dH",y+1,x+1);
      cout<<"\033["<<y+1<<";"<<x+1<<"H"<<flush; //equiv printf+flush ?
}
void effacer_ecran() {
//    printf("\033[2J");
      cout<<"\033[2J"<<flush; //equiv printf+flush ?
      _gotoxy(0,0);
}

/** from http://www.linuxquestions.org/questions/programming-9/get-cursor-position-in-c-947833/
with minor modifs
 Limited testing => appears to work (ubuntu, local) but perhaps not through network.
 Other possible sources: http://pipeviewer.googlecode.com/svn-history/r104/trunk/src/main/cursor.c
 or (perhaps better documented) http://unixwiz.net/techtips/wsize.c
*/
#include <fcntl.h>   // no newer c++ header
#include <termios.h> // idem
#include <cerrno>
#include <unistd.h>

#define   RD_EOF   -1
#define   RD_EIO   -2

static inline int rd(const int fd) {
    unsigned char   buffer[4];
    ssize_t         n;

    while (1) {
        n = read(fd, buffer, 1);
        if (n > (ssize_t)0) return buffer[0];
        else
        if (n == (ssize_t)0) return RD_EOF;
        else
        if (n != (ssize_t)-1) return RD_EIO;
        else
        if (errno != EINTR && errno != EAGAIN && errno != EWOULDBLOCK)
            return RD_EIO;
    }
}

static inline int wr(const int fd, const char *const data, const size_t bytes) {
    const char       *head = data;
    const char *const tail = data + bytes;
    ssize_t           n;

    while (head < tail) {
        n = write(fd, head, (size_t)(tail - head));
        if (n > (ssize_t)0) head += n;
        else
        if (n != (ssize_t)-1) return EIO;
        else
        if (errno != EINTR && errno != EAGAIN && errno != EWOULDBLOCK)
            return errno;
    }
    return 0;
}

/* Return a new file descriptor to the current TTY.
*/
int current_tty(void) {
    const char *dev;
    int         fd;
    dev = ttyname(STDIN_FILENO);
    if (!dev)
        dev = ttyname(STDOUT_FILENO);
    if (!dev)
        dev = ttyname(STDERR_FILENO);
    if (!dev) {
        errno = ENOTTY;
        return -1;
    }
    do {
        fd = open(dev, O_RDWR | O_NOCTTY);
    } while (fd == -1 && errno == EINTR);
    if (fd == -1) return -1;
    return fd;
}

/* As the tty for current cursor position.
 * This function returns 0 if success, errno code otherwise.
 * Actual errno will be unchanged.
*/
int cursor_position(const int tty, int *const rowptr, int *const colptr)
{
    struct termios  saved, temporary;
    int             retval, result, rows, cols, saved_errno;

    /* Bad tty? */
    if (tty == -1) return ENOTTY;

    saved_errno = errno;

    /* Save current terminal settings. */
    do {
        result = tcgetattr(tty, &saved);
    } while (result == -1 && errno == EINTR);
    if (result == -1) {
        retval = errno;
        errno = saved_errno;
        return retval;
    }
    /* Get current terminal settings for basis, too. */
    do {
        result = tcgetattr(tty, &temporary);
    } while (result == -1 && errno == EINTR);
    if (result == -1) {
        retval = errno;
        errno = saved_errno;
        return retval;
    }
    /* Disable ICANON, ECHO, and CREAD. */
    //temporary.c_lflag &= ~ICANON;
    //temporary.c_lflag &= ~ECHO;
    temporary.c_lflag &= ~(ICANON | ECHO);  /// FR alternative but does change the result...
    temporary.c_cflag &= ~CREAD;
    /* This loop is only executed once. When broken out,
     * the terminal settings will be restored, and the function
     * will return retval to caller. It's better than goto.
    */
    do {

        /* Set modified settings. */
        do {
            result = tcsetattr(tty, TCSANOW, &temporary);
        } while (result == -1 && errno == EINTR);
        if (result == -1) { // (FR) this happened over a network
            // (FR) then, reset to initial values is useful...
            tcsetattr(tty, TCSANOW, &saved);
            ///
            retval = errno;
            break;
        }

        /* Request cursor coordinates from the terminal. */
        retval = wr(tty, "\033[6n", 4);
        if (retval) break;
        /* Assume coordinate reponse parsing fails. */
        retval = EIO;
        /* Expect an ESC. */
        result = rd(tty);
        if (result != 27) break;
        /* Expect [ after the ESC. */
        result = rd(tty);
        if (result != '[') break;
        /* Parse rows. */
        rows = 0;
        result = rd(tty);
        /// il lit character par character et doit reconstruire l'entier...
        while (result >= '0' && result <= '9') {
            rows = 10 * rows + result - '0';
            result = rd(tty);
        }
        if (result != ';') break;
        /* Parse cols. */
        cols = 0;
        result = rd(tty);
        while (result >= '0' && result <= '9') {
            cols = 10 * cols + result - '0';
            result = rd(tty);
        }
        if (result != 'R') break;
        /* Success! */
        *rowptr = rows;
        *colptr = cols;
        retval = 0;
    } while (0);

    /* Restore saved terminal settings. */
    do {
        result = tcsetattr(tty, TCSANOW, &saved);
    } while (result == -1 && errno == EINTR);
    if (result == -1 && !retval)
        retval = errno;
    /* Done. */
    return retval;
}

int wherexy(int &row,int &col) {
    int         fd;
    char        buffer[64];
    char *const tail = buffer + sizeof(buffer);
    char       *head = buffer + sizeof(buffer);
    fd = current_tty();
    if (fd == -1) return 1;
    cursor_position(fd, &row, &col);
    return 0;
}

int wherex() {
    int row=0,col=0; /// default values if wherexy fails...
    wherexy(row,col);
    return(col);
}
int wherey() {
    int row=0,col=0;
    wherexy(row,col);
    return(row);
}
#endif


// utilitaire de comparaison de chaines insensible [`a] la casse
//attention sortie non intuitive! voir des usages prec['e]dents
int cmp_nocase(const string& s, const string& s2) {
	string::const_iterator p = s.begin();
	string::const_iterator p2 = s2.begin();

	while(p != s.end() && p2 != s2.end()) {
		if (toupper(*p) != toupper(*p2)) return((toupper(*p)<toupper(*p2)) ? -1 : 1);
		++p; ++p2;
	}
	return((s2.size()==s.size()) ? 0 : (s2.size()<s.size()) ? -1 : 1);
}
// vire les blancs [`a] droite
void rtrim(string *s) {
	while ((s->length()>0)  && (s->substr(s->length()-1,1)) == " ") {
			s->erase(s->length()-1,s->length());
	}
}

string ordonne(string geno) {// ordonne les genotypes diploides (STRING->STRING)
//int typ = atoi(geno.c_str());
string g1,g2;
	if (geno.length() == 4) { // cas [`a] 2*2 chiffres
	    g1 = geno.substr(0,2); // 2 chiffres de gauche
		g2 = geno.substr(2,2); // 2 chiffres de droite
	} else if (geno.length() == 6) {
		g1 = geno.substr(0,3); // 3 chiffres de gauche
		g2 = geno.substr(3,3); // 3 chiffres de droite
    } else {
        cerr <<"\nError: ordonne(string geno) called with geno= "<<geno<<endl;
        if (cinGetOnError) cin.get();
        exit(-1);
    }
    if (atoi(g2.c_str())>=atoi(g1.c_str())) return geno;
    else return (g2+g1);
}

int minAllele(int genotype,char coding) {
int g1,g2;
    if (coding<4) {cerr<<"useless call to CGenotypes::minAllele for haploid data";if (cinGetOnError) cin.get();exit(-1);}
    if (coding==4) {
        g1=genotype/100;
        g2=genotype%100;
        if (g1<g2) return g1; else return g2;
    }
    if (coding==6) {
        g1=genotype/1000;
        g2=genotype%1000;
        if (g1<g2) return g1; else return g2;
    }
	return(-1);
}

int maxAllele(int genotype,char coding) {
int g1,g2;
    if (coding<4) {cerr<<"useless call to CGenotypes::minAllele for haploid data";if (cinGetOnError) cin.get();exit(-1);}
    if (coding==4) {
        g1=genotype/100;
        g2=genotype%100;
        if (g1<g2) return g2; else return g1;
    }
    if (coding==6) {
        g1=genotype/1000;
        g2=genotype%1000;
        if (g1<g2) return g2; else return g1;
    }
	return(-1);
}




// ------------------------------------------------------------------------------------------
// -------------- impl['e]mentation ------------------------------------------------------------
// ------------------------------------------------------------------------------------------



// ------------------------------------------------------------------------------------------
// CAllele ----------------------------------------------------------------------------------

// constructeur
CAllele::CAllele(int identif, int effective) {
	_identif = identif;
	_effective = effective;
}
// incr['e]mente l'effectif
void CAllele::incEffective() {_effective ++;}

// renvoie l'effectif
int CAllele::getEffective() {return _effective;}


// ------------------------------------------------------------------------------------------
// CLocus -----------------------------------------------------------------------------------

// constructeur
CLocus::CLocus(int identifiant, string locName) {
	identif = identifiant;
	locusName = locName;
	alleleMax=0; //NOM d'allele max
	galleleMax=0; //NOM d'allele max
}


// ------------------------------------------------------------------------------------------
// CLocusGP ---------------------------------------------------------------------------------

// constructeur
CLocusGP::CLocusGP(int identifiant, string locName) : CLocus(identifiant, locName),alleles() {}

CLocusGP::~CLocusGP() {
//uses !=; not < for this iterator; also address the second element of the pair
    for (map<int,CAllele *>::iterator p=alleles.begin();p!=alleles.end();p++) {
        if (debug) cerr<<"destr called for CAllele* allele[...] in CLocusGP::~CLocusGP()\n";
        delete (*p).second;
    }
    for (map<int,CAllele *>::iterator p=galleles.begin();p!=galleles.end();p++) {
        if (debug) cerr<<"destr called for CAllele* gallele[...] in CLocusGP::~CLocusGP()\n";
        delete (*p).second;
    }
}


// teste l'existence d'un all[`e]le
int CLocusGP::alleleExists(int num) {
	map<int, CAllele *>::iterator p = alleles.find(num);
	return(p != alleles.end() ? -1 : 0);
}

// teste l'existence d'un all[`e]le
int CLocusGP::galleleExists(int num) {
	map<int, CAllele *>::iterator p = galleles.find(num);
	return(p != galleles.end() ? -1 : 0);
}

void CLocusGP::resetgIterator() { iter = galleles.begin(); }
void CLocusGP::resetIterator() { iter = alleles.begin(); }

int CLocusGP::getgNext() {
	int resu;
	if (iter == galleles.end()) {
		resu = -1;
	} else {
		resu = iter->first;
		iter ++;
	}
	return resu;
}

int CLocusGP::getNext() {
	int resu;
	if (iter == alleles.end()) {
		resu = -1;
	} else {
		resu = iter->first;
		iter ++;
	}
	return resu;
}

// renvoie le nombre d'all[`e]les
unsigned int CLocusGP::getNumber() {
	return alleles.size();
}


// cr['e]e ou incr['e]mente un all[`e]le
int CLocusGP::declareAllele(int num) {
	if (alleleExists(num)) {
	alleles[num]->incEffective();

	} else {
	alleles[num] = new CAllele(num,1);
	if(num>alleleMax) alleleMax=num; //NOM d'allele max
	}
	return(alleles.size());
}

// cr['e]e ou incr['e]mente un (g)all[`e]le
int CLocusGP::declaregAllele(int num) {
	if (galleleExists(num)) {
	   galleles[num]->incEffective();
	} else {
	   galleles[num] = new CAllele(num,1);
	   if(num>galleleMax) galleleMax=num; //NOM d'allele max
	}
	return(galleles.size());
}


// renvoie l'effectif de l'all[`e]le <num>
int CLocusGP::getEffective(int num) {
	if (alleleExists(num)) {
		return alleles[num]->getEffective();
	} else {
		return 0;
	}
}

int CLocusGP::getgEffective(int num) {
	if (alleleExists(num)) { // vs getggEffective below
		return galleles[num]->getEffective();
	} else {
		return 0;
	}
}

int CLocusGP::getggEffective(int num) {
	if (galleleExists(num)) {
		return galleles[num]->getEffective();
	} else {
		return 0;
	}
}

int CLocusGP::AlleleIdentif(int num) {
	return alleles[num]->_identif;
}



// ------------------------------------------------------------------------------------------
// CTypage ---------------------------------------------------------------------------------

// retourne le plus petit des 2 all[`e]les
int CTypage::minAllele(){return(gene1 > gene2 ? gene2 : gene1);}

// retourne le plus grand des 2 all[`e]les
int CTypage::maxAllele(){return(gene1 > gene2 ? gene1 : gene2);}


// ------------------------------------------------------------------------------------------
// CIndividual ------------------------------------------------------------------------------

// constructeur
CIndividual::CIndividual(string indName, int nbLoc){
	_indName = indName;
	_typages.reserve(nbLoc);
}

// destructeur
CIndividual::~CIndividual(){_typages.reserve(0);}

// ajoute un typage diploide
void CIndividual::addTypage(int allele1, int allele2) {
	_typages.resize(_typages.size() + 1);
	_typages[_typages.size()-1].gene1 = allele1;
	_typages[_typages.size()-1].gene2 = allele2;
	_typages[_typages.size()-1].valid = true;
}
// ajoute un typage diploide eventuellement incomplet aux *g*typages (genic)
void CIndividual::addgTypage(int allele1, int allele2) {
	_gtypages.resize(_gtypages.size() + 1);
	_gtypages[_gtypages.size()-1].gene1 = allele1;
	_gtypages[_gtypages.size()-1].gene2 = allele2;
	_gtypages[_gtypages.size()-1].valid = true;
}


// ajoute un typage haploide
void CIndividual::addTypage(int allele1) {
	_typages.resize(_typages.size() + 1);
	_typages[_typages.size()-1].gene1 = allele1;
	// attention .gene2 existe mais non valu['e]
	_typages[_typages.size()-1].valid = true;
}

// ajoute un typage vide
void CIndividual::addEmptyTypage() {
	_typages.resize(_typages.size() + 1);
	_typages[_typages.size()-1].valid = false;
}

int CIndividual::getTypage(int iLoc,int all) {
	if (all==1) return _typages[iLoc].gene1;
	else return _typages[iLoc].gene2;
}
;

// renvoie le plus petit des 2 all[`e]les du typage du locus <iLoc>
int CIndividual::getMinAllele(int iLoc){return(_typages[iLoc].minAllele());}

// renvoie le plus grand des 2 all[`e]les du typage du locus <iLoc>
int CIndividual::getMaxAllele(int iLoc){return(_typages[iLoc].maxAllele());}

// renvoie la validit['e] du typage au locus <iLoc>
bool CIndividual::isValid(int iLoc) {return(_typages[iLoc].valid);}

// renvoie le nom de l'individu
string CIndividual::indName() {return(_indName);}

// ------------------------------------------------------------------------------------------
// CPopulation ------------------------------------------------------------------------------

// constructeur, copie la liste des locus (sans les all[`e]les) depuis la liste globale
CPopulation::CPopulation(vector<CLocusGP *>lociSource)
: loci(lociSource.size()), inds() {
	vector<CLocusGP *>::iterator pSource;
	vector<CLocusGP *>::iterator p;
	p = loci.begin();
	for(pSource = lociSource.begin(); pSource != lociSource.end(); pSource ++, p++) {
		*p = new CLocusGP((*pSource)->identif, (*pSource)->locusName);
	}
}

// destructeur
CPopulation::~CPopulation() {
for (vector<CLocusGP *>::iterator p=loci.begin();p<loci.end();p++) {
        if (debug) cerr<<"destr called for CLocusGP * loci[...] in CPopulation::~CPopulation()\n";
        delete (*p);
    }
for (vector<CIndividual *>::iterator p=inds.begin();p<inds.end();p++) {
        if (debug) cerr<<"destr called for CIndividual * inds[...] in CPopulation::~CPopulation()\n";
        delete (*p);
    }
}

// alloue la place pour les genos d'un individu et retient son nom
// pour valeur des genos voir fillGenotypes ou semblable
int CPopulation::addIndividual(string indName) {
	CIndividual *ind = new CIndividual(indName, loci.size());
	// reserve ?
	if (inds.size() == inds.capacity()) {inds.reserve(inds.capacity() + 10);}
	inds.resize(inds.size()+1);
	inds[inds.size()-1] = ind;
	return(inds.size());
}

// renvoie le nom de la pop = nom du dernier individu (par convention genepop)
string CPopulation::popName() {
	return(inds.size() > 0 ? inds[inds.size() - 1]->indName() : strEMPTY);
}

// renvoie le dernier individu (= courant lors du parsing)
CIndividual *CPopulation::lastIndividual() {
	return(inds.size() > 0 ? inds[inds.size() - 1] : NULL);
}

// ------------------------------------------------------------------------------------------
// CFichier_genepop ----------------------------------------------------------------------

CFichier_genepop::CFichier_genepop(string fn) {
   fileName=fn;
}

int CFichier_genepop::checkName() { // should not be called from constructor because (1) int return value; (2) if return -1 exit in migraine, not exit in genepop
    string fileNameOri;
    fstream ORIfile(fileName.c_str(),ios::in|ios::out);
    if ( ! ORIfile.is_open()) {
	  fileNameOri=fileName;
	  fileName+=".txt";
	  ORIfile.clear();
	  ORIfile.open(fileName.c_str(),ios::in|ios::out);
	}
    if ( ! ORIfile.is_open()) {
        remove(fichierIn.c_str()); // otherwise, gets stuck if fichier.in contains an incorrect file name
        cerr<<"Cannot open "<<fileName.c_str()<<" or "<<fileNameOri.c_str()<<endl;
        cerr<<"Check input file name"<<endl;
        cin.ignore(); // vire le \n restant apr[`e]s une lecture dans cin
        return(-1); //
    }
    ///ELSE
    ORIfile.close();
    return(0); /// Normal exit
}


CFichier_genepop::~CFichier_genepop() {
for (vector<CLocusGP *>::iterator p=loci.begin();p<loci.end();p++) {
        if (debug) cerr<<"destr called for CLocusGP * loci[...] in CFichier_genepop::~CFichier_genepop()\n";
        delete (*p);
    }
for (vector<CPopulation *>::iterator p=pops.begin();p<pops.end();p++) {
        if (debug) cerr<<"destr called for CPopulation * pops[...] in CFichier_genepop::~CFichier_genepop()\n";
        delete (*p);
    }
}


// ajoute une population, transmet la liste de locus comme mod[`e]le
int CFichier_genepop::addPopulation() {
	CPopulation *pop = new CPopulation(loci);
	// reserve ?
	if (pops.size() == pops.capacity()) {pops.reserve(pops.capacity() + 10);}
	pops.resize(pops.size()+1);
	pops[pops.size()-1] = pop;
	return(pops.size());
}

// renvoie la derni[`e]re population (= courante lors du parsing)
CPopulation *CFichier_genepop::lastPop() {
	return(pops.size() > 0 ? pops[pops.size() - 1] : NULL);
}

// ajoute un locus [`a] la liste globale
int CFichier_genepop::addLocus(string locName) {
	CLocusGP *loc = new CLocusGP(loci.size(), locName);
	// reserve ?
	if (loci.size() == loci.capacity()) {loci.reserve(loci.capacity() + 10);}
	loci.resize(loci.size()+1);
	loci[loci.size()-1] = loc;
	return(loci.size());
}


int CFichier_genepop::parseFile() {
//cout<<"debut parseFile "<<EOLtype;getchar();
// A T T E N T I O N modifications paralleles dans read_bilocus
// lit le fichier / stocke des individus / ...
// POUR LES DONNEES HAPLOIDES
// la table genique est alleles
// POUR LES DONNEES DIPLOIDES:
//... constitue des tables geniques (alleles) avec les genos complets OU complets et incomplets (galleles)
// donc allele est util par CT pour les haploides, galleles pour les haploides
// les tables genotypiques doivent etre reconstitu['e]es (par ex [`a] partir de individus) dans les analyses qui les exploitent
// struc(), isolde, DL...
//***************
// dos/linux compatibility
// First implementation: (still there)
//   detecting EOL: \r\n under Windows, \n under linux. When a sample file written under Windows is read by Genepop
//   under Linux, the \r is present at the end of the string. Hence it must be read out before any buf.length()==0 test.
// BUT \r appears in output files as ^M and cause various troubles
// Second layer:
// Delete \r from input files. Found three methods
// (1) Implemented cat ... | tr -d '\r' > ... (tested)
// (2) vi "+set ff=unix" "+wq" (but buggy, not fixed)
// (3) GNU's recode command (not tested, not standard)
// (4) dos2unix fait la meme chose que (1) en plus simple d'['e]criture (impl['e]ment['e])
	int iLoc;
//	char toto[MAX_LINE_SIZE + 1];
//    string tempname;
//    char * tempnam;
	string buf,typst;
	stringstream stst(stringstream::in | stringstream::out);
	int nbLoc=0;
	int typ, g1, g2;
	string::size_type pos; // position de la virgule
	string::size_type lpos; // position des espaces entre genotypes
	//string::size_type possp; // position d'un eventuel espace
	string fileName2;
	ifstream inFile(fileName.c_str(), ios::in);
    if ( ! inFile.is_open()) {
        cerr<<"From parseFile(): Cannot open "<<fileName.c_str()<<" although it could previously be opened."<<endl;
        cerr<<"I exit."<<endl;
        if (cinGetOnError) cin.get();
        exit(-1); // vire le \n restant apr[`e]s une lecture dans cin
    }
    inFile.close();

    set_UNIX_EOLtype(fileName);

    inFile.open(fileName.c_str(), ios::in);
    getline(inFile,fileTitle);
	if (inputCheckBool) {cerr<<"|"<<fileTitle<<"|\n";}
  	if (inFile.eof()) {
  		cerr << "The file '" << fileName.c_str() << "' exists but does not even contain a title line." << endl;
  		if (cinGetOnError) cin.get();
  		exit(-1);
  	}
	do {    //recherche "pop" et lit les noms de locus au passage
//		inFile.get(toto, MAX_LINE_SIZE); buf = toto; inFile.ignore(1);
        	getline(inFile,buf);
	  	if (inFile.eof()) {
	  		cerr << "(!) From ParseFile(): The file '" << fileName << "' exists but no 'pop' keyword was found." << endl;
	  		if (cinGetOnError) cin.get();
	  		exit(-1);
	  	}
//		if (inputCheckBool) {cerr<<"|"<<buf<<"|\n";}
		while ((buf[0] == ' ') || (buf[0] == '\t') || (buf[0] == '\r')) {buf.erase(0, 1);} // vire les 'blancs' s.l. initiaux
#ifdef TESTBLANK
        if (buf.length()==0) goto next_line;
#else
		if (buf.length()==0) {
	  		cerr << "The file '" << fileName << "' contains a blank line. Check this file." << endl; //FER 08/2006
            if (cinGetOnError) cin.get();
	  		exit(-1);
	  	}
#endif
	    stst.clear();
        if (inputCheckBool) {cerr<<"(b)|"<<buf<<"|\n";}
        replace(buf.begin(), buf.end(), ',', '\b'); /// converts , to \b
        replace(buf.begin(), buf.end(), ' ', ','); /// nice hack converts all blanks to comma
        replace(buf.begin(), buf.end(), '\t', ','); /// likewise converts all tabs to comma
        replace(buf.begin(), buf.end(), '\b', ' '); /// back  converts \b to comma
		while (buf[buf.length()-1] == ',') {buf.erase(buf.length()-1, 1);} /// removes FINAL blanks and tabs
//if (inputCheckBool) {cerr<<"(b)|"<<buf<<"|\n";getchar();}
        if (cmp_nocase(buf, "POP") != 0) { // pas un "pop" (mais un nom de locus peut commencer par pop s'il n'y a pas un blanc au milieu....)
  	      stst.str(buf);
          while (!stst.eof()) {
            stst>>typst; // reads one string in the stst
///            if (inputCheckBool) {cerr<<"|"<<typst<<"|\n";}
            replace(typst.begin(), typst.end(), ',', ' '); /// whitespaces within names are put back
              while (typst[0] == ' ') {typst.erase(0, 1);} // vire les 'blancs'  initiaux
            if (inputCheckBool) {cerr<<"|"<<typst<<"|\n";}
            addLocus(typst);nbLoc++;
          }
        }
//next_line:;
	} while (cmp_nocase(buf, "POP") != 0);
// a ce point on a trouv['e] la premiere pop...
    coding.resize(nbLoc,0);
    if (inputCheckBool) {cerr<<nbLoc<<" loci declared; Press any key to resume."<<endl;getchar();}
	// pop
	addPopulation(); // nouvelle pop
	do {
//		inFile.get(toto, MAX_LINE_SIZE); buf = toto; inFile.ignore(1);
        getline(inFile,buf);
	    if (inputCheckBool) {cerr<<"|"<<buf<<"|\n";}
		while ((buf[0] == ' ') || (buf[0] == '\t')) {buf.erase(0, 1);} // vire les blancs initiaux
#ifdef TESTBLANK
        if (buf.length()==0) goto next_line;
#else
		if (buf.length()==0) {
	  		cerr << "The file '" << fileName << "' contains a blank line. Check this file." << endl; //FER 08/2006
            if (cinGetOnError) cin.get();
	  		exit(-1);
	  	}
#endif
		pos = buf.find(',');  // rajout FER pour verif qu'on n'est pas sur la ligne d'un genotype
//FR->FR ici il faudrait pouvoir virer des tabs precedents un nouveau POP, ou au moins les detecter, de la même facon que pour le premier pop...
		if (cmp_nocase(buf.substr(0,3), "POP") == 0 && pos == string::npos) { // d['e]but de nouvelle pop
			addPopulation();
		} else { // individu, pas pop
//			pos = buf.find(',');
//			if (pos == -1) {continue;} // ligne vide ou foireuse
			lastPop()->addIndividual(buf.substr(0,pos));
			buf.erase(0, pos+1);
//probleme pour genotypes sur plusieurs lignes:
//			for (iLoc=0; iLoc<nbLoc; iLoc ++) {
            iLoc=0;
			do { //while (iLoc<nbLoc) on va revenir la dessus
				while ((buf[0] == ' ') || (buf[0] == '\t') || (buf[0] == '\r')) {buf.erase(0, 1);} // vire les blancs/final \r
				while (buf.length()==0) { // si on est la c'est qu'on a pas tous les loci sur une ligne
//                    inFile.get(toto, MAX_LINE_SIZE); buf = toto; inFile.ignore(1);
                    getline(inFile,buf);
                    while ((buf[0] == ' ') || (buf[0] == '\t') || (buf[0] == '\r')) {buf.erase(0, 1);} // vire les blancs
                    // or /final \r since _next_ tested for buf.length()==0 [outer while()]
                } //[`a] ce point on doit avoir qq chose, soit un nouveau buf, soit la fin du buf prec
/*				lpos = min(buf.find(' '), min(buf.find('\t'), buf.length()));
				typ = atoi(buf.substr(0, lpos).c_str());
				buf.erase(0,lpos);*/
//cerr<<"\n"<<buf;
	    stst.clear();
	    stst.str(buf);
//cerr<<"\n"<< stst.str();
				// r['e]cup['e]ration des typages
		stst>>typst; // lit un mot dans la stst
		lpos=typst.length(); //en prend la lgr (typiquement 4)
		buf.erase(0,lpos); // le vire du buf
		typ = atoi(typst.c_str()); // converti en entier
//cerr<<endl<< stst.str()<<" > "<<typst<<" "<<lpos<<" @ "<<typ<<endl;getchar();
            if (coding[iLoc]==0) coding[iLoc]=lpos; //stocke le nombre de chiffre [`a] la premi[`e]re rencontre
            if (coding[iLoc]!=lpos) {//RL -> FR probleme comparaison unsigned / signed
               cout<<"(!) Problem in population "<<pops.size()<<", individual "<<pops.back()->inds.size()<<":"<<endl;
               cout<<"    inconsistent number of digits in genotypes of locus "<<iLoc+1<<" ("<<int(coding[iLoc])<<" expected, "<<lpos<<" found)"<<endl;
               cout<<"I must exit.";if (cinGetOnError) cin.get();exit(-1);
            }
				if (lpos == 4) { // cas [`a] 2*2 chiffres
					g1 = typ / 100; // 2 chiffres de gauche
					g2 = typ % 100; // 2 chiffres de droite
				} else if (lpos == 6) {
					g1 = typ / 1000; // 3 chiffres de gauche
					g2 = typ % 1000; // 3 chiffres de droite
				} else if (lpos==2 || lpos==3) { //haploid data, do nothing
                } else { // erreur
//cerr<<lpos;
					cerr<<"(!) Incorrect input '"<<typ <<"' in input file" << endl;
					cerr<<"    Check population "<<pops.size() << ", individual "<<lastPop()->inds.size()<< ", locus " << iLoc;
					if (cinGetOnError) cin.get();
					exit(-1);
				}
				// ATTENTION : ceci implique que pour etre valide un typage doit avoir les 2 g[`e]nes non nuls
// il faudra rajouter un cas et cr['e]er des variables pour des types alleliques (_typage etant un geno duplo)
                if  (lpos<4) {
                    if (typ>0) { //typage haploide valide
    					// ['e]criture du typage
    					lastPop()->lastIndividual()->addTypage(typ);
    					// ajout des all[`e]les pour la population
    					lastPop()->loci[iLoc]->declareAllele(typ);
    					// ajout des all[`e]les pour le fichier global
    					loci[iLoc]->declareAllele(typ);
                    } else { // typage haploide nul
					     lastPop()->lastIndividual()->addEmptyTypage();
				    }
                }
				else // typage diploide
                if ((g1 > 0) || (g2 > 0)) {
//if (iLoc==1) {cerr<<g1<<" "<<g2;getchar();}
					if ((g1 > 0) && (g2 > 0)) { // typage diploide complet
    					// ['e]criture du typage
    					lastPop()->lastIndividual()->addTypage(g1, g2);
    					// ajout des all[`e]les pour la population
    					lastPop()->loci[iLoc]->declareAllele(g1);
    					lastPop()->loci[iLoc]->declareAllele(g2);
    					// ajout des all[`e]les pour le fichier global
    					loci[iLoc]->declareAllele(g1);
    					loci[iLoc]->declareAllele(g2);
    				} else { // either g1 XOR g2 is >0
    				    lastPop()->lastIndividual()->addEmptyTypage(); //OUPS! 22/10/2006
    				    // more below
    				}
//cout<<"bof";getchar();
					// type complet ou incomplet mais non nul: declaration des alleles
					lastPop()->lastIndividual()->addgTypage(g1, g2);
					// ajout des all[`e]les pour la population
					if (g1>0) lastPop()->loci[iLoc]->declaregAllele(g1);
					if (g2>0) lastPop()->loci[iLoc]->declaregAllele(g2);
					// ajout des all[`e]les pour le fichier global
					if (g1>0) loci[iLoc]->declaregAllele(g1);
					if (g2>0) loci[iLoc]->declaregAllele(g2);
//cout<<"meuh"<<lastPop()->loci[iLoc]->getgEffective(g1);getchar();
				} else { // typage diploide nul
//if (iLoc==1) {cerr<<g1<<" nul "<<g2;getchar();}
					lastPop()->lastIndividual()->addEmptyTypage();
				}
//if (lastPop()->lastIndividual()->isValid(iLoc)) cerr<<"true"; else cerr<<"false";  getchar();
				iLoc ++;
			} while (iLoc<nbLoc);
		//check whether there is still something on the line
            while ((buf[0] == ' ') || (buf[0] == '\t') || (buf[0] == '\r')) {buf.erase(0, 1);} // vire les blancs
            if ( buf.length()>0 ) {
				cerr<<"(!) Invalid information for some individual: some information was found" << endl;
				cerr<<"   after all genotypes have been read for an individual." << endl;
				cerr<<"   Check individual "<<lastPop()->inds.size()<<" in population "<<pops.size();
				cerr<<"   Suspect information: '"<<buf<<"'";
				if (cinGetOnError) cin.get();
				exit(-1);
            }
//cout<<endl;
//nextline:;
		}
	} while (! inFile.eof());
  	inFile.close();
//	if (lpos==3 || lpos==6) digit3indic=false; else digit3indic=true;
	return(0);	 // TODO: pr['e]voir l'['e]chec de lecture
}

// derived from https://ccrma.stanford.edu/~craig/utility/flip/
void translateToUnix(const char* filename) {
#ifdef WIN32
 fstream infile(filename, ios::in | ios::binary);
#else
 fstream infile(filename, ios::in);
#endif
   if (!infile) {
      cout << "Error: cannot find file: " << filename << endl;
      return;
   }

   stringstream outstring;
   char ch, lastch;
   infile.get(ch);
   while (!infile.eof()) {
      if (ch == 0x0d) {
         outstring << (char)0x0a;
      } else if (ch == 0x0a) {
         if (lastch == 0x0d) {
            // do nothing: already converted MSDOS newline to Unix form
         } else {
            outstring << (char)0x0a;   // convert newline from Unix to Unix
         }
      } else {
         outstring << ch;
      }
      lastch = ch;
      infile.get(ch);
   }

   infile.close();
#ifdef MSWIN32
  fstream outfile(filename, ios::out | ios::binary);
#else
  fstream outfile(filename, ios::out);
#endif
   if (!outfile.is_open()) {
      cout << "Error: cannot write to file: " << filename << endl;
      return;
   }
   outstring << ends;
   outfile << outstring.str().c_str();
   outfile.close();
}

int set_eof_check_EOLtype(string EOLFileName, bool set_eof){
	long pos,lngth;
	char c, buf;
    EOLtype="";
    fstream ORIfile(EOLFileName.c_str(),ios::in|ios::out);
    if (!ORIfile.is_open()) {
        remove(fichierIn.c_str()); // otherwise, gets stuck if fichier.in contains an incorrect file name
        cerr << "From set_eof_check_EOLtype(): Cannot open file!" << endl;
		cerr << "Check input file name " << EOLFileName.c_str() << endl;
		cin.ignore(); // vire le \n restant apr[`e]s une lecture dans cin
//	cerr<<"(Return) to show working directory and list of `*.' files:"<<endl;if (cinGetOnError) cin.get();
#ifdef WIN32
//	system("dir *. /p");
#else
//	system("ls *. | more");
#endif
        return(-1);
    }
    while ( ! ORIfile.eof() ) {
         c=ORIfile.get();
         if (c=='\r' || c=='\n') break;
    }
    if (ORIfile.eof()) { //there was no \r nor \n
        cerr<<"No line terminator in the file!";
        cerr<<"I exit."<<endl;
        if (cinGetOnError) cin.get();
        exit(-1);
    } else if (c=='\n') {
        EOLtype="LinuxLF";
    } else /* we found \r first*/ if (ORIfile.get()=='\n'){
        EOLtype="WindowsCRLF"; //Windows \r\n
    } else { //\r not followed by \n
        cerr<<"(!) The file appears to contain a CR line terminator."<<endl;
        cerr<<"    This was the old MacIntosh end-of-line character."<<endl;
        cerr<<"    (i) If you are using Mac OS X, this should no longer occur (line terminator should be LF as in Unix)."<<endl;
        cerr<<"        However, some Microsoft editors for Mac OS X still use CR (sigh)."<<endl;
        cerr<<"        Genepop does not handle CR line terminators. Use a standard-compliant text editor."<<endl;
        cerr<<"    (ii) If you are using Windows/Linux, convert your input file to the correct format for your OS."<<endl;
        cerr<<"I exit."<<endl;
        if (cinGetOnError) cin.get();
        exit(-1);
    }

    if( ! set_eof )
    	return 0;

    /** ELSE **/

    ORIfile.clear();
    ORIfile.seekg(0,ios::end);
    lngth=ORIfile.tellg(); //longueur du fichier originel
// manitenant va lire charcater par char [`a] partir de la fin jusqu'[`a] trouver un chiffre
    pos=-1;
    ORIfile.seekg(pos,ios::end);
    ORIfile.read((char *) &buf, sizeof buf);
    while ( ! ( ((buf>47) && (buf<58)) || pos<= - lngth) ) { // stops if beginning of file reached [pos<= - lngth] or if  (numberfound && EOLtype.size()>0)
        pos--;
        ORIfile.seekg(pos,ios::end);
        ORIfile.read((char *) &buf, sizeof buf);
    }
    if (pos==-lngth) {
        cerr<<"No number, hence no genotype, in the file!";
        cerr<<"Exiting";
        if (cinGetOnError) cin.get();
        exit(-1);
    }
    // now we replace anything between the final number and the EOF by space chars
    ORIfile.seekp(lngth+pos+1);
    while (ORIfile.tellp()<lngth) ORIfile.put(' ');
    ORIfile.close();
return 1;
}

void set_UNIX_EOLtype(string EOLFileName) {
#ifdef WIN32
   // Unix-alike -> Windows
   // do nothing, because file 'as is' works anyway
#else
    /* EOLtype has been previously determined by
    set_eof_check_EOLtype(gp_file,true);
    previous to calling this function through ::parseFile() -> set_UNIX_EOLtype(fileName);
    */
    if (EOLtype=="WindowsCRLF") { // determined by ::set_eof_check_EOLtype !
       // Windows -> Unix
       translateToUnix(EOLFileName.c_str());
if (false) {
        cout<<endl<<"Input file seems to follow Windows end-of-line format. Attempting conversion to Unix format."<<endl;
        cerr<<"(if this fails, find some way to convert the file before running Genepop)\n";
  	    if (cinGetOnError) cin.get();
        // plus simple mais moins portable ?
        {
  	        stringstream cmdline;
            cmdline<<"dos2unix -k -q "<<EOLFileName.c_str()<<"\n";
    		int abyss=system(cmdline.str().c_str());
        }
		// don't use system() reurn call... may bear no relation to the cmdline
		// the return code may come from somthing else than the cmdline, and perror is equally useless
		 set_eof_check_EOLtype(EOLFileName, false); // ONLY recheck EOL
		 if (EOLtype=="WindowsCRLF") {
			 cerr<<"Problem calling 'dos2unix' on this OS. Trying 'fromdos' for conversion...\n";
			 // ubuntu's
			{
				 stringstream cmdline;
				 cmdline<<"fromdos "<<EOLFileName.c_str()<<"\n";
				 int abyss=system(cmdline.str().c_str());
			}
			 set_eof_check_EOLtype(EOLFileName, false); // ONLY recheck EOL
			 if (EOLtype=="WindowsCRLF") {
				 int bidon=0;
				 cerr<<"Problem calling 'fromdos' on this OS";
				 cerr<<"Trying 'tr' for conversion...\n";
				 char masque[] = "fileXXXXXX";
				 bidon=mkstemp(masque);
				 if (bidon==-1) {
					 cerr<<"(!) From parseFile(): Problem calling mkstemp() in this directory. Exiting...";
					 if (cinGetOnError) cin.get();
					 exit(-1);
				 }
				 string tempname=masque;
				{
					 stringstream cmdline;
					 cmdline<<"tr -d '\\r' <"<<EOLFileName.c_str()<<" > "<<tempname.c_str()<<endl;
					 int abyss=system(cmdline.str().c_str());
				}
				 set_eof_check_EOLtype(EOLFileName, false); // ONLY recheck EOL
				 if (EOLtype=="WindowsCRLF") {
					 cerr<<"(!) From parseFile(): Problem calling 'tr' on this OS. Exiting...";
					 if (cinGetOnError) cin.get();
					 exit(-1);
				 } else {
					 remove(EOLFileName.c_str());
					 rename(tempname.c_str(),EOLFileName.c_str());
				 }
			 }
			 else {
				 cout << "\nConversion from Windows EOL to Unix EOL successful! Press any key to continue...";
				 if (cinGetOnError)
					 cin.get();
			 }
		 }
		 else {
			 cout << "\nConversion from Windows EOL to Unix EOL successful! Press any key to continue...";
			 if (cinGetOnError)
				 cin.get();
		 }
}
    }
#endif
}





//======================================================================
//                      Lecture pour DL
//======================================================================

// it['e]rateur sur les cl['e]s
void CGenobilocus::resetIterator() { iter = mapmap.begin(); }

int CGenobilocus::getNext() {
	int resu;
	if (iter == mapmap.end()) {
		resu = -1;
	} else {
		resu = iter->first;
		iter ++;
	}
	return resu;
}

CGenobilocus CFichier_genepop::read_bilocus(int pop,int loc1,int loc2) { //loc2>loc1!
/*** haploid locus ->haploid info; diploid locus ->complete genotypic info retained only (10/12/2011)***/
// pour la pair loc1,loc2   attention loc1<loc2 => le loc1 est le premier lu sur la ligne
//    if (loc2<loc1) {cout<<"(!) From CGenobilocus CFichier_genepop::read_bilocus: loc2<loc1";getchar();}
//cout<<"debut CFichier_genepop::read_bilocus"<<flush;getchar();
	int nbLocMinusOne=loci.size()-1;
	stringstream stst(stringstream::in | stringstream::out);
    string buf,typst;
	int popit=0; // numero de la derni[`e]re pops rencontr['e]e
	int iLoc=0; // nb de loci deja rencontr['e]es
	size_t pos; // position de la virgule
	int lpos; // position des espaces entre genotypes
	//int possp; // position d'un eventuel espace
	int typ1=0,typ/*,bidonit*/,g1,g2;
	CGenobilocus counts;
	string fileName2;
	ifstream inFile(fileName.c_str(), ios::in);
    if ( ! inFile.is_open()) {
        cerr<<"From read_bilocus(): Cannot open "<<fileName.c_str()<<" although it could previously be opened."<<endl;
        cerr<<"I exit."<<endl;
        if (cinGetOnError) cin.get();
        exit(-1); // vire le \n restant apr[`e]s une lecture dans cin
    }
//	inFile.get(toto, MAX_LINE_SIZE); fileTitle = toto; inFile.ignore(1);
    getline(inFile,buf);
	while ((buf[0] == ' ') || (buf[0] == '\t') || (buf[0] == '\r')) buf.erase(0, 1); // vire les blancs/final \r
  	if (inFile.eof()) {
  		cerr << "(!) From read_bilocus(): The file '" << fileName << "' exists but does not contain data" << endl;
  		if (cinGetOnError) cin.get();
  		exit(-1);
  	}
	do {
//		inFile.get(toto, MAX_LINE_SIZE); buf = toto; inFile.ignore(1);
        getline(inFile,buf);
		while ((buf[0] == ' ') || (buf[0] == '\t') || (buf[0] == '\r')) {buf.erase(0, 1);} // vire les blancs/final \r
	  	if (inFile.eof()) {
	  		cerr << "(!) From read_bilocus(): The file '" << fileName << "' exists but but no 'pop' was found" << endl;
	  		if (cinGetOnError) cin.get();
	  		exit(-1);
	  	}
	    stst.clear();
	    stst.str(buf);
		stst>>typst; // lit un mot dans la stst
// the nice thing with the last line of code is that final tabulation chars are not pushed into typst
// so whatever text editor mess may be after "pop" is ignored, There is no need to check the length of the
// string beginning with "pop" and to write special code for spaces, tabs...
// in order that locus name of form popXX is recognized as not "pop"
// while pop\t is recognized as "pop"
	} while (cmp_nocase(typst, "POP") != 0);
// a ce point on a trouv['e] la premiere pop...
	// pop
//cout<<pos<<" "<<popit<<"! "<<buf<<flush;getchar();
	while ((!inFile.eof()) && !(popit==pop)) {
//		inFile.get(toto, MAX_LINE_SIZE); buf = toto; inFile.ignore(1);
        getline(inFile,buf);
		while ((buf[0] == ' ') || (buf[0] == '\t')) {buf.erase(0, 1);} // vire les blancs
		pos = buf.find(',');  // pour verif qu'on n'est pas sur la ligne d'un genotype
		if (cmp_nocase(buf.substr(0,3), "POP") == 0 && pos == string::npos) popit++;
//cout<<pos<<" "<<popit<<"@ "<<buf<<flush;getchar();
	}  //la on est arriv['e] sur la bonne pop
	do {
//		inFile.get(toto, MAX_LINE_SIZE); buf = toto; inFile.ignore(1);
        getline(inFile,buf);
		while ((buf[0] == ' ') || (buf[0] == '\t')) {buf.erase(0, 1);} // vire les blancs
		pos = buf.find(',');  // pour verif qu'on n'est pas sur la ligne d'un genotype
		if (cmp_nocase(buf.substr(0,3), "POP") == 0 && pos == string::npos) break; //on a rencontr['e] la pop suivante
        else { // individu, pas pop
			buf.erase(0, pos+1);
            iLoc=-1;
			do {
				while ((buf[0] == ' ') || (buf[0] == '\t') || (buf[0] == '\r')) {buf.erase(0, 1);} // vire les blancs/final \r
				while (buf.length()==0) { // si on est la c'est qu'on a pas tous les loci sur une ligne
                    getline(inFile,buf);
                    while ((buf[0] == ' ') || (buf[0] == '\t') || (buf[0] == '\r')) {buf.erase(0, 1);} // vire les blancs/final \r
                } //[`a] ce point on doit avoir qq chose
				lpos = min(buf.find(' '), min(buf.find('\t'), buf.length()));
				if (lpos>3) typ = atoi(ordonne(buf.substr(0, lpos)).c_str()); else typ=atoi(buf.substr(0, lpos).c_str());
				buf.erase(0,lpos);
				iLoc++;
				// r['e]cup['e]ration des typages
                if (iLoc==loc1) { //loc2>loc1!
                   if (counts.coding1==0) counts.coding1=lpos; //stocke le nombre de chiffre [`a] la premi[`e]re rencontre
                   if (counts.coding1!=lpos) {
                       cout<<"Inconsistent number of digits in genotypes of locus "<<loc1+1<<endl;
                       cout<<"I must exit.";exit(-1);
                   }
                   //tri geno ab, ba
                    // if lpos<4, haploide, 2 ou 3 chiffres: leave typ as is
                    if (lpos == 4) { // cas [`a] 2*2 chiffres   // fait la meme chose que fillGenotypes...
    					g1 = typ / 100; // 2 chiffres de gauche
    					g2 = typ % 100; // 2 chiffres de droite
                        if (min(g1,g2)==0) {
                            typ=0;
                        } else if (g1>g2) typ=100*g2+g1; //0201 ->0102
    				} else if (lpos == 6) {
    					g1 = typ / 1000; // 3 chiffres de gauche
    					g2 = typ % 1000; // 3 chiffres de droite
    					if (min(g1,g2)==0) {
                            typ=0;
                        } else if (g1>g2) typ=1000*g2+g1;
    				}
                   typ1=typ; //enregistre le typage au 1e locus
//cout<<typ1<<" / "<<typ<<flush;getchar();
//               	   testloc = counts.find(typ1);
//	               if (testloc==counts.end()) counts.insert(typ1);
                }
//cout<<"a";getchar();
                if (iLoc==loc2) { //loc2>loc1!
//cout<<typ1<<" // "<<typ<<flush;getchar();
//cerr<<counts.coding2<<endl;
                   if (counts.coding2==0) counts.coding2=lpos; //stocke le nombre de chiffre [`a] la premi[`e]re rencontre
                   if (counts.coding2!=lpos) {
                       cerr<<"Inconsistent number of digits in genotypes of locus "<<loc2+1<<endl;
                       cerr<<"I must exit.";if (cinGetOnError) cin.get();exit(-1);
                   }
                   if ((typ1!=0) && (typ!=0)) {
                    if (lpos == 4) { // cas [`a] 2*2 chiffres
    					g1 = typ / 100; // 2 chiffres de gauche
    					g2 = typ % 100; // 2 chiffres de droite
                        if (min(g1,g2)==0) {
                            typ=0;
                        } else if (g1>g2) typ=100*g2+g1; //0201 ->0102
    				} else if (lpos == 6) {
    					g1 = typ / 1000; // 3 chiffres de gauche
    					g2 = typ % 1000; // 3 chiffres de droite
                        if (min(g1,g2)==0) {
                            typ=0;
                        } else if (g1>g2) typ=1000*g2+g1;
    				}
//cout<<typ1<<" // "<<typ<<flush;getchar();
                        if (min(typ1,typ)>0) { // (genic): info at both loci; (genotypic) complete info at both loci
                           counts.marginal.declareGenotype(typ); //une facon de compter le nombre de colonnes totales du tableau
                           counts.mapmap[typ1].declareGenotype(typ); //ajoute 1 dans la case du geno bilocus (typ1,typ)
                        } //else ignoreincompleteinfo
                   }
//cout<<typ1<<" /// "<<typ<<counts.mapmap[typ1].getEffective(typ)<<flush;getchar();
                } //iLcc==loc2
//cout<<iLoc;getchar();
            //ici while (iLoc<loc2); planterait si genotypes sur 2 lignes et loc2 sur la premiere ligne
            // car iLoc serait r['e]initialis['e] [`a] la lecture de la 2e ligne
			} while (iLoc<nbLocMinusOne);
        }	// fin analyse individu
	} while (!inFile.eof()); //on a pu sortir au break....
  	inFile.close();
	return(counts);	 // TODO: pr['e]voir l'['e]chec de lecture
}

int CGenobilocus::getMinDim() {
    if (mapmap.size()<marginal.getNumber()) return mapmap.size(); else return marginal.getNumber();// RL->FR comparaison signed / unsigned
}

int CGenobilocus::getligNbr() {
    return mapmap.size();
}

int CGenobilocus::getcolNbr() {
    return marginal.getNumber();
}



vector<vector<int> >CGenobilocus::tabule() {
vector<vector<int> > sortie;
int ii,jj,typ/*,typ1*/;
  	sortie.resize(mapmap.size());
  	ii=0;
  	for (map<int,CGenotypes>::iterator iit=mapmap.begin();iit!=mapmap.end();iit++) {
        sortie[ii].resize(marginal.getNumber(),0);
        marginal.resetIterator();
        jj=0;
        while ((typ=marginal.getNext()) >= 0) { //puis parcoure la 2e dimension
//cout<<typ<<" "<<marginal.getEffective(typ);getchar();
             sortie[ii][jj]=iit->second.getEffective(typ);
//cout<<ii<<" "<<jj<<" "<<sortie[ii][jj];getchar();
             jj++;
        }
        ii++;
    }
//cout<<"ici"<<getchar();
return(sortie);
}




int CFichier_genepop::affiche_nb_alleles() {
//int gnourf=loci.size(); // loci.size() passe pas comme arg (by ref) de max
int /*borne=std::max(gnourf,49),*/ligne,colonne;
	for (int ii=8;ii<13;ii++) {
        _gotoxy(0,ii); cout<<"                                                                   ";
    }
	_gotoxy(3,9); cout<<"Largest allele index detected:\n";
	for (unsigned int ii=0;ii<loci.size();ii++) {
		ligne=11+ii-((ii+1)/10)*10;
		colonne=std::min((int) (5+((ii+1)/10)*15),65); //15 chars per col
		_gotoxy(colonne,ligne);
		if ((ligne==19)&&(colonne==65)&&(loci.size()>49)) cout<<"Etc.\n";
		else {
			cout<<loci[ii]->locusName.substr(0,8);
			_gotoxy(colonne+8,ligne);
			if (loci[ii]->galleles.size()>0)
			   cout<<": "<<loci[ii]->galleles.rbegin()->first;
			else if (loci[ii]->alleles.size()>0)
			   cout<<": "<<loci[ii]->alleles.rbegin()->first;
			else
			   cout<<": 0";
		}
	}
	_gotoxy(1,20);
	return 0;
}





// creation du fichier.in a partir du fichier genepop deja lu
int CFichier_genepop::createFichierIN(){
vector<CPopulation *>::iterator pP;
vector<CLocusGP *>::iterator pL;
time_t t1;//pour avoir la date et l'heure
struct tm *currentTime;


//ouverture du fichier de sortie FICHIER.IN
ofstream fichier_out("fichier.in", ios::out);

//remplissage du fichier

fichier_out << fileName << endl;//nom du fichier genepop
fichier_out << " " << pops.size() << "  " << loci.size() << endl;//nbre de pops et nbre de locus
for(pL = loci.begin();pL != loci.end();pL++){ //NOM d'allele max
fichier_out << " " << (*pL)->alleleMax << "  " << (*pL)->locusName << endl;
}
for(pP = pops.begin(); pP != pops.end(); pP++) {
			fichier_out << (*pP)->popName() <<endl;
		}

time(&t1);
currentTime=localtime(&t1);
//affiche la date
fichier_out << (currentTime->tm_mon) +1 << "-" << currentTime->tm_mday << "-" << currentTime->tm_year + 1900 << endl;
//affiche l'heure
fichier_out <<currentTime->tm_hour << ":" << currentTime->tm_min << ":" << currentTime->tm_sec << endl;

fichier_out.close();
return(0);

}
// ------------------------------------------------------------------------------------------
// remplit coord_pop;  ------------------------------------------------------------------------------
// n'itilise pas fichier.in
int CFichier_genepop::extract_coord_pop() {
vector<CPopulation *>::iterator pP;
vector<vector<double> >::iterator popu;
//int pos;
string buf,var;
char bidon;
stringstream stst(stringstream::in|stringstream::out);
     coord_pop.resize(pops.size());
//cerr<<pops.size()<<endl;
     popu=coord_pop.begin();
     for(pP = pops.begin(); pP != pops.end(); pP++) {
			(*popu).resize(2);
//cout<<(*pP)->popName();
//			stst.str((*pP)->popName()); // reads string (*pP)->popName()
			stst<<(*pP)->popName(); // reads string (*pP)->popName()
			stst>>(*popu)[0]>> (*popu)[1]; // writes to doubles
//cerr<<"\nCoordinates were read as '"<<(*popu)[0]<<"' and '"<< (*popu)[1]<<endl;getchar();
			if (isnan((*popu)[0]) || isnan((*popu)[1])) {
// unlikely to be useful. >> Assigns 0 when nan are read.
               cerr<<"Population coordinates not numeric for population "<<(*pP)->popName()<<":";
               cerr<<"Coordinates were read as '"<<(*popu)[0]<<"' and '"<< (*popu)[1]<<"'\n";
               cerr<<"Exiting...";
               if (cinGetOnError) cin.get();
               exit(-1);
            }
			while (stst>>bidon); // clear le contenu !
			stst.clear(); // clear les bits !
/*			ancienne version avant .clear()
            buf=(*pP)->popName();
			while (buf[0]==' ') buf.erase(0,1);
			pos=buf.find(' ');
			{stringstream ss(stringstream::in|stringstream::out);
			ss<<buf.substr(0,pos);
			ss>>(*popu)[0];
			}
			buf=buf.substr(pos+1,(*pP)->popName().length());
			while (buf[0]==' ') buf.erase(0,1);
			{stringstream ss(stringstream::in|stringstream::out);
			ss<<buf;
			ss>>(*popu)[1];
//cout<<(*popu)[0]<<" "<<(*popu)[1]<<endl;getchar();
			}*/
		    popu++;
     }
return(0);
}

/*int CFichier_genepop::renum3(char* Nnom) {
// structure correcte mais sortie incorrecte: 001110 deviendrait 1110
// pas d'utilit['e] pour l'instant mais garder la structure
ofstream Nin(Nnom);
	Nin<<"Renumbered from "<<Nnom<<endl;
	for (vector<CLocusGP *>::iterator p=loci.begin();p<loci.end();p++)
		Nin<<(*p)->locusName<<endl;
	for (vector<CPopulation *>::iterator p=pops.begin();p<pops.end();p++) {
		Nin<<"Pop\n";
		for (vector<CIndividual *>::iterator ii=(*p)->inds.begin();ii<(*p)->inds.end();ii++) {
			Nin<<(*ii)->indName()<<", ";
			for (int ll=0;ll<loci.size();ll++) {
				Nin<<(*ii)->getTypage(ll,1)<<(*ii)->getTypage(ll,2)<<" ";
			}
			Nin<<endl;
		}
	}
	Nin.close();
cout<<"****************";getchar();
return 0;
}
*/

/*********************************************************************/
int CFichier_genepop::computeCheckWriteDistMat(const char file[]) {
    /// computes distances
    /// checks values
    /// writes
#ifdef VERBOSE
	cout<<"debut computeCheckWriteDistMat"<<endl<<flush;
#endif
double maxpot=0;

vector<vector<double> >::iterator popi,popj;
  double pot,carre1,carre2;
  FILE *stre;;
  if((stre=fopen(file,"a"))==0)
    {cerr<<"computeCheckWriteDistMat cannot open file "<<file;if (cinGetOnError) cin.get();exit(1);}
  for(popj=coord_pop.begin()+1;popj!=coord_pop.end();popj++) {
     for(popi=coord_pop.begin();popi!=popj;popi++) {//cout<<(*popj)[0]<<" "<<(*popi)[0]<<" "<<(*popj)[1]<<" "<<(*popi)[1];getchar();
       carre1=pow((*popj)[0]-(*popi)[0],2);
       carre2=pow((*popj)[1]-(*popi)[1],2);
       pot=sqrt(carre1+carre2);
       if (pot>maxpot) maxpot=pot;
       //printf("distance ii jj =%f",pot);
       //getchar();
       fprintf(stre,"%.15le ",pot);
      }
      fprintf(stre,"\n");
   }
   fclose(stre);
#ifdef VERBOSE
	cout<<"fin computeCheckWriteDistMat"<<endl<<flush;
#endif
if (maxpot==0) return -1; else return 0;
}
/********************************************************************/




// ------------------------------------------------------------------------------------------
// CGenotypes -------------------------------------------------------------------------------
// tentative template version dans genotypes.h

// renvoie vrai si le g['e]notype <genotype> existe
bool CGenotypes::genotypeExists(int genotype) {
	map<int, int>::iterator p = genotypes.find(genotype);
	return(!(p == genotypes.end()));
}

// renvoie l'effectif du g['e]notype <genotype>
int CGenotypes::getEffective(int genotype) {
	if (genotypeExists(genotype)) { return genotypes[genotype]; } else { return 0; }
}


// cr['e]e ou incr['e]mente le g['e]notype <genotype> ; renvoie le nb de g['e]notypes
int CGenotypes::declareGenotype(int genotype) {
	if (genotypeExists(genotype)) {
		genotypes[genotype]++;
	} else {
		genotypes[genotype] = 1;
//if (alwaysAskBool) {cerr<<"from CGenotypes::declareGenotype "<<genotype<<" "<<genotypes.size();getchar();}
	}
	_sum ++;
	return(genotypes.size());
}

int CGenotypes::getMaxIdx() {
    return genotypes.rbegin()->first;
}


// remplit les g['e]notypes avec les informations de typages issues du fichier genepop
/*id CGenotypes::fillGenotypes(int locgeno, CPopulation *pop) {
	for (int iInd=0; iInd < pop->inds.size(); iInd ++) {
		if (pop->inds[iInd]->isValid(locgeno)) {
			declareGenotype(pop->inds[iInd]->getMaxAllele(locgeno) * 100 + pop->inds[iInd]->getMinAllele(locgeno));
		}
	}
}*/

// remplit les g['e]notypes avec les informations de typages issues du fichier genepop
void CGenotypes::fillGenotypes(int celocus, CPopulation *pop,char coding) {
//loc celocus doit etre un indice de locus, pas un typage
//getMax... recupere les alleles Maxet Min du typage (celocus pour l'indiv iInd), non du locus dans la pop
	for (unsigned int iInd=0; iInd < pop->inds.size(); iInd ++) {
/*cerr<<celocus<<" "<<pop->inds[iInd]->isValid(celocus)<<endl;
if (celocus==1) {
if (pop->inds[iInd]->isValid(celocus)) cerr<<"true"; else cerr<<"false";
cerr<<pop->inds[iInd]->getMaxAllele(celocus)<<" "<<pop->inds[iInd]->getMinAllele(celocus);
getchar();
}*/
		if (pop->inds[iInd]->isValid(celocus)) {
			if (coding==4) declareGenotype(pop->inds[iInd]->getMaxAllele(celocus) * 100 + pop->inds[iInd]->getMinAllele(celocus));
            else if (coding==6) declareGenotype(pop->inds[iInd]->getMaxAllele(celocus) * 1000 + pop->inds[iInd]->getMinAllele(celocus));
            else {// declareGenotype(celocus); incorrect [`a] cause de la remarque sur "celocus"
//                 cerr<<"call to CGenotypes::fillGenotypes() for haploid data.";getchar();
// not an error to be there when isolde_etc -> genotip2 calls it
//                 cout<<"\nLocus "<<celocus+1<<" will be ignored (haploid data)";
// this doesn't work because pauseGP is not known to genepop.cpp (needs change ?)
//                 if (pauseGP) { cout<<"(Return) to continue"<<endl; getchar();}
            }
		}
	}
}


/*void CGenotypes::fillAlleles(int geno, CPopulation *pop,char coding) {
	for (int iInd=0; iInd < pop->inds.size(); iInd ++) {
		if (pop->inds[iInd]->isValid(geno)) {
			if (coding==4) {
               declareGenotype(geno/100);
               declareGenotype(geno%100);
            }
            else if (coding==6) {
               declareGenotype(geno/1000);
               declareGenotype(geno%1000);
            }
            else declareGenotype(geno);
		}
	}
}*/


// purge
void CGenotypes::clear() {
	genotypes.clear();
	resetIterator();
	_sum = 0;
}

// it['e]rateur sur les cl['e]s
void CGenotypes::resetIterator() { iter = genotypes.begin(); }

int CGenotypes::getNext() {
	int resu;
	if (iter == genotypes.end()) {
//if (alwaysAskBool) {cerr<<"CGenotypes::getNext() = genotypes.end()";getchar(); }
		resu = -1;
	} else {
		resu = iter->first;
//if (alwaysAskBool) {cerr<<"CGenotypes::getNext() "<<iter->first<<" "<<resu;getchar(); }
		iter ++;
	}
	return resu;
}

// renvoie le nombre total de g['e]notypes (doit etre = au nb d'individus si pas de trous 0000)
int CGenotypes::getSum() {return _sum;}

// renvoie le nombre de g['e]notypes diff['e]rents
unsigned int CGenotypes::getNumber() {return genotypes.size();} //nb complete genotypes

// sortie d'une ligne (['e]ventuellement repli['e]e) au format "genotip2" des g['e]notypes,
//     le rang d['e]finit lequel des deux all[`e]les (faible ou fort) on ['e]crit.
void CGenotypes::printKeys(int rank, string *output, int keyWidth, int foldLines) {
	int genotype; // g['e]notype courant
	int colCount; // comptage des colonnes en sortie
	char buf[10];
	char fmt[10];

	resetIterator();
	sprintf(fmt, "%%-%dd", keyWidth);
	colCount = output->length();
	while ((genotype = getNext()) >= 0) {
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
// sortie d'une ligne (['e]ventuellement repli['e]e) des effectifs des g['e]notypes, le mod[`e]le donne la liste des g['e]notypes
//     dont on doit donner l'effectif
void CGenotypes::printValues(CGenotypes *model, string *output, int effectiveWidth, int foldLines) {
	int genotype; // g['e]notype courant
	int colCount; // comptage des colonnes en sortie
	char buf[10];

	colCount = output->length();
	model->resetIterator();
	while ((genotype = model->getNext()) >= 0) {
//cout<<genotype;getchar();
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


