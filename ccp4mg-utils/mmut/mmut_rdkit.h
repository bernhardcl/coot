#ifndef _MMUT_RDKIT_WRAP_
#define _MMUT_RDKIT_WRAP_
#include "mmdb2/mmdb_manager.h"
#include "mmdb2/mmdb_atom.h"
#include "ccp4srs/ccp4srs_monomer.h"
namespace RDKit {
class RWMol;
};
const char* SmilesToPDB(const char *smile,const char *fname,const char *cifname=NULL, int nconf=20, int maxIters=200);
int MolMinimize(RDKit::RWMol *mol, int nconf=20, int maxIters=200, ccp4srs::Monomer *monomer=NULL);
int MMDBMinimize(mmdb::Manager *molHnd, int nconf=20, int maxIters=200, ccp4srs::Monomer *monomer=NULL);
const char* RDKitLoadFileToPDB(const char *input_fname,const char *fname,const char *cifname=NULL);
RDKit::RWMol *RDKitAromaticityAnalyze(mmdb::PResidue res, bool aromaticize);
std::string RDKitResidueToSVG(RDKit::RWMol *);
#endif // _MMUT_RDKIT_WRAP_
