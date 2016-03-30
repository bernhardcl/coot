/*
     mmut/mmut_bonds.cc: CCP4MG Molecular Graphics Program
     Copyright (C) 2001-2008 University of York, CCLRC
     Copyright (C) 2009-2010 University of York
     Copyright (C) 2012 STFC

     This library is free software: you can redistribute it and/or
     modify it under the terms of the GNU Lesser General Public License
     version 3, modified in accordance with the provisions of the 
     license to address the requirements of UK law.
 
     You should have received a copy of the modified GNU Lesser General 
     Public License along with this library.  If not, copies may be 
     downloaded from http://www.ccp4.ac.uk/ccp4license.php
 
     This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU Lesser General Public License for more details.
*/


#include <math.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include "mmdb2/mmdb_manager.h"
#include "mmdb2/mmdb_math_graph.h"
#include "mmdb2/mmdb_mattype.h"
#include "mmut_bonds.h"
#include "mmut_manager.h"
#include "mman_manager.h"
#include "mmut_rdkit.h"
#include "mmdb2/mmdb_tables.h"
#include "ccp4srs/ccp4srs_monomer.h"
#include "mgutil.h"
#include "GraphMol/RWMol.h"

using namespace std;

typedef mmdb::Atom** PPCAtom;
typedef mmdb::Atom* PCAtom;
typedef mmdb::Atom CAtom;
typedef mmdb::Residue** PPCResidue;
typedef mmdb::Residue* PCResidue;
typedef mmdb::Residue CResidue;
typedef mmdb::Chain** PPCChain;
typedef mmdb::Chain* PCChain;
typedef mmdb::Chain CChain;
typedef mmdb::Model** PPCModel;
typedef mmdb::Model* PCModel;
typedef mmdb::Model CModel;
typedef mmdb::Contact* PSContact;
using namespace mmdb;

// Stores the trimmed input string into the given output buffer, which must be
// large enough to store the result.  If it is too small, the output is
// truncated.
size_t trimwhitespace(char *out, size_t len, const char *str)
{
  if(len == 0)
    return 0;

  const char *end;
  size_t out_size;

  // Trim leading space
  while(isspace(*str)) str++;

  if(*str == 0)  // All spaces?
  {
    *out = 0;
    return 1;
  }

  // Trim trailing space
  end = str + strlen(str) - 1;
  while(end > str && isspace(*end)) end--;
  end++;

  // Set output size to minimum of trimmed string length and buffer size minus 1
  out_size = (end - str) < len-1 ? (end - str) : len-1;

  // Copy trimmed string and add null terminator
  memcpy(out, str, out_size);
  out[out_size] = 0;

  return out_size;
}

realtype MGCovalentRadius[nElementNames] = {
    0.28, 0.31,
    1.28, 0.96, 0.84, 0.76, 0.71, 0.66, 0.57, 0.58,
    1.66, 1.41, 1.21, 1.11, 1.07, 1.05, 1.02, 1.06,
    2.03, 1.76,
          1.70, 1.60, 1.53, 1.39, 1.39, 1.32,
          1.26, 1.24, 1.32, 1.22,
                1.22, 1.20, 1.19, 1.20, 1.20, 1.16,
    2.20, 1.95,
          1.90, 1.75, 1.64, 1.54, 1.47, 1.46,
          1.42, 1.39, 1.45, 1.44,
                1.42, 1.39, 1.39, 1.38, 1.39, 1.40,
    2.44, 2.15,
          2.07, 2.04, 2.03, 2.01, 1.99, 1.98, 1.93,
          1.96, 1.94, 1.92, 1.92, 1.89, 1.90, 1.87,
                1.87, 1.75, 1.70, 1.62, 1.51, 1.44,
                1.41, 1.36, 1.36, 1.32,
                      1.45, 1.46, 1.48, 1.40, 1.50, 1.50,
    2.60, 2.21,
          2.15, 2.04, 2.00, 1.96, 1.90, 1.87, 1.80,
          1.65, 1.65, 1.65, 1.65, 1.65, 1.65, 1.65,
                1.69, /*SJM*/ 0.32, 0.10, /**/
                                  0.20, 0.20, 0.20,
                0.20, 0.20, 0.20, 0.20,
                      0.20, 0.20, 0.20,
    0.32, 0.32
  };

int GetGraphMatches(mmdb::Residue *R, ccp4srs::Monomer* monomer, ccp4srs::Manager *SRS, mmdb::ivector &anmatch, mmdb::cpstr altLoc=NULL){

//#define DEBUG_GRAPH_MATCH 1

  int rcSRS;
  mmdb::math::Graph *gSRS = monomer->getGraph(&rcSRS);
  //std::cout << "getGraphSRS rc " << rcSRS << "\n";

  gSRS->RemoveChirality();

  gSRS->ExcludeType(mmdb::getElementNo(mmdb::pstr("H")));
  rcSRS = gSRS->Build(false);
  //std::cout << "Build rc(SRS) " << rcSRS << "\n";

  mmdb::math::Graph *gPDB = new  mmdb::math::Graph();

  mmdb::Residue *altR = new mmdb::Residue();
  std::vector<int> altLocMatch;
  int iadd = 0;
  for(int iat=0;iat<R->GetNumberOfAtoms();iat++){
    if(!R->GetAtom(iat)->isTer()&&(!altLoc||strlen(altLoc)==0||!R->GetAtom(iat)->altLoc||strlen(R->GetAtom(iat)->altLoc)==0||strncmp(altLoc,R->GetAtom(iat)->altLoc,10)==0)){
      if(strncmp(R->GetAtom(iat)->element," H",2)!=0){
        mmdb::Atom *atom = new mmdb::Atom();
        char atID[1024];
        R->GetAtom(iat)->GetAtomID(atID);
        atom->Copy(R->GetAtom(iat));
        altR->AddAtom(atom);
        altLocMatch.push_back(iat);
        //std::cout << "Adding " << atID << "\n";
      }
    }
  }

  //std::cout << "altLocMatch.size(): " << altLocMatch.size() << "\n";

  int pdbMakeGraphErr = gPDB->MakeGraph(altR,altLoc);
  //std::cout << "MakeGraph rc " << pdbMakeGraphErr << "\n";

  if(pdbMakeGraphErr==mmdb::math::MKGRAPH_Ok){
#ifdef DEBUG_GRAPH_MATCH
    std::cout << "MakeGraph OK!\n";
#endif
  } else if(pdbMakeGraphErr==mmdb::math::MKGRAPH_NoAtoms){
    std::cout << "MakeGraph no atoms!\n";
  } else if(pdbMakeGraphErr==mmdb::math::MKGRAPH_ChangedAltLoc){
    std::cout << "MakeGraph changed altLoc!\n";
  } else if(pdbMakeGraphErr==mmdb::math::MKGRAPH_MaxOccupancy){
    std::cout << "MakeGraph max occupancy!\n";
  }


  gPDB->RemoveChirality();

  gPDB->ExcludeType(mmdb::getElementNo(mmdb::pstr("H")));
  int rc = gPDB->Build(false);
#ifdef DEBUG_GRAPH_MATCH
  std::cout << "Build rc " << rc << "\n";

  std::cout << "graph vertices from PDB " << gPDB->GetNofVertices() << "\n";
  std::cout << "graph edges from PDB " << gPDB->GetNofEdges() << "\n";
  std::cout << "graph vertices from SRS " << gSRS->GetNofVertices() << "\n";
  std::cout << "graph edges from SRS " << gSRS->GetNofEdges() << "\n";
#endif

  if(gPDB->GetNofVertices()>gSRS->GetNofVertices()){
    // This is harsh, but we are currently only trying to see if a ligand is what it claims to be.
    altR->DeleteAllAtoms();
    delete altR;
    delete gPDB;
    delete gSRS;
    return 1;
  }

  mmdb::math::GraphMatch *U = new mmdb::math::GraphMatch();

#ifdef DEBUG_GRAPH_MATCH
  std::cout << "PDB Graph:\n\n";
  gPDB->Print();
  std::cout << "SRS Graph:\n\n";
  gSRS->Print();
#endif

  std::cout.flush();
  fflush(stdout);

  int nvPDB, nvSRS;
  mmdb::math::PPVertex vPDB, vSRS; 
  gPDB->GetVertices(vPDB,nvPDB);
  gSRS->GetVertices(vSRS,nvSRS);
  
  mmdb::Residue *newR = new mmdb::Residue();

  // Aargh, aargh, aargh!!! altLocs will bite ....
  std::vector<int> vmatch(gSRS->GetNofVertices());
  for(int ivSRS=0;ivSRS<nvSRS;ivSRS++){
    vmatch[ivSRS] = -1;
  }

  int imatch = 0;
  char tmp[4];
  char tmp2[4];
  for(int ivPDB=0;ivPDB<nvPDB;ivPDB++){
    trimwhitespace(tmp,4,vPDB[ivPDB]->GetName());
    bool matched = false;
    for(int ivSRS=0;ivSRS<nvSRS;ivSRS++){
      trimwhitespace(tmp2,4,vSRS[ivSRS]->GetName());
      if(vSRS[ivSRS]->GetName()&&vPDB[ivPDB]->GetName()){
        if(strcmp(tmp2,tmp)==0){
#ifdef DEBUG_GRAPH_MATCH
          std::cout << vSRS[ivSRS]->GetName() << " " << vPDB[ivPDB]->GetName() << " " << ivSRS << " " << ivPDB << "\n";
#endif
          imatch++;
          matched = true;
          vmatch[ivSRS] = ivPDB;
          break;
        }
      }
    }
    if(!matched){
      std::cout << vPDB[ivPDB]->GetName() << " no match\n";
    }
  }

  //std::cout << imatch << "\n\n";

  //for(int ivSRS=0;ivSRS<nvSRS;ivSRS++){
    //std::cout << ivSRS << " " << vmatch[ivSRS] << "\n";
  //}
  for(int ivSRS=0;ivSRS<nvSRS;ivSRS++){
    if(vmatch[ivSRS]==-1){
      for(int ivPDB=0;ivPDB<nvPDB;ivPDB++){
        if(std::find(vmatch.begin(),vmatch.end(),ivPDB)==vmatch.end()){
          // We'll take the first one not accounted for. Gotta do something simple (for starters).
          char tmp[4];
          trimwhitespace(tmp,4,vPDB[ivPDB]->GetName());
          std::cout << ivPDB << " " << tmp << " not yet accounted for\n";
          vmatch[ivSRS] = ivPDB;
          break;
        }
      }
    }
  }

  std::vector<int> vmatch2;
  for(unsigned ivSRS=0;ivSRS<vmatch.size();ivSRS++){
    //if(vmatch[ivSRS]!=-1){
      vmatch2.push_back(vmatch[ivSRS]);
    //}
  }
  std::vector<int> rvmatch2(vmatch2.size());
  //std::cout << "\nAfter second pass and filter\n\n";
  for(int ivSRS=0;ivSRS<vmatch2.size();ivSRS++){
    rvmatch2[ivSRS] = -1;
  }
  std::vector<int> skipmatch;
  int iskip=0;
  for(int ivSRS=0;ivSRS<vmatch2.size();ivSRS++){
    //std::cout << ivSRS << " " << vmatch2[ivSRS] << "\n";
    if(vmatch2[ivSRS]!=-1){
      rvmatch2[vmatch2[ivSRS]] = ivSRS;
      skipmatch.push_back(iskip);
      iskip++;
    } else {
      skipmatch.push_back(-1);
    }
  }
  //std::cout << "\nSkip map\n\n";
  //for(int ivSRS=0;ivSRS<vmatch2.size();ivSRS++){
    //std::cout << ivSRS << " " << skipmatch[ivSRS] << "\n";
  //}

  //std::cout << "\nReverse map\n\n";
  //for(int ivSRS=0;ivSRS<rvmatch2.size();ivSRS++){
    //std::cout << ivSRS << " " << rvmatch2[ivSRS] << "\n";
  //}

  for(int ivSRS=0;ivSRS<vmatch2.size();ivSRS++){
    if(vmatch2[ivSRS]>=gPDB->GetNofVertices()){
      // We've gone badly wrong!
      std:: cout << "We've gone badly wrong!\n";
      delete U;
      newR->DeleteAllAtoms();
      delete newR;
      altR->DeleteAllAtoms();
      delete altR;
      delete gPDB;
      delete gSRS;
      return -1;
    }
    if(vmatch2[ivSRS]!=-1){
      mmdb::Atom *atom = new mmdb::Atom();
      atom->Copy(altR->GetAtom(vmatch2[ivSRS]));
      //char AtomID[1024];
      //altR->GetAtom(vmatch2[ivSRS])->GetAtomID(AtomID);
      //std::cout << "Adding " << AtomID << "\n";
      newR->AddAtom(atom);
    }
  }

  mmdb::math::Graph *gPDBNew = new  mmdb::math::Graph();
  pdbMakeGraphErr = gPDBNew->MakeGraph(newR,altLoc);
  //std::cout << "MakeGraph rc " << pdbMakeGraphErr << "\n";

  if(pdbMakeGraphErr==mmdb::math::MKGRAPH_Ok){
    //std::cout << "MakeGraph OK!\n";
  } else if(pdbMakeGraphErr==mmdb::math::MKGRAPH_NoAtoms){
    std::cout << "MakeGraph no atoms!\n";
  } else if(pdbMakeGraphErr==mmdb::math::MKGRAPH_ChangedAltLoc){
    std::cout << "MakeGraph changed altLoc!\n";
  } else if(pdbMakeGraphErr==mmdb::math::MKGRAPH_MaxOccupancy){
    std::cout << "MakeGraph max occupancy!\n";
  }


  gPDBNew->RemoveChirality();

  gPDBNew->ExcludeType(mmdb::getElementNo(mmdb::pstr("H")));
  rc = gPDBNew->Build(false);
  //std::cout << "Build rc(new) " << rc << "\n";

#ifdef DEBUG_GRAPH_MATCH
  std::cout << "New PDB Graph:\n\n";
  gPDBNew->Print();
#endif

  if(gPDB->GetNofVertices()!=gPDBNew->GetNofVertices()||gPDB->GetNofEdges()!=gPDBNew->GetNofEdges()){
    std::cout << "Error in creating permuted graph \n";
    delete U;
    newR->DeleteAllAtoms();
    delete newR;
    altR->DeleteAllAtoms();
    delete altR;
    delete gPDB;
    delete gPDBNew;
    delete gSRS;
    return -1;
  }

  if(abs(gSRS->GetNofVertices()-gPDBNew->GetNofVertices())<4&&abs(gSRS->GetNofEdges()-gPDBNew->GetNofEdges())<4&&gSRS->GetNofEdges()>0&&gPDBNew->GetNofEdges()>0){
    //std::cout << "Now matching graphs...\n";
    U->SetTimeLimit(10);
    U->MatchGraphs(gPDBNew,gSRS,1,true,mmdb::math::EXTTYPE_Ignore);
    if(U->GetStopSignal()){
      std::cout << "Timed out !!!!\n";
      delete U;
      newR->DeleteAllAtoms();
      delete newR;
      altR->DeleteAllAtoms();
      delete altR;
      delete gPDB;
      delete gPDBNew;
      delete gSRS;
      return 1;
    }
  
#ifdef DEBUG_GRAPH_MATCH
    std::cout << "number of input PDB vertices " << gPDBNew->GetNofVertices() << "\n";
    std::cout << "number of input SRS vertices " << gSRS->GetNofVertices() << "\n";
    std::cout << "number of graph matches " << U->GetNofMatches() << "\n";
    std::cout << "max graph match size " << U->GetMaxMatchSize() << "\n";
    U->PrintMatches();
#endif
  
    // We assume (blindly) that atoms are at least in usual order. There are probably
    // better ways.
    int mindist = std::numeric_limits<int>::max();
    int minmatch = -1;
    for(int imatch=0;imatch<U->GetNofMatches();imatch++){
      int dist = 0;
      int nv;
      mmdb::realtype p1,p2;
      mmdb::ivector FV1,FV2;
      U->GetMatch ( imatch, FV1, FV2, nv, p1, p2 );
      //std::cout << "Match " << imatch << " " << nv << " " << p1 << " " << p2 << "\n";
      for(int iv=1;iv<=nv;iv++){
#ifdef DEBUG_GRAPH_MATCH
        std::cout << FV1[iv] << " " << FV2[iv] << "\n";
#endif
        dist += abs(FV1[iv]-FV2[iv]);
      }
      //std::cout << "dist: " << dist << "\n\n";
      if(dist<mindist) {
        minmatch = imatch;
        mindist = dist;
      }
      if(dist==0) break;
    }
  
    //std::cout << "Best match " << minmatch << "\n";
    if(minmatch>-1){
      int nv;
      mmdb::realtype p1,p2;
      mmdb::ivector FV1,FV2;
      U->GetMatch ( minmatch, FV1, FV2, nv, p1, p2 );
      //std::cout << "Match " << minmatch << " " << nv << " " << p1 << " " << p2 << "\n"; std::cout.flush();
      for(int iv=1;iv<=nv;iv++){
        //std::cout << FV1[iv] << " " << FV2[iv] << "\n"; std::cout.flush();
        //std::cout << "   " << rvmatch2[FV1[iv]-1] << " " << FV2[iv]-1 << "\n"; std::cout.flush();
        if(true||(FV1[iv]-1<rvmatch2.size()&&FV2[iv]-1<rvmatch2.size())){
          if(rvmatch2[FV1[iv]-1]!=-1&&FV2[iv]-1<skipmatch.size()){
            anmatch[rvmatch2[FV1[iv]-1]] = altLocMatch[skipmatch[FV2[iv]-1]];
          }
        }
      }
      delete U;
      newR->DeleteAllAtoms();
      delete newR;
      altR->DeleteAllAtoms();
      delete altR;
      delete gPDB;
      delete gPDBNew;
      delete gSRS;
      return 0;
    }
  }
  delete U;
  newR->DeleteAllAtoms();
  delete newR;
  altR->DeleteAllAtoms();
  delete altR;
  delete gPDB;
  delete gPDBNew;
  delete gSRS;
  return 1;
}

//-----------------------------------------------------------------------
CMolBondParams::CMolBondParams() {
//-----------------------------------------------------------------------
   //printf ("Sbase maxAtomInRes %i\n",sbase->maxAtomInRes);
   interResCut = 2.4; 
   intraResCut = 2.4;
   maxBondRad = 2.4;
   // Mx bond length is ( at1->maxBondRad + at2->maxBondRad)
   // Make some allowance upwards
   maxBondRadFactor = 1.2; 
}

//-----------------------------------------------------------------------
CMolBondParams::~CMolBondParams() {
//-----------------------------------------------------------------------
 
}


//-----------------------------------------------------------------------
CMolBonds::CMolBonds(const PCMMUTManager molHndin, CMolBondParams *paramsin) :
  CMMANBase ( molHndin ) {
//-----------------------------------------------------------------------

  params = paramsin;
  
  sbaseCompoundID = NULL;
  sbaseAtomIndex = NULL;

  nAtoms = 0;
  nRes = 0;

}
//----------------------------------------------------------------------
CMolBonds::~CMolBonds(){
//----------------------------------------------------------------------
  //cout << "CMolBonds destructor" << endl;
}
 
//----------------------------------------------------------------------
std::string CMolBonds::FindBonds ( int udd_sbaseCompoundID,
        int udd_sbaseAtomOrdinal, int udd_atomEnergyType, const std::map<std::string,std::string> &customResCIFFiles, bool aromaticize, bool checkGraphs) {
//----------------------------------------------------------------------

  //mmdb::io::PFile graphFile = dynamic_cast<PCMMANManager>(molHnds[0])->GetSRS()->getGraphFile();
  mmdb::io::PFile structFile = dynamic_cast<PCMMANManager>(molHnds[0])->GetSRS()->getStructFile();

  std::ostringstream output;
  char AtomID[30];

  PPCAtom selAtom = NULL;
  PPCAtom selAtom0 = NULL;
  PPCResidue selRes = NULL;
  PCChain pCh;
  PCModel pMdl;
  ivector nMatchAtom = NULL;
  imatrix matchAtom = NULL;
  int nAtominRes,nAtoms,ir,i,j,na1,na2,ia1,ia2,ib,selH;
  int nModels,nm,nAtominModel; 

  int nr;
  PSContact contacts = NULL;
  int ic,nContacts;
  PCAtom pa1,pa2;
  PCResidue pr1;

  int nAlt;
  int nAltTot = 0 ;
  pstr sbaseCompoundID;
  bool doContacts;  
  bool modelBondsSame = true;
  int maxModels = 1;
  int firstModel = 1;
  int lastModel = 1;
  int maxMatch = 5;

  int restype1,restype2;

  // Pointer to atoms in residue for all models
  // Big bovver if >100 models
  /*
  PPCAtom modelSelAtom[100];
  int nSelAtom[100];
  for (nm=1;nm<=100;nm++) {  modelSelAtom[nm]=NULL; }
  */

  // For NMR structure with multiple models we need to 
  // find the bonds in one model but populate the bonds data
  // structure for all models - but first make sure all models
  // have the same composition 
  nModels = molHnds[0]->GetNumberOfModels();

  // Avoid big bovver. May use a lot of memory.
  PPCAtom *modelSelAtom = new PPCAtom[nModels+1];
  int* nSelAtom = new int[nModels+1];
  for (nm=0;nm<=nModels;nm++) {  modelSelAtom[nm]=NULL; }

  if (nModels>1) {
    lastModel = nModels;
    firstModel = 0;
    PCModel pMdl = NULL;
    while ( firstModel < lastModel && pMdl == NULL ) {
      firstModel++;
      pMdl = molHnds[0]->GetModel(firstModel);
    }
    //cout << "firstModel " << firstModel << endl;
  }
  nAtominModel = molHnds[0]->GetModel(firstModel)->GetNumberOfAtoms(0);
  if (lastModel>firstModel) {
    for (nm=firstModel+1;nm<=lastModel;nm++) {
      pMdl = molHnds[0]->GetModel(nm);
      if (pMdl != NULL && pMdl->GetNumberOfAtoms(0) != nAtominModel) {
        modelBondsSame = false;
        output << "Models do not have the same number of atoms" << endl;
      }
    }
  }
  if ( !modelBondsSame) {
    maxModels = nModels;
  } else {
    maxModels = firstModel;
  }
  //cout << endl << "modelBondsSame " << modelBondsSame << " " << maxModels << endl;

  // Get memory for   imatrix matchAtom; 
  // the pointers to atom in the Sbase structure
  GetMatrixMemory(matchAtom,maxMatch,2000,0,0);
  GetVectorMemory(nMatchAtom,2000,0);

  // Loop over the unique models (normally just one model)
  for (int iMod=firstModel;iMod<=maxModels;iMod++) {

  selH = molHnds[0]->NewSelection();
  molHnds[0]->SelectAtoms(selH,iMod,"*", ANY_RES,"*", ANY_RES,"*","*","*","*","*",SKEY_NEW);
  molHnds[0]->SelectAtoms(selH,iMod,"*", ANY_RES,"*", ANY_RES,"*","*","*","H","*",SKEY_CLR);
  selAtom0 = NULL;
  molHnds[0]->GetSelIndex ( selH,selAtom0,nAtoms );

  if (contacts) delete [] contacts;
  contacts = NULL;
  nContacts = 0;
  molHnds[0]->SeekContacts(selAtom0,nAtoms,selAtom0,nAtoms,
		 0.0,params->interResCut,1,contacts,nContacts,0,NULL,0);

  if  (molHnds[0]->GetModel(iMod)!= NULL) {

    //cout << "iMod " << iMod << endl;
    for (int ich = 0; ich< molHnds[0]->GetModel(iMod)->GetNumberOfChains();ich++) {
    //Get the residue selections
    pCh = NULL;
    pCh = molHnds[0]->GetModel(iMod)->GetChain(ich);
    selRes = NULL;
    molHnds[0]->GetModel(iMod)->GetResidueTable (pCh->GetChainID(), selRes, nr);


    // To find the INTRA-residue bonds
    // Loop over all residues finding the required Sbase residue
    nAtominRes=0;
    for (ir = 0; ir < nr; ir++ ) {
      bool haveHs = false;
      molHnds[0]->GetAtomTable1(iMod,pCh->GetChainID(),ir,selAtom,nAtominRes);
      //std::cout << endl << "ich,ir,nAtominRes " << ich << " " << ir << " " << selAtom[0]->serNum << " " << nAtominRes << "\n";
      if (modelBondsSame && nModels>1) {
        for (nm=firstModel+1;nm<=lastModel;nm++) {
          molHnds[0]->GetAtomTable1(nm,pCh->GetChainID(),ir,modelSelAtom[nm],nSelAtom[nm]);
        }
      }

      doContacts = false;
        //doContacts = true;

      // FIXME - makeBonds throws in a lots of bonds to the same atom.... Possiblt SRS bug, reported to Eugene 24/06/2014.
      //         So I reimplement partially for now.
      
      mmdb::PPAtom   A;
      mmdb::ivector  anmatch;
      ccp4srs::Monomer*       monomer=NULL;
      int natoms, mon_natoms;
      selRes[ir]->GetAtomTable ( A,natoms );

      // Find inter-res close contacts
/*
      bool freeLigand = false;
*/

      //std::cout << "All to this res contacts: " << nContactsAll << "\n";
      bool freeLigand = true;
      for(int icontact=0;icontact<nContacts;icontact++){
        double dist = contacts[icontact].dist;
        pa1 = selAtom0[contacts[icontact].id1];
        pa2 = selAtom0[contacts[icontact].id2];
        if((pa1->residue==selRes[ir]||pa2->residue==selRes[ir])&&(pa1->residue!=pa2->residue) && CMMUTSRS::CheckCovalentDistance(pa1->element,pa2->element,dist)==1){
          freeLigand = false;
          break;
        }
      }
      if(natoms<2){
        freeLigand = false;
      }

      for (i=0;i<natoms;i++)
        if (A[i]) A[i]->FreeBonds();

      bool doneCIF = false;
      bool deleteMonomer;
      if(customResCIFFiles.count(std::string( selRes[ir]->GetResName()))>0){
        // FIXME - dodgy assumes data loop before atoms before bonds. No guarantee!
        std::string cifname = customResCIFFiles.find(std::string( selRes[ir]->GetResName()))->second;
        const char *cifname_s = cifname.c_str();
        int RC;
        mmdb::mmcif::Loop *Loop1;
        mmdb::mmcif::Loop *Loop2;
        mmdb::mmcif::Loop *Loop3;
        mmdb::mmcif::File* CIFFile =  new mmdb::mmcif::File();
        mmdb::mmcif::Data* CIF;
        RC = CIFFile->ReadMMCIFFile ( cifname_s );
        //std::cout << CIFFile->GetNofData() << "\n";
        if ( RC ) {
          printf ("Error reading %s\n",cifname_s);
        } else {
          std::map<std::string,int> atomNameMap;
          for(int idata=0;idata<CIFFile->GetNofData();idata++){
          mmdb::ivector anmatchCIF=NULL;
          CIF = CIFFile->GetCIFData(idata);
          Loop1 = CIF->GetLoop ( "_chem_comp" );
          Loop2 = CIF->GetLoop ( "_chem_comp_atom" );
          Loop3 = CIF->GetLoop ( "_chem_comp_bond" );
          if(Loop1){
            int nDescr=0;
            while(nDescr<Loop1->GetLoopLength()){
              mmdb::pstr id=0;
              mmdb::pstr three_letter_code=0;
              mmdb::pstr name=0;
              mmdb::pstr group=0;
              int number_atoms_all;
              RC = Loop1->GetString ( id, "id",nDescr);
              RC = Loop1->GetString ( three_letter_code, "three_letter_code",nDescr);
              RC = Loop1->GetString ( name, "name",nDescr);
              RC = Loop1->GetString ( group, "group",nDescr);
              RC = Loop1->GetInteger ( number_atoms_all, "number_atoms_all",nDescr);
              nDescr++;
            }
          }
          if(Loop2){
            int nLoopAtoms=0;
            mmdb::GetVectorMemory ( anmatchCIF,Loop2->GetLoopLength(),0 );
            while(nLoopAtoms<Loop2->GetLoopLength()){
              anmatchCIF[nLoopAtoms] = -1;
              mmdb::pstr comp_id=0;
              mmdb::pstr atom_id=0;
              mmdb::pstr type_symbol=0;
              mmdb::pstr type_energy=0;
              mmdb::realtype x,y,z;
              RC = Loop2->GetString ( comp_id, "comp_id",nLoopAtoms);
              RC = Loop2->GetString ( atom_id, "atom_id",nLoopAtoms);
              //char *trimmed = new char[32];
              if (!RC &&atom_id ){ 
                for(int j=0;(j<natoms) && (anmatchCIF[nLoopAtoms]<0);j++){
                  std::string s(A[j]->name);
                  s = ccp4mg_trim(s);
                  const char* trimmed = s.c_str();
                  if (!strcmp(trimmed,atom_id)){
                    anmatchCIF[nLoopAtoms] = j;
                    break;
                  }
                }
                atomNameMap[std::string(atom_id)] = nLoopAtoms;
              }
              RC = Loop2->GetString ( type_symbol, "type_symbol",nLoopAtoms);
              RC = Loop2->GetString ( type_energy, "type_energy",nLoopAtoms);
              RC = Loop2->GetReal ( x, "x",nLoopAtoms);
              RC = Loop2->GetReal ( y, "y",nLoopAtoms);
              RC = Loop2->GetReal ( z, "z",nLoopAtoms);
              nLoopAtoms++;
            }
          }
          if(Loop3){
            int nLoopBonds=0;
            while(nLoopBonds<Loop3->GetLoopLength()){
              mmdb::pstr comp_id=0;
              mmdb::pstr atom_id_1=0;
              mmdb::pstr atom_id_2=0;
              mmdb::pstr type=0;
              mmdb::realtype value_dist;
              RC = Loop3->GetString ( comp_id, "comp_id",nLoopBonds);
              RC = Loop3->GetString ( atom_id_1, "atom_id_1",nLoopBonds);
              RC = Loop3->GetString ( atom_id_2, "atom_id_2",nLoopBonds);
              RC = Loop3->GetString ( type, "type",nLoopBonds);
              RC = Loop3->GetReal ( value_dist, "value_dist",nLoopBonds);
              int i1 = anmatchCIF[atomNameMap[std::string(atom_id_1)]];
              int i2 = anmatchCIF[atomNameMap[std::string(atom_id_2)]];
              if ((i1>=0) && (i2>=0))  {
                  if(strlen(type)>5&&(!strncmp(type,"triple",6))){
                    A[i1]->AddBond ( A[i2],ccp4srs::Bond::Triple,2 );
                    A[i2]->AddBond ( A[i1],ccp4srs::Bond::Triple,2 );
                  } else if(strlen(type)>7&&(!strncmp(type,"aromatic",8))){
                    A[i1]->AddBond ( A[i2],ccp4srs::Bond::Aromatic,2 );
                    A[i2]->AddBond ( A[i1],ccp4srs::Bond::Aromatic,2 );
                  } else if(strlen(type)>5&&(!strncmp(type,"double",6))){
                    A[i1]->AddBond ( A[i2],ccp4srs::Bond::Double,2 );
                    A[i2]->AddBond ( A[i1],ccp4srs::Bond::Double,2 );
                  } else {
                    A[i1]->AddBond ( A[i2],ccp4srs::Bond::Single,2 );
                    A[i2]->AddBond ( A[i1],ccp4srs::Bond::Single,2 );
                  }
              } 
              nLoopBonds++;
            }
          }
          if(anmatchCIF) mmdb::FreeVectorMemory ( anmatchCIF,0 );
          }
        }
        doneCIF = true;
      } else {
        if(CMMUTSRS::MonomerCache.find(selRes[ir]->GetResName())!=CMMUTSRS::MonomerCache.end()){
          monomer = CMMUTSRS::MonomerCache[selRes[ir]->GetResName()];
          deleteMonomer = false;
        } else {
          monomer = dynamic_cast<PCMMANManager>(molHnds[0])->GetSRS()->getMonomer ( selRes[ir]->GetResName(), structFile );
          deleteMonomer = true;
        }
      }
      if(!monomer){
        doContacts = true;
      } else {
      mon_natoms = monomer->n_atoms();
      mmdb::GetVectorMemory ( anmatch,mon_natoms,0 );

      // Loop over altLocs here
      std::set<std::string> altLocs;
      for(int j=0;j<natoms;j++){
        if(strlen(A[j]->altLoc)>0){
          altLocs.insert(std::string(A[j]->altLoc));
        }
      }
      altLocs.insert("");

      for(i=0;i<mon_natoms;i++){
        anmatch[i] = -1;
      }

      std::set<std::string>::const_iterator altIter;
      for(altIter = altLocs.begin();altIter != altLocs.end();altIter++){

      const char *altLocCStr = (*altIter).c_str();

      int rc=-1;

      if(checkGraphs&&freeLigand&&!doneCIF){
        doContacts = false;
        rc = GetGraphMatches(selRes[ir],monomer,dynamic_cast<PCMMANManager>(molHnds[0])->GetSRS(),anmatch,altLocCStr);
#ifdef DEBUG_GRAPH_MATCH
        std::cout << "rc from GetGraphMatches " << rc << "\n";
#endif
        if(rc!=0){
          doContacts = true;
        } else {
          mmdb::AtomName anamei;
          mmdb::ivector  nonHmatch;
          mmdb::GetVectorMemory ( nonHmatch,mon_natoms,0 );
          int nonHidx = 0;
          int nHs = 0;
          for(int im=0;im<mon_natoms;im++){
            nonHmatch[im] = -1;
            mmdb::cpstr ele;
            ele = monomer->atom(im)->element(); 
            if(ele&&strlen(ele)==1&&strncmp(ele,"H",1)==0){
              nHs++;
              continue;
            }
            nonHmatch[im] = nonHidx;
            nonHidx++;
          }
          if(nHs>0){
            mmdb::ivector  anmatch2;
            mmdb::GetVectorMemory ( anmatch2,mon_natoms,0 );
            for(int im=0;im<mon_natoms;im++){
              anmatch2[im] = -1;
              if(nonHmatch[im]>-1){
                anmatch2[im] = anmatch[nonHmatch[im]];
              }
            }
            for(int im=0;im<mon_natoms;im++){
              anmatch[im] = anmatch2[im];
            }
            mmdb::FreeVectorMemory ( anmatch2,0 );
          }
          mmdb::FreeVectorMemory ( nonHmatch,0 );
        }
      }

      if(rc!=0&&!doContacts&&!doneCIF){
      
      for(i=0;i<mon_natoms;i++){
        mmdb::AtomName anamei;
        mmdb::cpstr anameo;
        anmatch[i] = -1;
        anameo = monomer->atom(i)->name_pdb ( anamei ); 
        //std::cout <<  monomer->atom(i)->id() << "\n";
        for(int j=0;(j<natoms) && (anmatch[i]<0);j++){
          //std::cout << "Test: --" <<  A[j]->name << "--, --" << anameo << "--\n";
          //std::cout << strcmp(A[j]->name," W1 ") << "\n";
          //std::cout << strcmp(anameo,"W1  ") << "\n";
          //std::cout << "--" << A[j]->element << "--\n";
          // FIXME - Oh why is anameo not right all the time?
          if(strcmp(A[j]->element," H")==0){
            haveHs = true;
          }
          if(!altLocCStr||strlen(A[j]->altLoc)==0||(!strncmp(A[j]->altLoc,altLocCStr,1))){
          if (!strcmp(A[j]->name,anameo)){
            anmatch[i] = j;
            //std::cout << "Match: " << A[j]->name << " " << " " <<  anameo << "\n";
            break;
          } else if(!strncmp(A[j]->name,"CL",2)&&!strncmp(anameo," CL",3)&&!strcmp(A[j]->element,"CL")&&!strncmp(anameo+1,A[j]->name,3)) {
            anmatch[i] = j;
            //std::cout << "Match(CL): " << A[j]->name << " " << " " <<  anameo << "\n";
            break;
          } else if(!strncmp(A[j]->name,"BR",2)&&!strncmp(anameo," BR",3)&&!strcmp(A[j]->element,"BR")&&!strncmp(anameo+1,A[j]->name,3)) {
            anmatch[i] = j;
            //std::cout << "Match(BR): " << A[j]->name << " " << " " <<  anameo << "\n";
            break;
          } else if(!strcmp(A[j]->name," W1 ")&&!strcmp(anameo,"W1  ")&&!strcmp(A[j]->element," W")) {
            anmatch[i] = j;
            //std::cout << "Match(W): " << A[j]->name << " " << " " <<  anameo << "\n";
            break;
          } else if(!strcmp(A[j]->name," V  ")&&!strcmp(anameo,"V   ")&&!strcmp(A[j]->element," V")) {
            anmatch[i] = j;
            //std::cout << "Match(W): " << A[j]->name << " " << " " <<  anameo << "\n";
            break;
          } else if(!strcmp(A[j]->name,"SE  ")&&!strcmp(anameo," SE ")&&!strcmp(A[j]->element,"SE")) {
            anmatch[i] = j;
            std::cout << "Match(SE): " << A[j]->name << " " << " " <<  anameo << "\n";
            break;
          }
          }
        }
      }
      }

      //std::cout << monomer->chem_name() << " " << monomer->chem_type() << " " << monomer->chem_formula() << "\n";
      //std::cout << rc << " " << int(doContacts) << ", nbonds: " << monomer->n_bonds() << "\n";
      //std::cout << "Free ligand? " << int(freeLigand) << "\n";
      int iadd = 0;
      int ibnonh = 0;
      if(!doContacts&&monomer->n_bonds()>0&&!(strncmp(selRes[ir]->GetResName(),"HOH",3)==0)&&!(strcmp(monomer->chem_type(),"L-PEPTIDE LINKING")==0)&&!(strcmp(monomer->chem_type(),"D-PEPTIDE LINKING")==0)&&!(strcmp(monomer->chem_type(),"PEPTIDE LINKING")==0)&&!(strcmp(monomer->chem_type(),"RNA LINKING")==0)&&!(strcmp(monomer->chem_type(),"DNA LINKING")==0)){
        for (i=0;i<monomer->n_bonds();i++)  {
          ccp4srs::Bond* bondi = monomer->bond(i);
          int a1 = bondi->atom1();
          int a2 = bondi->atom2();
          ccp4srs::Atom* at1 = monomer->atom(a1);
          ccp4srs::Atom* at2 = monomer->atom(a2);
          mmdb::cpstr at1ele = at1->element();
          mmdb::cpstr at2ele = at2->element();
          if((strcmp(at1ele,"H")==0)||(strcmp(at2ele,"H")==0)){
            continue;
          }
          ibnonh++;
        }
      }
      if(!doContacts){
      for (i=0;i<monomer->n_bonds();i++)  {
        ccp4srs::Bond* bondi = monomer->bond(i);
        int i1 = anmatch[bondi->atom1()];
        int i2 = anmatch[bondi->atom2()];
        if ((i1>=0) && (i2>=0))  {
          if(!freeLigand&&(strcmp(monomer->chem_type(),"L-PEPTIDE LINKING")==0||strcmp(monomer->chem_type(),"D-PEPTIDE LINKING")==0||strcmp(monomer->chem_type(),"PEPTIDE LINKING")==0)&&!(strcmp(monomer->chem_type(),"RNA LINKING")==0)&&!(strcmp(monomer->chem_type(),"DNA LINKING")==0)){
            iadd++;
          } else {
            iadd++;
          }
        } 
      }
      }
      if(!doContacts&&monomer->n_bonds()>0&&!(strncmp(selRes[ir]->GetResName(),"HOH",3)==0)&&!(strcmp(monomer->chem_type(),"L-PEPTIDE LINKING")==0)&&!(strcmp(monomer->chem_type(),"D-PEPTIDE LINKING")==0)&&!(strcmp(monomer->chem_type(),"PEPTIDE LINKING")==0)&&!(strcmp(monomer->chem_type(),"RNA LINKING")==0)&&!(strcmp(monomer->chem_type(),"DNA LINKING")==0)){
        // This is an attempt to test for badly named ligands.
        if((float(iadd)/ibnonh)<0.85){
          std::cout << "Possibly badly named ligand " << selRes[ir]->GetResName() << ". Found only " << iadd << " bonds\n";
          doContacts = true;
        }
      }

      if(!doContacts){
        for (i=0;i<monomer->n_bonds();i++)  {
          ccp4srs::Bond* bondi = monomer->bond(i);
          int i1 = anmatch[bondi->atom1()];
          int i2 = anmatch[bondi->atom2()];
#ifdef DEBUG_GRAPH_MATCH
          std::cout << "Should bond " << bondi->atom1() << " to " << bondi->atom2() << std::endl;
          mmdb::AtomName anamei1;
          mmdb::cpstr anameo1;
          anameo1 = monomer->atom(bondi->atom1())->name_pdb ( anamei1 ); 
          mmdb::AtomName anamei2;
          mmdb::cpstr anameo2;
          anameo2 = monomer->atom(bondi->atom2())->name_pdb ( anamei2 ); 
          std::cout << "Should bond " << anameo1 << " to " << anameo2 << " " << bondi->atom1() << " to " << bondi->atom2() << " " << i1 << " " << i2 << std::endl;
#endif

          //char AtomID1[1024];
          //char AtomID2[1024];
          if ((i1>=0) && (i2>=0) && (i1<natoms) && (i2<natoms))  {
            //std::cout << "Adding i1,i2: " << i1 << ", " << i2 << ", natoms: " << natoms << "\n";
            //std::cout << "monomer numbering " << bondi->atom1() << ", " << bondi->atom2() << "\n";
            //A[i1]->GetAtomID(AtomID1);
            //A[i2]->GetAtomID(AtomID2);
            //std::cout << AtomID1 << " " << AtomID2 << std::endl;
            if(!freeLigand&&(strcmp(monomer->chem_type(),"L-PEPTIDE LINKING")==0||strcmp(monomer->chem_type(),"D-PEPTIDE LINKING")==0||strcmp(monomer->chem_type(),"PEPTIDE LINKING")==0)&&!(strcmp(monomer->chem_type(),"RNA LINKING")==0)&&!(strcmp(monomer->chem_type(),"DNA LINKING")==0)){
              A[i1]->AddBond ( A[i2],ccp4srs::Bond::Single,2 );
              A[i2]->AddBond ( A[i1],ccp4srs::Bond::Single,2 );
            } else {
              A[i1]->AddBond ( A[i2],bondi->order(),2 );
              A[i2]->AddBond ( A[i1],bondi->order(),2 );
            }
          } 
        }
        if (modelBondsSame && nModels>1) {
          for (nm=firstModel+1;nm<=lastModel;nm++) {
            PPCResidue selResMod = NULL;
            mmdb::PPAtom AMod;
            int nrMod;
            int natomsMod;
            molHnds[0]->GetModel(nm)->GetResidueTable (pCh->GetChainID(), selResMod, nrMod);
            selResMod[ir]->GetAtomTable ( AMod,natomsMod );

            for (i=0;i<monomer->n_bonds();i++)  {
              ccp4srs::Bond* bondi = monomer->bond(i);
              int i1 = anmatch[bondi->atom1()];
              int i2 = anmatch[bondi->atom2()];
              if ((i1>=0) && (i2>=0) && (i1<natoms) && (i2<natoms))  {
                if(!freeLigand&&(strcmp(monomer->chem_type(),"L-PEPTIDE LINKING")==0||strcmp(monomer->chem_type(),"D-PEPTIDE LINKING")==0||strcmp(monomer->chem_type(),"PEPTIDE LINKING")==0)&&!(strcmp(monomer->chem_type(),"RNA LINKING")==0)&&!(strcmp(monomer->chem_type(),"DNA LINKING")==0)){
                  AMod[i1]->AddBond ( AMod[i2],ccp4srs::Bond::Single,2 );
                  AMod[i2]->AddBond ( AMod[i1],ccp4srs::Bond::Single,2 );
                } else {
                  AMod[i1]->AddBond ( AMod[i2],bondi->order(),2 );
                  AMod[i2]->AddBond ( AMod[i1],bondi->order(),2 );
                }
              } 
            }

          }
        }
      }

      }
      mmdb::FreeVectorMemory ( anmatch,0 );
      if(deleteMonomer) delete monomer;
      }

      // THIS SHOULD SOMEHOW BE OPTIONAL!
      //std::cout << "GetBonds, RDKitAromaticityAnalyze checks: " << selRes[ir]->name << "\n";
      //std::cout << strncmp(selRes[ir]->name,"HOH",3) << "\n";
      //std::cout << strncmp(selRes[ir]->name,"WAT",3) << "\n";
        if(CMMUTSRS::MonomerCache.find(selRes[ir]->GetResName())!=CMMUTSRS::MonomerCache.end()){
          monomer = CMMUTSRS::MonomerCache[selRes[ir]->GetResName()];
          deleteMonomer = false;
        } else {
          monomer = dynamic_cast<PCMMANManager>(molHnds[0])->GetSRS()->getMonomer ( selRes[ir]->GetResName(), structFile );
          deleteMonomer = true;
        }
      bool possLig = true;
      if(monomer){
         //std::cout << monomer->chem_type() << "\n";
         if(strcmp(monomer->chem_type(),"L-PEPTIDE LINKING")!=0&&strcmp(monomer->chem_type(),"D-PEPTIDE LINKING")!=0&&strcmp(monomer->chem_type(),"PEPTIDE LINKING")!=0&&!(strcmp(monomer->chem_type(),"RNA LINKING")==0)&&!(strcmp(monomer->chem_type(),"DNA LINKING")==0)){
           possLig = true;
         } else {
           possLig = false;
         }
         if(deleteMonomer) delete monomer;
      }
      //std::cout << freeLigand << " " << possLig << " " << int(strlen(selRes[ir]->name)<3||(strncmp(selRes[ir]->name,"HOH",3)!=0&&strncmp(selRes[ir]->name,"WAT",3)!=0)) << "\n";
      if((freeLigand||possLig)&&(strlen(selRes[ir]->name)<3||(strncmp(selRes[ir]->name,"HOH",3)!=0&&strncmp(selRes[ir]->name,"WAT",3)!=0))){
        try {
          RDKit::RWMol *mol = RDKitAromaticityAnalyze(selRes[ir],aromaticize);
        
          std::string svg = RDKitResidueToSVG(mol);
          char ResID[1024];
          selRes[ir]->GetResidueID(ResID);
          monomerSVGS[std::string(ResID)] = svg;
          delete mol;
        } catch (...) {
          std::cerr << "RDKitAromaticityAnalyze: Some unhandled exception.\n";
          std::cerr << selRes[ir]->name << " ," << natoms << " atoms\n";
        }
      }

/*
      // FIXME, what do we do about altLocs (second arg)?
      if((dynamic_cast<PCMMANManager>(molHnds[0])->GetSRS()->makeBonds(selRes[ir],NULL,NULL,NULL,false))!=ccp4srs::CCP4SRS_Ok){
        doContacts = true;
      }
*/

      if(doneCIF){
        doContacts = false;
      }
      // There was no match in the Structure database so find
      // bonds on basis of close contacts
      if ( doContacts )  {
         if (modelBondsSame && nModels > 1 ) {
           IntraResContacts ( molHnds[0], selRes[ir], 1, modelSelAtom, nSelAtom, firstModel, lastModel );
         } else {
           IntraResContacts ( molHnds[0], selRes[ir], 1);
         }
      }
      else if ((checkGraphs&&freeLigand&&!doneCIF)||haveHs)  {
         if (modelBondsSame && nModels > 1 ) {
           IntraResContacts ( molHnds[0], selRes[ir], 1, modelSelAtom, nSelAtom, firstModel, lastModel, true );
         } else {
           IntraResContacts ( molHnds[0], selRes[ir], 1, NULL, NULL, 0, 0, true);
         }
      }
    }
  } }

  //INTER-residue bonds

  // Find inter-res close contacts
  if (contacts) delete [] contacts;
  contacts = NULL;
  nContacts = 0;
  molHnds[0]->SeekContacts(selAtom0,nAtoms,selAtom0,nAtoms,
		 0.0,params->interResCut,1,contacts,nContacts,0,NULL,0);
  //cout << "FindBonds INTER-residue bonds nContacts " << nAtoms << " " << nContacts << endl;

  // Check if each contact is between possible altLoc matches
  // and that this corresponds to recognised chemical link
  // if ( contacts[ic].id1 < contacts[ic].id2 ) {
  if ( contacts && nContacts > 0 ) { 
    for ( ic = 0; ic < nContacts; ic++) {
      if ( contacts[ic].id1 < contacts[ic].id2 ) {
        pa1 = selAtom0[contacts[ic].id1];
        pa2 = selAtom0[contacts[ic].id2];
        if(pa1->GetSeqNum()!=pa2->GetSeqNum()){
        if (  nAltTot <= 0 || molHnds[0]->doAltLocMatch( pa1, pa2 )  ) {
          //cout << "testing contact " << pa1->residue->seqNum << " " << pa2->residue->seqNum << endl;
          restype1 =  dynamic_cast<PCMMANManager>(molHnds[0])->GetRestypeCode (pa1->residue);
          restype2 = dynamic_cast<PCMMANManager>( molHnds[0])->GetRestypeCode (pa2->residue);
          //cout << "testing interres " <<  pa1->residue->seqNum << pa1->residue->name << " " << restype1 << " "  << pa2->residue->seqNum <<  pa2->residue->name  << " " << restype2 << endl;
          if (restype1 ==  RESTYPE_PEPTIDE && restype2 == RESTYPE_PEPTIDE ) {
            if (isInterResBond(pa1,pa2)) {
              AddConnection (contacts[ic].id1,contacts[ic].id2,selAtom0); 
            }else{
              double dist = contacts[ic].dist;
              if(CMMUTSRS::CheckCovalentDistance(pa1->element,pa2->element,dist)==1)
                 AddConnection (contacts[ic].id1,contacts[ic].id2,selAtom0); 
            }
          } else if ( (restype1==RESTYPE_NUCL || restype1==RESTYPE_RNA || restype1==RESTYPE_DNA) && 
                 (restype2 == RESTYPE_NUCL || restype2 == RESTYPE_RNA || restype2 == RESTYPE_DNA) ) {
 
            if (isInterResBond(pa1,pa2)){ 
               AddConnection (contacts[ic].id1,contacts[ic].id2,selAtom0); 
            } else {
               double dist = contacts[ic].dist;
               if(CMMUTSRS::CheckCovalentDistance(pa1->element,pa2->element,dist)==1){
                 AddConnection (contacts[ic].id1,contacts[ic].id2,selAtom0); 
               }
            }
          } else if ( (restype1==RESTYPE_SACH || restype1==RESTYPE_DSACH || restype1==RESTYPE_LSACH) && 
                 (restype2 == RESTYPE_SACH || restype2 == RESTYPE_DSACH || restype2 == RESTYPE_LSACH) ) {
 
            if (isInterResBond(pa1,pa2)) {
               AddConnection (contacts[ic].id1,contacts[ic].id2,selAtom0); 
            } else {
               double dist = contacts[ic].dist;
               if(CMMUTSRS::CheckCovalentDistance(pa1->element,pa2->element,dist)==1){
                 AddConnection (contacts[ic].id1,contacts[ic].id2,selAtom0); 
               }
            }

	  } else if ((restype1<RESTYPE_SOLVENT || restype1>RESTYPE_NONPOLY) &&
               (restype2<RESTYPE_SOLVENT || restype2>RESTYPE_NONPOLY ) ) {
            if (ltBondDistance(pa1,pa2,contacts[ic].dist) )  {
	      //cout << "non-peptide" << pa1->name << " " << pa2->name << endl;
              if((isMetal(pa1->name)&&!isMetal(pa2->name))||(isMetal(pa2->name)&&!isMetal(pa1->name))){
	        //cout << "metal coordination, ignoring" << endl;
              } else {
                AddConnection (contacts[ic].id1,contacts[ic].id2,selAtom0);
              }
            }
          } else {
            if (isInterResBond(pa1,pa2)) {
               AddConnection (contacts[ic].id1,contacts[ic].id2,selAtom0); 
            } else {
               double dist = contacts[ic].dist;
               if(CMMUTSRS::CheckCovalentDistance(pa1->element,pa2->element,dist)==1){
                 //std::cout << "Adding possible covalent contact\n";
                 //std::cout << pa1->element << " " <<  pa2->element << " " << dist << "\n";
                 AddConnection (contacts[ic].id1,contacts[ic].id2,selAtom0); 
               }
            }
          }
        }
        } // SeqNum check.
      }
    }
    if (modelBondsSame && nModels > 1 ) {
      //int nb;
      for (nm = firstModel+1; nm <= lastModel; nm++) {
      if (molHnds[0]->GetModel(nm) != NULL) {
        //nb = 0;
        molHnds[0]->SelectAtoms(selH,nm,"*", ANY_RES,"*", ANY_RES,"*","*","*","*","*",SKEY_NEW);
        molHnds[0]->SelectAtoms(selH,nm,"*", ANY_RES,"*", ANY_RES,"*","*","*","H","*",SKEY_CLR);
        selAtom0 = NULL;
        molHnds[0]->GetSelIndex ( selH,selAtom0,nAtoms );
        //cout << "selAtom0 " << nm << " nAtoms " << nAtoms << " nContacts " << nContacts << endl;
        for ( ic = 0; ic < nContacts; ic++) {
          if ( contacts[ic].id1 < contacts[ic].id2 && contacts[ic].id2 < nAtoms && contacts[ic].id2 > 0 ) {
            pa1 = selAtom0[contacts[ic].id1];
            pa2 = selAtom0[contacts[ic].id2];
            if ( nAltTot <= 0 || molHnds[0]->doAltLocMatch( pa1, pa2 )) { 
              restype1 =  dynamic_cast<PCMMANManager>(molHnds[0])->GetRestypeCode (pa1->residue);
              restype2 =  dynamic_cast<PCMMANManager>(molHnds[0])->GetRestypeCode (pa2->residue);
	      if (restype1 ==  RESTYPE_PEPTIDE && 
                      restype2 == RESTYPE_PEPTIDE) { 
                if (isInterResBond(pa1,pa2)) {
                  //nb++;
                  AddConnection (contacts[ic].id1,contacts[ic].id2,selAtom0);
                }
              } else if (ltBondDistance(pa1,pa2,contacts[ic].dist)) {
		//nb++; 
                  AddConnection (contacts[ic].id1,contacts[ic].id2,selAtom0);
	      }
            }
	  }
        }
        //cout << "nm " << nm << " nb " << nb << endl; 
      } }
    }                     
  }

  }  // End of loop over models

  molHnds[0]->DeleteSelection(selH);
  if (nMatchAtom) FreeVectorMemory(nMatchAtom, 0);
  nMatchAtom = NULL;
  if (matchAtom) FreeMatrixMemory(matchAtom,maxMatch,0,0);
  matchAtom = NULL;
  if (contacts) delete [] contacts;
  if (selAtom) delete [] selAtom;
  /*
  for (nm=2;nm<=nModels;nm++) {
    if (modelSelAtom[nm]) delete [] modelSelAtom[nm];
  }
  */
  //cout << "done FindBonds" << endl;
  if(modelSelAtom) delete [] modelSelAtom;
  if(nSelAtom) delete [] nSelAtom;


  delete structFile;
  return output.str();
}

std::string CMolBonds::GetMonomerSVG(const std::string &atomID){
  if(monomerSVGS.count(atomID)){
    return monomerSVGS[atomID];
  }
  return std::string("");
}

//-----------------------------------------------------------
void CMolBonds::AddConnection ( PCAtom pa1 , PCAtom pa2) {
//-----------------------------------------------------------
  
  //cout << "AddConnection " <<  GetMolHnd(0)->AtomLabel_atom(pa1) << " " 
  //	<< GetMolHnd(0)->AtomLabel_atom(pa2) << endl;

  if ( fabs(pa1->x - pa2->x) < 0.0001 &&
       fabs(pa1->y - pa2->y) < 0.0001 &&
       fabs(pa1->z - pa2->z)  < 0.0001 ) {
    cout << "Bonded atoms superposed " <<
        GetMolHnd(0)->AtomLabel_atom(pa1) << " " <<
        GetMolHnd(0)->AtomLabel_atom(pa2) << endl;
    return;
  }
      
   if ( pa1->GetNBonds() == 0 )
     pa1->AddBond(pa2,1,3);
   else 
     pa1->AddBond(pa2,1,1);
   if ( pa2->GetNBonds() == 0 )
     pa2->AddBond(pa1,1,3);
   else 
     pa2->AddBond(pa1,1,1);

}


//-----------------------------------------------------------
int CMolBonds::DeleteConnection ( PCAtom pa1 , PCAtom pa2) {
//-----------------------------------------------------------
// Delete the connection between two atom - need to zap
// from the bond list of both atoms

  mmdb::AtomBond *AtomBond;
  int nAtomBonds;
  PCAtom atoms[10];
  byte order[10];
  PCAtom pA,pB;

  //char AtomID1[30];
  //char AtomID2[30];
  
  int idel,j;

  for (int nn=0;nn<=1;nn++) {
    if (nn == 0) {
      pA = pa1;
      pB = pa2;
    } else {
      pA = pa2;
      pB = pa1;
    }
    idel = -1;
    j = 0;

    pA->GetBonds ( AtomBond, nAtomBonds);
    if(nAtomBonds>0) {
      //pA->GetAtomID ( AtomID1 );
      //pB->GetAtomID ( AtomID2 );
      //cout << AtomID1 << " " << AtomID2 << " nAtomBonds " << nAtomBonds;
      for (int i=0;i<nAtomBonds;i++) {
	//AtomBond[i].atom->GetAtomID(AtomID2);
        //cout << " " << AtomID2;
        if (AtomBond[i].atom == pB) {
          idel = i;
        } else {
          atoms[j] = AtomBond[i].atom;
          order[j] = AtomBond[i].order;
          j++;
        }
      }
      //cout << endl << "idel " << idel << endl;
      if (idel>=0) {
        pA->FreeBonds();
        for (int i=0;i<nAtomBonds-1;i++) {
          //cout << "adding " << i << " " <<  pA->GetNBonds() << " " <<  atoms[i]->name << endl;
          if ( pA->GetNBonds() == 0 )
            pA->AddBond( atoms[i],order[i],3);
          else 
            pA->AddBond( atoms[i],order[i],1);
        }
      }
    }
  }
  if (idel>=0)
    return 0;
  else
    return 1;
 
}

//-----------------------------------------------------------
void CMolBonds::AddConnection (
    int ia1 , int ia2, PPCAtom selAtom, int offset ) {
//-----------------------------------------------------------
  //cout << "ia1,ia2,offset " << ia1 << " " << ia2 << " " << offset << endl;
  //printf ( "Bond %s %s %i %s %s %i\n",selAtom[ia1]->residue->name,
  //selAtom[ia1]->name,selAtom[ia1]->serNum,selAtom[ia2]->residue->name,selAtom[ia2]->name),selAtom[ia2]->serNum;
  int ja1,ja2;
  ja1 = ia1 + offset;
  ja2 = ia2 + offset;
  //cout << "ja1,ja2 " << ja1 << " " << ja2 << endl;
  if ( fabs(selAtom[ja1]->x - selAtom[ja2]->x) < 0.0001 &&
       fabs(selAtom[ja1]->y - selAtom[ja2]->y) < 0.0001 &&
       fabs(selAtom[ja1]->z - selAtom[ja2]->z)  < 0.0001 ) {
    //cout << "Bonded atoms superposed " <<
    //  GetMolHnd(0)->AtomLabel_atom(selAtom[ja1]) <<
    // " " <<  GetMolHnd(0)->AtomLabel_atom(selAtom[ja2]) << endl;
    return;
  }
  int rv1,rv2;

   if ( selAtom[ja1]->GetNBonds() == 0 ) {
     rv1 = selAtom[ja1]->AddBond(selAtom[ja2],1,3);
   }
   else {
     rv1 = selAtom[ja1]->AddBond(selAtom[ja2],1,1);
   }
   if ( selAtom[ja2]->GetNBonds() == 0 ) {
     rv2 = selAtom[ja2]->AddBond(selAtom[ja1],1,3);
   }
   else {
     rv2 = selAtom[ja2]->AddBond(selAtom[ja1],1,1);
   }
   //cout << "rv " << rv1 << " " << rv2 << endl;
}
//-----------------------------------------------------------
void CMolBonds::AddConnection ( int ia1 , int ia2,
                    PPCAtom selAtom1, PPCAtom selAtom2,
                    int offset1 , int offset2 ) {
//-----------------------------------------------------------
  int ja1 = ia1 + offset1;
  int ja2 = ia2 + offset2;
  //printf ( "Bond %s %s %s %s\n",selAtom1[ja1]->residue->name,
  // selAtom1[ja1]->name,selAtom2[ja2]->residue->name,selAtom2[ja2]->name);

  if ( fabs(selAtom1[ja1]->x - selAtom2[ja2]->x) < 0.0001 &&
       fabs(selAtom1[ja1]->y - selAtom2[ja2]->y) < 0.0001 &&
       fabs(selAtom1[ja1]->z - selAtom2[ja2]->z)  < 0.0001 ) {
    //cout << "Bonded atoms superposed " <<
    //  GetMolHnd(0)->AtomLabel_atom(selAtom1[ja1]) <<
    // " " <<  GetMolHnd(0)->AtomLabel_atom(selAtom2[ja2]) << endl;
    return;
  }


   if ( selAtom1[ja1]->GetNBonds() == 0 ) {
     selAtom1[ja1]->AddBond(selAtom2[ja2],1,3);
   }
   else {
     selAtom1[ja1]->AddBond(selAtom2[ja2],1,1);
   }
   if ( selAtom2[ja2]->GetNBonds() == 0 ) {
     selAtom2[ja2]->AddBond(selAtom1[ja1],1,3);
   }
   else {
     selAtom2[ja2]->AddBond(selAtom1[ja1],1,1);
   }
}

//---------------------------------------------------------------------
int CMolBonds::IntraResContacts ( mmdb::Manager *molHnd, PCResidue pRes, int nAlt, 
     PPCAtom modelSelAtom[], int nSelAtom[], int firstModel, int lastModel, bool HsOnly ) {
//---------------------------------------------------------------------
  PPCAtom pAtom1=0;
  PPCAtom pAtom2=0;
  int nAtom1,nAtom2,nContacts,ic;
  PSContact contacts = NULL;
  PCAtom pa1,pa2;

  pAtom1 = NULL;
  pAtom2 = NULL;

  int selHndH,selHndNoH;

  if(HsOnly) {
    //std::cout << "CMolBonds::IntraResContacts Now think about how to add Hs\n";
    selHndH =  molHnd->NewSelection();
    selHndNoH =  molHnd->NewSelection();
    molHnd->SelectResidue(selHndH,pRes,mmdb::STYPE_ATOM,mmdb::SKEY_NEW,true);
    molHnd->SelectResidue(selHndNoH,pRes,mmdb::STYPE_ATOM,mmdb::SKEY_NEW,true);
    molHnd->SelectAtoms(selHndH,0,"*", ANY_RES,"*", ANY_RES,"*","*","*","H","*",SKEY_AND);
    molHnd->SelectAtoms(selHndNoH,0,"*", ANY_RES,"*", ANY_RES,"*","*","*","H","*",SKEY_XOR);
    molHnd->GetSelIndex(selHndH,pAtom1,nAtom1);
    molHnd->GetSelIndex(selHndNoH,pAtom2,nAtom2);
    //std::cout << "CMolBonds::IntraResContacts Number of Hs " << nAtom1 << "\n";
    //std::cout << "CMolBonds::IntraResContacts Number of non-Hs " << nAtom2 << "\n";
    if(nAtom1==0||nAtom2==0){
      molHnd->DeleteSelection(selHndH);
      molHnd->DeleteSelection(selHndNoH);
      return 0;
    }
  } else {
    pRes->GetAtomTable1(pAtom1,nAtom1);
    pRes->GetAtomTable1(pAtom2,nAtom2);
  }
 
  //cout << "IntraResContacts nAtom1 " << nAtom1 << " " << firstModel << " " << lastModel << endl;

  char AtomID1[1024];
  char AtomID2[1024];

  // One atom in residue - probably water
  if ( nAtom1 == 1 ){
    if(!HsOnly){
      if(pAtom1) delete [] pAtom1;
      if(pAtom2) delete [] pAtom2;
    } else {
      molHnd->DeleteSelection(selHndH);
      molHnd->DeleteSelection(selHndNoH);
    }
    return 0;
  }

  GetMolHnd(0)->SeekContacts(pAtom1,nAtom1,pAtom2,nAtom2,
		 0.0,params->intraResCut,0,
		 contacts,nContacts,0,NULL,0);

  //printf ( "nContacts %i\n",nContacts);
  if ( contacts && nContacts > 0 ) { 
    for ( ic = 0; ic < nContacts; ic++) {
      if ( HsOnly||contacts[ic].id1 < contacts[ic].id2 ) {
        pa1 = pAtom1[contacts[ic].id1];
        pa2 = pAtom2[contacts[ic].id2];
        if(pa1&&pa2){
        pa1->GetAtomID ( AtomID1 );
        pa2->GetAtomID ( AtomID2 );
        //std::cout << "IntraResContacts " << AtomID1 << " " << AtomID2 << " " << contacts[ic].dist << std::endl;
        if ( ( nAlt <= 0 || (molHnds[0]->doAltLocMatch( pa1, pa2 )) ) 
            && (ltBondDistance(pa1,pa2,contacts[ic].dist)) ) {
          AddConnection (contacts[ic].id1,contacts[ic].id2,pAtom1,pAtom2);
          if (lastModel>0) {
            for (int nm=firstModel+1;nm<=lastModel;nm++) {
            if (GetMolHnd(0)->GetModel(nm) != NULL) {
              if (nSelAtom[nm] == nAtom1) AddConnection (contacts[ic].id1,contacts[ic].id2,modelSelAtom[nm],modelSelAtom[nm]);
            }}
          }
        } 
        }
      }
    }
  }
  
  if(HsOnly){
    if(molHnd){
      molHnd->DeleteSelection(selHndH);
      molHnd->DeleteSelection(selHndNoH);
    }
  } else {
    if(pAtom1) delete [] pAtom1;
    if(pAtom2) delete [] pAtom2;
  }
  if ( contacts ) delete [] contacts;
  return 0;
}

//---------------------------------------------------------------------
bool CMolBonds::ltBondDistance ( PCAtom pa1, PCAtom pa2, realtype dist ) {
//---------------------------------------------------------------------
  return CMMUTSRS::CheckCovalentDistance(pa1->element,pa2->element,dist)==1;
  /*
  PCLibElement la1,la2;
  realtype rmax;
  if ( (la1 = params->sbase->LibElement(pa1->element)) &&
       (la2 = params->sbase->LibElement(pa2->element)) ) {
    rmax = (la1->maxBondRad + la2->maxBondRad) * params->maxBondRadFactor;
    //printf( "%s %s dist %f rmax %f\n",pa1->name,pa2->name,dist,rmax);
    return dist < rmax;
  }
  else
    return dist < params->maxBondRad;
  */
}

//----------------------------------------------------------------------
bool CMolBonds::isInterResBond ( PCAtom pa1, PCAtom pa2 ) {
//----------------------------------------------------------------------
//FIXME - Is this really necessary. Is not checking covalent distance good enough?
//        Perhaps we'll find out.
  return false;
  /*
  int type1,type2,n;
  AtomName a1,a2;

  type1 =  dynamic_cast<PCMMANManager>(molHnds[0])->GetRestypeCode (pa1->residue);
  type2 =  dynamic_cast<PCMMANManager>(molHnds[0])->GetRestypeCode (pa2->residue);

  //cout << "CMolBonds::isInterResBond " << first << " " <<  second << " " << type1 << " " << type2 << " " << params->sbase->nLinks << endl;
  
  strcpy_css(a1,pa1->name);
  strcpy_css(a2,pa2->name);
  
  for ( n = 0; n <= params->sbase->nLinks; n++ ) {
    if (strcmp(params->sbase->link[n]->id,"symmetry") != 0 &&
        strcmp(params->sbase->link[n]->id,"gap") != 0 ) {
      //cout << "CMolBonds::isInterResBond " << n << " " << params->sbase->link[n]->id << endl;
      if ( params->sbase->link[n]->lg1.Match(type1, pa1->residue->name, a1)&&
	  params->sbase->link[n]->lg2.Match(type2, pa2->residue->name, a2)) {
	//printf ("Type of link: %i\n",n);
        return true;
      } else if ( params->sbase->link[n]->lg2.Match(type1, pa1->residue->name, a1)&&
	  params->sbase->link[n]->lg1.Match(type2, pa2->residue->name, a2)) {
	//printf ("Type of link: %i\n",n);
          return true;
      }
    }
  }
  return false;
*/
}



