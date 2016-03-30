/*
     mmut/mmut_basepairs.cc: CCP4MG Molecular Graphics Program
     Copyright (C) 2001-2008 University of York, CCLRC
     Copyright (C) 2009-2010 University of York

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

#include <string.h>
#include <vector>
#include <map>
#include <utility>
#include <stdlib.h>
#include <math.h>
#include "mman_manager.h"
#include "mmut_basepairs.h"
#include "mmut_hbond.h"
#include "plane.h"

CNABasePairs::CNABasePairs(CMMANManager *molHnd, int selHnd, mmdb::Atom** selAtoms, int nSelAtoms, const AtomColourVector &atom_colour_vector){
  Calculate(molHnd,selHnd,selAtoms,nSelAtoms,atom_colour_vector);
}

void CNABasePairs::Calculate(CMMANManager *molHnd, int selHnd, mmdb::Atom** selAtoms, int nSelAtoms, const AtomColourVector &atom_colour_vector, double *hb_params_array){
  
  int C5sel = molHnd->NewSelection();
  
  molHnd->Select(C5sel,mmdb::STYPE_ATOM,0,"*",mmdb::ANY_RES,"*",mmdb::ANY_RES,"*","*","C5*","C","*",mmdb::SKEY_NEW);
  molHnd->Select(C5sel,mmdb::STYPE_ATOM,0,"*",mmdb::ANY_RES,"*",mmdb::ANY_RES,"*","*","C5'","C","*",mmdb::SKEY_OR);

  CHBond hbonds(molHnd, selHnd);
  if(hb_params_array) hbonds.SetParams(8,hb_params_array);
  hbonds.Calculate();

  std::vector<std::pair<mmdb::Residue*,Plane> > plane_res_pairs;
  std::vector<double*> plane_res_colours;

  //std::cout << "In base pairs calc " << nSelAtoms << " atoms\n"; std::cout.flush();
  //std::cout << "In base pairs calc " << hbonds.hbonds.GetNofConnections() << " hbonds\n"; std::cout.flush();
  
  char ID[50];
  char ID1[50];
  char ID2[50];
  char AtomID1[50];
  char AtomID2[50];
  std::map<std::string,std::map<std::string,int> > alreadyDone;
  for(int i=0;i<nSelAtoms;i++){
    if(selAtoms[i]->isInSelection(C5sel)){
      mmdb::Residue* res = selAtoms[i]->GetResidue();
      if(res){
        if(res->isNucleotide()||res->isDNARNA()||molHnd->isDNARNA(res)){
          //std::cout << "Considering atom " << selAtoms[i]->name << " in base " << selAtoms[i]->GetChainID() << "/" << selAtoms[i]->GetResidueNo() << "\n";
	  for(size_t j=0;j<hbonds.hbonds.GetNofConnections();j++){
	     const char* thisID = res->GetResidueID(ID);
             mmdb::Atom* atom1 = hbonds.hbonds.pAtom1[j];
             mmdb::Atom* atom2 = hbonds.hbonds.pAtom2[j];
             mmdb::Residue* res1 = atom1->GetResidue();
             mmdb::Residue* res2 = atom2->GetResidue();
	     const char* id1 = res1->GetResidueID(ID1);
	     const char* id2 = res2->GetResidueID(ID2);
	     atom1->GetAtomID(AtomID1);
	     atom2->GetAtomID(AtomID2);
             //std::cout << "   Testing atom " << atom1->name << " in base " << atom1->GetChainID() << "/" << atom1->GetResidueNo() << "\n";
             //std::cout << "      with atom " << atom2->name << " in base " << atom2->GetChainID() << "/" << atom2->GetResidueNo() << "\n";
             //std::cout << "      strncmp(thisID,id1,50): " << strncmp(thisID,id1,50) << "\n";
             //std::cout << "      strncmp(thisID,id2,50): " << strncmp(thisID,id2,50) << "\n";
             //std::cout << "      strncmp(id1,id2,50): " << strncmp(id1,id2,50) << "\n";
             //std::cout << "      res1->GetSeqNum()-res2->GetSeqNum(): " << res1->GetSeqNum()-res2->GetSeqNum() << "\n";
             //std::cout << "      res1->isNucleotide(): " << res1->isNucleotide() << " " << res1->name << "\n";
             //std::cout << "      res1->isDNARNA(): " << res1->isDNARNA() << "\n";
             //std::cout << "      molHnd->isDNARNA(res1): " << molHnd->isDNARNA(res1) << "\n";
             //std::cout << "      molHnd->isDNARNA(res2): " << molHnd->isDNARNA(res2) << "\n";
             bool adjacent = (abs(res1->GetSeqNum()-res2->GetSeqNum())==1)&&(strncmp(res1->GetChainID(),res2->GetChainID(),10)==0);
             
	     if((strncmp(thisID,id2,50)==0&&strncmp(id1,id2,50)!=0&&(!adjacent)&&(res1->isNucleotide()||res1->isDNARNA()||molHnd->isDNARNA(res1)))||(strncmp(thisID,id1,50)==0&&strncmp(id1,id2,50)!=0&&(!adjacent)&&(res2->isNucleotide()||res2->isDNARNA()||molHnd->isDNARNA(res2)))){
                     //std::cout << "Passed first lot of tests\n";
                     //std::cout << "molHnd->BondLength(atom1,atom2) " << molHnd->BondLength(atom1,atom2) << "\n";
                     //std::cout << "atom1->name " << atom1->name << "\n";
                     //std::cout << "atom2->name " << atom2->name << "\n";
		     if(
			   /*
			   ((
			   std::string(res2->name)==std::string("C")||
			   std::string(res2->name)==std::string("DC")||
			   std::string(res2->name)==std::string("Cd")||
			   std::string(res2->name)==std::string("CYT")
			   )&&
			   (
			   std::string(res1->name)==std::string("T")||
			   std::string(res1->name)==std::string("DT")||
			   std::string(res1->name)==std::string("Td")||
			   std::string(res1->name)==std::string("THY")
			   ))||
			   ((
			   std::string(res2->name)==std::string("T")||
			   std::string(res2->name)==std::string("DT")||
			   std::string(res2->name)==std::string("Td")||
			   std::string(res2->name)==std::string("THY")
			   )&&
			   (
			   std::string(res1->name)==std::string("C")||
			   std::string(res1->name)==std::string("DC")||
			   std::string(res1->name)==std::string("Cd")||
			   std::string(res1->name)==std::string("CYT")
			   ))||
			   ((
			   std::string(res2->name)==std::string("A")||
			   std::string(res2->name)==std::string("DA")||
			   std::string(res2->name)==std::string("Ad")||
			   std::string(res2->name)==std::string("ADE")
			   )&&
			   (
			   std::string(res1->name)==std::string("G")||
			   std::string(res1->name)==std::string("DG")||
			   std::string(res1->name)==std::string("Gd")||
			   std::string(res1->name)==std::string("GUA")
			   ))||
			   ((
			   std::string(res2->name)==std::string("G")||
			   std::string(res2->name)==std::string("DG")||
			   std::string(res2->name)==std::string("Gd")||
			   std::string(res2->name)==std::string("GUA")
			   )&&
			   (
			   std::string(res1->name)==std::string("A")||
			   std::string(res1->name)==std::string("DA")||
			   std::string(res1->name)==std::string("Ad")||
			   std::string(res1->name)==std::string("ADE")
			   ))||
                           */
			   (
			    molHnd->BondLength(atom1,atom2)>3.3
			   )||
			   (
			    std::string(atom1->name).find('\'')!=std::string::npos
			   )||
			   (
			    std::string(atom1->name).find('*')!=std::string::npos
			   )||
			   (
			    std::string(atom2->name).find('\'')!=std::string::npos
			   )||
			   (
			    std::string(atom2->name).find('*')!=std::string::npos
			   )){
                              continue;
                 }
                 if(strncmp(thisID,id2,50)==0){
                   if(alreadyDone[std::string(thisID)].find(std::string(id1))==alreadyDone[thisID].end()){
                      alreadyDone[std::string(thisID)][std::string(id1)] = 1;
                   } else {
                      //std::cout << "Rejecting duplicate " << ID << " " << id1 << ", with length: " << molHnd->BondLength(atom1,atom2) << "\n"; std::cout.flush();
                      continue;
                   }
                 }
                 if(strncmp(thisID,id1,50)==0){
                   if(alreadyDone[thisID].find(std::string(id2))==alreadyDone[thisID].end()){
                      alreadyDone[thisID][id2] = 1;
                   } else {
                      //std::cout << "Rejecting duplicate " << ID << " " << id2 << ", with length: " << molHnd->BondLength(atom1,atom2) << "\n"; std::cout.flush();
                      continue;
                   }
                 }
		 const double *col = atom_colour_vector.GetRGB(i); 
                 //std::cout << "Accepting base pair " << AtomID1 << " " << AtomID2 << ", with length: " << molHnd->BondLength(atom1,atom2) << "\n"; std::cout.flush();
	         if(strncmp(thisID,id2,50)==0&&strncmp(id1,id2,50)!=0&&(res1->isNucleotide()||res1->isDNARNA()||molHnd->isDNARNA(res1))){
		     base_pairs.push_back(std::pair<mmdb::Residue*,mmdb::Residue*>(res,res1));
		     colours.push_back(std::pair<const double*,const double*>(col,col));
	         }
	         if(strncmp(thisID,id1,50)==0&&strncmp(id1,id2,50)!=0&&(res2->isNucleotide()||res2->isDNARNA()||molHnd->isDNARNA(res2))){
		     base_pairs.push_back(std::pair<mmdb::Residue*,mmdb::Residue*>(res,res2));
		     colours.push_back(std::pair<const double*,const double*>(col,col));
	         }
	     }
	  }
	  //std::cout.flush();
        }
      }
    }
  }

  //std::cout << "BP size, BP colours size " << base_pairs.size() << " " << colours.size() << "\n"; std::cout.flush();

}

const mmdb::Residue* CNABasePairs::GetPairedResidue(const mmdb::Residue* res_in) const {
  return res_in;
}

int CNABasePairs::GetPairedResidueIndex(const int res_in) const {
  return res_in;
}
