/*
     mmut/mmut_bonds.h: CCP4MG Molecular Graphics Program
     Copyright (C) 2001-2008 University of York, CCLRC

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


#ifndef __CCP4MolBonds__ 
#define __CCP4MolBonds__
#include <map>
#include "mmdb2/mmdb_manager.h"
#include "mman_base.h"

class CMGSBase;

DefineClass(CMolBondParams);

class CMolBondParams {

  friend class CMolBonds;
  friend class CMMANManager;
public:
  CMolBondParams();
  ~CMolBondParams();
protected:
  mmdb::realtype interResCut;
  mmdb::realtype intraResCut;
  mmdb::realtype maxBondRad;
  mmdb::realtype maxBondRadFactor;

};

DefineClass(CMolBonds);
 
class CMolBonds : public CMMANBase {
public :
 
  // Constructors
  CMolBonds( const PCMMUTManager molHndin, CMolBondParams *params  );
  // , int selHndin, CMolBondParams *params);
 
 // Destructor
  ~CMolBonds();
  std::string FindBonds (int udd_sbaseCompoundID,
		 int udd_sbaseAtomOrdinal, int udd_atomEnergyType, const std::map<std::string,std::string> &customResCIFFiles, bool aromaticize, bool checkGraphs );
  std::string GetMonomerSVG(const std::string &atomID);
  std::map<std::string,std::string> GetMonomerSVGs(){return monomerSVGS;};

  void AddConnection (int ia1, int ia2, mmdb::Atom** selAtom1, mmdb::Atom** selAtom2,
                       int offset1 =0 , int offset2 = 0);
  void AddConnection (int ia1, int ia2, mmdb::Atom** selAtom,int offset=0);
  void AddConnection (mmdb::Atom* pa1, mmdb::Atom* pa2);
  int DeleteConnection ( mmdb::Atom* pa1 , mmdb::Atom* pa2);

  bool isInterResBond ( mmdb::Atom* p1, mmdb::Atom* p2);
  int IntraResContacts ( mmdb::Manager *molHnd, mmdb::Residue* p1, int nAlt,  mmdb::Atom** modelSelAtom[]=NULL, int nSelAtom[]=NULL, int firstModel=0, int lastModel=0, bool HsOnly=false);
  bool ltBondDistance ( mmdb::Atom* pa1, mmdb::Atom* pa2, mmdb::realtype dist);

 private:
 
  CMolBondParams *params;

  int nAtoms;
  int nRes;
  mmdb::psvector  sbaseCompoundID;
  mmdb::ivector sbaseAtomIndex;

  // Structures for holding the selected bonds
  int nB;
  mmdb::imatrix bonds;
  std::map<std::string,std::string> monomerSVGS;

};
#endif
