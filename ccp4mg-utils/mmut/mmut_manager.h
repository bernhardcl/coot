/*
     mmut/mmut_manager.h: CCP4MG Molecular Graphics Program
     Copyright (C) 2001-2008 University of York, CCLRC
     Copyright (C) 2009-2011 University of York
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


#ifndef __MMUT_Manager__
#define __MMUT_Manager__

#include <string>
#include <vector>
#include "mmdb2/mmdb_manager.h"
#include "cartesian.h"

#define PI 3.141592653589793238462643

// =======================  CMMUTManager  ===========================

class PisaAssemblyTransformation;

DefineClass(CMMUTManager);
DefineStreamFunctions(CMMUTManager) ;

class CMMUTManager : public mmdb::Manager  {

    std::vector<matrix> biomatrices;
    bool doneBiomolecule;
  public :

    CMMUTManager();
    ~CMMUTManager();

    const mmdb::Cryst& get_cell() { return cryst; } 
    mmdb::Cryst* get_cell_p() { return & cryst; } 
  
    void ListAtomInfo(int selHnd);
    void PrintAtomicComposition(int selHnd);
    void PrintResidueComposition(int selHnd);
    void PrintSequence(int selHnd);
    void PrintBValues(int selHnd);
    //void PrintLengthsAndAngles(int selHnd);
    //void PrintLengthsAndAngles(realtype min, realtype max);
 
    int NumberOfHydrogens(int selHnd);
    int ResNoLookup(mmdb::pstr resname);
    int TotalNumRes(int selHnd);
    mmdb::pstr GetSequence(int selHnd);
    mmdb::pstr GetSequenceFromResidues(int selHnd);
    void SelectOneAtomForNonAminoNucleotide(int selHnd, int atomSelHnd);
    int *AtomicComposition(int selHnd);
    int *ResidueComposition(int selHnd);
    mmdb::Residue* NextResidue( mmdb::Residue* pRes ,int increment=1);

    mmdb::realtype  BondLength(mmdb::Atom* A, mmdb::Atom* B);
    mmdb::realtype  BondAngle(mmdb::Atom* A, mmdb::Atom* B, mmdb::Atom* C);
    mmdb::realtype  TorsionAngle(mmdb::Atom* A, mmdb::Atom* B, mmdb::Atom* C, mmdb::Atom* D);
    mmdb::realtype  MolWeight(int selHnd);
    mmdb::realtype  MolWeightWithH(int selHnd);
    mmdb::realtype *GetBValues(int selHnd);
    std::vector<double> GetBValuesDoubleVector(int selHnd);
    mmdb::realtype *CentreOfMass(int selHnd);
    Cartesian CentreOfMassAsCartesian(int selHnd);
    mmdb::realtype  Mass(int selHnd);
    mmdb::realtype *CentreOfCoordinates(int selHnd);
    Cartesian CentreOfCoordinatesAsCartesian(int selHnd);
    int  NumberOfAtoms(int selHnd);
    std::vector<double> Extent(int selHnd);
    mmdb::realtype ExtentSize(int selHnd);
    std::vector<Cartesian> GetPrincipalComponents(int selHnd);

    bool isMainChain(mmdb::Atom* p_atom);
    bool doAltLocMatch ( mmdb::Atom* pa1, mmdb::Atom* pa2 ); 
    int NameComparison ( const char *name , int ntypes , const char *types[] );
    std::string  TrimString(mmdb::pstr inp);
    std::string AtomLabel(mmdb::Atom* p_atom, int mask[]);
    bool ChainIDisDigit(mmdb::Chain* p_ch);
    const char* AtomLabel_atom1(mmdb::Atom* p_atom);
    const char* AtomLabel_atom(mmdb::Atom* p_atom);
    const char* AtomLabel_residue(mmdb::Atom* p_atom);
    const char* AtomLabel_chain(mmdb::Atom* p_atom);
    const char* AtomLabel_residue1(mmdb::Residue* p_res);
    const char* AtomLabel_mask(mmdb::Atom* p_atom, int mask[]);
    //int ApplyTransform(int selHnd,double rotmat[],double transv[]);

    //Editor
    int WriteSelection (int selHnd,char *file, const char *format="PDB");
    int PutSelectedAtoms (int selHnd , mmdb::Manager* mmdb2);
    int CopySelection (int selHnd,mmdb::Manager* mmdb2);
    int FindCloseAtomPairs ( int selHnd, double min_distance, 
			     double max_distance);
   int FixElementNames();

   std::string Source();
   std::string Unknowns();
   std::string GetRemarksString();
   std::string SiteInfo();
   double Resolution();
   std::string StructureTitle();

   int ApplyPDBSecStructure(int model);

   std::vector<double> GetCellInfo();
   std::string MMUTGetSpaceGroup();

   std::string SelectionToSCOP(int selHnd);

   std::vector<matrix> GetBiomoleculeAsMatrices(int nBiomol,int nModel=1);
   
   static mmdb::Manager* GetCAModel(mmdb::Manager *molHnd);
   int GenerateTransformedChain(const char *chainID, mmdb::realtype *vmat, mmdb::Manager *molHnd2);

   int RemoveSmallHelices(int model,int minHelices=4);

   void  SelectAminoNotHet (
             int   selHnd,   // must be obtained from NewSelection()
             mmdb::SELECTION_TYPE selType,  // selection type STYPE_XXXXX
             int   iModel,   // model number; iModel=0 means
                             // 'any model'
             mmdb::cpstr Chains,   // may be several chains "A,B,W"; "*"
                             // means 'any chain' (in selected
                             // model(s))
             int   ResNo1,   // starting residue sequence number
             mmdb::cpstr Ins1,     // starting residue insertion code; "*"
                             // means 'any code'
             int   ResNo2,   // ending residue sequence number.
             mmdb::cpstr Ins2,     // ending residue insertion code; "*"
                             // means 'any code'. Combination of
                             // ResNo1=ResNo2=ANY_RES and
                             // Ins1=Ins2="*" means 'any residue'
                             // (in selected chain(s))
             mmdb::cpstr RNames,   // may be several residue names
                             // "ALA,GLU,CIS"; "*" means
                             // 'any residue name'
             mmdb::cpstr ANames,   // may be several names "CA,CB"; "*"
                             // means 'any atom' (in selected
                             // residue(s))
             mmdb::cpstr Elements, // may be several element types
                             // 'H,C,O,CU'; "*" means 'any element'
             mmdb::cpstr altLocs,  // may be several alternative
                             // locations 'A,B'; "*" means
                             // 'any alternative location'
             mmdb::SELECTION_KEY selKey=mmdb::SKEY_OR // selection key
           );
   bool isNTerminusBound(mmdb::Residue* res);
   bool isCTerminusBound(mmdb::Residue* res);
   bool isPeptideBound(mmdb::Residue* res);

   int GetNonTerAllSelHnd();

   int CopySelectedAtomsToChain(int selHnd, mmdb::Chain* newChain);

   private: 

     //Analysis
    int *iatom_types;
    int *iatom_type_lookup;

 
};

class PisaAssemblyTransformation {
  std::string selection;
  std::string visual_id;
  matrix transformation;
 public:
  PisaAssemblyTransformation(){};
  PisaAssemblyTransformation(const std::string &selection_in, const std::string &visual_id_in, const matrix& transformation_in){
    selection = selection_in;
    visual_id = visual_id_in;
    transformation = transformation_in;
  };
  std::string Selection() const { return selection; };
  std::string VisualID() const { return visual_id; };
  matrix Transformation() const { return transformation; };
};

#endif
