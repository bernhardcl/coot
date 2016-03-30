/*
     mmut/mman_manager.h: CCP4MG Molecular Graphics Program
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


#ifndef __MMAN_Manager__
#define __MMAN_Manager__

#include <string>
#include <vector>
#include <map>
#include "ccp4srs/ccp4srs_manager.h"
#include "cartesian.h"
#include "matrix.h"
#include "mmut_manager.h"
#include "mmut_sasarea.h"
#include "mmut_srs.h"
#include "mmut_bonds.h"
//FIXME SRS
//#include "mmut_bonds.h"

enum enum_atomTypeData { HBTYPE, VDWRADIUS, VDWHRADIUS, IONRADIUS };
enum enum_property { PROPERTY_B, PROPERTY_OCC, PROPERTY_CHARGE,
                     PROPERTY_X, PROPERTY_Y, PROPERTY_Z,
                     PROPERTY_SEC, PROPERTY_ATOM_SAS, PROPERTY_RES_SAS,
                     PROPERTY_ATOM_CONTACT, PROPERTY_RES_CONTACT,
                     PROPERTY_SERIAL }; 

enum enum_editMode { MMAN_COORDINATES, MMAN_XTLDATA };

enum SRS_coord_types { MMUT_SRS_MLIB, MMUT_SRS_IDEAL, MMUT_SRS_RCSB, MMUT_SRS_DEFAULT };

// =======================  CMMANManager  ===========================

DefineClass(CMMANManager);
DefineStreamFunctions(CMMANManager) ;

// Function to return the parent PCMMANManager of an atom
CMMANManager* GetMMANManager(mmdb::Atom* pAtom);
std::string GetMMANManagerAddress(mmdb::Atom* pAtom);

// Get entry from ccp4srs
mmdb::Manager *GetPDBFromSRSEntryName(ccp4srs::Manager *SRS, const char *entry_name);
  
class CMMANManager : public CMMUTManager  {
    mmdb::realtype MetalCoordinationDistance[92];
 
  ccp4srs::Manager* SRS;
  int udd_vdwRadius;
  int udd_vdwHRadius;
  int udd_ionRadius;
  int udd_hbType;
  int udd_atomEnergyType;
  std::map<std::string,MMUTSRSAtomInfo> srsAtoms;
  std::map<std::string,std::string> customResCIFFiles;
  public :

    CMMANManager();
    CMMANManager(mmdb::Manager*);
//FIXME SRS
    //CMMANManager(PCMGSBase p_sbase_in, PCMolBondParams p_bond_params_in);
    CMMANManager(ccp4srs::Manager* SRS_in,PCMolBondParams p_bond_params_in);
    ~CMMANManager();
    std::string GetMonomerSVG(const std::string &atomID){return p_bonds->GetMonomerSVG(atomID);};
    std::map<std::string,std::string> GetMonomerSVGs(){return p_bonds->GetMonomerSVGs();};

    void setCustomResCIFFiles(const std::map<std::string,std::string> &customResCIFFiles_in){customResCIFFiles = customResCIFFiles_in;};
    ccp4srs::Manager *GetSRS();
    const std::map<std::string,MMUTSRSAtomInfo> &GetSrsAtoms() {return srsAtoms;};
    // SAS
    CSASArea *GetSASArea(int selHndin=-1);

    //Bonds
    std::string GetMolBonds (bool aromaticize=true, bool checkGraphs=true);
    int EditBonds (int mode, mmdb::Atom* p_atom1, mmdb::Atom* p_atom2);

    std::vector<double> GetAtomRadii ( int selHnd, int type, double scale );
    mmdb::pstr GetAtomEnergyType(mmdb::Atom* p_atom);

    mmdb::realtype GetMetalCoordinationDistance(mmdb::Atom* p_atom);
    mmdb::realtype GetAtomVDWRadius(mmdb::Atom* p_atom);
    mmdb::realtype GetAtomIonRadius(mmdb::Atom* p_atom);
    int GetAtomHBondType1(mmdb::Atom* p_atom);
    int LoadCharge(std::string loadfrom);
    std::string PrintCharges(void);

    void SetLabelMask(int i, int value);
    std::string AtomLabel(mmdb::Atom* p_atom);
    std::string AtomLabel(mmdb::Atom* p_atom,int mask[]);
    void ListBonds(int selHnd,int natoms,mmdb::Atom** selAtom);
    std::string ListSecStructure (int mask_in[], mmdb::Atom* pAtom=NULL );
    int TestBonding ( mmdb::Atom* patom1, mmdb::Atom* patom2, int max=5);
    int RestoreData (mmdb::Manager* restore_molHnd, int mode=MMAN_COORDINATES);
    int LoadUDDData( const int property=PROPERTY_B );

    int CopyModel(int model);
    int GenerateSymmetryModel(int model,int nsym,int i,int j,int k);
    int GenerateTransformedModel(int model,mmdb::realtype *vmat);
    int ApplyTransformtoModel(int model,mmdb::realtype *vmat,bool undo=false);
    std::string GetSymOpTitle(int nsym,int i,int j,int k);
    int ApplySymmetrytoModel(int model,int nsym,int i,int j,int k,bool undo=false);
    int IfSymmetryNeighbours(int selHnd, int model, int nsym, 
			     int i, int j, int k, double dist );
    int MoveFragment(int nMove, mmdb::Atom** moveAtoms, Cartesian dxyz); 
    int SelectChainTermini( void );
    int SelectSSETermini( int selHnd=-1 );

    bool isAminoacid (mmdb::Residue* pres);
    bool isDNARNA (mmdb::Residue* pres);
    int GetRestypeCode ( mmdb::Residue* pres);
//FIXME SRS
    //std::map<std::string,PCSBStructure> monlib; 

    int SetCustomRestype ( const std::string &resname , const int &restype , bool clear=false);
    int SetCustomResSynonym ( const std::string &resname , const std::string &alias , bool clear=false);


    int ExcludeOverlappedAtoms ( const int selHnd ,const mmdb::realtype cutoff, int theModel=0 );
    int SetTransform ( const std::vector<float>& transf , const std::vector<float>& transl , const bool reset);
    int SetTransform ( const mmdb::realtype rot ,const std::vector<float>& axis, const int selHnd = -1 );
    int SetTransform ( const matrix tMat, const bool reset);
    int SetTransform ( mmdb::mat44 &TMatrix , const bool reset);
    int SetTransform ( const std::vector<float>&, const bool reset);
    int ReApplyTransform( const bool reset=0);
    void UnSetTransform(bool apply_inverse=1);
    std::vector<float> GetTransform();
    std::string GetTransformString();
    bool GetIsTransformed() { return isTransformed; }
    double AtomicRMSDistance( mmdb::Atom** A1, int nA, mmdb::Atom** A2);
    int TransformToSuperposeAtoms (  mmdb::Atom** A1, int nA, mmdb::Atom** A2 );
    int TransformToSuperposeCloseAtoms(PCMMANManager fxMolHnd, int fxSelHnd , mmdb::realtype central_cutoff, mmdb::realtype cutoff, int mvSuperposeHnd,int fxSuperposeHnd);
    double DeltaResidueOrientation (mmdb::Residue* pRes, mmdb::Residue* pResFx);
    int CopyCoordinates(mmdb::Manager* fromMolHnd,int fromModel=1);
    int LoadSerial(mmdb::Manager* fromMolHnd );
    int LoadSerialFromDifferentModel(mmdb::Manager* fromMolHnd , int uddSerial);
    bool GetUnremediated() { return unremediated; };
    int GetNumberOfSecStructure(int type);
    std::string PrintSecStructure (void);
    int GetLibTMatrix(mmdb::mat44 &TMatrix,int nsym,int i,int j,int k);
    int ApplyCartesiansDeltas(const std::vector<Cartesian> &dxyz, int selHnd, double scale=1.0); 
    std::string GetAddress();
    mmdb::Manager *GeneratePISAAssembly(const std::vector<PisaAssemblyTransformation>&);
 private:

    // SAS
    PCSASArea p_sas;
    int udd_atomSAS;
    int udd_resSAS;

    PCMolBonds p_bonds;
    PCMolBondParams p_bond_params;

    // Bonds
//FIXME SRS
/*
    int SetSBaseAndBondParams(PCMGSBase p_sbase_in,PCMolBondParams p_bond_params_in);
    PCMGSBase p_sbase;
    //Atom energy types
    int udd_sbaseCompoundID;
    int udd_sbaseAtomOrdinal;

*/

    std::string loaded_charge;
    int label_mask[20];

    std::map<std::string,int> customResTypes;
    std::map<std::string,std::string> customResSynonym;

    bool isTransformed,transform_com_set ;
    mmdb::mat44 current_transform;
    mmdb::realtype transform_com[3];

    bool unremediated;


};

#endif
