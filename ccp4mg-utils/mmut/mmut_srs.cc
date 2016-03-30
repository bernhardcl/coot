#include <algorithm>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string.h>
#include "mmdb2/mmdb_manager.h"
#include "mmdb2/mmdb_mattype.h"
#include "mmdb2/mmdb_mmcif_.h"
#include "mmdb2/mmdb_tables.h"
#include "mmut_srs.h"

//------------------------------------------------------------------
int CMMUTSRS::CheckCovalentDistance (mmdb::Element a1, mmdb::Element a2,mmdb::realtype d) {
//------------------------------------------------------------------
  std::string s1 = std::string(a1);
  std::string s2 = std::string(a2);
  s1.erase(std::remove(s1.begin(), s1.end(),' '), s1.end() );
  s2.erase(std::remove(s2.begin(), s2.end(),' '), s2.end() );
  bool found = false;
  for(unsigned i=0;i<covalentDistances.size();i++){
     CMGCovalentDistance cd = covalentDistances[i];
     if((s1==cd.GetFirstAtom()&&s2==cd.GetSecondAtom())||(s2==cd.GetFirstAtom()&&s1==cd.GetSecondAtom())){
       found = true;
       if(d>=cd.GetMinLength()&&d<=cd.GetMaxLength()){
         return 1;
       }
     }
  }
  if(!found) return -1; // We have to trust dictionaries, or else we need complete covalent_distances file
  return 0;
}

//------------------------------------------------------------------
int CMMUTSRS::LoadCovalentDistances(mmdb::pstr filename ) {
//------------------------------------------------------------------
  covalentDistances.clear();
  std::ifstream cdFile(filename);
  char line[1024];
  std::string a1,a2;
  mmdb::realtype mind,maxd;
  if(cdFile){
    while(cdFile.getline(line,1024)){
      std::istringstream sl(std::string(line),std::istringstream::in);
      sl >> a1;
      sl >> a2;
      sl >> mind;
      sl >> maxd;
      CMGCovalentDistance cd(a1,a2,mind,maxd);
      covalentDistances.push_back(cd);
    }
    cdFile.close();
  }
  return 0;
}

//------------------------------------------------------------------
int CMMUTSRS::AddCovalentDistance(mmdb::Element a1, mmdb::Element a2,mmdb::realtype mind,mmdb::realtype maxd) {
//------------------------------------------------------------------
  CMGCovalentDistance cd(std::string(a1),std::string(a2),mind,maxd);
  covalentDistances.push_back(cd);
  return 0;
}

//---------------------------------------------------------------
int CMMUTSRS::GetCharge( mmdb::mmcif::Loop *Loop, int N ) {
//----------------------------------------------------------------
  int RC;
  mmdb::pstr type=0;
  mmdb::pstr hbtype=0;
  mmdb::pstr el=0;
  mmdb::realtype vdwRadius,vdwHRadius,ionRadius,charge;
  if (!Loop){
    return -1;
  }
  if (N >= Loop->GetLoopLength()) {
   return -1;
  }

  RC = Loop->GetString ( type, CIFTAG_LIBATOM_TYPE,N);
  if ( RC || (!type) ) {
    return -1;
  }
  RC = Loop->GetString(hbtype, CIFTAG_LIBATOM_HBTYPE,N);
  if(strlen(hbtype)>0){
    typeHBonds[std::string(type)] = hbtype[0];
  } else {
    typeHBonds[std::string(type)] = 'N';
  }
  RC = Loop->GetReal(vdwRadius, CIFTAG_LIBATOM_VDWRAD, N);
  RC = Loop->GetReal(vdwHRadius, CIFTAG_LIBATOM_VDWHRAD, N);
  RC = Loop->GetReal(ionRadius,CIFTAG_LIBATOM_IONRAD, N );
  RC = Loop->GetString (el,CIFTAG_LIBATOM_ELEMENT,N);
  RC = Loop->GetReal(charge,CIFTAG_LIBATOM_CHARGE, N );
  typeCharges[std::string(type)] = charge;
  //std::cout << type << " " << charge << "\n";
  return 0;
}


int CMMUTSRS::LoadEnerLib(const std::string &elib) {
  int RC;
  mmdb::mmcif::Loop *Loop1;
  mmdb::mmcif::Data* CIF =  new mmdb::mmcif::Data();
  const char *filename = elib.c_str();
  RC = CIF->ReadMMCIFData ( filename );
  if ( RC ) {
    printf ("Error reading %s\n",filename);
    delete CIF;
    return -1;
  }
  Loop1 = CIF->GetLoop ( CIFCAT_LIBATOM );
  if ( !Loop1 ) {
    delete CIF;
    return -1;
  }

  RC =0;
  int nLibAtoms = 0;
  do {
    RC =  GetCharge ( Loop1 , nLibAtoms );
    if ( !RC) {
      nLibAtoms++;
    }
  }while (!RC );

  delete CIF;
  return 0;

}

//---------------------------------------------------------------
int CMMUTSRS::GetElement( mmdb::mmcif::Loop * Loop, int N ) {
//----------------------------------------------------------------
  int RC;
  mmdb::pstr el=0;
  mmdb::pstr type=0;
  
  if (!Loop) return -1;
  if (N >= Loop->GetLoopLength()) return -1;
   
  RC = Loop->GetString (el,"name",N);
  if (RC||!el) return RC;
  RC = Loop->GetString ( type, "lib_atom_type",N);
  if (RC||!type) return RC;
  //std::cout << "Adding default atom type--" << el << "--  %%" << type << "%%\n";
  defaultAtomType[std::string(el)] = std::string(type);
  return RC;
}

//------------------------------------------------------------------
int CMMUTSRS::LoadEleLib (const std::string &elelib) {
//------------------------------------------------------------------
  // Load the monomer link info from CIF file
  int RC;
  mmdb::mmcif::Loop *Loop2;
  mmdb::mmcif::Data* CIF =  new mmdb::mmcif::Data();
  const char *filename = elelib.c_str();
  RC = CIF->ReadMMCIFData ( filename );
  if ( RC ) {
    printf ("Error reading %s\n",filename);
    delete CIF;
    return -1;
  }

  Loop2 = CIF->GetLoop ( "_lib_element" );
  if ( !Loop2 ) {
    printf("Error reading element info from ener_lib.cif\n");
    delete CIF;
    return -1;
  }
  
  int nLibElements = 0;
  RC =0;
  do {
    RC =  GetElement ( Loop2 , nLibElements );
    if ( !RC) {
      nLibElements++;
    }
  }while (!RC);


  delete CIF;

  return 0;
}

//------------------------------------------------------------------
int CMMUTSRS::LoadMonomerCache (ccp4srs::Manager *SRS) {
//------------------------------------------------------------------
  mmdb::io::PFile structFile = SRS->getStructFile();
  for(int iaa=0;iaa<mmdb::nResNames;iaa++){
     ccp4srs::Monomer* monomer = SRS->getMonomer(mmdb::ResidueName[iaa],structFile);
     MonomerCache[mmdb::ResidueName[iaa]] = monomer;
  }
  for(int iaa=0;iaa<mmdb::nNucleotideNames;iaa++){
     ccp4srs::Monomer* monomer = SRS->getMonomer(mmdb::NucleotideName[iaa],structFile);
     MonomerCache[mmdb::NucleotideName[iaa]] = monomer;
  }
  for(int iaa=0;iaa<mmdb::nSolventNames;iaa++){
     ccp4srs::Monomer* monomer = SRS->getMonomer(mmdb::StdSolventName[iaa],structFile);
     MonomerCache[mmdb::StdSolventName[iaa]] = monomer;
  }
  delete structFile;
  return 0;
}

std::vector<CMGCovalentDistance> CMMUTSRS::covalentDistances;
std::map<std::string,mmdb::realtype> CMMUTSRS::typeCharges;
std::map<std::string,char> CMMUTSRS::typeHBonds;
std::map<std::string,std::string> CMMUTSRS::defaultAtomType;
std::map<std::string,ccp4srs::Monomer*> CMMUTSRS::MonomerCache;
