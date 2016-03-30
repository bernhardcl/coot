#ifndef __MMUT_SRS__
#define __MMUT_SRS__
#include "ccp4srs/ccp4srs_manager.h"
#include "ccp4srs/ccp4srs_monomer.h"
#include "mmdb2/mmdb_mattype.h"
#include "mmdb2/mmdb_cifdefs.h"
#include <string>
#include <map>
#include <vector>

enum { RESTYPE_PEPTIDE, RESTYPE_DPEPTIDE, RESTYPE_LPEPTIDE,
       RESTYPE_NUCL,RESTYPE_DNA, RESTYPE_RNA,
       RESTYPE_SACH, RESTYPE_DSACH, RESTYPE_LSACH,
       RESTYPE_SOLVENT, RESTYPE_SOLUTE,RESTYPE_NONPOLY, RESTYPE_PSEUDO,
       RESTYPE_METAL,
       RESTYPE_UNKNOWN };

// CIF tags for ener_lib.cif file the lib_atom table
#define CIFCAT_LIBATOM "_lib_atom"
#define CIFTAG_LIBATOM_TYPE "type"
#define CIFTAG_LIBATOM_WEIGHT "weight"
#define CIFTAG_LIBATOM_HBTYPE "hb_type"
#define CIFTAG_LIBATOM_VDWRAD "vdw_radius"
#define CIFTAG_LIBATOM_VDWHRAD "vdwh_radius"
#define CIFTAG_LIBATOM_IONRAD "ion_radius"
#define CIFTAG_LIBATOM_ELEMENT "element"
#define CIFTAG_LIBATOM_VALENCY "valency"
#define CIFTAG_LIBATOM_CHARGE "surface_potential_charge"

class MMUTSRSAtomInfo {
  std::string _id;
  std::string energyType;
  mmdb::realtype vdwRadius;
  mmdb::realtype vdwHRadius;
  mmdb::realtype ionRadius;
  char hbType;
 public:
  std::string energy_type() const {return energyType;};
  mmdb::realtype vdw_radius() const {return vdwRadius;};
  mmdb::realtype vdwh_radius() const {return vdwHRadius;};
  mmdb::realtype ion_radius() const {return ionRadius;};
  char hb_type() const {return hbType;};
  void set_energy_type(const std::string &energyType_in) {energyType = energyType_in;};
  void set_vdw_radius(const mmdb::realtype vdwRadius_in) {vdwRadius = vdwRadius_in;};
  void set_vdwh_radius(const mmdb::realtype vdwHRadius_in) {vdwHRadius = vdwHRadius_in;};
  void set_ion_radius(const mmdb::realtype ionRadius_in) {ionRadius = ionRadius_in;};
  void set_hb_type(const char hbType_in) {hbType = hbType_in;};
};

class CMGCovalentDistance {
  std::string a1;
  std::string a2;
  mmdb::realtype mind;
  mmdb::realtype maxd;
  public:
    CMGCovalentDistance(std::string a1_in,std::string a2_in,mmdb::realtype mind_in,mmdb::realtype maxd_in):a1(a1_in),a2(a2_in),mind(mind_in),maxd(maxd_in){};
    std::string GetFirstAtom() const {return a1;};
    std::string GetSecondAtom() const {return a2;};
    mmdb::realtype GetMinLength() const {return mind;};
    mmdb::realtype GetMaxLength() const {return maxd;};
};

class CMMUTSRS {
  static std::vector<CMGCovalentDistance> covalentDistances;
 public:
  CMMUTSRS(){};
  static int CheckCovalentDistance (mmdb::Element a1, mmdb::Element a2,mmdb::realtype d);
  static int LoadCovalentDistances(mmdb::pstr filename ) ;
  static int AddCovalentDistance(mmdb::Element a1, mmdb::Element a2,mmdb::realtype mind,mmdb::realtype maxd) ;
  static int GetCharge( mmdb::mmcif::Loop *Loop, int N );
  static int LoadEnerLib(const std::string &elib);
  static int LoadEleLib(const std::string &elib);
  static int GetElement( mmdb::mmcif::Loop * Loop, int N );
  static int LoadMonomerCache (ccp4srs::Manager *SRS);
  //FIXME Should have accessor;
  static std::map<std::string,mmdb::realtype> typeCharges;
  static std::map<std::string,char> typeHBonds;
  static std::map<std::string,std::string> defaultAtomType;
  static std::map<std::string,ccp4srs::Monomer*> MonomerCache;
};

#endif //__MMUT_SRS__
