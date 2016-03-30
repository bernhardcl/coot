#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <exception>
#include <string.h>
#include <math.h>
#include "GraphMol/SmilesParse/SmilesParse.h"
#include "GraphMol/FileParsers/MolSupplier.h"
#include "GraphMol/FileParsers/MolWriters.h"
#include "GraphMol/FileParsers/FileParsers.h"
#include "GraphMol/DistGeomHelpers/Embedder.h"
#include "GraphMol/ROMol.h"
#include "GraphMol/RWMol.h"
#include "GraphMol/MolOps.h"
#include "ForceField/ForceField.h"
#include "GraphMol/ForceFieldHelpers/UFF/AtomTyper.h"
#include "GraphMol/ForceFieldHelpers/UFF/Builder.h"
#include "GraphMol/AtomIterators.h"
#include "GraphMol/BondIterators.h"
#include "GraphMol/Conformer.h"
#include "GraphMol/MonomerInfo.h"
#include "GraphMol/SanitException.h"
#include "GraphMol/PeriodicTable.h"
#include "mmdb2/mmdb_atom.h"
#include "ccp4srs/ccp4srs_bond.h"
#ifdef _DO_GLYCO_TESTING_
#include "clipper-glyco.h"
#endif

#include "GraphMol/Depictor/RDDepictor.h"
#include "GraphMol/MolDrawing/MolDrawing.h"
#include "GraphMol/MolDrawing/DrawingToSVG.h"

#include "mmut_rdkit.h"

int MapAtomToResidueAtomNum(mmdb::PResidue res,mmdb::PAtom at){
  int iat = 0;
  for(int i=0;i<res->GetNumberOfAtoms();i++){
    if(strlen(res->GetAtom(i)->altLoc)==0||(!strncmp(res->GetAtom(i)->altLoc,"A",1))){
      if(res->GetAtom(i)==at){
        return iat;
      }
      iat++;
    }
  }
  return -1;
}

std::string RDKitResidueToSVG(RDKit::RWMol *mol){
  /* 
    Return SVG from RDKit::RWMol *mol;
  */

  if(!mol) return std::string("");

  RDKit::RWMol *cp = new RDKit::RWMol(*mol);
  RDKit::MolOps::removeHs(*cp);
  RDDepict::compute2DCoords(*cp);
  std::vector<int> drawing=RDKit::Drawing::MolToDrawing(*cp);
  std::string svg=RDKit::Drawing::DrawingToSVG(drawing);
  delete cp;

  return svg;

}

RDKit::RWMol *RDKitAromaticityAnalyze(mmdb::PResidue res, bool aromaticize){
  // FIXME - Needs to consider altLocs!
  /* 
    Uses rdkit to determine if atoms are involved in aromatic bonds or not.
    Returns 0 if errors trapped, the RDKit::RWMol otherwise.
  */
  //std::cout << "RDKitAromaticityAnalyze, residue: " << res->name << "\n";

  if(!res) return 0;
  if(res->GetNumberOfAtoms()<1) return 0;
  RDKit::RWMol *mol = 0;

  mmdb::io::File f;
  int buffer_size = res->GetNumberOfAtoms()*162+2;
  char *fbuffer = new char[buffer_size];
  for(int i=0;i<buffer_size;i++) fbuffer[i] ='\0';
  f.assign(buffer_size,1,fbuffer);

  res->PDBASCIIAtomDump(f);

  f.flush();
  f.shut();
  //std::cout << fbuffer << "\n"; std::cout.flush();

  try {
    mol = RDKit::PDBBlockToMol(fbuffer,true,false);
  } catch (RDKit::MolSanitizeException &e) {
    std::cerr << "RDKitAromaticityAnalyze: Some problem reading my own pdb output.\n";
    std::cerr << res->name << "\n";
    std::cerr << fbuffer << "\n";
  } catch (...) {
    std::cerr << "RDKitAromaticityAnalyze: Some problem reading my own pdb output.\n";
    std::cerr << res->name << "\n";
    std::cerr << fbuffer << "\n";
  }
  delete [] fbuffer;
  if(!mol) return 0;

  mmdb::PPAtom atomTable = 0;
  int nAtoms;
  res->GetAtomTable(atomTable,nAtoms);

  //std::cout << "mmdb nAtoms " << nAtoms << "\n";
  //std::cout << "rdKit nAtoms " << mol->getNumAtoms() << "\n";

  int iat = 0;
  std::vector<int> halfPlusCharges;
  std::vector<int> halfMinusCharges;
  for(int i=0;i<nAtoms;i++){
    if(strlen(atomTable[i]->altLoc)==0||(!strncmp(atomTable[i]->altLoc,"A",1))){
      mmdb::AtomBond *atomBond;
      mmdb::AtomBondI *atomBondI=0;
      int nAtomBond;
      atomTable[i]->GetBonds(atomBondI,nAtomBond);
      atomTable[i]->GetBonds(atomBond,nAtomBond);
      //std::cout  << "RDKitAromaticityAnalyze: mmdb charge: " << atomTable[i]->charge << "\n"; std::cout.flush();
      mol->getAtomWithIdx(iat)->setFormalCharge(int(atomTable[i]->charge));
      if(fabs(atomTable[i]->charge-0.5)<1e-4){
        halfPlusCharges.push_back(iat);
        //std::cout << "atom " << iat << "(" << mol->getAtomWithIdx(iat)->getSymbol() << ") has +.5 charge\n";
        //std::cout << atomTable[i]->name << "\n";
      }
      if(fabs(atomTable[i]->charge+0.5)<1e-4) {
       halfMinusCharges.push_back(iat);
        //std::cout << "atom " << iat << "(" << mol->getAtomWithIdx(iat)->getSymbol() << ") has -.5 charge\n";
        //std::cout << atomTable[i]->name << "\n";
      }
      //std::cout << "RDKitAromaticityAnalyze: charge set to: " << mol->getAtomWithIdx(iat)->getFormalCharge() << "\n"; std::cout.flush();
      for(int j=0;j<nAtomBond;j++){
        if(strlen(atomBond[j].atom->altLoc)==0||(!strncmp(atomBond[j].atom->altLoc,"A",1))){
          int bondedAtom = MapAtomToResidueAtomNum(res,atomBond[j].atom);
          //std::cout << "RDKitAromaticityAnalyze: check bond indices " << iat << " " << bondedAtom << " " << mol->getNumAtoms() << "\n"; std::cout.flush();
          if(bondedAtom>-1){
            RDKit::Bond *bond = mol->getBondBetweenAtoms(iat,bondedAtom);
            if(bond) {
              if(atomBondI[j].order==ccp4srs::Bond::Double){
                //std::cout << "Trying to add double " << iat << " " << bondedAtom << "\n"; std::cout.flush();
                mol->getBondBetweenAtoms(iat,bondedAtom)->setBondType(RDKit::Bond::DOUBLE);
              }else if(atomBondI[j].order==ccp4srs::Bond::Triple){
                //std::cout << "Trying to add triple " << iat << " " << bondedAtom << "\n"; std::cout.flush();
                mol->getBondBetweenAtoms(iat,bondedAtom)->setBondType(RDKit::Bond::TRIPLE);
              }else if(atomBondI[j].order==ccp4srs::Bond::Aromatic){
                //std::cout << "Trying to add aromatic " << iat << " " << bondedAtom << "\n";
                // Do nothing.
              } else {
                //std::cout << "Trying to add single " << iat << " " << bondedAtom << "\n"; std::cout.flush();
                mol->getBondBetweenAtoms(iat,bondedAtom)->setBondType(RDKit::Bond::SINGLE);
              }
            }
          }
        }
      }
      delete atomBondI;
      iat++;
    }
  }

  RDKit::ROMol::BondIterator bondIter = mol->beginBonds();
  
  std::map<int,std::vector<int> > halfPlusBondCharges;
  std::map<int,std::vector<int> > halfMinusBondCharges;
  while(bondIter!=mol->endBonds()){
    int idx1 = (*bondIter)->getBeginAtomIdx();
    int idx2 = (*bondIter)->getEndAtomIdx();
    if(std::count(halfPlusCharges.begin(),halfPlusCharges.end(),idx2)){
      halfPlusBondCharges[idx1].push_back(idx2);
    }
    if(std::count(halfPlusCharges.begin(),halfPlusCharges.end(),idx1)){
      halfPlusBondCharges[idx2].push_back(idx1);
    }
    if(std::count(halfMinusCharges.begin(),halfMinusCharges.end(),idx2)){
      halfMinusBondCharges[idx1].push_back(idx2);
    }
    if(std::count(halfMinusCharges.begin(),halfMinusCharges.end(),idx1)){
      halfMinusBondCharges[idx2].push_back(idx1);
    }
    bondIter++;
  }
  
  std::map<int,std::vector<int> >::iterator plusIter;
  std::map<int,std::vector<int> >::iterator minusIter;

  for(plusIter=halfPlusBondCharges.begin();plusIter!=halfPlusBondCharges.end();plusIter++){
    //std::cout << ".5+ " << plusIter->first << "\n";
    std::vector<int>::iterator connIter;
    for(connIter=plusIter->second.begin();connIter!=plusIter->second.end();connIter++){
      //std::cout << *connIter << " ";
    }
    //std::cout << "\n";
  }
 
  std::vector<int> simpleMinus;
  for(minusIter=halfMinusBondCharges.begin();minusIter!=halfMinusBondCharges.end();minusIter++){
    //std::cout << ".5- " << minusIter->first << "\n";
    std::vector<int>::iterator connIter;
    for(connIter=minusIter->second.begin();connIter!=minusIter->second.end();connIter++){
      //std::cout << *connIter << " ";
      simpleMinus.push_back(*connIter);
    }
    //std::cout << "\n";
  }
 
  std::vector<int>::iterator simpleMinusIter;
  for(simpleMinusIter=simpleMinus.begin();simpleMinusIter!=simpleMinus.end();simpleMinusIter++){
    //std::cout << *simpleMinusIter << " " << std::count(simpleMinus.begin(),simpleMinus.end(),*simpleMinusIter)<< "\n";
    if(std::count(simpleMinus.begin(),simpleMinus.end(),*simpleMinusIter)==1){
      //mol->getAtomWithIdx(*simpleMinusIter)->setFormalCharge(-1);
      bondIter = mol->beginBonds();
      while(bondIter!=mol->endBonds()){
        int idx1 = (*bondIter)->getBeginAtomIdx();
        int idx2 = (*bondIter)->getEndAtomIdx();
        if(*simpleMinusIter==idx1||*simpleMinusIter==idx2){
          //std::cout << "Would like to add double bond between " << idx1 << " " << idx2 << " " << mol->getAtomWithIdx(idx1)->getSymbol() << " " << mol->getAtomWithIdx(idx2)->getSymbol() << "\n";
          //std::cout << "Valencies before : " << mol->getAtomWithIdx(idx1)->getExplicitValence() << " " << mol->getAtomWithIdx(idx2)->getExplicitValence() << "\n";
          const std::vector<int> &valens1 = RDKit::PeriodicTable::getTable()->getValenceList(mol->getAtomWithIdx(idx1)->getAtomicNum());
          int maxValence1=*(valens1.rbegin());
          const std::vector<int> &valens2 = RDKit::PeriodicTable::getTable()->getValenceList(mol->getAtomWithIdx(idx2)->getAtomicNum());
          int maxValence2=*(valens2.rbegin());
          if(mol->getAtomWithIdx(idx1)->getExplicitValence()<maxValence1&&mol->getAtomWithIdx(idx2)->getExplicitValence()<maxValence2){

            (*bondIter)->setBondType(RDKit::Bond::DOUBLE);
          }
          try {
            RDKit::MolOps::sanitizeMol(*mol);
          } catch (...) {
            std::cerr << "RDKitAromaticityAnalyze: Some problem adding bonds (whilst looping) for aromaticity test.\n";
            std::cerr << res->name << "\n";
            return 0;
          } 
          //std::cout << "Valencies after : " << mol->getAtomWithIdx(idx1)->getExplicitValence() << " " << mol->getAtomWithIdx(idx2)->getExplicitValence() << "\n";
        }
        bondIter++;
      }
    }
  }

  try {
    RDKit::MolOps::sanitizeMol(*mol);
  } catch (...) {
    std::cerr << "RDKitAromaticityAnalyze: Some problem adding bonds for aromaticity test.\n";
    std::cerr << res->name << "\n";
    return 0;
  }
 
  if(!aromaticize) return mol;

  try {
  char AtomID1[1024];
  char AtomID2[1024];
  iat = 0;
  for(int i=0;i<nAtoms;i++){
    if(strlen(atomTable[i]->altLoc)==0||(!strncmp(atomTable[i]->altLoc,"A",1))){
      mmdb::AtomBond *atomBond;
      int nAtomBond;
      atomTable[i]->GetAtomID(AtomID1);
      atomTable[i]->GetBonds(atomBond,nAtomBond);
      for(int j=0;j<nAtomBond;j++){
        atomBond[j].atom->GetAtomID(AtomID2);
        //std::cout << "Consider " << AtomID1 << " " << AtomID2 << "\n";
        if(strlen(atomBond[j].atom->altLoc)==0||(!strncmp(atomBond[j].atom->altLoc,"A",1))){
          int bondedAtom = MapAtomToResidueAtomNum(res,atomBond[j].atom);
          if(bondedAtom>-1){
            RDKit::Bond *bond = mol->getBondBetweenAtoms(iat,bondedAtom);
            if(bond) {
              if(bond->getIsAromatic()) {
                //std::cout << "Setting " << AtomID1 << " " << AtomID2 << " aromatic\n";
                atomBond[j].order=ccp4srs::Bond::Aromatic;
              } else {
                //std::cout << AtomID1 << " " << AtomID2 << " is not aromatic\n";
              }
            } else {
              std::cout << "RDKitAromaticityAnalyze: No matching bond\n";
            }
          } else {
            std::cout << "RDKitAromaticityAnalyze: Fail to get bonded atom\n";
          }
        } else {
          //std::cout << "Fail altLoc test\n";
        }
      }
      iat++;
    }
  }
  } catch (...) {
    std::cerr << "RDKitAromaticityAnalyze: Some matching bonds for aromaticity test.\n";
    std::cerr << res->name << "\n";
    return 0;
  }

  return mol;

}

void writeCIF(RDKit::RWMol *mol,const char *cifname, int confId=0){
  RDKit::AtomMonomerInfo *info = mol->getAtomWithIdx(0)->getMonomerInfo();
  std::map<std::string,int> elemMap;
  std::map<int,std::string> atomMap;
    bool outFileRequested;
    std::ofstream realOutFile;
    if(strlen(cifname)==1&&strncmp(cifname,"-",1)==0){
      outFileRequested = false;
    } else {
      outFileRequested = true;
      realOutFile.open(cifname, std::ios::out);
    }
    std::ostream &output = (outFileRequested ? realOutFile : std::cout);
    output << "global_\n";
    output << "_lib_name         ?\n";
    output << "_lib_version      ?\n";
    output << "_lib_update       ?\n";
    output << "# ------------------------------------------------\n";
    output << "#\n";
    output << "# ---   LIST OF MONOMERS ---\n";
    output << "#\n";
    output << "data_comp_list\n";
    output << "loop_\n";
    output << "_chem_comp.id\n";
    output << "_chem_comp.three_letter_code\n";
    output << "_chem_comp.name\n";
    output << "_chem_comp.group\n";
    output << "_chem_comp.number_atoms_all\n";
    output << "_chem_comp.number_atoms_nh\n";
    output << "_chem_comp.desc_level\n";
    int nAtAll = 0;
    int nAtNoH = 0;
    RDKit::ROMol::AtomIterator atomIter = mol->beginAtoms();
    while(atomIter!=mol->endAtoms()){
      std::string symbol = (*atomIter)->getSymbol();
      nAtAll++;
      if(symbol!="H"){
        nAtNoH++;
      }
      atomIter++;
    }
    output << "UNL      UNL 'UNKNOWN LIGAND                      ' non-polymer        " << nAtAll << " " <<  nAtNoH << " .\n";
    output << "# ------------------------------------------------------\n";
    output << "# ------------------------------------------------------\n";
    output << "#\n";
    output << "# --- DESCRIPTION OF MONOMERS ---\n";
    output << "#\n";
    output << "data_comp_UNL\n";
    output << "#\n";
    output << "loop_\n";
    output << "_chem_comp_atom.comp_id\n";
    output << "_chem_comp_atom.atom_id\n";
    output << "_chem_comp_atom.type_symbol\n";
    output << "_chem_comp_atom.type_energy\n";
    output << "_chem_comp_atom.partial_charge\n";
    output << "_chem_comp_atom.x\n";
    output << "_chem_comp_atom.y\n";
    output << "_chem_comp_atom.z\n";
    // FIXME Loop over atoms.
    // FAD           O2P    O    OP       -0.500      0.000    0.000    0.000
    atomIter = mol->beginAtoms();
    const RDKit::Conformer *conf=&(mol->getConformer(confId));
    while(atomIter!=mol->endAtoms()){
      double x, y, z;
      const RDGeom::Point3D pos = conf->getAtomPos((*atomIter)->getIdx());
      x = pos.x; y = pos.y; z = pos.z;
  
      std::string symbol = (*atomIter)->getSymbol();
      int charge = (*atomIter)->getFormalCharge();
  
      if(info){
        std::string name = info->getName();
        output << "UNL           " << name << " " << symbol  << " " << name << " " << charge << " " << x << " " << y << " " << z << "\n";
      } else {
        if(elemMap.count(symbol)){
           elemMap[symbol]++;
        } else {
           elemMap[symbol] = 1;
        }
        std::stringstream s;
        s << elemMap[symbol];
        std::string iStr = s.str();
        std::string name = symbol+iStr;
        atomMap[(*atomIter)->getIdx()] = name;
        output << "UNL           " << name << " " << symbol  << " " << name << " " << charge << " " << x << " " << y << " " << z << "\n";
      }
  
      atomIter++;
    }
    output << "loop_\n";
    output << "_chem_comp_bond.comp_id\n";
    output << "_chem_comp_bond.atom_id_1\n";
    output << "_chem_comp_bond.atom_id_2\n";
    output << "_chem_comp_bond.type\n";
    output << "_chem_comp_bond.value_dist\n";
    output << "_chem_comp_bond.value_dist_esd\n";
    RDKit::ROMol::BondIterator bondIter = mol->beginBonds();
    while(bondIter!=mol->endBonds()){
      std::string beginAtom =  atomMap[(*bondIter)->getBeginAtomIdx()];
      std::string endAtom =  atomMap[(*bondIter)->getEndAtomIdx()];
      const RDGeom::Point3D pos1 = conf->getAtomPos((*bondIter)->getBeginAtomIdx());
      const RDGeom::Point3D pos2 = conf->getAtomPos((*bondIter)->getEndAtomIdx());
      double bondLength = sqrt((pos1.x-pos2.x)*(pos1.x-pos2.x) + (pos1.y-pos2.y)*(pos1.y-pos2.y) + (pos1.z-pos2.z)*(pos1.z-pos2.z));
      std::string bondType;
      if((*bondIter)->getBondType()==RDKit::Bond::DOUBLE){
         bondType = "double";
      } else if((*bondIter)->getBondType()==RDKit::Bond::TRIPLE){
         bondType = "triple";
      } else if((*bondIter)->getBondType()==RDKit::Bond::AROMATIC){
         bondType = "aromatic";
      } else {
         bondType = "single";
      }
      output << "UNL  " << "    " << beginAtom <<   " "  << endAtom  << " " << bondType  << " " << bondLength << "  0.020\n";
      bondIter++;
    }
    output << "# ------------------------------------------------------\n";
    if(outFileRequested){
      realOutFile.close();
    }

}

const char* SmilesToPDB(const char *smile,const char *fname, const char *cifname, int nconf, int maxIters){
  if(!smile||strlen(smile)<1) return "Zero length SMILES string.";
  RDKit::RWMol *mol = 0;
  try {
    mol =  RDKit::SmilesToMol(smile);
  } catch (RDKit::MolSanitizeException &e) {
    return e.message();
  } catch (...) {
    return "SMILES parse error?";
  }
  if(!mol) return "SMILES parse error?";
  try {
    RDKit::MolOps::addHs(*mol);
  } catch (...) {
    return "Error adding hydrogens.";
  }

  try {
    RDKit::DGeomHelpers::EmbedMolecule(*mol);
  } catch (...) {
    return "Error embedding molecule.";
  }

  int minCid = MolMinimize(mol, nconf, maxIters);
  if(minCid==-1) return "Minimize from SMILES error";

  try {
    RDKit::PDBWriter writer(fname);
    try {
      writer.write(*mol,minCid);
    } catch (...) {
      return "Error writing mol.";
    }
    try {
        writer.flush();
    } catch (...) {
      return "Error flushing writer.";
    }
    try {
      writer.close();
    } catch (...) {
      return "Error closing writer.";
    }
  } catch (...) {
    return "Error creating writer.";
  }

  try {
    if(cifname&&strlen(cifname)>0){
      writeCIF(mol,cifname,minCid);
    }
  } catch (...) {
    return "Error writing CIF";
  }

  return 0;

}

int MMDBMinimize(mmdb::Manager *molHnd, int nconf, int maxIters, ccp4srs::Monomer *monomer){
  std::cout << "Now minimize .... ?\n";

  RDKit::RWMol *mol = 0;

  int nAtoms;
  mmdb::PPAtom atomTable = 0;
  molHnd->GetAtomTable(atomTable,nAtoms);

  std::cout << nAtoms << std::endl;

  mmdb::io::File f;
  int buffer_size = nAtoms*162+2;
  char *fbuffer = new char[buffer_size];
  for(int i=0;i<buffer_size;i++) fbuffer[i] ='\0';
  f.assign(buffer_size,1,fbuffer);

  std::cout << buffer_size << std::endl;

  molHnd->GetModel(1)->PDBASCIIDump(f);

  f.flush();
  f.shut();

  std::cout << fbuffer << "\n"; std::cout.flush();

  try {
    mol = RDKit::PDBBlockToMol(fbuffer,true,true);
  } catch (RDKit::MolSanitizeException &e) {
    std::cerr << "MMDBMinimize: Some problem reading my own pdb output.\n";
  }

  RDKit::ROMol::AtomIterator atomIter = mol->beginAtoms();
  while(atomIter!=mol->endAtoms()){
    std::string symbol = (*atomIter)->getSymbol();
    std::cout << symbol << std::endl;
    atomIter++;
  }

  std::cout << MolMinimize(mol,nconf,maxIters,monomer);

  return 0;
}

int MolMinimize(RDKit::RWMol *mol, int nconf, int maxIters, ccp4srs::Monomer *monomer){
  //RDKit::MolOps::Kekulize(*mol); // This futzes with bond lengths, which I do not want!

#ifdef _DO_GLYCO_TESTING_
  bool force4C1chair = false;
#endif

  double vdwThresh=10.0;
  int confId=-1;
  bool ignoreInterfragInteractions=true;

  double minE = 1e+32;
  int minCid = -1;
  RDKit::INT_VECT cids=RDKit::DGeomHelpers::EmbedMultipleConfs(*mol, nconf);
  for(unsigned icid=0;icid<cids.size();icid++){
   
#ifdef _DO_GLYCO_TESTING_
    if(force4C1chair) {
      clipper::MiniMol minimol;
      clipper::MPolymer mp;
      clipper::MMonomer mm;
      mm.set_type("BGC");
      mm.set_seqnum(1);
      mm.set_id(1);
      minimol.init(clipper::Spacegroup::p1(),clipper::Cell(clipper::Cell_descr(300,300,300,90,90,90)));
      const RDKit::Conformer *conf=&(mol->getConformer(icid));
      RDKit::ROMol::AtomIterator atomIter = mol->beginAtoms();
      std::map<std::string,int> elemMap;
      std::map<int,std::string> atomMap;
      while(atomIter!=mol->endAtoms()){
        double x, y, z;
        const RDGeom::Point3D pos = conf->getAtomPos((*atomIter)->getIdx());
        x = pos.x; y = pos.y; z = pos.z;
        std::string symbol = (*atomIter)->getSymbol();
        clipper::MAtom mat;
        mat.set_coord_orth(clipper::Coord_orth(x,y,z));
        mat.set_element(symbol);
        if(elemMap.count(symbol)){
           elemMap[symbol]++;
        } else {
           elemMap[symbol] = 1;
        }
          std::stringstream s;
          s << elemMap[symbol];
          std::string iStr = s.str();
          std::string name = symbol+iStr;
        mat.set_id(name);
        mat.set_occupancy(1.0);
        mat.set_u_iso(20.0);
        mm.insert(mat);
        atomIter++;
      }
      //std::cout << "Monomer size " << mm.size() << std::endl;
      mp.insert(mm);
      //std::cout << "Polymer size " << mp.size() << std::endl;
      minimol.insert(mp);
      //std::cout << "Polymer atom list " << mp.atom_list().size() << std::endl;
      clipper::MMonomer theCopy;
      theCopy.copy(mm,clipper::MM::COPY_MPC);
      //std::cout << "Copy size " << theCopy.size() << std::endl;
      const clipper::MAtomNonBond nb = clipper::MAtomNonBond(minimol, 5.0);
      clipper::MSugar sugar(minimol,mm,nb);
      std::vector<clipper::ftype> cpParams = sugar.cremer_pople_params();
      if(cpParams.size()>2){
        //std::cout << "Q = " << cpParams[0] << std::endl;
        //std::cout << "PHI = " << cpParams[1] << std::endl;
        //std::cout << "THETA = " << cpParams[2] << std::endl;
        clipper::ftype theta = cpParams[2];
        if(theta<0||theta>20) continue;
      } else {
        // Not a sugar!
        continue;
      }
      std::cout << sugar.type_of_sugar() << std::endl;
    }
#endif

    ForceFields::ForceField *ff;
    try {
      ff=RDKit::UFF::constructForceField(*mol,vdwThresh, cids[icid],ignoreInterfragInteractions);
    } catch (...) {
      std::cout << "Error constructing forcefield.\n";
      return -1;
    }
    try {
      ff->initialize();
    } catch (...) {
      std::cout << "Error initializing forcefield.\n";
      return -1;
    }
    int res;
    try {
      res=ff->minimize(maxIters);
    } catch (...) {
      std::cout << "Error minimizing forcefield.\n";
      return -1;
    }

    double E;
    try {
      E = ff->calcEnergy();
    } catch (...) {
      std::cout << "Error calculating energy.\n";
      return -1;
    }
    if(E<minE){
      minE = E;
      minCid = icid;
    }
    delete ff;
  }
  std::cout << minE << "\n";

  return minCid;
}

const char* RDKitLoadFileToPDB(const char *input_fname,const char *fname,const char *cifname){
  std::cout << "RDKitLoadFileToPDB\n"; std::cout.flush();
  
  RDKit::RWMol *mol = 0;
  try {
    std::cout << "Try mol2\n"; std::cout.flush();
    mol =  RDKit::Mol2FileToMol(input_fname);
    std::cout << "Done mol2\n"; std::cout.flush();
  } catch (std::exception &e) {
    std::cout << e.what() << "\n";
    try {
      std::cout << "Try mol\n"; std::cout.flush();
      mol =  RDKit::MolFileToMol(input_fname,true,false);
      std::cout << "Done mol\n"; std::cout.flush();
    } catch (std::exception &e) {
      std::cout << e.what() << "\n";
      try {
        std::cout << "Try tpl\n"; std::cout.flush();
        mol =  RDKit::TPLFileToMol(input_fname);
        std::cout << "Done tpl\n"; std::cout.flush();
      } catch (std::exception &e) {
        std::cout << "Last exception\n"; std::cout.flush();
        return e.what();
      }
    }
  }
  std::cout << "Made mol?\n"; std::cout.flush();
  if(!mol) return "RDKit::mol from file error?";

  RDKit::PDBWriter writer(fname);
  writer.write(*mol);
  writer.flush();
  writer.close();

  if(cifname&&strlen(cifname)>0){
    writeCIF(mol,cifname);
  }

  return 0;
}


