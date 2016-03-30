/*
     util/CParamsManager.cc: CCP4MG Molecular Graphics Program
     Copyright (C) 2001-2008 University of York, CCLRC
     Copyright (C) 2009 University of York

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

#include "CParamsManager.h"
#include <map>
#include <string>
#include <iostream>

CParamsManager::CParamsManager() {}
CParamsManager::~CParamsManager() {}

void CParamsManager::SetString ( const std::string &key, const std::string &value ){
  Strings[key] = value;
}

void CParamsManager::SetFloat ( const std::string &key, const float value ){
  Floats[key] = value;
}

void CParamsManager::SetInt ( const std::string &key, const int value ) {
  Ints[key]=value; 
}

void CParamsManager::PrintElements () const {
  PrintElements(std::cout);
}

void CParamsManager::PrintElements ( std::ostream &c ) const {
	std::map<std::string,int>::const_iterator ielem;
	for(ielem=Ints.begin();ielem!=Ints.end();++ielem){
		c << "int: " << ielem->first << " " << ielem->second << "\n";
	}
	std::map<std::string,float>::const_iterator felem;
	for(felem=Floats.begin();felem!=Floats.end();++felem){
		c << "float: " << felem->first << " " << felem->second << "\n";
	}
	std::map<std::string,std::string>::const_iterator selem;
	for(selem=Strings.begin();selem!=Strings.end();++selem){
		c << "string: " << selem->first << " " << selem->second << "\n";
	}
}

std::string CParamsManager::GetString (  const std::string &key ) const {
  std::map<std::string,std::string>::const_iterator elem;
  elem = Strings.find( key );
  if ( elem == Strings.end() ) {
    std::cout << "GetInt no value " << key << "\n";
    return std::string("");
  } else {
      return elem->second ;
  }
  
}

int CParamsManager::GetInt (  const std::string &key ) const {
  std::map<std::string,int>::const_iterator elem;
  elem = Ints.find( key );
  if ( elem == Ints.end() ) {
    std::cout << "GetInt no value " << key << "\n";
    return 0;
  } else {
      return elem->second ;
  }
  
}

float CParamsManager::GetFloat (  const std::string &key ) const {
  std::map<std::string,float>::const_iterator elem = Floats.find( key );
  if ( elem == Floats.end() ) {
    std::cout << "GetFloat no value " << key << "\n";
    return 0;
  } else {
    return elem->second ;
  }

}
