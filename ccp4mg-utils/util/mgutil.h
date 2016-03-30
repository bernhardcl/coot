/*
     util/mgutil.h: CCP4MG Molecular Graphics Program
     Copyright (C) 2001-2005 University of York, CCLRC

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


#ifndef _CCP4MG_UTIL_
#define _CCP4MG_UTIL_
#include <string>
#include <algorithm>
std::string IntToString(int fn_isize,const char *format="%d");
//std::string IntToString(int fn_isize);
std::string FloatToString(float fn_isize,const char *format="%f");
//std::string FloatToString(float fn_isize);
std::string ToUpper(const std::string &in);
// trim from start
static inline std::string &ccp4mg_ltrim(std::string &s) {
        s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
        return s;
}

// trim from end
static inline std::string &ccp4mg_rtrim(std::string &s) {
        s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
        return s;
}

// trim from both ends
static inline std::string &ccp4mg_trim(std::string &s) {
        return ccp4mg_ltrim(ccp4mg_rtrim(s));
}

#endif
