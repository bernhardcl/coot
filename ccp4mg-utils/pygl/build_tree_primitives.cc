/*
     pygl/build_tree_primitives.cc: CCP4MG Molecular Graphics Program
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

#ifdef _WIN32
#include <windows.h>
#endif

#include <iostream>
#include <stdlib.h>
#include "mgtree.h"
#include "mgutil.h"
#include "cdisplayobject.h"
#include "cprimitive.h"
#include "CParamsManager.h"
#include "cartesian.h"
#include <mman_manager.h>
#include "cbuild.h"
#include "rgbreps.h"
#include "help_globals.h"
#include "catmull.h"
#include <vector>
#include <utility>
#include <algorithm>
#include <map>
#include <math.h>
#include <string>
#include <string.h>
#include "cartesian.h"
#include "texture.h"
#include "matrix.h"
#include <stdlib.h>
#include "connect.h"
#include "splineinfo.h"
#include <mmut_connectivity.h>
#include <mmut_basepairs.h>
#include <mmut_lipids.h>
#include <mmut_util.h>
#include <mmdb2/mmdb_manager.h>
#include <atom_util.h>
#include <sstream>
#include <iostream>
#include <iomanip>
#include "atom_texture_coords.h"

#include "mmut/mmut_hbond.h"

#ifndef M_PI
#define M_PI 3.141592653589793238462643
#endif
#define PIBY2 (M_PI * 2)

using namespace mmdb;

enum enum_SecStr { NOSECSTR, BETA, BULGE, TURN3, TURN4, TURN5, ALPHA }; // We want a #include enum_secstr from mmut_sec ....

ClickedLine FindLine(Displayobject &obj, const std::vector<Cartesian> &primorigin, const std::vector<Cartesian> &xyzbox, const Connectivity &conn, int symmetry){
  std::vector<std::vector<std::vector<int> > > conn_order_lists = conn.GetConnectivityLists();
  std::vector<std::vector<int> > conn_lists = conn_order_lists[0];
  return FindLine(obj, primorigin, xyzbox, conn_lists, symmetry);
}

ClickedLine FindLine(Displayobject &obj, const std::vector<Cartesian> &primorigin, const std::vector<Cartesian> &xyzbox, const std::vector<std::vector<int> > &conn_lists, int symmetry){
  /* Find lines nearest to, eg., a mouse click. */

  int clicked_symm = -1;

  std::vector<int> lines;
  lines.push_back(-1);
  lines.push_back(-1);

  matrix objrotmat = obj.quat.getInvMatrix();

  Cartesian front = Cartesian::MidPoint(objrotmat*xyzbox[2],objrotmat*xyzbox[6]);
  Cartesian back  = Cartesian::MidPoint(objrotmat*xyzbox[3],objrotmat*xyzbox[7]);

  std::vector<Plane> planes;
  std::vector<Cartesian> points;

  //planes.push_back(Plane(objrotmat*xyzbox[0],objrotmat*xyzbox[2],objrotmat*xyzbox[6])); // Front clipping plane
  //planes.push_back(Plane(objrotmat*xyzbox[1],objrotmat*xyzbox[7],objrotmat*xyzbox[3])); // Back clipping plane

  Volume v = GetClippingPlanes();
  planes.push_back(v.GetPlane(5));
  planes.push_back(v.GetPlane(4));
  points.push_back(planes[0].find_points_on_plane()[0]);
  points.push_back(planes[1].find_points_on_plane()[0]);

  //points.push_back(objrotmat*xyzbox[0]);
  //points.push_back(objrotmat*xyzbox[0]);

  double mindist = 1.0e+8;
  Cartesian prim0;
  Cartesian prim1;

  int i = 0;
  std::vector<std::vector<int> >::const_iterator conn_iter = conn_lists.begin();

  std::vector<Cartesian> unconnected;
  std::vector<int> unconnected_map;

  int nsym = obj.GetNumSymmetryMatrices();
  if(!symmetry) nsym = 0;
  int isym = 0;

  while(conn_iter!=conn_lists.end()){
    std::vector<int>::const_iterator k=conn_iter->begin();
    /* If k has no connections then we have a problem */
    if(k==conn_iter->end()){
      unconnected.push_back(primorigin[i]);
      unconnected_map.push_back(i);
    }
    while(k!=conn_iter->end()){
      isym = -1;
      do{
      prim0 = primorigin[i];
      prim1 = primorigin[*k];
      if(nsym>0&&isym>-1){
	matrix T = obj.GetSymmetryMatrix(isym);
	prim0 = T*prim0;
	prim1 = T*prim1;
      }
      std::vector<Plane>::iterator plane = planes.begin();
      std::vector<Cartesian>::const_iterator point = points.begin();
      int in_clip_planes0 = 1;
      int in_clip_planes1 = 1;
      while(plane!=planes.end()){
        Cartesian n = plane->get_normal();
        n.normalize();
        Cartesian p2prim = *point-prim0;
        p2prim.normalize();
        if(n.DotProduct(n,p2prim)>1e-3)
	  in_clip_planes0 = 0;
        Cartesian p2prim2 = *point-prim1;
        p2prim2.normalize();
        if(n.DotProduct(n,p2prim2)>1e-3)
	  in_clip_planes1 = 0;
        point++;
        plane++;
      }
      if(in_clip_planes0||in_clip_planes1){
	std::vector<double> linedist = DistanceBetweenTwoLines(front,back,prim0,prim1);
	double dist = fabs(linedist[0]);
	double u = linedist[2];
        if(dist<1.5){
        }
	if(u<0.25&&u>-0.25&&dist<0.5&&dist<mindist&&in_clip_planes0){
          mindist = dist;
	  lines[0] = i;
	  lines[1] = -1;
	  clicked_symm = isym;
	}
	if(u>0.75&&u<1.25&&dist<0.5&&dist<mindist&&in_clip_planes1){
          mindist = dist;
	  lines[0] = -1;
	  lines[1] = *k;
	  clicked_symm = isym;
	}
        if(u>=0.25&&u<=0.75&&dist<0.5&&dist<mindist&&in_clip_planes0&&in_clip_planes1) {
          mindist = dist;
	  lines[0] = i;
	  lines[1] = *k;
	  clicked_symm = isym;
        }
      }
      isym++;
      }while(isym<nsym);
      k++;
    }
    conn_iter++; i++;
  }

  ClickedLine nearprim = FindPoint(obj,unconnected,xyzbox,symmetry);
  if(nearprim.first>-1){
    Cartesian prim = primorigin[unconnected_map[nearprim.first]];
    if(nearprim.symm>-1){
      matrix T = obj.GetSymmetryMatrix(nearprim.symm);
      prim = T*prim;
    }
    std::vector<double> linedisttmp = DistanceBetweenTwoLines(front,back,prim,prim);
    double dist = fabs(linedisttmp[0]);
    if(dist<mindist){
      lines[0] = unconnected_map[nearprim.first];
      lines[1] = -1;
      mindist  = dist;
      clicked_symm = nearprim.symm;
    }
  }

  //cout << "Minimum distance: " << mindist << ", between: " << lines[0] << ", " << lines[1] << "\n";
  ClickedLine line;
  line.first  = lines[0];
  line.second = lines[1];
  line.dist   = mindist;
  line.symm = clicked_symm;
  return line;
}

ClickedLine FindLine(Displayobject &obj, const std::vector<SimpleConnection> &conn, const std::vector<Cartesian> &xyzbox, int symmetry){

  std::vector<int> lines;
  int clicked_symm = -1;

  lines.push_back(-1);
  lines.push_back(-1);

  matrix objrotmat = obj.quat.getInvMatrix();
  Cartesian dum;
  Cartesian front = dum.MidPoint(objrotmat*xyzbox[2],objrotmat*xyzbox[6]);
  Cartesian back  = dum.MidPoint(objrotmat*xyzbox[3],objrotmat*xyzbox[7]);

  std::vector<Plane> planes;
  std::vector<Cartesian> points;

  //planes.push_back(Plane(objrotmat*xyzbox[0],objrotmat*xyzbox[2],objrotmat*xyzbox[6])); // Front clipping plane
  //planes.push_back(Plane(objrotmat*xyzbox[1],objrotmat*xyzbox[7],objrotmat*xyzbox[3])); // Back clipping plane
  points.push_back(objrotmat*xyzbox[0]);
  points.push_back(objrotmat*xyzbox[1]);

  Volume v = GetClippingPlanes();
  planes.push_back(v.GetPlane(5));
  planes.push_back(v.GetPlane(4));

  double mindist = 1.0e+8;
  Cartesian prim0;
  Cartesian prim1;

  int i = 0;
  std::vector<SimpleConnection>::const_iterator conn_iter = conn.begin();

  int nsym = obj.GetNumSymmetryMatrices();
  if(!symmetry) nsym = 0;
  int isym = 0;

  while(conn_iter!=conn.end()){
      isym = -1;
      do{
      prim0 = conn_iter->first;
      prim1 = conn_iter->second;
      if(nsym>0&&isym>-1){
	matrix T = obj.GetSymmetryMatrix(isym);
	prim0 = T*prim0;
	prim1 = T*prim1;
      }
      std::vector<Plane>::iterator plane = planes.begin();
      std::vector<Cartesian>::const_iterator point = points.begin();
      int in_clip_planes = 1;
      while(plane!=planes.end()){
        Cartesian n = plane->get_normal();
        n.normalize();
        Cartesian p2prim = *point-prim0;
        p2prim.normalize();
        if(n.DotProduct(n,p2prim)>0.0)
	  in_clip_planes = 0;
        p2prim = *point-prim1;
        p2prim.normalize();
        if(n.DotProduct(n,p2prim)>0.0)
	  in_clip_planes = 0;
        point++;
        plane++;
      }
      if(in_clip_planes){
	std::vector<double> linedist = DistanceBetweenTwoLines(front,back,prim0,prim1);
	double dist = linedist[0];
	double u = linedist[2];
	if(u<0.25&&u>-0.25&&dist<0.5&&dist<mindist){
          mindist = dist;
	  lines[0] = i;
	  lines[1] = -1;
	  clicked_symm = isym;
	}
	if(u>0.75&&u<1.25&&dist<0.5&&dist<mindist){
          mindist = dist;
	  lines[0] = -1;
	  lines[1] = i;
	  clicked_symm = isym;
	}
        if(u>=0.25&&u<=0.75&&dist<0.5&&dist<mindist) {
          mindist = dist;
	  lines[0] = i;
	  lines[1] = i;
	  clicked_symm = isym;
        }
      }
      isym++;
      }while(isym<nsym);
    conn_iter++; i++;
  }

  ClickedLine line;
  line.first  = lines[0];
  line.second = lines[1];
  line.dist   = mindist;
  line.symm = clicked_symm;
  return line;
}



void DrawSimpleConnection(Displayobject &obj, const std::vector<SimpleConnection> &conn, const std::vector<double> &colour, int style, int width, int labelstyle, const std::string &labelcolour,
  const std::string &family,  const std::string &weight, 
  const std::string &slant, const std::string &size){
  std::vector<double> defcol;
  if(labelcolour!=""&&labelcolour!="default"&&labelcolour!="complement"){
       defcol = RGBReps::GetColour(labelcolour);
  }
  DrawSimpleConnectionColourVec(obj,conn,colour,style,width,labelstyle,defcol,
		  family,weight,slant,size);
}

void DrawSimpleConnectionColourVec(Displayobject &obj, const std::vector<SimpleConnection> &conn, const std::vector<double> &colour, int style, int width, int labelstyle, const std::vector<double>&labelcolf,
  const std::string &family,  const std::string &weight, 
  const std::string &slant, const std::string &size){
  double col[3] ={double(colour[0]),double(colour[1]),double(colour[2])};
  double alpha = double(colour[3]);
  int cylinders_accu = 8;
  double cylinder_radius_scale = 0.02;

  std::vector<Cartesian> carts(2);
  Text *text;

  LineCollection *lines = new LineCollection();
  PolyCollection *polys = new PolyCollection();
  DashLinesCollection *dash_lines_collection = new DashLinesCollection();
  CylinderCollection *cylinders = new CylinderCollection();
  DashCylinderCollection *dash_cylinders = new DashCylinderCollection();

  bool have_lines = false;
  bool have_polys = false;
  std::vector<SimpleConnection>::const_iterator i = conn.begin();
  while(i!=conn.end()){
    carts[0] = i->first;
    carts[1] = i->second;
    // style == NOLINE => do nothing
    if(style==DASHLINE){
      dash_lines_collection->add(carts,col,0);
      //DashLineElement *line;
      //line = new DashLineElement(carts,col,carts[0],double(width),alpha);
      //lines->add_primitive(line);
      have_lines=true;
    }
    if(style==LINE){
      LineElement *line;
      line = new LineElement(carts,col,carts[0],double(width),alpha);
      lines->add_primitive(line);
      have_lines=true;
    }
    if(style==ARROW){
      Arrow *line;
      line = new Arrow(carts,col,carts[0],double(width),alpha);
      lines->add_primitive(line);
      have_lines=true;
    }
    if(style==DASHARROW){
      DashArrow *line;
      line = new DashArrow(carts,col,carts[0],double(width),alpha);
      lines->add_primitive(line);
      have_lines=true;
    }
    if(style==DASHCYLINDER){
      DashCylinderElement *line;
      line = new DashCylinderElement(carts,col,carts[0],double(width)*cylinder_radius_scale,
                                  alpha,cylinders_accu);
      line->SetDashLength(0.2);
      line->SetDashEnd(1);
      dash_cylinders->add_primitive(line);
      have_polys=true;
    }
    if(style==CYLINDER){
      Cylinder *line;
      line = new Cylinder(carts,col,carts[0],double(width)*cylinder_radius_scale,cylinders_accu);
      cylinders->add_primitive(line);
      have_polys=true;
    }
    if(style==CYLINDERARROW){
      Cartesian p = 0.7 * carts[1] + 0.3 * carts[0];
      Cartesian p0 = carts[0];
      carts[0] = p;
      Cone *cone;
      cone = new Cone(carts,col,carts[0],double(width)*cylinder_radius_scale*2.0,alpha,cylinders_accu);
      polys->add_primitive(cone);
      carts[0] = p0;
      carts[1] = p;
      Cylinder *line;
      line = new Cylinder(carts,col,carts[0],double(width)*cylinder_radius_scale,alpha,cylinders_accu);
      polys->add_primitive(line);
      have_polys=true;
    }
    i++;
  }


  lines->SetSize((double)width);
  dash_lines_collection->SetSize((double)width);
  dash_lines_collection->SetDashLength(0.2);
  dash_cylinders->SetDashLength(0.2);
  if(have_lines)obj.add_primitive(lines);
  if(have_polys)obj.add_primitive(polys);
  if(have_polys)obj.add_primitive(cylinders);
  if(have_polys)obj.add_primitive(dash_cylinders);
  if(have_lines)obj.add_primitive(dash_lines_collection);
  // Add text label
  if (labelstyle == NOTLABELLED) return;

  i = conn.begin();

  while(i!=conn.end()){
    /*
    if (labelcolf.size() != 0) 
     label = "<colour=\""+labelcolour+"\">"+ i->label+"</colour>";
    else
    */
    if ( labelstyle == LABELLEDCENTRE )
      text = new Text ( i->first.MidPoint(i->first,i->second) , i->label,
              i->first.MidPoint(i->first,i->second));
    else if ( labelstyle == LABELLEDSTART)
      text = new Text ( i->first , i->label,i->first);
    else if ( labelstyle == LABELLEDEND)
      text = new Text ( i->second ,i->label,i->second);
    else
      text = new Text ( i->first , i->label,i->first);

    text->SetFontSize(18);
    if(family!="") text->SetFontFamily(family);
    if(weight!="") text->SetFontWeight(weight);
    if(slant!="") text->SetFontSlant(slant);
    if(size!="") text->SetFontSize(atoi(size.c_str()));
    if(labelcolf.size()!=0){
       text->SetColour(labelcolf[0],labelcolf[1],labelcolf[2],1.0);
    } else {
       text->SetDefaultColour();
    }
    //text->initialize();
    obj.add_text_primitive(text);
    i++;
  }
}


void DrawSimpleConnectionTags(Displayobject &obj,
  const std::vector<SimpleConnection> &conn,
  const std::vector<double> &colour, int style, int width,
  int labelstyle, const std::string &labelcolour,
  const std::string &family,  const std::string &weight, 
  const std::string &slant, const std::string &size,
  const std::vector<int> &tags,const std::vector<int> &selTags){

  std::vector<double> labelcolf = RGBReps::GetColour(labelcolour);
  DrawSimpleConnectionColourVecTags(obj,conn,colour,style,width,labelstyle,labelcolf,
		  family,weight,slant,size,tags,selTags);

}

void DrawSimpleConnectionColourVecTags(Displayobject &obj,
  const std::vector<SimpleConnection> &conn,
  const std::vector<double> &colour, int style, int width,
  int labelstyle, const  std::vector<double> &labelcolf,
  const std::string &family,  const std::string &weight, 
  const std::string &slant, const std::string &size,
  const std::vector<int> &tags,const std::vector<int> &selTags){

  double col[3] ={double(colour[0]),double(colour[1]),double(colour[2])};
  double alpha = double(colour[3]);
  int cylinders_accu = 8;
  double cylinder_radius_scale = 0.02;
  bool apply_selection = false;
  std::string label;
  std::vector<Cartesian> carts(2);
  Text *text;

  std::string labelcolr = FloatToString(floor(labelcolf[0])*255,"%x");
  std::string labelcolg = FloatToString(floor(labelcolf[1])*255,"%x");
  std::string labelcolb = FloatToString(floor(labelcolf[2])*255,"%x");
  std::string labelcolour = "#" + labelcolr + labelcolg + labelcolb;

  if (selTags.size()>=1)apply_selection = true; 
  std::vector<SimpleConnection>::const_iterator i = conn.begin();
  std::vector<int>::const_iterator j = tags.begin();
  while(i!=conn.end()){
    if (!apply_selection || 
         std::find(selTags.begin(),selTags.end(),*j)!=selTags.end()) {

      carts[0] = i->first;
      carts[1] = i->second;
      if(style==DASHLINE){
        DashLine *line;
        line = new DashLine(carts,col,carts[0],double(width),alpha);
        obj.add_primitive(line);
      }
      if(style==LINE){
        Line *line;
        line = new Line(carts,col,carts[0],double(width),alpha);
        obj.add_primitive(line);
      }
      if(style==ARROW){
        Arrow *line;
        line = new Arrow(carts,col,carts[0],double(width),alpha);
        obj.add_primitive(line);
      }
      if(style==DASHARROW){
        DashArrow *line;
        line = new DashArrow(carts,col,carts[0],double(width),alpha);
        obj.add_primitive(line);
      }

      if(style==DASHCYLINDER){
        DashCylinderElement *line;
        line = new DashCylinderElement(carts,col,carts[0],double(width)*cylinder_radius_scale,
                                  alpha,cylinders_accu);
        line->SetDashLength(0.2);
        line->SetDashEnd(1);
        obj.add_primitive(line);
      }
 

      if(style==CYLINDER){
        Cylinder *line;
        line = new Cylinder(carts,col,carts[0],double(width)*cylinder_radius_scale,alpha,cylinders_accu);
        obj.add_primitive(line);
      }
      if(style==CYLINDERARROW){
        Cartesian p = 0.7 * carts[1] + 0.3 * carts[0];
        Cartesian p0 = carts[0];
        carts[0] = p;
        Cone *cone;
        cone = new Cone(carts,col,carts[0],double(width)*cylinder_radius_scale*2.0,alpha,cylinders_accu);
        obj.add_primitive(cone);
        carts[0] = p0;
        carts[1] = p;
        Cylinder *line;
        line = new Cylinder(carts,col,carts[0],double(width)*cylinder_radius_scale,alpha,cylinders_accu);
        obj.add_primitive(line);
      }
    }
    i++;
    j++;
  }

  // Add text label
  if (labelstyle == NOTLABELLED) return;

  i = conn.begin();
  j = tags.begin();
  while(i!=conn.end()){
    if (!apply_selection || 
           std::find(selTags.begin(),selTags.end(),*j)!=selTags.end()) {
      if (labelcolour != "") 
        label = "<colour=\""+labelcolour+"\">"+ i->label+"</colour>";
      else
        label =  i->label;
      if ( labelstyle == LABELLEDCENTRE )
        text = new Text ( i->first.MidPoint(i->first,i->second) , label, i->first.MidPoint(i->first,i->second));
      else if ( labelstyle == LABELLEDSTART)
        text = new Text ( i->first , label,i->first);
      else if ( labelstyle == LABELLEDEND)
        text = new Text ( i->second ,label,i->second);
      else 
        text = new Text ( i->first , label,i->first);
  
      obj.add_text_primitive(text);
      text->SetFontFamily(family);
      text->SetFontWeight(weight);
      text->SetFontSlant(slant);
      text->SetFontSize(atoi(size.c_str()));
    }
    i++;
    j++;
  }

}

std::vector<int> GetPointsInVolume(Displayobject &obj, const std::vector<Cartesian> &atoms, const Volume &volume){

  Cartesian atom_origin;
  int clicked;
  Plane plane;
  std::vector <Cartesian> points;
  Cartesian n;
  Cartesian p2atom;
  std::vector<int> clicked_atoms;

  for(unsigned int j=0;j<atoms.size();j++){
     clicked = 1;
     atom_origin = atoms[j];
     matrix mat = obj.quat.getMatrix();
     atom_origin = mat*atom_origin;
     atom_origin += obj.origin;
     for(int ii=0;ii<volume.GetNumberOfPlanes();ii++){
       plane = volume.GetPlane(ii);
       n = plane.get_normal();
       points = plane.find_points_on_plane();
       p2atom = points[0] - atom_origin;
       n.normalize();
       p2atom.normalize();
       if(n.DotProduct(n,p2atom)<0.0)
         clicked = 0;
     }
     if(clicked) clicked_atoms.push_back(j);
  }

  return clicked_atoms;

}

ClickedLine FindPoint(Displayobject &obj, const std::vector<Cartesian> &primorigin, const std::vector<Cartesian> &xyzbox, int symmetry){

  int clicked_symm = -1;
  int nearprim = findprimc(xyzbox,primorigin,obj.origin,obj.quat.getMatrix());
  std::vector<Cartesian> primorigin2 = primorigin;

  unsigned int nsym = obj.GetNumSymmetryMatrices();
  if(!symmetry) nsym = 0;

  if(nearprim==-1){
    for(unsigned int j=0;j<nsym;j++){
      matrix T = obj.GetSymmetryMatrix(j);
      primorigin2.clear();
      for(unsigned int i=0;i<primorigin.size();i++){
        Cartesian cart = primorigin[i];
        cart = T*cart;
        primorigin2.push_back(cart);
      }
      nearprim = findprimc(xyzbox,primorigin2,obj.origin,obj.quat.getMatrix());
      //std::cout << "For symmetry " << j << " found " << nearprim << std::endl;
      if(nearprim>-1){
	 clicked_symm = j;
         break;
      }
    }
  }

  ClickedLine line;
  line.first  = nearprim;
  line.second = -1;
  line.dist   = -1.0; // Calculate this later...
  line.symm = clicked_symm;

  return line;

}          

int *GetTextIDS(Displayobject &obj){
  return obj.GetTextIDS();
}

int GetNumberOfTextIDS(Displayobject &obj){
  return obj.GetNumberOfTextIDS();
}

void SetTextString(Displayobject &obj,int text_id, const char* new_string){
  obj.SetTextString(text_id,new_string);
}

void SetTextString(Displayobject &obj,int text_id, const std::string &new_string){
  obj.SetTextString(text_id,new_string);
}

const char* GetTextString(Displayobject &obj,int text_id){
  return obj.GetTextString(text_id);
}

void DeleteTextLabel(Displayobject &obj, int text_id){
  obj.DeleteTextPrimitive(text_id);
}

int AddTextLabel(Displayobject &obj, double x, double y, double z, const std::string &label){
  Text *text;
  Cartesian primorigin = Cartesian(x,y,z,1.0);
  text = new Text(primorigin,label,primorigin);
  obj.add_text_primitive(text);
  return text->GetID();
}

int AddTextLabel(Displayobject &obj, double x, double y, double z, const char *label){
  int newtextid = AddTextLabel(obj,x,y,z,std::string(label));
  return newtextid;
}

void AddBillBoardTextLabel(Displayobject &obj, double x, double y, const std::string &label){
  BillBoardText *text;
  Cartesian primorigin = Cartesian(x,y,0);
  text = new BillBoardText(primorigin,label,primorigin);
  obj.add_text_primitive(text);
}

void AddBillBoardTextLabel(Displayobject &obj, double x, double y, const char *label){
  AddBillBoardTextLabel(obj,x,y,std::string(label));
}

void FitToPolynomial(std::vector<Cartesian> &carts, int pass){
  std::vector <Cartesian> spline = SplineCurve(carts,(carts.size()-1)*4,2,pass);
  for(unsigned ii=1;ii<carts.size()-2;ii++){
    carts[ii] = spline[ii*(spline.size()+1)/(carts.size()-1)];
   }
}


void build_beta_surface(CMMANManager *molH, int atom_selHnd_in, Displayobject &obj, const CParamsManager &params, const AtomColourVector &atom_colour_vector){


  std::cout << "build_beta_surface\n";
  int sec_str_mask[] = {1,1,1,0,0,0,0,0};
  std::string sec_strucs = molH->ListSecStructure(sec_str_mask);

  std::cout << sec_strucs << "\n";
  std::cout << "Done build_beta_surface\n";
  return;

  int CAselHnd;
  mmdb::Atom** atomTable;
  int nAtoms;

  int atom_selHnd = molH->NewSelection();
  molH->Select(atom_selHnd,STYPE_ATOM,0,"*",ANY_RES,"*",ANY_RES,"*","*","*","*","*",SKEY_NEW);
  molH->Select(atom_selHnd,STYPE_ATOM,atom_selHnd_in,SKEY_AND);
  molH->ExcludeOverlappedAtoms(atom_selHnd,0.8);

  // Find all CA - to use as quick check if atom is in this set
  CAselHnd = molH->NewSelection();
  molH->Select(CAselHnd,STYPE_ATOM,0,"*",ANY_RES,"*",ANY_RES,"*","*","CA","C","*",SKEY_NEW);
  molH->Select(CAselHnd,STYPE_ATOM,atom_selHnd_in,SKEY_AND);
  molH->ExcludeOverlappedAtoms(CAselHnd,0.8);
  molH->GetSelIndex ( atom_selHnd, atomTable, nAtoms );
  std::vector<Cartesian>  cavertices;
  const double *colour_array=0;
  double red[] = {1.0,0.0,0.0,1.0};

  /* We won't try anything fancy with colours just yet. */
  //if(!atom_colour_vector)
     colour_array = red;

  double min_x = 1e+8;
  double min_y = 1e+8;
  double max_x = 1e-8;
  double max_y = 1e-8;

  std::cout << "\n";
  //PolyCollection *polys = new PolyCollection();
  for(int j=0;j<nAtoms;j++){
    if(atomTable[j]->isInSelection(CAselHnd) && molH->isAminoacid(atomTable[j]->residue)){
      mmdb::Atom* pCA = atomTable[j];
      mmdb::Residue* pRes = pCA->residue;
      // Save the CA pointer
      if(int(pRes->SSE)== SSE_Strand || int(pRes->SSE)== SSE_Bulge){
        cavertices.push_back(Cartesian(pCA->x,pCA->y,pCA->z));
        //SphereElement *sphere = new SphereElement(cavertices.back(),colour_array,cavertices.back(),0.4,1.0,2);
        //polys->add_primitive(sphere);
        if(pCA->x<min_x) min_x = pCA->x;
        if(pCA->y<min_y) min_y = pCA->y;
        if(pCA->x>max_x) max_x = pCA->x;
        if(pCA->y>max_y) max_y = pCA->y;
        std::cout << cavertices.back() << "\n";
        //if( atom_colour_vector ){
          //colour_array = atom_colour_vector->GetRGB(j);
        //}
      }
    }
  }
  //obj.add_primitive(polys);

  std::cout << "\n";
  std::cout << min_x << " " << max_x << "\n";
  std::cout << min_y << " " << max_y << "\n";
  std::cout << "\n";

  //std::cout << cavertices.size() << "\n";
  if(cavertices.size()>5){//Need at least 6 sets of coords to satisfy our 6 unknown coeffs.
    std::vector<Cartesian> carts(2);
    double width = 2.0;
    LineCollection *lines = new LineCollection();
    min_x -= 3;
    min_y -= 3;
    max_x += 3;
    max_y += 3;
    std::vector<double> poly_params = LeastSquaresQuadraticFit3D(cavertices);
    //std::cout << "Draw function " << poly_params[0] << "x^2 + " << poly_params[1] << "y^2 + " << poly_params[2] << "xy + " << poly_params[3] << "x + " << poly_params[4] << "y + " << poly_params[5] << "\n";
    //std::cout << "In range: " << min_x << " -> " << max_x << ", " << min_y << " -> " << max_y << "\n";
    double x = min_x;
    double delta_x = (max_x-min_x)/30.0;
    double delta_y = (max_y-min_y)/30.0;
    while(x<max_x){
      double y = min_y;
      while(y<max_y){
        double z = poly_params[0]*x*x + poly_params[1]*y*y + poly_params[2]*x*y + poly_params[3]*x + poly_params[4]*y + poly_params[5];
        carts[0] = Cartesian(x,y,z);
        z = poly_params[0]*x*x + poly_params[1]*(y+delta_y)*(y+delta_y) + poly_params[2]*x*(y+delta_y) + poly_params[3]*x + poly_params[4]*(y+delta_y) + poly_params[5];
        carts[1] = Cartesian(x,y+delta_y,z);
        //std::cout << carts[1].get_x()-carts[0].get_x() << " " << carts[1].get_y()-carts[0].get_y() << "\n";
        LineElement *line = new LineElement(carts,colour_array,carts[0],width,1.0);
        lines->add_primitive(line);
        z = poly_params[0]*(x+delta_x)*(x+delta_x) + poly_params[1]*y*y + poly_params[2]*(x+delta_x)*y + poly_params[3]*(x+delta_x) + poly_params[4]*y + poly_params[5];
        carts[1] = Cartesian(x+delta_x,y,z);
        line = new LineElement(carts,colour_array,carts[0],width,1.0);
        lines->add_primitive(line);
        y += delta_y;
      }
      x += delta_x;
    }
    lines->SetSize(width);
    obj.add_primitive(lines);
  }

}

void build_spline(const SplineInfo &splineinfo, Displayobject &obj, int mode, const CParamsManager &params,  const CParamsManager &global_params, const std::string &texture, const std::string &bumpmap){

#ifdef _DO_TIMINGS_
  clock_t t1 = clock();
#endif
  unsigned int i;
  int ribbon_accus[] = {6, 9, 12, 18, 30, 36};

  //std::cout << "into build_spline" << std::endl;

  int multicolour = 1;
  int spline_accu = 4+4*global_params.GetInt("solid_quality");
  

  if(bumpmap!=""&&mode!=BONDS&&mode!=FATBONDS&&mode!=THINBONDS){
    image_info iinfo = image_info(bumpmap.c_str());
    load_texture(iinfo,MIPMAP);
    glEnable(GL_TEXTURE_2D);
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
    glEnable(GL_TEXTURE_GEN_S);
    glEnable(GL_TEXTURE_GEN_T);
    glTexGeni(GL_S, GL_TEXTURE_GEN_MODE, GL_SPHERE_MAP);
    glTexGeni(GL_T, GL_TEXTURE_GEN_MODE, GL_SPHERE_MAP);
  }


  std::vector<std::vector<Cartesian> >exp_atomColourVector;
  if(multicolour&&(int)splineinfo.colours.size()>0) {
    for(i=0;i<splineinfo.colours.size();++i){
      exp_atomColourVector.push_back(std::vector<Cartesian>(0));
      for(int j=0;j<int(splineinfo.colours[i].size())-1;++j){
	for(int k=0;k<spline_accu;++k){
          /* Need to be more clever since multicolour is too much of a catchall */
          double frac = 0.0;// double(k)/spline_accu;
          Cartesian frac_col = (1.0-frac) * splineinfo.colours[i][j] + frac * splineinfo.colours[i][j+1];
	  exp_atomColourVector[i].push_back(frac_col);
        }
      }
      for(int k=0;k<spline_accu;++k){
        double frac = 0.0;// double(k)/spline_accu;
        Cartesian frac_col = frac * splineinfo.colours[i][splineinfo.colours[i].size()-1] + (1.0-frac) * splineinfo.colours[i][splineinfo.colours[i].size()-2];
        exp_atomColourVector[i].push_back(frac_col);
      } 
    }
  }

  std::vector<Cartesian> cartesians;
  std::vector<Cartesian> n1_cartesians;
  std::vector<Cartesian> n2_cartesians;


  std::vector<std::vector<Cartesian> >::const_iterator splines_iter=splineinfo.splines.begin();
  std::vector<Cartesian>::const_iterator spline_iter;
  std::vector<std::vector<Cartesian> >::const_iterator n1_splines_iter=splineinfo.n1_splines.begin();
  std::vector<Cartesian>::const_iterator n1_spline_iter;
  std::vector<std::vector<Cartesian> >::const_iterator n2_splines_iter=splineinfo.n2_splines.begin();
  std::vector<Cartesian>::const_iterator n2_spline_iter;
  std::vector<std::vector<Cartesian> >::const_iterator colour_vecs=exp_atomColourVector.begin();
  std::vector<Cartesian>::const_iterator colour_vec;
  std::vector<Cartesian> colours_vec;

  std::vector<std::vector<Cartesian> >::const_iterator nasplines_iter=splineinfo.nasplines.begin();
  std::vector<std::vector<Cartesian> >::const_iterator n1_nasplines_iter=splineinfo.n1_nasplines.begin();
  std::vector<std::vector<Cartesian> >::const_iterator n2_nasplines_iter=splineinfo.n2_nasplines.begin();
  std::vector<std::vector<Cartesian> > naexp_atomColourVector;
  if(multicolour&&(int)splineinfo.nacolours.size()>0) {
    for(i=0;i<splineinfo.nacolours.size();++i){
      naexp_atomColourVector.push_back(std::vector<Cartesian>(0));
      for(unsigned int j=0;j<splineinfo.nacolours[i].size();++j)
	for(int k=0;k<spline_accu;++k)
	   naexp_atomColourVector[i].push_back(splineinfo.nacolours[i][j]);
    }
  }
  std::vector<std::vector<Cartesian> >::const_iterator nacolours_iter=naexp_atomColourVector.begin();
  while(nasplines_iter!=splineinfo.nasplines.end()){
    if(nasplines_iter->size()>2){
    float worm_width = params.GetFloat("worm_width");
    float ribbon_width = params.GetFloat("ribbon_width");
    int ribbon_style = params.GetInt("ribbon_style");
    int ribbon_accu = ribbon_accus[global_params.GetInt("solid_quality")];
    int spline_accu = 4+4*global_params.GetInt("solid_quality");
    if (mode == SPLINE|| mode == FATWORM){
      double *col_tmp = RGBReps::GetColourP(1);
      Ribbon *ribbon = new Ribbon(*nasplines_iter,*n1_nasplines_iter,*n2_nasplines_iter,col_tmp,(*nasplines_iter)[0],*nacolours_iter,worm_width,ribbon_width,worm_width,1.0,ribbon_accu*2,spline_accu,0,ribbon_style);
      obj.add_primitive(ribbon);
      delete [] col_tmp;
    }else{
      double *col_tmp = RGBReps::GetColourP(1);
      Ribbon *ribbon = new Ribbon(*nasplines_iter,*n1_nasplines_iter,*n2_nasplines_iter,col_tmp,(*nasplines_iter)[0],*nacolours_iter,worm_width,worm_width,worm_width,1.0,ribbon_accu*2,spline_accu);
      obj.add_primitive(ribbon);
      delete [] col_tmp;
    }
    }
    ++nasplines_iter;
    ++n1_nasplines_iter;
    ++n2_nasplines_iter;
    if(multicolour) ++nacolours_iter;
  }

  int totalpoints = 0;
  int nchain = 0;
  while(splines_iter!=splineinfo.splines.end()){
    std::vector<std::vector<int> > secstr_indices = splineinfo.secstr_indices[nchain];
    std::vector<std::vector<int> >::const_iterator secstr_iter=secstr_indices.begin();
    ++nchain;
    if(splines_iter->size()<2){
      ++splines_iter;
      ++n1_splines_iter;
      ++n2_splines_iter;
      if(multicolour) ++colour_vecs;
      secstr_indices = splineinfo.secstr_indices[nchain];
      secstr_iter=secstr_indices.begin();
      ++nchain;
    }
    spline_iter=(*splines_iter).begin();
    n1_spline_iter=(*n1_splines_iter).begin();
    n2_spline_iter=(*n2_splines_iter).begin();

    totalpoints += (*splines_iter).size();
    if(multicolour) colour_vec=(*colour_vecs).begin();
    //std::cout << "New chain, size: " <<  (*splines_iter).size() << "\n";
    while(secstr_iter!=secstr_indices.end()){
      if(secstr_iter<secstr_indices.end()-1) {
        //std::cout << (*secstr_iter)[0] << " to " << (*(secstr_iter+1))[0] << "\n";
        int begin = (*secstr_iter)[0]*spline_accu;
        int end = (*(secstr_iter+1))[0]*spline_accu;
        for(int i=begin;i<=end&&(*splines_iter).size()>0&&spline_iter!=(*splines_iter).end();++i){ 
	  cartesians.push_back(*spline_iter);
	  n1_cartesians.push_back(*n1_spline_iter);
	  n2_cartesians.push_back(*n2_spline_iter);
	  if(multicolour) colours_vec.push_back(*colour_vec);
	  ++n1_spline_iter;
	  ++n2_spline_iter;
	  ++spline_iter;
	  if(multicolour) ++colour_vec;
        }
        bool doOneMore = false;
	if(spline_iter+1==(*splines_iter).end()) doOneMore = true;
	if((*splines_iter).size()>0&&spline_iter!=(*splines_iter).end()&&spline_iter+1<(*splines_iter).end() && spline_iter+2<(*splines_iter).end() && (*(spline_iter+1)-*(spline_iter+2)).length()>15./spline_accu) doOneMore = true;
	if((*splines_iter).size()>0&&doOneMore&&spline_iter!=(*splines_iter).end()&&n1_spline_iter!=(*n1_splines_iter).end()&&n2_spline_iter!=(*n2_splines_iter).end()){
	  cartesians.push_back(*spline_iter);
	  n1_cartesians.push_back(*n1_spline_iter);
	  n2_cartesians.push_back(*n2_spline_iter);
	  if(multicolour&&(colour_vec!=(*colour_vecs).end())) colours_vec.push_back(*colour_vec);
        }
        //std::cout << "cartesians.size() " << cartesians.size() << "\n"; std::cout.flush();
        ++totalpoints;
        --spline_iter;
        --n1_spline_iter;
        --n2_spline_iter;
	if(multicolour) --colour_vec;
      }else{
	while((*splines_iter).size()>0&&spline_iter<(*splines_iter).end()){
	  cartesians.push_back(*spline_iter);
	  n1_cartesians.push_back(*n1_spline_iter);
	  n2_cartesians.push_back(*n2_spline_iter);
	  if(multicolour) colours_vec.push_back(*colour_vec);
	  ++n1_spline_iter;
	  ++n2_spline_iter;
	  ++spline_iter;
	  if(multicolour) ++colour_vec;
	}
        --spline_iter;
        --n1_spline_iter;
        --n2_spline_iter;
	if(multicolour) --colour_vec;
      }
      if(cartesians.size()>2) {
        if((cartesians[cartesians.size()-1]-cartesians[cartesians.size()-2]).length()<1e-4){
           cartesians.pop_back();
           n1_cartesians.pop_back();
           n2_cartesians.pop_back();
        }
        //std::cout << cartesians.size() << " " << n1_cartesians.size() << " " << n2_cartesians.size() << "\n";
        //std::cout << cartesians[0] << " " << cartesians[cartesians.size()-1] << "\n";
        //if(cartesians.size()>1) cartesians.pop_back();
        float worm_width = params.GetFloat("worm_width");
        float arrow_length = params.GetFloat("arrow_width");
        float arrow_width = params.GetFloat("arrow_width");
        int ribbon_accu = ribbon_accus[global_params.GetInt("solid_quality")];
        int ribbon_style = params.GetInt("ribbon_style");
        int helix_style = params.GetInt("helix_style");
        int spline_accu = 4+4*global_params.GetInt("solid_quality");
        bool two_colour_ribbon = params.GetInt("two_colour_ribbon");
        bool grey_ribbon_edge = params.GetInt("grey_ribbon_edge");
        float alpha_helix_width;
        float beta_sheet_width;
	float helix_tube_diameter;
        if (mode == SPLINE || mode == FATWORM) {
          alpha_helix_width = params.GetFloat("alpha_helix_width");
          helix_tube_diameter =  params.GetFloat("helix_tube_diameter");
          beta_sheet_width = params.GetFloat("alpha_helix_width");
          //beta_sheet_width = params.GetFloat("beta_sheet_width");
        } else {
          alpha_helix_width = params.GetFloat("worm_width");
          beta_sheet_width = params.GetFloat("worm_width");
        }
        if(fabs(arrow_length)<1e-2) arrow_length = 1.0;
        if(arrow_width<beta_sheet_width) arrow_width = beta_sheet_width;
        if(colours_vec.size()>1){
          colours_vec.pop_back();
          colours_vec.push_back(colours_vec.back());
        }
        if((*secstr_iter)[1]==ALPHA&&cartesians.size()>4){
          if (mode == SPLINE||mode==WORM||mode==VARIABLEWORM) {
	    Ribbon *ribbon;
            double *col_tmp = RGBReps::GetColourP((*secstr_iter)[1]+1);
            if(mode==WORM){
	      ribbon = new Ribbon(cartesians,n1_cartesians,n2_cartesians,col_tmp,cartesians[0],colours_vec,worm_width,alpha_helix_width,worm_width,1.0,ribbon_accu*2,spline_accu);
            }else if(mode==VARIABLEWORM){
	      ribbon = new VariableWorm(cartesians,n1_cartesians,n2_cartesians,col_tmp,cartesians[0],colours_vec,worm_width,alpha_helix_width,worm_width,1.0,ribbon_accu*2,spline_accu);
            }else{
	      ribbon = new Ribbon(cartesians,n1_cartesians,n2_cartesians,col_tmp,cartesians[0],colours_vec,worm_width,alpha_helix_width,worm_width,1.0,ribbon_accu*2,spline_accu,0,helix_style,two_colour_ribbon,grey_ribbon_edge);
            }
	    obj.add_primitive(ribbon);
            delete [] col_tmp;
          }else{
            double *col_tmp = RGBReps::GetColourP((*secstr_iter)[1]+1);
            col_tmp = RGBReps::GetColourP(RGBReps::GetColourNumber("black"));
            /* This is cylinder rep of alpha helices !! */
            std::vector<Cartesian> lead_in(cartesians.begin(),cartesians.begin()+1*spline_accu);
            std::vector<Cartesian> n1_lead_in(n1_cartesians.begin(),n1_cartesians.begin()+1*spline_accu);
            std::vector<Cartesian> n2_lead_in(n2_cartesians.begin(),n2_cartesians.begin()+1*spline_accu);
            std::vector<Cartesian> cols_lead_in(colours_vec.begin(),colours_vec.begin()+1*spline_accu);
	    Worm *ribbon = new Worm(lead_in,n1_lead_in,n2_lead_in,col_tmp,cartesians[0],cols_lead_in,worm_width,worm_width,worm_width,1.0,ribbon_accu*2,spline_accu);
	    obj.add_primitive(ribbon);
            if(cartesians.size()>2*spline_accu+1){
            std::vector<Cartesian> main(cartesians.begin()+1*spline_accu-1,cartesians.end()-1*spline_accu+1);
            std::vector<Cartesian> n1_main(n1_cartesians.begin()+1*spline_accu-1,n1_cartesians.end()-1*spline_accu+1);
            std::vector<Cartesian> n2_main(n2_cartesians.begin()+1*spline_accu-1,n2_cartesians.end()-1*spline_accu+1);
            std::vector<Cartesian> cols_main(colours_vec.begin()+1*spline_accu-1,colours_vec.end()-1*spline_accu+1);
	    ribbon = new Worm(main,n1_main,n2_main,col_tmp,cartesians[0],cols_main,helix_tube_diameter,helix_tube_diameter,helix_tube_diameter,1.0,ribbon_accu*2,spline_accu);
	    obj.add_primitive(ribbon);
            std::vector<Cartesian> lead_out(cartesians.end()-1*spline_accu,cartesians.end());
            std::vector<Cartesian> n1_lead_out(n1_cartesians.end()-1*spline_accu,n1_cartesians.end());
            std::vector<Cartesian> n2_lead_out(n2_cartesians.end()-1*spline_accu,n2_cartesians.end());
            std::vector<Cartesian> cols_lead_out(colours_vec.end()-1*spline_accu,colours_vec.end());
	    ribbon = new Worm(lead_out,n1_lead_out,n2_lead_out,col_tmp,cartesians[0],cols_lead_out,worm_width,worm_width,worm_width,1.0,ribbon_accu*2,spline_accu);
	    obj.add_primitive(ribbon);
            }
            delete [] col_tmp;
          }
        }else if((*secstr_iter)[1]==BETA&&cartesians.size()>4){
          double *col_tmp = RGBReps::GetColourP((*secstr_iter)[1]+1);
	  Ribbon *ribbon;
          if (mode == SPLINE|| mode == FATWORM) {
	    ribbon = new ArrowHeadRibbon(cartesians,n1_cartesians,n2_cartesians,col_tmp,cartesians[0],colours_vec,worm_width,beta_sheet_width,worm_width,arrow_length,arrow_width,1.0,ribbon_accu*2,spline_accu,0,ribbon_style,false,grey_ribbon_edge);
            obj.add_primitive(ribbon);
          } else {
            if(mode==VARIABLEWORM){
	      ribbon = new VariableWorm(cartesians,n1_cartesians,n2_cartesians,col_tmp,cartesians[0],colours_vec,worm_width,beta_sheet_width,worm_width,1.0,ribbon_accu*2,spline_accu);
            }else{
	      ribbon = new Ribbon(cartesians,n1_cartesians,n2_cartesians,col_tmp,cartesians[0],colours_vec,worm_width,beta_sheet_width,worm_width,1.0,ribbon_accu*2,spline_accu);
            }
   	   obj.add_primitive(ribbon);
	  }
          delete [] col_tmp;
        }else{
          double *col_tmp = RGBReps::GetColourP((*secstr_iter)[1]+1);
          if(mode==VARIABLEWORM){
	    VariableWorm *ribbon = new VariableWorm(cartesians,n1_cartesians,n2_cartesians,col_tmp,cartesians[0],colours_vec,worm_width,worm_width,worm_width,1.0,ribbon_accu*2,spline_accu);
	    obj.add_primitive(ribbon);
          }else{
	    Worm *ribbon = new Worm(cartesians,n1_cartesians,n2_cartesians,col_tmp,cartesians[0],colours_vec,worm_width,worm_width,worm_width,1.0,ribbon_accu*2,spline_accu);
	    obj.add_primitive(ribbon);
          }
          delete [] col_tmp;
        }
      }
      cartesians.clear();
      n1_cartesians.clear();
      n2_cartesians.clear();
      colours_vec.clear();
      //std::cout << (*secstr_iter)[0]*params.spline_accu << " " << totalpoints << "\n";
      ++secstr_iter;
    }
    ++splines_iter;
    ++n1_splines_iter;
    ++n2_splines_iter;
    if(multicolour) ++colour_vecs;
  }
 //std::cout << "done build_spline" << std::endl;
#ifdef _DO_TIMINGS_
  clock_t t2 = clock();
  std::cout << "Time for build_spline" << ((t2-t1)*1000.0/CLOCKS_PER_SEC)/1000.0<< "\n";
#endif
}

ConnectivityDraw::ConnectivityDraw(){
}

void ConnectivityDraw::SetParametersAndCalculate(const Connectivity &connectivity_in, PCMMANManager molhnd_in, Displayobject &obj, int mode, const CParamsManager &params, const CParamsManager &global_params, int nSelAtoms, const AtomColourVector &atom_colour_vector, const std::vector<double> &atomRadii, const SplineInfo &splineinfo_ribbon, const SplineInfo &splineinfo_worm, const std::string &texture, const std::string &bumpmap, int stick_colour, int side_to_ribbon, int side_to_worm, int bonds_mode){
  molhnd = molhnd_in;
  RedrawPrimitives(obj,connectivity_in,mode,params,global_params,nSelAtoms,atom_colour_vector,atomRadii,splineinfo_ribbon,splineinfo_worm,texture,bumpmap,stick_colour,side_to_ribbon,side_to_worm,bonds_mode);
}

ConnectivityDraw::ConnectivityDraw(const Connectivity &connectivity_in, PCMMANManager molhnd_in, Displayobject &obj, int mode, const CParamsManager &params, const CParamsManager &global_params, int nSelAtoms, const AtomColourVector &atom_colour_vector, const std::vector<double> &atomRadii, const SplineInfo &splineinfo_ribbon, const SplineInfo &splineinfo_worm, const std::string &texture, const std::string &bumpmap, int stick_colour, int side_to_ribbon, int side_to_worm,int bonds_mode ){

  molhnd = molhnd_in;
  RedrawPrimitives(obj,connectivity_in,mode,params,global_params,nSelAtoms,atom_colour_vector,atomRadii,splineinfo_ribbon,splineinfo_worm,texture,bumpmap,stick_colour,side_to_ribbon,side_to_worm,bonds_mode);
}

bool GetOtherChildren(int original_bonded, int original_root, int root, int bonded, const std::vector<std::vector<int> > &conn_lists, int depth, int &max_depth,std::vector<int> &current_branch,bool storeNow,int &ring_count, std::vector<std::vector<int> > &rings){
  if(bonded==-1) return false;
  std::vector<int> adjacent_children = conn_lists[bonded];
  adjacent_children.erase(std::remove(adjacent_children.begin(), adjacent_children.end(), root), adjacent_children.end());
  std::vector<int>::iterator child = adjacent_children.begin();
  if(bonded==original_root){
      //std::cout << "success ...." << depth << "\n";
      //std::cout << "original_root ...." << original_root << "\n";
      //std::cout << "original_bonded ...." << original_bonded << "\n";
      //for(int i=0;i<current_branch.size();i++) {std::cout << current_branch[i] << " ";}; std::cout << "\n";
    if(depth==max_depth&&!storeNow) {
      //std::cout << "Another ring of depth " << depth << "\n";
      ring_count++;
      rings.push_back(current_branch);
    } else if(depth<max_depth&&!storeNow) {
      ring_count=0;
      rings.clear();
      rings.push_back(current_branch);
    }
    max_depth = depth;
    return true;
  }
  if(depth>=max_depth) return false;
  while(child!=adjacent_children.end()){
    /*
    if(storeNow){
    }
    */
    //for(int i=0;i<depth;i++) std::cout << " ";
    //std::cout << *child << " (depth " << depth << ", bonded: " << bonded << ")\n";
    if((std::count(current_branch.begin(),current_branch.end(),*child)>0)||*child==original_bonded){
      //std::cout << "bombing ....\n";
      //for(int i=0;i<current_branch.size();i++) {std::cout << current_branch[i] << " ";}; std::cout << "\n";
      child++;
      continue;
    }
    current_branch[depth] = *child;
    for(int ib=depth+1;ib<current_branch.size();ib++){current_branch[ib]=-1;}
    bool cycle = GetOtherChildren(original_bonded, original_root, bonded,*child,conn_lists,depth+1,max_depth,current_branch,storeNow,ring_count,rings);
    if(cycle&&storeNow){
      //std::cout << "ring_count " << ring_count << "\n";
      //FIXME - Hmm (original_root>bonded) makes bond be on same side for both atoms when in fused rings,
      //        but causes some rings not to be identified!
      //if(ring_count==0||(original_root>bonded&&std::count(current_branch.begin(),current_branch.end(),original_bonded)<1)){
      if(ring_count==0||(std::count(current_branch.begin(),current_branch.end(),original_bonded)<1)) {
        return true;
      }
    }
    child++;
  }
  for(int ib=depth+1;ib<current_branch.size();ib++){current_branch[ib]=-1;}
  return false;
}

std::vector<int> GetFuncGroups(int root, int bonded, const std::vector<std::vector<int> > &conn_lists){
  std::vector<int> adjacent_children = conn_lists[bonded];
  adjacent_children.erase(std::remove(adjacent_children.begin(), adjacent_children.end(), root), adjacent_children.end());
  return adjacent_children;
}

std::vector<std::vector<int> > CheckRing(int root, int bonded, const std::vector<std::vector<int> > &conn_lists){
  if(bonded==-1)  return std::vector<std::vector<int> >();
  int max_depth = 7;
  std::vector<int> current_branch (max_depth);
  for(int id=0;id<max_depth;id++)
     current_branch[id] = -1;
  std::vector<int> dum;
  std::vector<int> adjacent_children = conn_lists[bonded];
  adjacent_children.erase(std::remove(adjacent_children.begin(), adjacent_children.end(), root), adjacent_children.end());

  std::vector<std::vector<int> > rings;
  int ring_count = 0;
  std::vector<int>::iterator child = adjacent_children.begin();
  while(child!=adjacent_children.end()){
    current_branch[0] = *child; for(int ib=1;ib<current_branch.size();ib++){current_branch[ib]=-1;}
    GetOtherChildren(bonded,root,bonded,*child,conn_lists,1,max_depth,current_branch,false,ring_count,rings);
    child++;
  }


  if(rings.size()>0){
    std::vector<std::vector<int> > theRings;
    std::vector<int> theRing;
    for(unsigned iring=0;iring<rings.size();iring++){
      theRing = rings[iring];
      theRing.insert(theRing.begin(),bonded);
      std::vector<int>::iterator it;
      it=std::find(theRing.begin(),theRing.end(),-1);
      theRing = std::vector<int>(theRing.begin(),it);
      theRings.push_back(theRing);
    }
    return theRings;
  } else {
    return std::vector<std::vector<int> >();
  }


  current_branch = std::vector<int> (max_depth);
  for(int id=0;id<max_depth;id++)
     current_branch[id] = -1;
  child = adjacent_children.begin();
  bool cycle = false;
  while(child!=adjacent_children.end()){
    //std::cout << *child << " (depth " << 0 << ", bonded: " << bonded << ")\n";
    current_branch[0] = *child;
    cycle = GetOtherChildren(bonded,root,bonded,*child,conn_lists,1,max_depth,current_branch,true,ring_count,rings);
    if(cycle){
      break;
    }
    child++;
  }
  if(cycle){
     current_branch.insert(current_branch.begin(),bonded);
     std::cout << "Returning with " << rings.size() << " rings\n";
     std::vector<std::vector<int> > theRings;
     theRings.push_back(current_branch);
     return theRings;
  }
  return std::vector<std::vector<int> >();
}

void DrawDashedCylinder(int stick_colour,int mode, const std::vector<Cartesian> &carts, const double *colour_array, const double *stick_colour_array, double cylinders_size, int spheres_accu, SphereCollection *spheres, int catom_index){
  Cartesian b = carts[1]-carts[0];
  double bl = b.length();
  Cartesian sp = carts[1]-0.1/bl*b;
  SphereElement *sphere = new SphereElement(sp,colour_array,carts[0],cylinders_size,1.0,spheres_accu);
  sphere->SetColourOverride(catom_index);
  spheres->add_primitive(sphere);
  double dist = 0.3;
  while(dist+.2<bl){
    sp = carts[1]-dist/bl*b;
    sphere = new SphereElement(sp,colour_array,carts[0],cylinders_size,1.0,spheres_accu);
    sphere->SetColourOverride(catom_index);
    spheres->add_primitive(sphere);
    dist += 0.2;
  }
}

void DrawCylinder(int stick_colour,int mode, const std::vector<Cartesian> &carts, const double *colour_array, const double *stick_colour_array, double cylinders_size, int cylinders_accu, CylinderCollection *cylinders, int catom_index){
          if (stick_colour > 0 && mode==BALLSTICK ) {
            CylinderElement *line = new CylinderElement(carts,stick_colour_array,carts[0],cylinders_size,1.0,cylinders_accu);
            cylinders->add_primitive(line);
          } else {
            CylinderElement *line = new CylinderElement(carts,colour_array,carts[0],cylinders_size,1.0,cylinders_accu);
            line->SetColourOverride(catom_index);
            cylinders->add_primitive(line);
          }
}

void ConnectivityDraw::RedrawPrimitives(Displayobject &obj, const Connectivity &connectivity, int mode, const CParamsManager &params, const CParamsManager &global_params, int nSelAtoms,  const AtomColourVector &atom_colour_vector, const std::vector<double> &atomRadii, const SplineInfo &splineinfo_ribbon, const SplineInfo &splineinfo_worm, const std::string &texture, const std::string &bumpmap, int stick_colour, int side_to_ribbon, int side_to_worm, int bonds_mode){

  double width=1.0;
  bool warning = false;

  static int cylinders_accus[] = {8,18,24};

  double spheres_size=1.0;
  double cylinders_size=params.GetFloat("cylinder_width");
  int spheres_accu=1,cylinders_accu=4;
  double midpoint_frac0 = 0.5;
  double midpoint_frac1 = 0.5;
  std::vector <Cartesian> pyramid_carts;

  bool dashed = params.GetInt("dashed_bonds"); 
  double dash_length = params.GetFloat("dashed_bond_length"); 

  std::vector<Cartesian> all_ring_centres;

  //std::cout << "RedrawPrimitives " << bonds_mode << " " << side_to_ribbon << std::endl;

  std::vector <int> catom_indices = connectivity.GetCAtomIndex();

  double trace_cutoff = params.GetFloat("trace_cutoff"); 
  double trace_cutoff_sq = trace_cutoff*trace_cutoff;

  if(mode==SPHERES){
    spheres_size = 1.2;
    spheres_accu = global_params.GetInt("solid_quality"); 
    cylinders_accu = 8;
  }
  if(mode==BALLSTICK){
    spheres_size = 0.6;
    cylinders_size = params.GetFloat("ballstick_stick");
    spheres_accu = global_params.GetInt("solid_quality");
    if(global_params.GetInt("solid_quality")<3)
      cylinders_accu = cylinders_accus[global_params.GetInt("solid_quality")];
    else
      cylinders_accu = 8+8*global_params.GetInt("solid_quality");
  }
  if(mode==CYLINDERS){
    spheres_size = cylinders_size;
    spheres_accu = global_params.GetInt("solid_quality");
    if(global_params.GetInt("solid_quality")<3)
      cylinders_accu = cylinders_accus[global_params.GetInt("solid_quality")];
    else
      cylinders_accu = 8+8*global_params.GetInt("solid_quality");
  }
  if(mode==BONDS)
    width = params.GetInt("bond_width");
  else if (mode==FATBONDS)
    width = params.GetInt("fat_bond_width");
  else if (mode==PYRAMIDS) {
    float pyra_size = params.GetFloat("pyramid_size");
    pyramid_carts.push_back(Cartesian(2*pyra_size,pyra_size+0.001,0.001));
    pyramid_carts.push_back(Cartesian(-2*pyra_size+0.001,pyra_size,0.001));
    pyramid_carts.push_back(Cartesian(0.001,-pyra_size,2*pyra_size));
    pyramid_carts.push_back(Cartesian(0.001,-pyra_size,-2*pyra_size+0.001));
    pyramid_carts.push_back(Cartesian(2*pyra_size,pyra_size,0.001));
    pyramid_carts.push_back(Cartesian(0.001,-pyra_size+0.001,2*pyra_size));
  } else
    width = params.GetInt("thin_bond_width");

  bool showMultipleBonds = params.GetInt("show_multiple_bonds");

  bool drawDelocRing = params.GetInt("deloc_ring");
  bool drawDoubleRing;
  if(drawDelocRing){
    drawDoubleRing = false;
  } else {
    drawDoubleRing = true;
  }

  std::vector<std::vector<std::vector<int> > > conn_order_lists = connectivity.GetConnectivityLists();
  std::vector<std::vector<int> > conn_lists = conn_order_lists[0];
  std::vector<std::vector<int> > conn_orders = conn_order_lists[1];
  std::vector<std::vector<int> > ext_conn_lists_spline = connectivity.GetExternalSplineConnectivityLists();
  std::vector<std::vector<int> > ext_conn_lists = connectivity.GetExternalConnectivityLists();
  mmdb::Atom** atoms = connectivity.GetAtoms();
  int natoms = connectivity.GetNumberOfAtoms();
  std::vector <Cartesian> int_carts = CartesiansFromAtoms(atoms,natoms);
  if(atoms) delete [] atoms;

  std::vector<std::vector <Cartesian> > ext_carts = GetExternalCartesiansWithSplineInfo(molhnd,ext_conn_lists,splineinfo_ribbon,splineinfo_worm,side_to_ribbon,side_to_worm,params.GetFloat("trace_cutoff"));
  std::vector<std::vector <Cartesian> > ext_carts_spline = GetExternalCartesiansWithSplineInfo(molhnd,ext_conn_lists_spline,splineinfo_ribbon,splineinfo_worm,side_to_ribbon,side_to_worm,params.GetFloat("trace_cutoff"));

  LineCollection *lines = new LineCollection();
  PolyCollection *polys = new PolyCollection();
  SphereCollection *spheres = new SphereCollection();
  CylinderCollection *cylinders = new CylinderCollection();
  PointSphereCollection *point_spheres = new PointSphereCollection();

  LinesCollection *lines_collection = new LinesCollection();
  DashLinesCollection *dash_lines_collection = new DashLinesCollection();
  lines_collection->SetSize(width);
  dash_lines_collection->SetSize(width);
  dash_lines_collection->SetDashLength(dash_length);

  spheres->SetAccu(spheres_accu);
  cylinders->SetAccu(cylinders_accu);

  std::vector<Cartesian> carts(2);

  Cartesian tmp_v;
  Cartesian midpoint;

  const double *stick_colour_array=0;
  const double *colour_array=0;
  //double colour_array[4];
  if ( stick_colour > 0 ) stick_colour_array = RGBReps::GetColourP(stick_colour);
  if (side_to_ribbon>0 || side_to_worm>0 ) {
    midpoint_frac0 = 0.0;
    midpoint_frac1 = 1.0;
    //std::cout << "Setting midpoint to external cart\n";
  }

  lines_collection->reserve((ext_conn_lists_spline.size()+ext_conn_lists.size()+conn_lists.size())*3);
  if(dashed)
    dash_lines_collection->reserve((ext_conn_lists_spline.size()+ext_conn_lists.size()+conn_lists.size())*3);
  if (bonds_mode != DRAW_INTERNAL_BONDS) {

  //if(!atom_colour_vector) std::cout << "BIG PROBLEM: atom_colour_vector not valid\n";
  std::cout.flush();
  for(unsigned i=0;i<ext_conn_lists_spline.size()&&(side_to_ribbon>0 || side_to_worm>0);i++){
    colour_array = atom_colour_vector.GetRGB(i);
    
    for(unsigned j=0;j<ext_conn_lists_spline[i].size();j++){
      if(dashed==true&&(mode==BONDS||mode==FATBONDS||mode==THINBONDS)){
        //Cartesian midpoint =midpoint_frac0 *int_carts[i]+midpoint_frac1*ext_carts[i][j];
        midpoint.setMultiplyAndAdd(midpoint_frac0 ,int_carts[i],midpoint_frac1,ext_carts_spline[i][j]);
        carts[1] = int_carts[i];
        carts[0] = midpoint;
        tmp_v = carts[0] -  carts[1];
        if (tmp_v.length() < 0.001||tmp_v.length() > trace_cutoff_sq  ) {
          if (!warning) {
	    warning = true;
            std::cout << "Error drawing BONDS ext\n" ;
          }
          continue;
        }
        dash_lines_collection->add(carts,colour_array,catom_indices[i]);
        //DashLineElement *line = new DashLineElement(carts,colour_array,carts[1],width,1.0);
        //line->SetColourOverride(catom_indices[i]);
        //line->SetDashLength(dash_length);
        //lines->add_primitive(line);
      }
      if(dashed==false&&(mode==BONDS||mode==FATBONDS||mode==THINBONDS)){
        //Cartesian midpoint =midpoint_frac0 *int_carts[i]+midpoint_frac1*ext_carts[i][j];
        midpoint.setMultiplyAndAdd(midpoint_frac0,int_carts[i],midpoint_frac1,ext_carts_spline[i][j]);
        carts[0] = int_carts[i];
        carts[1] = midpoint;
        bool distOK = Cartesian::CheckDistanceRange(carts[0],carts[1],0.001,6);
        //tmp_v = carts[0] -  carts[1];
        //if (tmp_v.length() < 0.001 ||tmp_v.length() > 6  ) {
        if (!distOK) {
          if (!warning) {
	    warning = true;
            std::cout << "Error drawing BONDS ext\n" ;
          }
          continue;
        }
        lines_collection->add(carts,colour_array,catom_indices[i]);
        //LineElement *line = new LineElement(carts,colour_array,carts[0],width,1.0);
        //line->SetColourOverride(catom_indices[i]);
        //lines->add_primitive(line);
      }
      if(mode==CYLINDERS||mode==BALLSTICK){
        //Cartesian midpoint =midpoint_frac0 *int_carts[i]+midpoint_frac1*ext_carts[i][j];
        midpoint.setMultiplyAndAdd(midpoint_frac0 ,int_carts[i],midpoint_frac1,ext_carts_spline[i][j]);
        carts[0] = int_carts[i];
        carts[1] = midpoint;
        tmp_v = carts[0] -  carts[1];
        if (tmp_v.length() < 0.001 ) {
          if ( !warning ) {
            warning = true;
            std::cout << "Error drawing CYLINDER ext\n";
          }
          continue;
        }
        if (stick_colour > 0 && mode==BALLSTICK && stick_colour_array ) {
          CylinderElement *line = new CylinderElement(carts,stick_colour_array,carts[0],cylinders_size,1.0,cylinders_accu);
          line->SetColourOverride(catom_indices[i]);
          cylinders->add_primitive(line);
        } else {
          CylinderElement *line = new CylinderElement(carts,colour_array,carts[0],cylinders_size,1.0,cylinders_accu);
          line->SetColourOverride(catom_indices[i]);
          cylinders->add_primitive(line);
        }
      }
    }
    if(colour_array) delete [] colour_array;
  }

  //if(!atom_colour_vector) std::cout << "BIG PROBLEM: atom_colour_vector not valid\n";
  std::cout.flush();
  for(unsigned i=0;i<ext_conn_lists.size();i++){
    colour_array = atom_colour_vector.GetRGB(i);
    
    for(unsigned j=0;j<ext_conn_lists[i].size();j++){
      if(dashed==true&&(mode==BONDS||mode==FATBONDS||mode==THINBONDS)){
        //Cartesian midpoint =midpoint_frac0 *int_carts[i]+midpoint_frac1*ext_carts[i][j];
        midpoint.setMultiplyAndAdd(0.5 ,int_carts[i],0.5,ext_carts[i][j]);
        carts[1] = int_carts[i];
        carts[0] = midpoint;
        tmp_v = carts[0] -  carts[1];
        if (tmp_v.length() < 0.001||tmp_v.length() > trace_cutoff_sq  ) {
          if (!warning) {
	    warning = true;
            std::cout << "Error drawing BONDS ext\n" ;
          }
          continue;
        }
        dash_lines_collection->add(carts,colour_array,catom_indices[i]);
        //DashLineElement *line = new DashLineElement(carts,colour_array,carts[1],width,1.0);
        //line->SetColourOverride(catom_indices[i]);
        //line->SetDashLength(dash_length);
        //lines->add_primitive(line);
      }
      if(dashed==false&&(mode==BONDS||mode==FATBONDS||mode==THINBONDS)){
        //Cartesian midpoint =midpoint_frac0 *int_carts[i]+midpoint_frac1*ext_carts[i][j];
        midpoint.setMultiplyAndAdd(0.5,int_carts[i],0.5,ext_carts[i][j]);
        carts[0] = int_carts[i];
        carts[1] = midpoint;
        bool distOK = Cartesian::CheckDistanceRange(carts[0],carts[1],0.001,6);
        //tmp_v = carts[0] -  carts[1];
        //if (tmp_v.length() < 0.001 ||tmp_v.length() > 6  ) {
        if (!distOK) {
          if (!warning) {
	    warning = true;
            std::cout << "Error drawing BONDS ext\n" ;
          }
          continue;
        }
        lines_collection->add(carts,colour_array,catom_indices[i]);
        //LineElement *line = new LineElement(carts,colour_array,carts[0],width,1.0);
        //line->SetColourOverride(catom_indices[i]);
        //lines->add_primitive(line);
      }
      if(mode==CYLINDERS||mode==BALLSTICK){
        //Cartesian midpoint =midpoint_frac0 *int_carts[i]+midpoint_frac1*ext_carts[i][j];
        midpoint.setMultiplyAndAdd(0.5 ,int_carts[i],0.5,ext_carts[i][j]);
        carts[0] = int_carts[i];
        carts[1] = midpoint;
        tmp_v = carts[0] -  carts[1];
        if (tmp_v.length() < 0.001 ) {
          if ( !warning ) {
            warning = true;
            std::cout << "Error drawing CYLINDER ext\n";
          }
          continue;
        }
        if (stick_colour > 0 && mode==BALLSTICK && stick_colour_array ) {
          CylinderElement *line = new CylinderElement(carts,stick_colour_array,carts[0],cylinders_size,1.0,cylinders_accu);
          line->SetColourOverride(catom_indices[i]);
          cylinders->add_primitive(line);
        } else {
          CylinderElement *line = new CylinderElement(carts,colour_array,carts[0],cylinders_size,1.0,cylinders_accu);
          line->SetColourOverride(catom_indices[i]);
          cylinders->add_primitive(line);
        }
      }
    }
    if(colour_array) delete [] colour_array;
  }
 

  }

  if (bonds_mode != DRAW_EXTERNAL_BONDS) {

  for(unsigned i=0;i<conn_lists.size();i++){
    colour_array = atom_colour_vector.GetRGB(i);
    if(mode==POINTS_SPHERES){
      Point *sphere = new Point(int_carts[i],colour_array,int_carts[i],atomRadii[i],1.0);
      sphere->SetColourOverride(catom_indices[i]);
      point_spheres->add_primitive(sphere);
    }
    if(mode==SPHERES||mode==BALLSTICK){
      if (atomRadii[i]<0.01) {
	//std::cout << "SPHERE radii";
      } else {
      SphereElement *sphere = new SphereElement(int_carts[i],colour_array,int_carts[i],atomRadii[i],1.0,spheres_accu);
      sphere->SetColourOverride(catom_indices[i]);
      spheres->add_primitive(sphere);
      }
    }
    if(mode==CYLINDERS){
        //if (spheres_size<0.01) std::cout << "CYLINDER radii";
        SphereElement *sphere = new SphereElement(int_carts[i],colour_array,int_carts[i],spheres_size,1.0,spheres_accu);
        sphere->SetColourOverride(catom_indices[i]);
        spheres->add_primitive(sphere);
    }
    if(mode==PYRAMIDS){
      for (unsigned j=0;j<5;j++) {     
        carts[0] = int_carts[i]+pyramid_carts[j];
        carts[1] = int_carts[i]+pyramid_carts[j+1];
        LineElement *line = new LineElement(carts,colour_array,carts[1],width,1.0);
        line->SetColourOverride(catom_indices[i]);
        lines->add_primitive(line);
      }
      carts[0] = int_carts[i]+pyramid_carts[1];
      carts[1] = int_carts[i]+pyramid_carts[3];
      LineElement *line = new LineElement(carts,colour_array,carts[1],width,1.0);
      line->SetColourOverride(catom_indices[i]);
      lines->add_primitive(line);
    }

    for(unsigned j=0;j<conn_lists[i].size();j++){
      if(dashed==true&&(mode==BONDS||mode==FATBONDS||mode==THINBONDS)){
        //Cartesian midpoint = Cartesian::MidPoint(int_carts[i],int_carts[conn_lists[i][j]]);
        Cartesian midpoint = 0.50*int_carts[i]+0.50*int_carts[conn_lists[i][j]];
        carts[1] = int_carts[i];
        carts[0] = midpoint;
        tmp_v = carts[0] -  carts[1];
        if (tmp_v.length() < 0.001||tmp_v.length() > 6   ) {
          if (!warning) {
            warning = true;
            std::cout << "Error drawing BONDS \n";
          }
          continue;
        }
        dash_lines_collection->add(carts,colour_array,catom_indices[i]);
        //DashLineElement *line = new DashLineElement(carts,colour_array,carts[1],width,1.0);
        //line->SetColourOverride(catom_indices[i]);
        //line->SetDashLength(dash_length);
        //lines->add_primitive(line);
      }
      if(dashed==false&&(mode==BONDS||mode==FATBONDS||mode==THINBONDS)){
        //Cartesian midpoint = Cartesian::MidPoint(int_carts[i],int_carts[conn_lists[i][j]]);
        //Cartesian midpoint = 0.50*int_carts[i]+0.50*int_carts[conn_lists[i][j]];
        midpoint.setMultiplyAndAdd(0.50,int_carts[i],0.50,int_carts[conn_lists[i][j]]);
        carts[0] = int_carts[i];
        carts[1] = midpoint;
        bool distOK = Cartesian::CheckDistanceRange(carts[0],carts[1],0.001,6);
        if (!distOK) {
          if (!warning) {
            warning = true;
            std::cout << "Error drawing BONDS \n";
          }
          continue;
        }
        if(showMultipleBonds&&conn_orders[i][j]==ccp4srs::Bond::Triple){
          std::vector<int> thisFuncGroups = GetFuncGroups(conn_lists[i][j],i,conn_lists);
          std::vector<int> otherFuncGroups = GetFuncGroups(i,conn_lists[i][j],conn_lists);
          std::vector<int> adjFuncGroups;
          std::vector<Cartesian> planeAtoms;
          if(thisFuncGroups.size()>0&&otherFuncGroups.size()>0){
            //std::cout << "Both bigger than 0 " << thisFuncGroups.size() << " " << otherFuncGroups.size() << "\n";
            if(i>conn_lists[i][j]){
              adjFuncGroups = GetFuncGroups(i,thisFuncGroups[0],conn_lists);
              planeAtoms.push_back(int_carts[i]);
              planeAtoms.push_back(int_carts[thisFuncGroups[0]]);
            } else {
              adjFuncGroups = GetFuncGroups(conn_lists[i][j],otherFuncGroups[0],conn_lists);
              planeAtoms.push_back(int_carts[conn_lists[i][j]]);
              planeAtoms.push_back(int_carts[otherFuncGroups[0]]);
            }
          } else if(thisFuncGroups.size()>0){
            adjFuncGroups = GetFuncGroups(i,thisFuncGroups[0],conn_lists);
            planeAtoms.push_back(int_carts[i]);
            planeAtoms.push_back(int_carts[thisFuncGroups[0]]);
          } else if(otherFuncGroups.size()>0){
            adjFuncGroups = GetFuncGroups(conn_lists[i][j],otherFuncGroups[0],conn_lists);
            planeAtoms.push_back(int_carts[conn_lists[i][j]]);
            planeAtoms.push_back(int_carts[otherFuncGroups[0]]);
          }
          for(unsigned iadj=0;iadj<adjFuncGroups.size();iadj++){
            planeAtoms.push_back(int_carts[adjFuncGroups[iadj]]);
          }

          lines_collection->add(carts,colour_array,catom_indices[i]);

          Cartesian atmp;
          Cartesian l = carts[0] - carts[1];
          l.normalize();

          if(planeAtoms.size()>2){
             Cartesian l1 = planeAtoms[1]-planeAtoms[0];
             Cartesian l2 = planeAtoms[2]-planeAtoms[1];
             l1.normalize();
             l2.normalize();
             Cartesian cross1 = Cartesian::CrossProduct(l1,l2);
             cross1.normalize();
             atmp = Cartesian::CrossProduct(cross1,l);
          } else {
            Cartesian xaxis(1,0,0);
            Cartesian yaxis(0,1,0);
            Cartesian zaxis(0,0,1);
            atmp = Cartesian::CrossProduct(l,zaxis);
            if(atmp.length() < 0.000000001){
              atmp = Cartesian::CrossProduct(l,yaxis);
              if(atmp.length() < 0.000000001){
                atmp = Cartesian::CrossProduct(l,xaxis);
              }
            }
          }

          if(atmp.length() > 0.000000001){
            atmp.normalize(0.1);
            std::vector<Cartesian> double_bond_carts;
            double_bond_carts.push_back(carts[0]+atmp);
            double_bond_carts.push_back(carts[1]+atmp);
            lines_collection->add(double_bond_carts,colour_array,catom_indices[i]);
            double_bond_carts[0] = carts[0]-atmp;
            double_bond_carts[1] = carts[1]-atmp;
            lines_collection->add(double_bond_carts,colour_array,catom_indices[i]);
          }
        } else if(showMultipleBonds&&(conn_orders[i][j]==ccp4srs::Bond::Double||conn_orders[i][j]==ccp4srs::Bond::Deloc||conn_orders[i][j]==ccp4srs::Bond::Aromatic)){
          //std::cout << "Drawing a multiple bond " << i << " " << conn_lists[i][j] << " " << conn_orders[i][j] << "\n";
          //std::cout << "                          " << int_carts[i] << " " << int_carts[conn_lists[i][j]] << "\n";
          // Ring detection
          bool skipBond = false;
          std::vector<std::vector<int> > rings = CheckRing(i,conn_lists[i][j],conn_lists);
          if(rings.size()>0&&rings[0].size()>2){
            //std::cout << "It's an aromatic ring hopefully " << rings[0].size() << "\n";
            std::vector<int> ring = rings[0];
            for(int irings=1;irings<rings.size();irings++){
              // FIXME, we should favour aromatic rings ...
              if((conn_lists[i][j]<i&&rings[irings][1]<ring[1])||(conn_lists[i][j]>i&&rings[irings][rings[irings].size()-2]<ring[ring.size()-2])){
                ring = rings[irings];
              }
            }

            // FIXME - Need to check bond lengths also as in literature
            bool drawThisDoubleRing;
            if(drawDelocRing&&(conn_orders[i][j]==ccp4srs::Bond::Aromatic||conn_orders[i][j]==ccp4srs::Bond::Deloc)){
              Cartesian zaxis = Cartesian(0,0,1);
              std::vector<int> ring2 = rings[0];
              for(int irings=0;irings<rings.size();irings++){
                ring2 = rings[irings];
                std::vector<Cartesian> ring_carts2;
                std::vector<Cartesian> ringNorms;
                //std::vector<double> ringLengths;
                bool allAromatic = true;
                for(int iring=0;iring<ring2.size();iring++){
                  if(iring<ring2.size()-1){
                    for(unsigned irc=0;irc<conn_lists[ring2[iring]].size();irc++){
                      if(conn_lists[ring2[iring]][irc]==ring2[iring+1]){
                        if((conn_orders[ring2[iring]][irc]!=ccp4srs::Bond::Aromatic)&&(conn_orders[ring2[iring]][irc]!=ccp4srs::Bond::Deloc)&&(conn_orders[ring2[iring]][irc]!=ccp4srs::Bond::Double))
                          allAromatic = false;
                          break;
                      }
                    }
                  } else {
                    for(unsigned irc=0;irc<conn_lists[ring2[iring]].size();irc++){
                      if(conn_lists[ring2[iring]][irc]==ring2[0]){
                        if((conn_orders[ring2[iring]][irc]!=ccp4srs::Bond::Aromatic)&&(conn_orders[ring2[iring]][irc]!=ccp4srs::Bond::Deloc)&&(conn_orders[ring2[iring]][irc]!=ccp4srs::Bond::Double))
                          allAromatic = false;
                          break;
                      }
                    }
                  }
                  /*
                  if((conn_orders[ring2[iring]]!=ccp4srs::Bond::Deloc)&&(conn_orders[ring2[iring]]!=ccp4srs::Bond::Aromatic)){
                   allAromatic = false;
                  }
                  */
                  ring_carts2.push_back(int_carts[ring2[iring]]);
                  if(iring>0){
                    //ringLengths.push_back(LineLength(int_carts[ring2[iring]],int_carts[ring2[iring-1]]));
                  }
                  if(iring>1){
                    Cartesian crossRing = Cartesian::CrossProduct(int_carts[ring2[iring]]-int_carts[ring2[iring-1]],int_carts[ring2[iring-1]]-int_carts[ring2[iring-2]]);
                    crossRing.normalize();
                    ringNorms.push_back(crossRing);
                  }
                } 
                Cartesian crossRing = Cartesian::CrossProduct(int_carts[ring2[0]]-int_carts[ring2[ring2.size()-1]],int_carts[ring2[ring2.size()-1]]-int_carts[ring2[ring2.size()-2]]);
                crossRing.normalize();
                ringNorms.push_back(crossRing);
                crossRing = Cartesian::CrossProduct(int_carts[ring2[1]]-int_carts[ring2[0]],int_carts[ring2[0]]-int_carts[ring2[ring2.size()-1]]);
                crossRing.normalize();
                ringNorms.push_back(crossRing);
                Cartesian ringNorm = Cartesian::MidPoint(ringNorms);
                ringNorm.normalize();
                //ringLengths.push_back(LineLength(int_carts[ring2[0]],int_carts[ring2[ring2.size()-1]]));
                //double lstddev = VectorDoubleStandardDeviation(ringLengths);
                //std::cout << "Ring length stddev " << lstddev << "\n";

                Cartesian ring_centre = Cartesian::MidPoint(ring_carts2);
                double np = Cartesian::RMSDistance(ringNorms);
                bool drawThisRing = true;
                for(unsigned int i=0;i<all_ring_centres.size();i++){
                  if((all_ring_centres[i]-ring_centre).length()<1e-4){
                    //std::cout << "Already drawn this ring\n";
                    skipBond = true;
                    drawThisRing = false;
                    break;
                  }
                }
                if(skipBond&&drawDelocRing){
                  drawThisDoubleRing = false;
                  break;
                }

                if(!allAromatic){
                  drawThisRing = false;
                }else if(np>0.3){
                  drawThisRing = false;
                  drawThisDoubleRing = true;
                } else {
                  drawThisDoubleRing = false;
                }

                if(drawThisRing){
                  all_ring_centres.push_back(ring_centre);
                  MGArc *circ = new MGArc(ring_centre,colour_array,ring_centre,.5*(ring_carts2[0]-ring_centre).length());
                  double angle = acos(Cartesian::DotProduct(ringNorm,zaxis))*180.0/M_PI;
                  if(fabs(angle)>0.00001){
                    Cartesian rotax = Cartesian::CrossProduct(ringNorm,zaxis);
                    Quat q1 = Quat(rotax,1,-angle);
                    matrix m = q1.getInvMatrix();
                    circ->setMatrix(m);
                  }
                  circ->setWidth(width);
                  obj.add_primitive(circ);
                }
              }
            } else {
              drawThisDoubleRing = true;
            } 
            if(drawDoubleRing||drawThisDoubleRing){
              std::vector<Cartesian> ring_carts;
              for(int iring=0;iring<ring.size();iring++){
                ring_carts.push_back(int_carts[ring[iring]]);
              } 
              Cartesian ring_centre = Cartesian::MidPoint(ring_carts);
              Cartesian mc = ring_centre - midpoint;
              mc.normalize(0.1);
              std::vector<Cartesian> double_bond_carts;
              double_bond_carts.push_back(carts[0]+mc);
              double_bond_carts.push_back(carts[1]+mc);
              Cartesian dbadjust = double_bond_carts[0]-double_bond_carts[1];
              dbadjust.normalize(0.05);
              double_bond_carts[0] -= dbadjust;
              if(conn_orders[i][j]==ccp4srs::Bond::Double){
                lines_collection->add(double_bond_carts,colour_array,catom_indices[i]);
              } else {
                dash_lines_collection->add(double_bond_carts,colour_array,catom_indices[i]);
              }
            }
            lines_collection->add(carts,colour_array,catom_indices[i]);
          } else {
            std::vector<int> thisFuncGroups = GetFuncGroups(conn_lists[i][j],i,conn_lists);
            std::vector<int> otherFuncGroups = GetFuncGroups(i,conn_lists[i][j],conn_lists);
            //std::cout << thisFuncGroups.size() << " " << otherFuncGroups.size() << "\n";
            if(thisFuncGroups.size()==3||otherFuncGroups.size()==3){
                std::vector<Cartesian> planeCarts;
                std::vector<Cartesian> doublePlaneCarts;
                std::vector<Cartesian> singlePlaneCarts;
                int nDouble = 0;
              if(thisFuncGroups.size()==3){
                for(unsigned k=0;k<conn_lists[i].size();k++){
                  if(j!=k){
                    if(conn_orders[i][k]==ccp4srs::Bond::Double||conn_orders[i][k]==ccp4srs::Bond::Deloc||conn_orders[i][k]==ccp4srs::Bond::Aromatic){
                      doublePlaneCarts.push_back(int_carts[conn_lists[i][k]]);
                      nDouble++;
                    } else {
                      singlePlaneCarts.push_back(int_carts[conn_lists[i][k]]);
                    }
                  }
                }
              } else {
                for(unsigned k=0;k<conn_lists[conn_lists[i][j]].size();k++){
                  if(i!=conn_lists[conn_lists[i][j]][k]){
                    if(conn_orders[conn_lists[i][j]][k]==ccp4srs::Bond::Double||conn_orders[conn_lists[i][j]][k]==ccp4srs::Bond::Deloc||conn_orders[conn_lists[i][j]][k]==ccp4srs::Bond::Aromatic){
                      doublePlaneCarts.push_back(int_carts[conn_lists[conn_lists[i][j]][k]]);
                      nDouble++;
                    } else {
                      singlePlaneCarts.push_back(int_carts[conn_lists[conn_lists[i][j]][k]]);
                    }
                  }
                }
              }
                if(nDouble==0){
                  Cartesian mc = singlePlaneCarts[0] - singlePlaneCarts[1];
                  mc.normalize(0.05);
                  std::vector<Cartesian> double_bond_carts;
                  double_bond_carts.push_back(carts[0]+mc);
                  double_bond_carts.push_back(carts[1]+mc);
                  lines_collection->add(double_bond_carts,colour_array,catom_indices[i]);
                  double_bond_carts[0] = carts[0]-mc;
                  double_bond_carts[1] = carts[1]-mc;
                  if(conn_orders[i][j]==ccp4srs::Bond::Double){
                    lines_collection->add(double_bond_carts,colour_array,catom_indices[i]);
                  } else {
                    dash_lines_collection->add(double_bond_carts,colour_array,catom_indices[i]);
                  }
                }
                if(nDouble==1){
                  Cartesian mc = singlePlaneCarts[0] - singlePlaneCarts[1];
                  mc.normalize(0.05);
                  std::vector<Cartesian> double_bond_carts;
                  double_bond_carts.push_back(carts[0]+mc);
                  double_bond_carts.push_back(carts[1]+mc);
                  lines_collection->add(double_bond_carts,colour_array,catom_indices[i]);
                  double_bond_carts[0] = carts[0]-mc;
                  double_bond_carts[1] = carts[1]-mc;
                  if(conn_orders[i][j]==ccp4srs::Bond::Double){
                    lines_collection->add(double_bond_carts,colour_array,catom_indices[i]);
                  } else {
                    dash_lines_collection->add(double_bond_carts,colour_array,catom_indices[i]);
                  }
                }
                if(nDouble==2){
                  Cartesian cross1;
                  Cartesian cross2;
                  if(thisFuncGroups.size()==3){
                    cross1 = singlePlaneCarts[0]-carts[0];
                    cross2 = carts[1]-carts[0];
                  } else {
                    cross1 = singlePlaneCarts[0]-int_carts[conn_lists[i][j]];
                    cross2 = carts[0]-carts[1];
                  }
                  cross1.normalize();
                  cross2.normalize();
                  Cartesian mc = Cartesian::CrossProduct(cross1,cross2);
                  mc.normalize(0.05);
                  std::vector<Cartesian> double_bond_carts;
                  double_bond_carts.push_back(carts[0]+mc);
                  double_bond_carts.push_back(carts[1]+mc);
                  Cartesian dbadjust = double_bond_carts[0]-double_bond_carts[1];
                  if(thisFuncGroups.size()==3){
                    dbadjust.normalize(0.05);
                    double_bond_carts[0] -= dbadjust;
                  }
                  lines_collection->add(double_bond_carts,colour_array,catom_indices[i]);
                  double_bond_carts[0] = carts[0]-mc;
                  double_bond_carts[1] = carts[1]-mc;
                  if(thisFuncGroups.size()==3){
                    double_bond_carts[0] -= dbadjust;
                  }
                  if(conn_orders[i][j]==ccp4srs::Bond::Double){
                    lines_collection->add(double_bond_carts,colour_array,catom_indices[i]);
                  } else {
                    dash_lines_collection->add(double_bond_carts,colour_array,catom_indices[i]);
                  }
                }
                if(nDouble==3){
                  Cartesian cross1;
                  Cartesian cross2;
                  if(thisFuncGroups.size()==3){
                    if(j<2){
                      cross1 = doublePlaneCarts[1]-doublePlaneCarts[2];
                      cross2 = carts[1]-carts[0];
                    } else {
                      cross1 = doublePlaneCarts[0]-doublePlaneCarts[1];
                      cross2 = carts[1]-carts[0];
                    }
                  } else {
                    if(i<conn_lists[conn_lists[i][j]][2]){
                      cross1 = doublePlaneCarts[1]-doublePlaneCarts[2];
                      cross2 = carts[0]-carts[1];
                    } else {
                      cross1 = doublePlaneCarts[0]-doublePlaneCarts[1];
                      cross2 = carts[0]-carts[1];
                    }
                  }
                  cross1.normalize();
                  cross2.normalize();
                  Cartesian mc = cross1;
                  mc.normalize(0.05);
                  std::vector<Cartesian> double_bond_carts;
                  double_bond_carts.push_back(carts[0]+mc);
                  double_bond_carts.push_back(carts[1]+mc);
                  Cartesian dbadjust = double_bond_carts[0]-double_bond_carts[1];
                  if(thisFuncGroups.size()==3){
                    dbadjust.normalize(0.05);
                    double_bond_carts[0] -= dbadjust;
                  }
                  lines_collection->add(double_bond_carts,colour_array,catom_indices[i]);
                  double_bond_carts[0] = carts[0]-mc;
                  double_bond_carts[1] = carts[1]-mc;
                  if(thisFuncGroups.size()==3){
                    double_bond_carts[0] -= dbadjust;
                  }
                  if(conn_orders[i][j]==ccp4srs::Bond::Double){
                    lines_collection->add(double_bond_carts,colour_array,catom_indices[i]);
                  } else {
                    dash_lines_collection->add(double_bond_carts,colour_array,catom_indices[i]);
                  }
                }
            }
            if(thisFuncGroups.size()==2&&(otherFuncGroups.size()<2||(otherFuncGroups.size()==2&&(i<conn_lists[i][j])))){
              Cartesian mc = int_carts[thisFuncGroups[0]]-int_carts[thisFuncGroups[1]];
              bool groupZeroDouble = false;
              bool groupOneDouble = false;
              bool otherDouble = false;
              for(int k=0;k<conn_orders[i].size();k++){
                 if(conn_lists[i][k] == thisFuncGroups[0]){
                   if(conn_orders[i][k]==ccp4srs::Bond::Double||conn_orders[i][k]==ccp4srs::Bond::Deloc){
                     groupZeroDouble = true;
                   }
                 }
                 if(conn_lists[i][k] == thisFuncGroups[1]){
                   if(conn_orders[i][k]==ccp4srs::Bond::Double||conn_orders[i][k]==ccp4srs::Bond::Deloc){
                     groupOneDouble = true;
                   }
                 }
              }
              if(groupZeroDouble&&!groupOneDouble){
                mc = -mc;
                otherDouble = true;
              }
              if(otherDouble){
                mc.normalize(0.1);
                std::vector<Cartesian> double_bond_carts;
                double_bond_carts.push_back(carts[0]);
                double_bond_carts.push_back(carts[1]);
                lines_collection->add(double_bond_carts,colour_array,catom_indices[i]);
                double_bond_carts[0] = carts[0]-mc;
                double_bond_carts[1] = carts[1]-mc;
                Cartesian dbadjust = double_bond_carts[0]-double_bond_carts[1];
                dbadjust.normalize(0.05);
                double_bond_carts[0] -= dbadjust;
                if(conn_orders[i][j]==ccp4srs::Bond::Double){
                  lines_collection->add(double_bond_carts,colour_array,catom_indices[i]);
                } else {
                  dash_lines_collection->add(double_bond_carts,colour_array,catom_indices[i]);
                }
              } else {
                mc.normalize(0.05);
                std::vector<Cartesian> double_bond_carts;
                double_bond_carts.push_back(carts[0]+mc);
                double_bond_carts.push_back(carts[1]+mc);
                lines_collection->add(double_bond_carts,colour_array,catom_indices[i]);
                double_bond_carts[0] = carts[0]-mc;
                double_bond_carts[1] = carts[1]-mc;
                if(conn_orders[i][j]==ccp4srs::Bond::Double){
                  lines_collection->add(double_bond_carts,colour_array,catom_indices[i]);
                } else {
                  dash_lines_collection->add(double_bond_carts,colour_array,catom_indices[i]);
                }
              }
            }
            if((thisFuncGroups.size()<2||(thisFuncGroups.size()==2&&(i>conn_lists[i][j])))&&otherFuncGroups.size()==2){
              Cartesian mc = int_carts[otherFuncGroups[0]]-int_carts[otherFuncGroups[1]];
              bool groupZeroDouble = false;
              bool groupOneDouble = false;
              bool otherDouble = false;
              for(int k=0;k<conn_orders[conn_lists[i][j]].size();k++){
                 if(conn_lists[conn_lists[i][j]][k] == otherFuncGroups[0]){
                   if(conn_orders[conn_lists[i][j]][k]==ccp4srs::Bond::Double||conn_orders[conn_lists[i][j]][k]==ccp4srs::Bond::Deloc){
                     groupZeroDouble = true;
                   }
                 }
                 if(conn_lists[conn_lists[i][j]][k] == otherFuncGroups[1]){
                   if(conn_orders[conn_lists[i][j]][k]==ccp4srs::Bond::Double||conn_orders[conn_lists[i][j]][k]==ccp4srs::Bond::Deloc){
                     groupOneDouble = true;
                   }
                 }
              }
              if(groupZeroDouble&&!groupOneDouble){
                mc = -mc;
                otherDouble = true;
              }
              if(otherDouble){
                mc.normalize(0.1);
                std::vector<Cartesian> double_bond_carts;
                double_bond_carts.push_back(carts[0]);
                double_bond_carts.push_back(carts[1]);
                lines_collection->add(double_bond_carts,colour_array,catom_indices[i]);
                double_bond_carts[0] = carts[0]-mc;
                double_bond_carts[1] = carts[1]-mc;
                if(conn_orders[i][j]==ccp4srs::Bond::Double){
                  lines_collection->add(double_bond_carts,colour_array,catom_indices[i]);
                } else {
                  dash_lines_collection->add(double_bond_carts,colour_array,catom_indices[i]);
                }
              } else {
                mc.normalize(0.05);
                std::vector<Cartesian> double_bond_carts;
                double_bond_carts.push_back(carts[0]+mc);
                double_bond_carts.push_back(carts[1]+mc);
                lines_collection->add(double_bond_carts,colour_array,catom_indices[i]);
                double_bond_carts[0] = carts[0]-mc;
                double_bond_carts[1] = carts[1]-mc;
                if(conn_orders[i][j]==ccp4srs::Bond::Double){
                  lines_collection->add(double_bond_carts,colour_array,catom_indices[i]);
                } else {
                  dash_lines_collection->add(double_bond_carts,colour_array,catom_indices[i]);
                }
              }
            }
            if((thisFuncGroups.size()==1&&otherFuncGroups.size()==0)||(thisFuncGroups.size()==0&&otherFuncGroups.size()==1)){
              Cartesian l;
              if((thisFuncGroups.size()==1&&otherFuncGroups.size()==0)){
                l = int_carts[i]-int_carts[conn_lists[i][j]];
              } else {
                l = int_carts[conn_lists[i][j]]-int_carts[i];
              }
              l.normalize();
              Cartesian xaxis(1,0,0);
              Cartesian yaxis(0,1,0);
              Cartesian zaxis(0,0,1);
              Cartesian mc = Cartesian::CrossProduct(l,zaxis);
              if(mc.length() < 0.000000001){
                mc = Cartesian::CrossProduct(l,yaxis);
                if(mc.length() < 0.000000001){
                   mc = Cartesian::CrossProduct(l,xaxis);
                }
              }
              mc.normalize(0.05);
              std::vector<Cartesian> double_bond_carts;
              double_bond_carts.push_back(carts[0]+mc);
              double_bond_carts.push_back(carts[1]+mc);
              lines_collection->add(double_bond_carts,colour_array,catom_indices[i]);
              double_bond_carts[0] = carts[0]-mc;
              double_bond_carts[1] = carts[1]-mc;
              if(conn_orders[i][j]==ccp4srs::Bond::Double){
                lines_collection->add(double_bond_carts,colour_array,catom_indices[i]);
              } else {
                dash_lines_collection->add(double_bond_carts,colour_array,catom_indices[i]);
              }
            }
            if((thisFuncGroups.size()==0&&otherFuncGroups.size()==0)&&conn_lists[i].size()>0){
              Cartesian l;
              l = int_carts[i]-int_carts[conn_lists[i][0]];
              l.normalize();
              Cartesian xaxis(1,0,0);
              Cartesian yaxis(0,1,0);
              Cartesian zaxis(0,0,1);
              Cartesian mc = Cartesian::CrossProduct(l,zaxis);
              if(mc.length() < 0.000000001){
                mc = Cartesian::CrossProduct(l,yaxis);
                if(mc.length() < 0.000000001){
                   mc = Cartesian::CrossProduct(l,xaxis);
                }
              }
              mc.normalize(0.05);
              std::vector<Cartesian> double_bond_carts;
              double_bond_carts.push_back(carts[0]+mc);
              double_bond_carts.push_back(carts[1]+mc);
              lines_collection->add(double_bond_carts,colour_array,catom_indices[i]);
              double_bond_carts[0] = carts[0]-mc;
              double_bond_carts[1] = carts[1]-mc;
              if(conn_orders[i][j]==ccp4srs::Bond::Double){
                lines_collection->add(double_bond_carts,colour_array,catom_indices[i]);
              } else {
                dash_lines_collection->add(double_bond_carts,colour_array,catom_indices[i]);
              }
            }
            if(thisFuncGroups.size()==1&&otherFuncGroups.size()==1){
              Cartesian c1 = int_carts[otherFuncGroups[0]];
              Cartesian c2 = int_carts[conn_lists[i][j]];
              Cartesian c3 = int_carts[i];
              Cartesian c4 = int_carts[thisFuncGroups[0]];
              //std::cout << DihedralAngle(c1,c2,c3,c4) << "\n";
              if(fabs(DihedralAngle(c1,c2,c3,c4))>M_PI/2){
                // TRANS
                Cartesian cross1 = Cartesian::CrossProduct(c4-c1,c2-c3);
                cross1.normalize();
                Cartesian cross2 = c2-c3;
                cross2.normalize();
                Cartesian mc = Cartesian::CrossProduct(cross1,cross2);
                mc.normalize(0.05);
                std::vector<Cartesian> double_bond_carts;
                double_bond_carts.push_back(carts[0]+mc);
                double_bond_carts.push_back(carts[1]+mc);
                lines_collection->add(double_bond_carts,colour_array,catom_indices[i]);
                double_bond_carts[0] = carts[0]-mc;
                double_bond_carts[1] = carts[1]-mc;
                if(conn_orders[i][j]==ccp4srs::Bond::Double){
                  lines_collection->add(double_bond_carts,colour_array,catom_indices[i]);
                } else {
                  dash_lines_collection->add(double_bond_carts,colour_array,catom_indices[i]);
                }
              } else {
                // CIS
                Cartesian mc = Cartesian::MidPoint(c4,c1) - Cartesian::MidPoint(c2,c3);
                mc.normalize(0.1);
                std::vector<Cartesian> double_bond_carts;
                double_bond_carts.push_back(carts[0]+mc);
                double_bond_carts.push_back(carts[1]+mc);
                Cartesian dbadjust = double_bond_carts[0]-double_bond_carts[1];
                dbadjust.normalize(0.05);
                double_bond_carts[0] -= dbadjust;
                if(conn_orders[i][j]==ccp4srs::Bond::Double){
                  lines_collection->add(double_bond_carts,colour_array,catom_indices[i]);
                } else {
                  dash_lines_collection->add(double_bond_carts,colour_array,catom_indices[i]);
                }
                lines_collection->add(carts,colour_array,catom_indices[i]);
              }
            }
          }
        } else {
          /*
          if(conn_orders[i][j]==ccp4srs::Bond::Metal){
            std::cout << "METAL !!!!!!!!!!!!!!!\n";
          }
          */
          lines_collection->add(carts,colour_array,catom_indices[i]);
        }
      }
      if(mode==CYLINDERS||mode==BALLSTICK){
        //Cartesian midpoint = Cartesian::MidPoint(int_carts[i],int_carts[conn_lists[i][j]]);
        midpoint.setMultiplyAndAdd(0.50,int_carts[i],0.50,int_carts[conn_lists[i][j]]);
        carts[0] = int_carts[i];
        carts[1] = midpoint;
        tmp_v = carts[0] -  carts[1];
        if (tmp_v.length() < 0.001||tmp_v.length() > 6   ) {
	  if (!warning) {
	    warning = true;
            std::cout << "Error drawing CYLINDER \n";
          }
          continue;
        }

        double tbOffset;
        double tbSize;
        double dbOffset;
        double dbSize;
        if(mode==CYLINDERS){
          tbOffset = cylinders_size*2./3.;
          tbSize = cylinders_size*.3;
          dbOffset = cylinders_size*.5;
          dbSize = cylinders_size*.4;
        } else {
          tbOffset = cylinders_size*.9;
          tbSize = cylinders_size*.42;
          dbOffset = cylinders_size*.675;
          dbSize = cylinders_size*.54;
        }
        if(showMultipleBonds&&conn_orders[i][j]==ccp4srs::Bond::Triple){
          std::vector<int> thisFuncGroups = GetFuncGroups(conn_lists[i][j],i,conn_lists);
          std::vector<int> otherFuncGroups = GetFuncGroups(i,conn_lists[i][j],conn_lists);
          std::vector<int> adjFuncGroups;
          std::vector<Cartesian> planeAtoms;
          if(thisFuncGroups.size()>0&&otherFuncGroups.size()>0){
            //std::cout << "Both bigger than 0 " << thisFuncGroups.size() << " " << otherFuncGroups.size() << "\n";
            if(i>conn_lists[i][j]){
              adjFuncGroups = GetFuncGroups(i,thisFuncGroups[0],conn_lists);
              planeAtoms.push_back(int_carts[i]);
              planeAtoms.push_back(int_carts[thisFuncGroups[0]]);
            } else {
              adjFuncGroups = GetFuncGroups(conn_lists[i][j],otherFuncGroups[0],conn_lists);
              planeAtoms.push_back(int_carts[conn_lists[i][j]]);
              planeAtoms.push_back(int_carts[otherFuncGroups[0]]);
            }
          } else if(thisFuncGroups.size()>0){
            adjFuncGroups = GetFuncGroups(i,thisFuncGroups[0],conn_lists);
            planeAtoms.push_back(int_carts[i]);
            planeAtoms.push_back(int_carts[thisFuncGroups[0]]);
          } else if(otherFuncGroups.size()>0){
            adjFuncGroups = GetFuncGroups(conn_lists[i][j],otherFuncGroups[0],conn_lists);
            planeAtoms.push_back(int_carts[conn_lists[i][j]]);
            planeAtoms.push_back(int_carts[otherFuncGroups[0]]);
          }
          for(unsigned iadj=0;iadj<adjFuncGroups.size();iadj++){
            planeAtoms.push_back(int_carts[adjFuncGroups[iadj]]);
          }

          DrawCylinder(stick_colour,mode,carts,colour_array,stick_colour_array,tbSize,cylinders_accu,cylinders,catom_indices[i]);

          Cartesian atmp;
          Cartesian l = carts[0] - carts[1];
          l.normalize();

          if(planeAtoms.size()>2){
             Cartesian l1 = planeAtoms[1]-planeAtoms[0];
             Cartesian l2 = planeAtoms[2]-planeAtoms[1];
             l1.normalize();
             l2.normalize();
             Cartesian cross1 = Cartesian::CrossProduct(l1,l2);
             cross1.normalize();
             atmp = Cartesian::CrossProduct(cross1,l);
          } else {
            Cartesian xaxis(1,0,0);
            Cartesian yaxis(0,1,0);
            Cartesian zaxis(0,0,1);
            atmp = Cartesian::CrossProduct(l,zaxis);
            if(atmp.length() < 0.000000001){
              atmp = Cartesian::CrossProduct(l,yaxis);
              if(atmp.length() < 0.000000001){
                atmp = Cartesian::CrossProduct(l,xaxis);
              }
            }
          }

          if(atmp.length() > 0.000000001){
            atmp.normalize(tbOffset);
            std::vector<Cartesian> double_bond_carts;
            double_bond_carts.push_back(carts[0]+atmp);
            double_bond_carts.push_back(carts[1]+atmp);
            DrawCylinder(stick_colour,mode,double_bond_carts,colour_array,stick_colour_array,tbSize,cylinders_accu,cylinders,catom_indices[i]);
            double_bond_carts[0] = carts[0]-atmp;
            double_bond_carts[1] = carts[1]-atmp;
            DrawCylinder(stick_colour,mode,double_bond_carts,colour_array,stick_colour_array,tbSize,cylinders_accu,cylinders,catom_indices[i]);
          }
        } else if(showMultipleBonds&&(conn_orders[i][j]==ccp4srs::Bond::Double||conn_orders[i][j]==ccp4srs::Bond::Deloc||conn_orders[i][j]==ccp4srs::Bond::Aromatic)){
          //std::cout << "Drawing a multiple bond " << i << " " << conn_lists[i][j] << " " << conn_orders[i][j] << "\n";
          //std::cout << "                          " << int_carts[i] << " " << int_carts[conn_lists[i][j]] << "\n";
          // Ring detection
          bool skipBond = false;
          std::vector<std::vector<int> > rings = CheckRing(i,conn_lists[i][j],conn_lists);
          if(rings.size()>0&&rings[0].size()>2&&conn_orders[i][j]==ccp4srs::Bond::Aromatic){
            //std::cout << "It's a ring\n";
            std::vector<int> ring = rings[0];
            int theRing = 0;
            for(int irings=1;irings<rings.size();irings++){
              //std::cout << i << " " << conn_lists[i][j] << " " << rings[irings][1] << " " << ring[1] << " " << rings[irings][rings[irings].size()-2] << " " << ring[ring.size()-2] << "\n";
              if((conn_lists[i][j]<i&&rings[irings][1]<ring[1])||(conn_lists[i][j]>i&&rings[irings][rings[irings].size()-2]<ring[ring.size()-2])){
                //std::cout << "swap rings\n";
                ring = rings[irings];
                theRing = irings;
              }
            }

            // FIXME - Need to check bond lengths also as in literature
            bool drawThisDoubleRing;
            if(drawDelocRing){
              Cartesian zaxis = Cartesian(0,0,1);
              std::vector<int> ring2 = rings[0];
              for(int irings=0;irings<rings.size();irings++){
                //if(irings==theRing){
                ring2 = rings[irings];
                std::vector<Cartesian> ring_carts2;
                std::vector<Cartesian> ringNorms;
                bool allAromatic = true;
                for(int iring=0;iring<ring2.size();iring++){
                  if(iring<ring2.size()-1){
                    for(unsigned irc=0;irc<conn_lists[ring2[iring]].size();irc++){
                      if(conn_lists[ring2[iring]][irc]==ring2[iring+1]){
                        if((conn_orders[ring2[iring]][irc]!=ccp4srs::Bond::Aromatic)&&(conn_orders[ring2[iring]][irc]!=ccp4srs::Bond::Deloc)&&(conn_orders[ring2[iring]][irc]!=ccp4srs::Bond::Double))
                          allAromatic = false;
                          break;
                      }
                    }
                  } else {
                    for(unsigned irc=0;irc<conn_lists[ring2[iring]].size();irc++){
                      if(conn_lists[ring2[iring]][irc]==ring2[0]){
                        if((conn_orders[ring2[iring]][irc]!=ccp4srs::Bond::Aromatic)&&(conn_orders[ring2[iring]][irc]!=ccp4srs::Bond::Deloc)&&(conn_orders[ring2[iring]][irc]!=ccp4srs::Bond::Double))
                          allAromatic = false;
                          break;
                      }
                    }
                  }
                  ring_carts2.push_back(int_carts[ring2[iring]]);
                  if(iring>1){
                    Cartesian crossRing = Cartesian::CrossProduct(int_carts[ring2[iring]]-int_carts[ring2[iring-1]],int_carts[ring2[iring-1]]-int_carts[ring2[iring-2]]);
                    crossRing.normalize();
                    ringNorms.push_back(crossRing);
                  }
                } 
                Cartesian crossRing = Cartesian::CrossProduct(int_carts[ring2[0]]-int_carts[ring2[ring2.size()-1]],int_carts[ring2[ring2.size()-1]]-int_carts[ring2[ring2.size()-2]]);
                crossRing.normalize();
                ringNorms.push_back(crossRing);
                crossRing = Cartesian::CrossProduct(int_carts[ring2[1]]-int_carts[ring2[0]],int_carts[ring2[0]]-int_carts[ring2[ring2.size()-1]]);
                crossRing.normalize();
                ringNorms.push_back(crossRing);
                Cartesian ringNorm = Cartesian::MidPoint(ringNorms);
                ringNorm.normalize();

                Cartesian ring_centre = Cartesian::MidPoint(ring_carts2);
                double np = Cartesian::RMSDistance(ringNorms);
                bool drawThisRing = true;
                for(unsigned int iallr=0;iallr<all_ring_centres.size();iallr++){
                  if((all_ring_centres[iallr]-ring_centre).length()<1e-4){
                    drawThisRing = false;
                    //std::cout << "Already drawn this ring " << i << " " << conn_lists[i][j] << "\n";
                    skipBond = true;
                    break;
                  }
                }
                if(skipBond&&drawDelocRing){
                  drawThisDoubleRing = false;
                  break;
                }

                //std::cout << np << "\n";
                if(!allAromatic){
                  drawThisRing = false;
                }else if(np>0.3){
                  drawThisRing = false;
                  drawThisDoubleRing = true;
                } else {
                  drawThisDoubleRing = false;
                }

                if(drawThisRing){
            //std::cout << "aromatic ring " << " " << i << " " << conn_lists[i][j] << "\n";
                  all_ring_centres.push_back(ring_centre);
                  Torus *circ= new Torus(ring_centre,colour_array,ring_centre,.5*(ring_carts2[0]-ring_centre).length()-.5*cylinders_size,.1);
                  double angle = acos(Cartesian::DotProduct(ringNorm,zaxis))*180.0/M_PI;
                  if(fabs(angle)>0.00001){
                    Cartesian rotax = Cartesian::CrossProduct(ringNorm,zaxis);
                    Quat q1 = Quat(rotax,1,-angle);
                    matrix m = q1.getInvMatrix();
                    circ->setMatrix(m);
                  }
                  obj.add_primitive(circ);
                }
              //}
              }
            } 

            //std::cout << "drawThisDoubleRing " << drawThisDoubleRing << " " << i << " " << conn_lists[i][j] << "\n";
            if(drawDoubleRing||drawThisDoubleRing){
            std::vector<Cartesian> ring_carts;
            for(int iring=0;iring<ring.size();iring++){
              ring_carts.push_back(int_carts[ring[iring]]);
            } 
            Cartesian ring_centre = Cartesian::MidPoint(ring_carts);
            Cartesian mc = ring_centre - midpoint;
            mc.normalize(dbOffset);
            std::vector<Cartesian> double_bond_carts;
            double_bond_carts.push_back(carts[0]+mc);
            double_bond_carts.push_back(carts[1]+mc);
            if(conn_orders[i][j]==ccp4srs::Bond::Double){
              DrawCylinder(stick_colour,mode,double_bond_carts,colour_array,stick_colour_array,dbSize,cylinders_accu,cylinders,catom_indices[i]);
            } else {
              DrawDashedCylinder(stick_colour,mode,double_bond_carts,colour_array,stick_colour_array,dbSize,spheres_accu,spheres,catom_indices[i]);
            }
            double_bond_carts[0] -= 2.*mc;
            double_bond_carts[1] -= 2.*mc;
            DrawCylinder(stick_colour,mode,double_bond_carts,colour_array,stick_colour_array,dbSize,cylinders_accu,cylinders,catom_indices[i]);
            } else {
              DrawCylinder(stick_colour,mode,carts,colour_array,stick_colour_array,cylinders_size,cylinders_accu,cylinders,catom_indices[i]);
            }
          } else {
            std::vector<int> thisFuncGroups = GetFuncGroups(conn_lists[i][j],i,conn_lists);
            std::vector<int> otherFuncGroups = GetFuncGroups(i,conn_lists[i][j],conn_lists);
            //std::cout << thisFuncGroups.size() << " " << otherFuncGroups.size() << "\n";
            if(thisFuncGroups.size()==3||otherFuncGroups.size()==3){
                std::vector<Cartesian> planeCarts;
                std::vector<Cartesian> doublePlaneCarts;
                std::vector<Cartesian> singlePlaneCarts;
                int nDouble = 0;
              if(thisFuncGroups.size()==3){
                for(unsigned k=0;k<conn_lists[i].size();k++){
                  if(j!=k){
                    if(conn_orders[i][k]==ccp4srs::Bond::Double||conn_orders[i][k]==ccp4srs::Bond::Deloc||conn_orders[i][k]==ccp4srs::Bond::Aromatic){
                      doublePlaneCarts.push_back(int_carts[conn_lists[i][k]]);
                      nDouble++;
                    } else {
                      singlePlaneCarts.push_back(int_carts[conn_lists[i][k]]);
                    }
                  }
                }
              } else {
                for(unsigned k=0;k<conn_lists[conn_lists[i][j]].size();k++){
                  if(i!=conn_lists[conn_lists[i][j]][k]){
                    if(conn_orders[conn_lists[i][j]][k]==ccp4srs::Bond::Double||conn_orders[conn_lists[i][j]][k]==ccp4srs::Bond::Deloc||conn_orders[conn_lists[i][j]][k]==ccp4srs::Bond::Aromatic){
                      doublePlaneCarts.push_back(int_carts[conn_lists[conn_lists[i][j]][k]]);
                      nDouble++;
                    } else {
                      singlePlaneCarts.push_back(int_carts[conn_lists[conn_lists[i][j]][k]]);
                    }
                  }
                }
              }
                if(nDouble==0){
                  Cartesian mc = singlePlaneCarts[0] - singlePlaneCarts[1];
                  mc.normalize(dbOffset);
                  std::vector<Cartesian> double_bond_carts;
                  double_bond_carts.push_back(carts[0]+mc);
                  double_bond_carts.push_back(carts[1]+mc);
                  DrawCylinder(stick_colour,mode,double_bond_carts,colour_array,stick_colour_array,dbSize,cylinders_accu,cylinders,catom_indices[i]);
                  double_bond_carts[0] = carts[0]-mc;
                  double_bond_carts[1] = carts[1]-mc;
                  if(conn_orders[i][j]==ccp4srs::Bond::Double){
                    DrawCylinder(stick_colour,mode,double_bond_carts,colour_array,stick_colour_array,dbSize,cylinders_accu,cylinders,catom_indices[i]);
                  } else {
                    DrawDashedCylinder(stick_colour,mode,double_bond_carts,colour_array,stick_colour_array,dbSize,spheres_accu,spheres,catom_indices[i]);
                  }
                }
                if(nDouble==1){
                  Cartesian mc = singlePlaneCarts[0] - singlePlaneCarts[1];
                  mc.normalize(dbOffset);
                  std::vector<Cartesian> double_bond_carts;
                  double_bond_carts.push_back(carts[0]+mc);
                  double_bond_carts.push_back(carts[1]+mc);
                  DrawCylinder(stick_colour,mode,double_bond_carts,colour_array,stick_colour_array,dbSize,cylinders_accu,cylinders,catom_indices[i]);
                  double_bond_carts[0] = carts[0]-mc;
                  double_bond_carts[1] = carts[1]-mc;
                  if(conn_orders[i][j]==ccp4srs::Bond::Double){
                    DrawCylinder(stick_colour,mode,double_bond_carts,colour_array,stick_colour_array,dbSize,cylinders_accu,cylinders,catom_indices[i]);
                  } else {
                    DrawDashedCylinder(stick_colour,mode,double_bond_carts,colour_array,stick_colour_array,dbSize,spheres_accu,spheres,catom_indices[i]);
                  }
                }
                if(nDouble==2){
                  Cartesian cross1;
                  Cartesian cross2;
                  if(thisFuncGroups.size()==3){
                    cross1 = singlePlaneCarts[0]-carts[0];
                    cross2 = carts[1]-carts[0];
                  } else {
                    cross1 = singlePlaneCarts[0]-int_carts[conn_lists[i][j]];
                    cross2 = carts[0]-carts[1];
                  }
                  cross1.normalize();
                  cross2.normalize();
                  Cartesian mc = Cartesian::CrossProduct(cross1,cross2);
                  mc.normalize(dbOffset);
                  std::vector<Cartesian> double_bond_carts;
                  double_bond_carts.push_back(carts[0]+mc);
                  double_bond_carts.push_back(carts[1]+mc);
                  Cartesian dbadjust = double_bond_carts[0]-double_bond_carts[1];
                  if(thisFuncGroups.size()==3){
                    dbadjust.normalize(0.05);
                    double_bond_carts[0] -= dbadjust;
                  }
                  DrawCylinder(stick_colour,mode,double_bond_carts,colour_array,stick_colour_array,dbSize,cylinders_accu,cylinders,catom_indices[i]);
                  double_bond_carts[0] = carts[0]-mc;
                  double_bond_carts[1] = carts[1]-mc;
                  if(thisFuncGroups.size()==3){
                    double_bond_carts[0] -= dbadjust;
                  }
                  if(conn_orders[i][j]==ccp4srs::Bond::Double){
                    DrawCylinder(stick_colour,mode,double_bond_carts,colour_array,stick_colour_array,dbSize,cylinders_accu,cylinders,catom_indices[i]);
                  } else {
                    DrawDashedCylinder(stick_colour,mode,double_bond_carts,colour_array,stick_colour_array,dbSize,spheres_accu,spheres,catom_indices[i]);
                  }
                }
                if(nDouble==3){
                  Cartesian cross1;
                  Cartesian cross2;
                  if(thisFuncGroups.size()==3){
                    if(j<2){
                      cross1 = doublePlaneCarts[1]-doublePlaneCarts[2];
                      cross2 = carts[1]-carts[0];
                    } else {
                      cross1 = doublePlaneCarts[0]-doublePlaneCarts[1];
                      cross2 = carts[1]-carts[0];
                    }
                  } else {
                    if(i<conn_lists[conn_lists[i][j]][2]){
                      cross1 = doublePlaneCarts[1]-doublePlaneCarts[2];
                      cross2 = carts[0]-carts[1];
                    } else {
                      cross1 = doublePlaneCarts[0]-doublePlaneCarts[1];
                      cross2 = carts[0]-carts[1];
                    }
                  }
                  cross1.normalize();
                  cross2.normalize();
                  Cartesian mc = cross1;
                  mc.normalize(dbOffset);
                  std::vector<Cartesian> double_bond_carts;
                  double_bond_carts.push_back(carts[0]+mc);
                  double_bond_carts.push_back(carts[1]+mc);
                  DrawCylinder(stick_colour,mode,double_bond_carts,colour_array,stick_colour_array,dbSize,cylinders_accu,cylinders,catom_indices[i]);
                  double_bond_carts[0] = carts[0]-mc;
                  double_bond_carts[1] = carts[1]-mc;
                  if(conn_orders[i][j]==ccp4srs::Bond::Double){
                    DrawCylinder(stick_colour,mode,double_bond_carts,colour_array,stick_colour_array,dbSize,cylinders_accu,cylinders,catom_indices[i]);
                  } else {
                    DrawDashedCylinder(stick_colour,mode,double_bond_carts,colour_array,stick_colour_array,dbSize,spheres_accu,spheres,catom_indices[i]);
                  }
                }
            }
            if(thisFuncGroups.size()==2&&(otherFuncGroups.size()<2||(otherFuncGroups.size()==2&&(i<conn_lists[i][j])))){
              Cartesian mc = int_carts[thisFuncGroups[0]]-int_carts[thisFuncGroups[1]];
              mc.normalize(dbOffset);
              std::vector<Cartesian> double_bond_carts;
              double_bond_carts.push_back(carts[0]+mc);
              double_bond_carts.push_back(carts[1]+mc);
              DrawCylinder(stick_colour,mode,double_bond_carts,colour_array,stick_colour_array,dbSize,cylinders_accu,cylinders,catom_indices[i]);
              double_bond_carts[0] = carts[0]-mc;
              double_bond_carts[1] = carts[1]-mc;
              if(conn_orders[i][j]==ccp4srs::Bond::Double){
                DrawCylinder(stick_colour,mode,double_bond_carts,colour_array,stick_colour_array,dbSize,cylinders_accu,cylinders,catom_indices[i]);
              } else {
                DrawDashedCylinder(stick_colour,mode,double_bond_carts,colour_array,stick_colour_array,dbSize,spheres_accu,spheres,catom_indices[i]);
              }
            }
            if((thisFuncGroups.size()<2||(thisFuncGroups.size()==2&&(i>conn_lists[i][j])))&&otherFuncGroups.size()==2){
              Cartesian mc = int_carts[otherFuncGroups[0]]-int_carts[otherFuncGroups[1]];
              mc.normalize(dbOffset);
              std::vector<Cartesian> double_bond_carts;
              double_bond_carts.push_back(carts[0]+mc);
              double_bond_carts.push_back(carts[1]+mc);
              DrawCylinder(stick_colour,mode,double_bond_carts,colour_array,stick_colour_array,dbSize,cylinders_accu,cylinders,catom_indices[i]);
              double_bond_carts[0] = carts[0]-mc;
              double_bond_carts[1] = carts[1]-mc;
              if(conn_orders[i][j]==ccp4srs::Bond::Double){
                DrawCylinder(stick_colour,mode,double_bond_carts,colour_array,stick_colour_array,dbSize,cylinders_accu,cylinders,catom_indices[i]);
              } else {
                DrawDashedCylinder(stick_colour,mode,double_bond_carts,colour_array,stick_colour_array,dbSize,spheres_accu,spheres,catom_indices[i]);
              }
            }
            if(thisFuncGroups.size()==1&&otherFuncGroups.size()==1){
              Cartesian c1 = int_carts[otherFuncGroups[0]];
              Cartesian c2 = int_carts[conn_lists[i][j]];
              Cartesian c3 = int_carts[i];
              Cartesian c4 = int_carts[thisFuncGroups[0]];
              //std::cout << DihedralAngle(c1,c2,c3,c4) << "\n";
              if(fabs(DihedralAngle(c1,c2,c3,c4))>M_PI/2){
                // TRANS
                Cartesian cross1 = Cartesian::CrossProduct(c4-c1,c2-c3);
                cross1.normalize();
                Cartesian cross2 = c2-c3;
                cross2.normalize();
                Cartesian mc = Cartesian::CrossProduct(cross1,cross2);
                mc.normalize(dbOffset);
                std::vector<Cartesian> double_bond_carts;
                double_bond_carts.push_back(carts[0]+mc);
                double_bond_carts.push_back(carts[1]+mc);
                DrawCylinder(stick_colour,mode,double_bond_carts,colour_array,stick_colour_array,dbSize,cylinders_accu,cylinders,catom_indices[i]);
                double_bond_carts[0] = carts[0]-mc;
                double_bond_carts[1] = carts[1]-mc;
                if(conn_orders[i][j]==ccp4srs::Bond::Double){
                  DrawCylinder(stick_colour,mode,double_bond_carts,colour_array,stick_colour_array,dbSize,cylinders_accu,cylinders,catom_indices[i]);
                } else {
                  DrawDashedCylinder(stick_colour,mode,double_bond_carts,colour_array,stick_colour_array,dbSize,spheres_accu,spheres,catom_indices[i]);
                }
              } else {
                // CIS
                Cartesian mc = Cartesian::MidPoint(c4,c1) - Cartesian::MidPoint(c2,c3);
                mc.normalize(dbOffset);
                std::vector<Cartesian> double_bond_carts;
                double_bond_carts.push_back(carts[0]+mc);
                double_bond_carts.push_back(carts[1]+mc);
                Cartesian dbadjust = double_bond_carts[0]-double_bond_carts[1];
                dbadjust.normalize(dbOffset);
                double_bond_carts[0] -= dbadjust;
                if(conn_orders[i][j]==ccp4srs::Bond::Double){
                  DrawCylinder(stick_colour,mode,double_bond_carts,colour_array,stick_colour_array,dbSize,cylinders_accu,cylinders,catom_indices[i]);
                } else {
                  DrawDashedCylinder(stick_colour,mode,double_bond_carts,colour_array,stick_colour_array,dbSize,spheres_accu,spheres,catom_indices[i]);
                }
                double_bond_carts[0] = carts[0]-mc;
                double_bond_carts[1] = carts[1]-mc;
                DrawCylinder(stick_colour,mode,double_bond_carts,colour_array,stick_colour_array,dbSize,cylinders_accu,cylinders,catom_indices[i]);
              }
            }
            if((thisFuncGroups.size()==1&&otherFuncGroups.size()==0)||(thisFuncGroups.size()==0&&otherFuncGroups.size()==1)){
              Cartesian l;
              if((thisFuncGroups.size()==1&&otherFuncGroups.size()==0)){
                l = int_carts[i]-int_carts[conn_lists[i][j]];
              } else {
                l = int_carts[conn_lists[i][j]]-int_carts[i];
              }
              l.normalize();
              Cartesian xaxis(1,0,0);
              Cartesian yaxis(0,1,0);
              Cartesian zaxis(0,0,1);
              Cartesian mc = Cartesian::CrossProduct(l,zaxis);
              if(mc.length() < 0.000000001){
                mc = Cartesian::CrossProduct(l,yaxis);
                if(mc.length() < 0.000000001){
                   mc = Cartesian::CrossProduct(l,xaxis);
                }
              }
              mc.normalize(dbOffset);
              std::vector<Cartesian> double_bond_carts;
              double_bond_carts.push_back(carts[0]+mc);
              double_bond_carts.push_back(carts[1]+mc);
              DrawCylinder(stick_colour,mode,double_bond_carts,colour_array,stick_colour_array,dbSize,cylinders_accu,cylinders,catom_indices[i]);
              double_bond_carts[0] = carts[0]-mc;
              double_bond_carts[1] = carts[1]-mc;
              if(conn_orders[i][j]==ccp4srs::Bond::Double){
                  DrawCylinder(stick_colour,mode,double_bond_carts,colour_array,stick_colour_array,dbSize,cylinders_accu,cylinders,catom_indices[i]);
              } else {
                  DrawDashedCylinder(stick_colour,mode,double_bond_carts,colour_array,stick_colour_array,dbSize,spheres_accu,spheres,catom_indices[i]);
              }
            }
            if((thisFuncGroups.size()==0&&otherFuncGroups.size()==0)&&conn_lists[i].size()>0) {
              Cartesian l;
              l = int_carts[i]-int_carts[conn_lists[i][0]];
              l.normalize();
              Cartesian xaxis(1,0,0);
              Cartesian yaxis(0,1,0);
              Cartesian zaxis(0,0,1);
              Cartesian mc = Cartesian::CrossProduct(l,zaxis);
              if(mc.length() < 0.000000001){
                mc = Cartesian::CrossProduct(l,yaxis);
                if(mc.length() < 0.000000001){
                   mc = Cartesian::CrossProduct(l,xaxis);
                }
              }
              mc.normalize(dbOffset);
              std::vector<Cartesian> double_bond_carts;
              double_bond_carts.push_back(carts[0]+mc);
              double_bond_carts.push_back(carts[1]+mc);
              DrawCylinder(stick_colour,mode,double_bond_carts,colour_array,stick_colour_array,dbSize,cylinders_accu,cylinders,catom_indices[i]);
              double_bond_carts[0] = carts[0]-mc;
              double_bond_carts[1] = carts[1]-mc;
              if(conn_orders[i][j]==ccp4srs::Bond::Double){
                  DrawCylinder(stick_colour,mode,double_bond_carts,colour_array,stick_colour_array,dbSize,cylinders_accu,cylinders,catom_indices[i]);
              } else {
                  DrawDashedCylinder(stick_colour,mode,double_bond_carts,colour_array,stick_colour_array,dbSize,spheres_accu,spheres,catom_indices[i]);
              }
            }
          }
        } else {
           DrawCylinder(stick_colour,mode,carts,colour_array,stick_colour_array,cylinders_size,cylinders_accu,cylinders,catom_indices[i]);
        }
      }
    }
    if(conn_lists[i].size()==0 && ext_conn_lists[i].size()==0 && ext_conn_lists_spline[i].size()==0){
      if(mode==BONDS||mode==FATBONDS||mode==THINBONDS){
        Cartesian p1 = int_carts[i];
        p1.set_x(p1.get_x()-0.2);
        Cartesian p2 = int_carts[i];
        p2.set_x(p2.get_x()+0.2);
        Cartesian origin = int_carts[i];
        carts[0] = p1; 
        carts[1] = p2;
        lines_collection->add(carts,colour_array,0);
        //LineElement* xaxis = new LineElement(carts,colour_array,origin,width,1.0);
        p1 = int_carts[i];
        p1.set_y(p1.get_y()-0.2);
        p2 = int_carts[i];
        p2.set_y(p2.get_y()+0.2);
        carts[0] = p1; 
        carts[1] = p2;
        lines_collection->add(carts,colour_array,0);
        //LineElement* yaxis = new LineElement(carts,colour_array,origin,width,1.0);
        p1 = int_carts[i];
        p1.set_z(p1.get_z()-0.2);
        p2 = int_carts[i];
        p2.set_z(p2.get_z()+0.2);
        carts[0] = p1; 
        carts[1] = p2;
        lines_collection->add(carts,colour_array,0);
        //LineElement* zaxis = new LineElement(carts,colour_array,origin,width,1.0);
        //lines->add_primitive(xaxis);
        //lines->add_primitive(yaxis);
        //lines->add_primitive(zaxis);
      }
    }
    if(colour_array) delete [] colour_array;
  }
  }

  lines->SetSize(width);
  obj.add_primitive(lines);
  obj.add_primitive(lines_collection);
  obj.add_primitive(dash_lines_collection);
  obj.add_primitive(polys);
  obj.add_primitive(spheres);
  obj.add_primitive(cylinders);
  obj.add_primitive(point_spheres);
  if ( stick_colour > 0 ) delete [] stick_colour_array;

}

std::vector<Cartesian> GetLineThroughBasePairs(mmdb::Residue* res1, mmdb::Residue* res2){
  int natoms1;
  mmdb::Atom** atoms1=0;
  res1->GetAtomTable1(atoms1,natoms1);
  int natoms2;
  mmdb::Atom** atoms2=0;
  res2->GetAtomTable1(atoms2,natoms2);

  std::vector<Cartesian> N1;
  std::vector<Cartesian> O1;
  std::vector<Cartesian> N2;
  std::vector<Cartesian> O2;

  for(int i=0;i<natoms1;i++){
    if(!strncmp(atoms1[i]->element," N",2)) N1.push_back(Cartesian(atoms1[i]->x,atoms1[i]->y,atoms1[i]->z));
    if(!strncmp(atoms1[i]->element," O",2)) O1.push_back(Cartesian(atoms1[i]->x,atoms1[i]->y,atoms1[i]->z));
  }
  for(int i=0;i<natoms2;i++){
    if(!strncmp(atoms2[i]->element," N",2)) N2.push_back(Cartesian(atoms2[i]->x,atoms2[i]->y,atoms2[i]->z));
    if(!strncmp(atoms2[i]->element," O",2)) O2.push_back(Cartesian(atoms2[i]->x,atoms2[i]->y,atoms2[i]->z));
  }

  //std::cout << N1.size() << " " << O1.size() << "\n";
  //std::cout << N2.size() << " " << O2.size() << "\n";

  std::vector<Cartesian> vecs;
  Cartesian midpoints(0,0,0);
  //std::cout << res1->GetResName() << "\n";
  std::vector<Cartesian>::iterator O1iter = O1.begin();
  while(O1iter!=O1.end()){
    std::vector<Cartesian>::iterator N2iter = N2.begin();
    while(N2iter!=N2.end()){
      double ONlength = ((*O1iter)-(*N2iter)).length();
      if(ONlength>2.6&&ONlength<3.2){
        //std::cout << "ON? " << ONlength << "\n";
        vecs.push_back((*O1iter)-(*N2iter));
        midpoints += Cartesian::MidPoint((*O1iter),(*N2iter));
      }
      N2iter++;
    }
    O1iter++;
  }
  std::vector<Cartesian>::iterator N1iter = N1.begin();
  while(N1iter!=N1.end()){
    std::vector<Cartesian>::iterator N2iter = N2.begin();
    while(N2iter!=N2.end()){
      double NNlength = ((*N1iter)-(*N2iter)).length();
      if(NNlength>2.6&&NNlength<3.2){
        //std::cout << "NN? " << NNlength << "\n";
        vecs.push_back((*N1iter)-(*N2iter));
	midpoints += Cartesian::MidPoint((*N1iter),(*N2iter));
      }
      N2iter++;
    }
    std::vector<Cartesian>::iterator O2iter = O2.begin();
    while(O2iter!=O2.end()){
      double NOlength = ((*N1iter)-(*O2iter)).length();
      if(NOlength>2.6&&NOlength<3.2){
        //std::cout << "NO? " << NOlength << "\n";
        vecs.push_back((*N1iter)-(*O2iter));
	midpoints += Cartesian::MidPoint((*N1iter),(*O2iter));
      }
      O2iter++;
    }
    N1iter++;
  }

  Cartesian v(0,0,0); 
  Cartesian m(0,0,0); 
  std::vector<Cartesian> vm;
  if(vecs.size()>1){
    for(unsigned ii=0;ii<vecs.size();ii++){
      vecs[ii].normalize();
      v += vecs[ii];
      //std::cout << vecs[ii] << "\n";
    }
    v /= vecs.size();
    m = midpoints/vecs.size();
    vm.push_back(v);
    vm.push_back(m);
  }
  //std::cout << v << ", " << m << "\n";
  return vm;
}

Cartesian GetClosestSplinePoint(const std::vector<Cartesian> &carts, const SplineInfo &splineinfo);
Cartesian GetClosestSplinePoint(const Cartesian &cart, const SplineInfo &splineinfo);

std::vector<Cartesian> GetBasePairEnds(mmdb::Residue* res1, mmdb::Residue* res2, const SplineInfo &splineinfo){
  std::vector<Cartesian> carts(2);
  mmdb::Atom* c11 = res1->GetAtom("C1\'");
  mmdb::Atom* c21 = res1->GetAtom("C2\'");
  mmdb::Atom* c31 = res1->GetAtom("C3\'");
  mmdb::Atom* c41 = res1->GetAtom("C4\'");
  mmdb::Atom* o41 = res1->GetAtom("O4\'");

  if(!c11) c11 = res1->GetAtom("C1*");
  if(!c21) c21 = res1->GetAtom("C2*");
  if(!c31) c31 = res1->GetAtom("C3*");
  if(!c41) c41 = res1->GetAtom("C4*");
  if(!o41) o41 = res1->GetAtom("O4*");

  mmdb::Atom* c12 = res2->GetAtom("C1\'");
  mmdb::Atom* c22 = res2->GetAtom("C2\'");
  mmdb::Atom* c32 = res2->GetAtom("C3\'");
  mmdb::Atom* c42 = res2->GetAtom("C4\'");
  mmdb::Atom* o42 = res2->GetAtom("O4\'");

  if(!c12) c12 = res2->GetAtom("C1*");
  if(!c22) c22 = res2->GetAtom("C2*");
  if(!c32) c32 = res2->GetAtom("C3*");
  if(!c42) c42 = res2->GetAtom("C4*");
  if(!o42) o42 = res2->GetAtom("O4*");

  mmdb::Atom* c51 = res1->GetAtom("C5\'");
  if(!c51) c51 = res1->GetAtom("C5*");
  mmdb::Atom* c52 = res2->GetAtom("C5\'");
  if(!c52) c52 = res2->GetAtom("C5*");
  mmdb::Atom* p1 = res1->GetAtom("P");
  mmdb::Atom* p2 = res2->GetAtom("P");
  Cartesian backbone1, backbone2;
  if(c51&&c52){
    backbone1 = Cartesian(c51->x,c51->y,c51->z);
    backbone2 = Cartesian(c52->x,c52->y,c52->z);
  } else if(p1&&p2){
    backbone1 = Cartesian(p1->x,p1->y,p1->z);
    backbone2 = Cartesian(p2->x,p2->y,p2->z);
  } else {
    carts.clear();
    return carts;
  }

  if((c11||c21||c31||c41||o41)&&(c12||c22||c32||c42||o42)){
    std::vector<Cartesian> carts1;
    std::vector<Cartesian> carts2;
    if(c11) carts1.push_back(Cartesian(c11->x,c11->y,c11->z));
    if(c21) carts1.push_back(Cartesian(c21->x,c21->y,c21->z));
    if(c31) carts1.push_back(Cartesian(c31->x,c31->y,c31->z));
    if(c41) carts1.push_back(Cartesian(c41->x,c41->y,c41->z));
    if(o41) carts1.push_back(Cartesian(o41->x,o41->y,o41->z));
    if(c12) carts2.push_back(Cartesian(c12->x,c12->y,c12->z));
    if(c22) carts2.push_back(Cartesian(c22->x,c22->y,c22->z));
    if(c32) carts2.push_back(Cartesian(c32->x,c32->y,c32->z));
    if(c42) carts2.push_back(Cartesian(c42->x,c42->y,c42->z));
    if(o42) carts2.push_back(Cartesian(o42->x,o42->y,o42->z));
    carts[0] = Cartesian::MidPoint(carts1);
    carts[1] = Cartesian::MidPoint(carts2);

    Cartesian cart1 = carts[0];
    Cartesian cart2 = carts[1];
    if((cart1-cart2).length()>9.0){
       std::vector<Cartesian> carts_tmp(2);
       Cartesian M = Cartesian::MidPoint(cart1,cart2);

       Cartesian MP = cart1 - M;
       Cartesian MC = backbone1 - M;

       Cartesian MPnorm = MP;
       MPnorm.normalize();
       double l = Cartesian::DotProduct(MP,MC)/MP.length();
        
       carts_tmp[0] = M;
       carts_tmp[1] = M+l*MPnorm;
       //carts[0] = GetClosestSplinePoint(carts_tmp,splineinfo);
       carts[0] = GetClosestSplinePoint(carts_tmp[1],splineinfo);

       MP = cart2 - M;
       MC = backbone2 - M;

       MPnorm = MP;
       MPnorm.normalize();
       l = Cartesian::DotProduct(MP,MC)/MP.length();

       carts_tmp[1] = M+l*MPnorm;
       //carts[1] = GetClosestSplinePoint(carts_tmp,splineinfo);
       carts[1] = GetClosestSplinePoint(carts_tmp[1],splineinfo);

    }
  } else {
    /*
    std::cout << "Do not have all atoms\n";
    if(c11) std::cout << "C11\n"; std::cout.flush();
    if(c21) std::cout << "C21\n"; std::cout.flush();
    if(c31) std::cout << "C31\n"; std::cout.flush();
    if(c41) std::cout << "C41\n"; std::cout.flush();
    if(o41) std::cout << "O41\n"; std::cout.flush();
    if(c12) std::cout << "C12\n"; std::cout.flush();
    if(c22) std::cout << "C22\n"; std::cout.flush();
    if(c32) std::cout << "C32\n"; std::cout.flush();
    if(c42) std::cout << "C42\n"; std::cout.flush();
    if(o42) std::cout << "O42\n"; std::cout.flush();
    */
    carts.clear();
  }

  return carts;
}

void DrawBaseBlockInt(PolyCollection *polys, mmdb::Residue* res1, const double *col1, const CParamsManager &params, int selHnd, mmdb::Atom* n1,mmdb::Atom* c2,mmdb::Atom* n3,mmdb::Atom* c4,mmdb::Atom* c5,mmdb::Atom* c6,mmdb::Atom* n7,mmdb::Atom* c8,mmdb::Atom* n9){

    //float thickness = params.GetFloat("cylinder_width")-0.02;
    float thickness = params.GetFloat("base_block_thickness")-0.02;
    if (thickness<0.02)thickness=0.1;
 
    TriangleElement* tri;
    if(n1&&c2&&n3&&c4&&c5&&c6){
      Cartesian n1cart(n1->x,n1->y,n1->z);
      Cartesian c2cart(c2->x,c2->y,c2->z);
      Cartesian n3cart(n3->x,n3->y,n3->z);
      Cartesian c4cart(c4->x,c4->y,c4->z);
      Cartesian c5cart(c5->x,c5->y,c5->z);
      Cartesian c6cart(c6->x,c6->y,c6->z);

      std::vector <Cartesian> carts;
      Cartesian up = Cartesian::CrossProduct(n1cart-c2cart,c2cart-n3cart);
      up.normalize(thickness);

      carts.push_back(n1cart+up);
      carts.push_back(c2cart+up);
      carts.push_back(n3cart+up);
      carts.push_back(c4cart+up);
      carts.push_back(c5cart+up);
      carts.push_back(c6cart+up);
      Cartesian midpoint = Cartesian::MidPoint(carts);

      carts.clear();
      carts.push_back(midpoint);
      carts.push_back(n1cart+up);
      carts.push_back(c2cart+up);
      tri = new TriangleElement(carts,col1,Cartesian::MidPoint(carts),1.0);
      polys->add_primitive(tri);

      carts.clear();
      carts.push_back(midpoint);
      carts.push_back(c2cart+up);
      carts.push_back(n3cart+up);
      tri = new TriangleElement(carts,col1,Cartesian::MidPoint(carts),1.0);
      polys->add_primitive(tri);

      carts.clear();
      carts.push_back(midpoint);
      carts.push_back(n3cart+up);
      carts.push_back(c4cart+up);
      tri = new TriangleElement(carts,col1,Cartesian::MidPoint(carts),1.0);
      polys->add_primitive(tri);

      carts.clear();
      carts.push_back(midpoint);
      carts.push_back(c4cart+up);
      carts.push_back(c5cart+up);
      tri = new TriangleElement(carts,col1,Cartesian::MidPoint(carts),1.0);
      polys->add_primitive(tri);

      carts.clear();
      carts.push_back(midpoint);
      carts.push_back(c5cart+up);
      carts.push_back(c6cart+up);
      tri = new TriangleElement(carts,col1,Cartesian::MidPoint(carts),1.0);
      polys->add_primitive(tri);

      carts.clear();
      carts.push_back(midpoint);
      carts.push_back(c6cart+up);
      carts.push_back(n1cart+up);
      tri = new TriangleElement(carts,col1,Cartesian::MidPoint(carts),1.0);
      polys->add_primitive(tri);

      carts.clear();
      carts.push_back(c6cart-up);
      carts.push_back(c5cart-up);
      carts.push_back(c4cart-up);
      carts.push_back(n3cart-up);
      carts.push_back(c2cart-up);
      carts.push_back(n1cart-up);
      carts.push_back(c6cart-up);
      midpoint = Cartesian::MidPoint(carts);

      carts.clear();
      carts.push_back(midpoint);
      carts.push_back(c2cart-up);
      carts.push_back(n1cart-up);
      tri = new TriangleElement(carts,col1,Cartesian::MidPoint(carts),1.0);
      polys->add_primitive(tri);

      carts.clear();
      carts.push_back(midpoint);
      carts.push_back(n3cart-up);
      carts.push_back(c2cart-up);
      tri = new TriangleElement(carts,col1,Cartesian::MidPoint(carts),1.0);
      polys->add_primitive(tri);

      carts.clear();
      carts.push_back(midpoint);
      carts.push_back(c4cart-up);
      carts.push_back(n3cart-up);
      tri = new TriangleElement(carts,col1,Cartesian::MidPoint(carts),1.0);
      polys->add_primitive(tri);

      carts.clear();
      carts.push_back(midpoint);
      carts.push_back(c5cart-up);
      carts.push_back(c4cart-up);
      tri = new TriangleElement(carts,col1,Cartesian::MidPoint(carts),1.0);
      polys->add_primitive(tri);

      carts.clear();
      carts.push_back(midpoint);
      carts.push_back(c6cart-up);
      carts.push_back(c5cart-up);
      tri = new TriangleElement(carts,col1,Cartesian::MidPoint(carts),1.0);
      polys->add_primitive(tri);

      carts.clear();
      carts.push_back(midpoint);
      carts.push_back(n1cart-up);
      carts.push_back(c6cart-up);
      tri = new TriangleElement(carts,col1,Cartesian::MidPoint(carts),1.0);
      polys->add_primitive(tri);

      carts.clear();
      carts.push_back(n1cart-up);
      carts.push_back(c2cart-up);
      carts.push_back(c2cart+up);
      //carts.push_back(n1cart+up);
      tri = new TriangleElement(carts,col1,Cartesian::MidPoint(carts),1.0);
      polys->add_primitive(tri);

      carts.clear();
      carts.push_back(n1cart-up);
      //carts.push_back(c2cart-up);
      carts.push_back(c2cart+up);
      carts.push_back(n1cart+up);
      tri = new TriangleElement(carts,col1,Cartesian::MidPoint(carts),1.0);
      polys->add_primitive(tri);

      carts.clear();
      carts.push_back(c2cart-up);
      carts.push_back(n3cart-up);
      carts.push_back(n3cart+up);
      //carts.push_back(c2cart+up);
      tri = new TriangleElement(carts,col1,Cartesian::MidPoint(carts),1.0);
      polys->add_primitive(tri);

      carts.clear();
      carts.push_back(c2cart-up);
      //carts.push_back(n3cart-up);
      carts.push_back(n3cart+up);
      carts.push_back(c2cart+up);
      tri = new TriangleElement(carts,col1,Cartesian::MidPoint(carts),1.0);
      polys->add_primitive(tri);

      carts.clear();
      carts.push_back(c5cart-up);
      carts.push_back(c6cart-up);
      carts.push_back(c6cart+up);
      //carts.push_back(c5cart+up);
      tri = new TriangleElement(carts,col1,Cartesian::MidPoint(carts),1.0);
      polys->add_primitive(tri);

      carts.clear();
      carts.push_back(c5cart-up);
      //carts.push_back(c6cart-up);
      carts.push_back(c6cart+up);
      carts.push_back(c5cart+up);
      tri = new TriangleElement(carts,col1,Cartesian::MidPoint(carts),1.0);
      polys->add_primitive(tri);

      carts.clear();
      carts.push_back(c6cart-up);
      carts.push_back(n1cart-up);
      carts.push_back(n1cart+up);
      //carts.push_back(c6cart+up);
      tri = new TriangleElement(carts,col1,Cartesian::MidPoint(carts),1.0);
      polys->add_primitive(tri);

      carts.clear();
      carts.push_back(c6cart-up);
      //carts.push_back(n1cart-up);
      carts.push_back(n1cart+up);
      carts.push_back(c6cart+up);
      tri = new TriangleElement(carts,col1,Cartesian::MidPoint(carts),1.0);
      polys->add_primitive(tri);

      carts.clear();
      carts.push_back(n3cart-up);
      carts.push_back(c4cart-up);
      carts.push_back(c4cart+up);
      //carts.push_back(n3cart+up);
      tri = new TriangleElement(carts,col1,Cartesian::MidPoint(carts),1.0);
      polys->add_primitive(tri);

      carts.clear();
      carts.push_back(n3cart-up);
      //carts.push_back(c4cart-up);
      carts.push_back(c4cart+up);
      carts.push_back(n3cart+up);
      tri = new TriangleElement(carts,col1,Cartesian::MidPoint(carts),1.0);
      polys->add_primitive(tri);

      if(n9&&c8&&n7){
        Cartesian n9cart(n9->x,n9->y,n9->z);
        Cartesian c8cart(c8->x,c8->y,c8->z);
        Cartesian n7cart(n7->x,n7->y,n7->z);

        carts.clear();
        carts.push_back(c5cart+up);
        carts.push_back(c4cart+up);
        carts.push_back(n9cart+up);
        carts.push_back(c8cart+up);
        carts.push_back(n7cart+up);
        carts.push_back(c5cart+up);
        midpoint = Cartesian::MidPoint(carts);

        carts.clear();
        carts.push_back(midpoint);
        carts.push_back(c5cart+up);
        carts.push_back(c4cart+up);
        tri = new TriangleElement(carts,col1,Cartesian::MidPoint(carts),1.0);
        polys->add_primitive(tri);

        carts.clear();
        carts.push_back(midpoint);
        carts.push_back(c4cart+up);
        carts.push_back(n9cart+up);
        tri = new TriangleElement(carts,col1,Cartesian::MidPoint(carts),1.0);
        polys->add_primitive(tri);

        carts.clear();
        carts.push_back(midpoint);
        carts.push_back(n9cart+up);
        carts.push_back(c8cart+up);
        tri = new TriangleElement(carts,col1,Cartesian::MidPoint(carts),1.0);
        polys->add_primitive(tri);

        carts.clear();
        carts.push_back(midpoint);
        carts.push_back(c8cart+up);
        carts.push_back(n7cart+up);
        tri = new TriangleElement(carts,col1,Cartesian::MidPoint(carts),1.0);
        polys->add_primitive(tri);

        carts.clear();
        carts.push_back(midpoint);
        carts.push_back(n7cart+up);
        carts.push_back(c5cart+up);
        tri = new TriangleElement(carts,col1,Cartesian::MidPoint(carts),1.0);
        polys->add_primitive(tri);

        carts.clear();
        carts.push_back(n7cart-up);
        carts.push_back(c8cart-up);
        carts.push_back(n9cart-up);
        carts.push_back(c4cart-up);
        carts.push_back(c5cart-up);
        carts.push_back(n7cart-up);
        midpoint = Cartesian::MidPoint(carts);

        carts.clear();
        carts.push_back(midpoint);
        carts.push_back(c4cart-up);
        carts.push_back(c5cart-up);
        tri = new TriangleElement(carts,col1,Cartesian::MidPoint(carts),1.0);
        polys->add_primitive(tri);

        carts.clear();
        carts.push_back(midpoint);
        carts.push_back(n9cart-up);
        carts.push_back(c4cart-up);
        tri = new TriangleElement(carts,col1,Cartesian::MidPoint(carts),1.0);
        polys->add_primitive(tri);

        carts.clear();
        carts.push_back(midpoint);
        carts.push_back(c8cart-up);
        carts.push_back(n9cart-up);
        tri = new TriangleElement(carts,col1,Cartesian::MidPoint(carts),1.0);
        polys->add_primitive(tri);

        carts.clear();
        carts.push_back(midpoint);
        carts.push_back(n7cart-up);
        carts.push_back(c8cart-up);
        tri = new TriangleElement(carts,col1,Cartesian::MidPoint(carts),1.0);
        polys->add_primitive(tri);

        carts.clear();
        carts.push_back(midpoint);
        carts.push_back(c5cart-up);
        carts.push_back(n7cart-up);
        tri = new TriangleElement(carts,col1,Cartesian::MidPoint(carts),1.0);
        polys->add_primitive(tri);

        carts.clear();
        carts.push_back(c4cart-up);
        carts.push_back(n9cart-up);
        carts.push_back(n9cart+up);
        //carts.push_back(c4cart+up);
        tri = new TriangleElement(carts,col1,Cartesian::MidPoint(carts),1.0);
        polys->add_primitive(tri);

        carts.clear();
        carts.push_back(c4cart-up);
        //carts.push_back(n9cart-up);
        carts.push_back(n9cart+up);
        carts.push_back(c4cart+up);
        tri = new TriangleElement(carts,col1,Cartesian::MidPoint(carts),1.0);
        polys->add_primitive(tri);

        carts.clear();
        carts.push_back(n9cart-up);
        carts.push_back(c8cart-up);
        carts.push_back(c8cart+up);
        //carts.push_back(n9cart+up);
        tri = new TriangleElement(carts,col1,Cartesian::MidPoint(carts),1.0);
        polys->add_primitive(tri);

        carts.clear();
        carts.push_back(n9cart-up);
        //carts.push_back(c8cart-up);
        carts.push_back(c8cart+up);
        carts.push_back(n9cart+up);
        tri = new TriangleElement(carts,col1,Cartesian::MidPoint(carts),1.0);
        polys->add_primitive(tri);

        carts.clear();
        carts.push_back(c8cart-up);
        carts.push_back(n7cart-up);
        carts.push_back(n7cart+up);
        //carts.push_back(c8cart+up);
        tri = new TriangleElement(carts,col1,Cartesian::MidPoint(carts),1.0);
        polys->add_primitive(tri);

        carts.clear();
        carts.push_back(c8cart-up);
        //carts.push_back(n7cart-up);
        carts.push_back(n7cart+up);
        carts.push_back(c8cart+up);
        tri = new TriangleElement(carts,col1,Cartesian::MidPoint(carts),1.0);
        polys->add_primitive(tri);

        carts.clear();
        carts.push_back(n7cart-up);
        carts.push_back(c5cart-up);
        carts.push_back(c5cart+up);
        //carts.push_back(n7cart+up);
        tri = new TriangleElement(carts,col1,Cartesian::MidPoint(carts),1.0);
        polys->add_primitive(tri);

        carts.clear();
        carts.push_back(n7cart-up);
        //carts.push_back(c5cart-up);
        carts.push_back(c5cart+up);
        carts.push_back(n7cart+up);
        tri = new TriangleElement(carts,col1,Cartesian::MidPoint(carts),1.0);
        polys->add_primitive(tri);
      }else{
        carts.clear();
        carts.push_back(c4cart-up);
        carts.push_back(c5cart-up);
        carts.push_back(c5cart+up);
        //carts.push_back(c4cart+up);
        tri = new TriangleElement(carts,col1,Cartesian::MidPoint(carts),1.0);
        polys->add_primitive(tri);

        carts.clear();
        carts.push_back(c4cart-up);
        //carts.push_back(c5cart-up);
        carts.push_back(c5cart+up);
        carts.push_back(c4cart+up);
        tri = new TriangleElement(carts,col1,Cartesian::MidPoint(carts),1.0);
        polys->add_primitive(tri);
      }
    }
}

void DrawBaseBlock(PolyCollection *polys, mmdb::Residue* res1, const double *col1, const CParamsManager &params, int selHnd ){

    mmdb::Atom* n1 = res1->GetAtom("N1");
    mmdb::Atom* c2 = res1->GetAtom("C2");
    mmdb::Atom* n3 = res1->GetAtom("N3");
    mmdb::Atom* c4 = res1->GetAtom("C4");
    mmdb::Atom* c5 = res1->GetAtom("C5");
    mmdb::Atom* c6 = res1->GetAtom("C6");
    mmdb::Atom* n7 = res1->GetAtom("N7");
    mmdb::Atom* c8 = res1->GetAtom("C8");
    mmdb::Atom* n9 = res1->GetAtom("N9");
    if(n1&&c2&&n3&&c4&&c5&&c6){
      DrawBaseBlockInt(polys, res1, col1, params, selHnd, n1,c2,n3,c4,c5,c6,n7,c8,n9 );
      return;
    }

    n1 = res1->GetAtom("N1",NULL,"A");
    c2 = res1->GetAtom("C2",NULL,"A");
    n3 = res1->GetAtom("N3",NULL,"A");
    c4 = res1->GetAtom("C4",NULL,"A");
    c5 = res1->GetAtom("C5",NULL,"A");
    c6 = res1->GetAtom("C6",NULL,"A");
    n7 = res1->GetAtom("N7",NULL,"A");
    c8 = res1->GetAtom("C8",NULL,"A");
    n9 = res1->GetAtom("N9",NULL,"A");
    if(n1&&c2&&n3&&c4&&c5&&c6){
      if(n1->isInSelection(selHnd)&&c2->isInSelection(selHnd)&&n3->isInSelection(selHnd)&&c4->isInSelection(selHnd)&&c5->isInSelection(selHnd)&&c6->isInSelection(selHnd)){
        DrawBaseBlockInt(polys, res1, col1, params, selHnd, n1,c2,n3,c4,c5,c6,n7,c8,n9 );
      }
    }
    n1 = res1->GetAtom("N1",NULL,"B");
    c2 = res1->GetAtom("C2",NULL,"B");
    n3 = res1->GetAtom("N3",NULL,"B");
    c4 = res1->GetAtom("C4",NULL,"B");
    c5 = res1->GetAtom("C5",NULL,"B");
    c6 = res1->GetAtom("C6",NULL,"B");
    n7 = res1->GetAtom("N7",NULL,"B");
    c8 = res1->GetAtom("C8",NULL,"B");
    n9 = res1->GetAtom("N9",NULL,"B");
    if(n1&&c2&&n3&&c4&&c5&&c6){
      if(n1->isInSelection(selHnd)&&c2->isInSelection(selHnd)&&n3->isInSelection(selHnd)&&c4->isInSelection(selHnd)&&c5->isInSelection(selHnd)&&c6->isInSelection(selHnd)){
        DrawBaseBlockInt(polys, res1, col1, params, selHnd, n1,c2,n3,c4,c5,c6,n7,c8,n9 );
      }
    }
    n1 = res1->GetAtom("N1",NULL,"C");
    c2 = res1->GetAtom("C2",NULL,"C");
    n3 = res1->GetAtom("N3",NULL,"C");
    c4 = res1->GetAtom("C4",NULL,"C");
    c5 = res1->GetAtom("C5",NULL,"C");
    c6 = res1->GetAtom("C6",NULL,"C");
    n7 = res1->GetAtom("N7",NULL,"C");
    c8 = res1->GetAtom("C8",NULL,"C");
    n9 = res1->GetAtom("N9",NULL,"C");
    if(n1&&c2&&n3&&c4&&c5&&c6){
      if(n1->isInSelection(selHnd)&&c2->isInSelection(selHnd)&&n3->isInSelection(selHnd)&&c4->isInSelection(selHnd)&&c5->isInSelection(selHnd)&&c6->isInSelection(selHnd)){
        DrawBaseBlockInt(polys, res1, col1, params, selHnd, n1,c2,n3,c4,c5,c6,n7,c8,n9 );
      }
    }
}

std::vector<Cartesian> DrawSugarBlockInt(TriangleCollection *tris, CylinderCollection *cyls, mmdb::Residue* res1, const double *col1, const CParamsManager &params, const CParamsManager &global_params, int selHnd, mmdb::Atom* at1,mmdb::Atom* at2,mmdb::Atom* at3,mmdb::Atom* at4,mmdb::Atom* at5,mmdb::Atom* at6, mmdb::Residue* res, mmdb::Manager *molHnd, bool two_colour, float sugar_block_thickness, float sugar_block_scale = 1.0){
  // FIXME, likely broken for F6P,RIB and friends.
  //std::cout << "Drawing " << res->GetResName() << "\n";
  std::vector<Cartesian> retval;

  // Should not be able to get here with these not already registered.
  int udd_C1X = molHnd->GetUDDHandle ( UDR_RESIDUE,"GLYCO_BLOCK_C1X" );
  int udd_C1Y = molHnd->GetUDDHandle ( UDR_RESIDUE,"GLYCO_BLOCK_C1Y" );
  int udd_C1Z = molHnd->GetUDDHandle ( UDR_RESIDUE,"GLYCO_BLOCK_C1Z" );
  int udd_C3X = molHnd->GetUDDHandle ( UDR_RESIDUE,"GLYCO_BLOCK_C3X" );
  int udd_C3Y = molHnd->GetUDDHandle ( UDR_RESIDUE,"GLYCO_BLOCK_C3Y" );
  int udd_C3Z = molHnd->GetUDDHandle ( UDR_RESIDUE,"GLYCO_BLOCK_C3Z" );
  int udd_C4X = molHnd->GetUDDHandle ( UDR_RESIDUE,"GLYCO_BLOCK_C4X" );
  int udd_C4Y = molHnd->GetUDDHandle ( UDR_RESIDUE,"GLYCO_BLOCK_C4Y" );
  int udd_C4Z = molHnd->GetUDDHandle ( UDR_RESIDUE,"GLYCO_BLOCK_C4Z" );
  int udd_C5X = molHnd->GetUDDHandle ( UDR_RESIDUE,"GLYCO_BLOCK_C5X" );
  int udd_C5Y = molHnd->GetUDDHandle ( UDR_RESIDUE,"GLYCO_BLOCK_C5Y" );
  int udd_C5Z = molHnd->GetUDDHandle ( UDR_RESIDUE,"GLYCO_BLOCK_C5Z" );
  int udd_C2X = molHnd->GetUDDHandle ( UDR_RESIDUE,"GLYCO_BLOCK_C2X" );
  int udd_C2Y = molHnd->GetUDDHandle ( UDR_RESIDUE,"GLYCO_BLOCK_C2Y" );
  int udd_C2Z = molHnd->GetUDDHandle ( UDR_RESIDUE,"GLYCO_BLOCK_C2Z" );

  if(at1&&at2&&at3&&at4&&at5){
    //std::cout << "DrawSugarBlockInt have at least 5 atoms.\n";
    Cartesian xaxis(1,0,0);
    Cartesian yaxis(0,1,0);
    Cartesian zaxis(0,0,1);
    std::vector<Cartesian> norm_vec;
    std::vector<Cartesian> mid_vec;
    Cartesian cat1 = Cartesian(at1->x,at1->y,at1->z);
    Cartesian cat2 = Cartesian(at2->x,at2->y,at2->z);
    Cartesian cat3 = Cartesian(at3->x,at3->y,at3->z);
    Cartesian cat4 = Cartesian(at4->x,at4->y,at4->z);
    Cartesian cat5 = Cartesian(at5->x,at5->y,at5->z);
    double d12 = (cat2-cat1).length();
    double d23 = (cat3-cat2).length();
    double d34 = (cat4-cat3).length();
    double d45 = (cat5-cat4).length();
    Cartesian c123 = Cartesian::CrossProduct((cat2-cat1),(cat3-cat2));
    Cartesian c234 = Cartesian::CrossProduct((cat3-cat2),(cat4-cat3));
    Cartesian c345 = Cartesian::CrossProduct((cat4-cat3),(cat5-cat4));
    if((d12>1.0&&d12<3.2)&&(d23>1.0&&d23<3.2)&&(d34>1.0&&d34<3.2)&&(d45>1.0&&d45<3.2)){
      if(at6){
        //std::cout << "Six membered test.\n";
        Cartesian cat6 = Cartesian(at6->x,at6->y,at6->z);
        double d56 = (cat6-cat5).length();
        double d61 = (cat6-cat1).length();
        if((d61>1.0&&d61<3.2)&&(d56>1.0&&d56<3.2)){
          norm_vec.push_back(c123);
          norm_vec.push_back(c234);
          norm_vec.push_back(c345);
          mid_vec.push_back(cat1);
          mid_vec.push_back(cat2);
          mid_vec.push_back(cat3);
          mid_vec.push_back(cat4);
          mid_vec.push_back(cat5);
          Cartesian c456 = Cartesian::CrossProduct((cat5-cat4),(cat6-cat5));
          Cartesian c561 = Cartesian::CrossProduct((cat6-cat5),(cat1-cat6));
          norm_vec.push_back(c456);
          norm_vec.push_back(c561);
          mid_vec.push_back(cat6);
          //std::cout << "We have 6 membered ring\n";
        } else if(d61>4.0&&d61<7.0&&strncmp(res->GetResName(),"XLS",3)==0) {
          norm_vec.push_back(c123);
          norm_vec.push_back(c234);
          norm_vec.push_back(c345);
          mid_vec.push_back(cat1);
          mid_vec.push_back(cat2);
          mid_vec.push_back(cat3);
          mid_vec.push_back(cat4);
          mid_vec.push_back(cat5);
          Cartesian c456 = Cartesian::CrossProduct((cat5-cat4),(cat6-cat5));
          Cartesian c561 = Cartesian::CrossProduct((cat6-cat5),(cat1-cat6));
          norm_vec.push_back(c456);
          norm_vec.push_back(c561);
          mid_vec.push_back(cat6);
          std::cout << "Linear?\n";
        }
      } else {
        //std::cout << "Five membered test.\n";
        double d51 = (cat5-cat1).length();
        if((d51>1.0&&d51<3.2)){
          norm_vec.push_back(c123);
          norm_vec.push_back(c234);
          norm_vec.push_back(c345);
          mid_vec.push_back(cat1);
          mid_vec.push_back(cat2);
          mid_vec.push_back(cat3);
          mid_vec.push_back(cat4);
          mid_vec.push_back(cat5);
          Cartesian c451 = Cartesian::CrossProduct((cat5-cat4),(cat1-cat5));
          norm_vec.push_back(c451);
          //std::cout << "We have 5 membered ring\n";
        }
      }
    } else {
      std::cout << "Fail first test\n";
    }
   
    int nsectors = 20;
    int quality = global_params.GetInt("solid_quality");
    if(quality==0){
      nsectors = 20;
    } else if(quality==1){
      nsectors = 40;
    } else if(quality==2){
      nsectors = 60;
    }

    float sugar_block_radius = 1.4 * sugar_block_scale;

    if(norm_vec.size()>0){
      Cartesian normal = Cartesian::MidPoint(norm_vec);
      normal.normalize();
      Cartesian centre = Cartesian::MidPoint(mid_vec);
      retval.push_back(centre);
      retval.push_back(normal);
      std::vector<Cartesian> carts = std::vector<Cartesian>(2);
      carts[0] = centre - sugar_block_thickness*normal;
      carts[1] = centre + sugar_block_thickness*normal;
      if(strncmp(res->GetResName(),"BMA",3)==0||strncmp(res->GetResName(),"MAN",3)==0||strncmp(res->GetResName(),"GAL",3)==0||strncmp(res->GetResName(),"GLC",3)==0||strncmp(res->GetResName(),"BGC",3)==0){
        CylinderElement *cylinder = new CylinderElement(carts,col1,centre,sugar_block_radius,1.0,nsectors);
        cyls->add_primitive(cylinder);
      }
      Cartesian p = normal;
      p.normalize();
      Cartesian up;
      Cartesian right;
      if(fabs(Cartesian::DotProduct(xaxis,p))<0.95){
        up = Cartesian::CrossProduct(xaxis,p);
        up.normalize();
      }else if(fabs(Cartesian::DotProduct(yaxis,p))<0.95){
        up = Cartesian::CrossProduct(yaxis,p);
        up.normalize();
      } else {
        up = Cartesian::CrossProduct(zaxis,p);
        up.normalize();
      }
      right = Cartesian::CrossProduct(up,p);
      up.normalize(sugar_block_radius);
      right.normalize(sugar_block_radius);
      if(strncmp(res->GetResName(),"NAG",3)==0||strncmp(res->GetResName(),"NBG",3)==0||strncmp(res->GetResName(),"NGA",3)==0||strncmp(res->GetResName(),"NG1",3)==0||strncmp(res->GetResName(),"NG6",3)==0||strncmp(res->GetResName(),"GCS",3)==0||strncmp(res->GetResName(),"6MN",3)==0||strncmp(res->GetResName(),"GLP",3)==0||strncmp(res->GetResName(),"GP4",3)==0){

        bool diagonal = false;
        if(strncmp(res->GetResName(),"GCS",3)==0) {
         diagonal = true;
        }
        two_colour = two_colour && diagonal;
        double white[] = {1.0,1.0,1.0,1.0};
        white[3] = col1[3];

        Cartesian c1c4 = cat1-cat4;
        c1c4.normalize();
        Cartesian c1c4Perp = Cartesian::CrossProduct(normal,c1c4);
        c1c4Perp.normalize();
        Cartesian norm2 = Cartesian::CrossProduct(c1c4,c1c4Perp);
        c1c4.normalize(sugar_block_radius);
        c1c4Perp.normalize(sugar_block_radius);
        // FIXME - First approximation. Need to change based on shape drawn.
        Cartesian ls = centre+c1c4+c1c4Perp;
        Cartesian le = centre+c1c4-c1c4Perp;
        double t = DistanceBetweenPointAndLine(ls,le,cat1)[1];
        Cartesian pp = ls + t*(le-ls);
        res->PutUDData(udd_C1X,pp.get_x()); res->PutUDData(udd_C1Y,pp.get_y()); res->PutUDData(udd_C1Z,pp.get_z());

        Cartesian ls3 = centre+c1c4+c1c4Perp;
        Cartesian le3 = centre-c1c4+c1c4Perp;
        double t3 = DistanceBetweenPointAndLine(ls3,le3,cat3)[1];
        Cartesian pp3 = ls3 + t3*(le3-ls3);
        res->PutUDData(udd_C3X,pp3.get_x()); res->PutUDData(udd_C3Y,pp3.get_y()); res->PutUDData(udd_C3Z,pp3.get_z());

        ls = centre-c1c4+c1c4Perp;
        le = centre-c1c4-c1c4Perp;
        t = DistanceBetweenPointAndLine(ls,le,cat4)[1];
        pp = ls + t*(le-ls);
        res->PutUDData(udd_C4X,pp.get_x()); res->PutUDData(udd_C4Y,pp.get_y()); res->PutUDData(udd_C4Z,pp.get_z());

        Cartesian ls5 = centre+c1c4-c1c4Perp;
        Cartesian le5 = centre-c1c4-c1c4Perp;
        double t5 = DistanceBetweenPointAndLine(ls5,le5,cat5)[1];
        Cartesian pp5 = ls5 + t5*(le5-ls5);
        res->PutUDData(udd_C5X,pp5.get_x()); res->PutUDData(udd_C5Y,pp5.get_y()); res->PutUDData(udd_C5Z,pp5.get_z());

        std::vector<Cartesian> cartsTri = std::vector<Cartesian>(3);

        cartsTri[0] = centre - c1c4 + c1c4Perp - sugar_block_thickness*norm2;
        cartsTri[1] = centre + c1c4 + c1c4Perp - sugar_block_thickness*norm2;
        cartsTri[2] = centre - c1c4 - c1c4Perp - sugar_block_thickness*norm2;
        TriangleElement *tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);

        cartsTri[1] = centre + c1c4 - c1c4Perp - sugar_block_thickness*norm2;
        cartsTri[0] = centre + c1c4 + c1c4Perp - sugar_block_thickness*norm2;
        cartsTri[2] = centre - c1c4 - c1c4Perp - sugar_block_thickness*norm2;
        if(two_colour)
          tri = new TriangleElement(cartsTri,white,centre,1.0);
        else
          tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);

        cartsTri[1] = centre - c1c4 + c1c4Perp + sugar_block_thickness*norm2;
        cartsTri[0] = centre + c1c4 + c1c4Perp + sugar_block_thickness*norm2;
        cartsTri[2] = centre - c1c4 - c1c4Perp + sugar_block_thickness*norm2;
        tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);

        cartsTri[0] = centre + c1c4 - c1c4Perp + sugar_block_thickness*norm2;
        cartsTri[1] = centre + c1c4 + c1c4Perp + sugar_block_thickness*norm2;
        cartsTri[2] = centre - c1c4 - c1c4Perp + sugar_block_thickness*norm2;
        if(two_colour)
          tri = new TriangleElement(cartsTri,white,centre,1.0);
        else
          tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);

        cartsTri[0] = centre - c1c4 - c1c4Perp + sugar_block_thickness*norm2;
        cartsTri[1] = centre - c1c4 - c1c4Perp - sugar_block_thickness*norm2;
        cartsTri[2] = centre + c1c4 - c1c4Perp - sugar_block_thickness*norm2;
        tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);
        cartsTri[0] = centre - c1c4 - c1c4Perp + sugar_block_thickness*norm2;
        cartsTri[1] = centre + c1c4 - c1c4Perp - sugar_block_thickness*norm2;
        cartsTri[2] = centre + c1c4 - c1c4Perp + sugar_block_thickness*norm2;
        tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);

        cartsTri[1] = centre - c1c4 + c1c4Perp + sugar_block_thickness*norm2;
        cartsTri[0] = centre - c1c4 + c1c4Perp - sugar_block_thickness*norm2;
        cartsTri[2] = centre + c1c4 + c1c4Perp - sugar_block_thickness*norm2;
        tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);
        cartsTri[1] = centre - c1c4 + c1c4Perp + sugar_block_thickness*norm2;
        cartsTri[0] = centre + c1c4 + c1c4Perp - sugar_block_thickness*norm2;
        cartsTri[2] = centre + c1c4 + c1c4Perp + sugar_block_thickness*norm2;
        tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);

        cartsTri[0] = centre - c1c4 + c1c4Perp + sugar_block_thickness*norm2;
        cartsTri[1] = centre - c1c4 + c1c4Perp - sugar_block_thickness*norm2;
        cartsTri[2] = centre - c1c4 - c1c4Perp - sugar_block_thickness*norm2;
        tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);
        cartsTri[0] = centre - c1c4 + c1c4Perp + sugar_block_thickness*norm2;
        cartsTri[1] = centre - c1c4 - c1c4Perp - sugar_block_thickness*norm2;
        cartsTri[2] = centre - c1c4 - c1c4Perp + sugar_block_thickness*norm2;
        tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);

        cartsTri[1] = centre + c1c4 + c1c4Perp + sugar_block_thickness*norm2;
        cartsTri[0] = centre + c1c4 + c1c4Perp - sugar_block_thickness*norm2;
        cartsTri[2] = centre + c1c4 - c1c4Perp - sugar_block_thickness*norm2;
        tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);
        cartsTri[1] = centre + c1c4 + c1c4Perp + sugar_block_thickness*norm2;
        cartsTri[0] = centre + c1c4 - c1c4Perp - sugar_block_thickness*norm2;
        cartsTri[2] = centre + c1c4 - c1c4Perp + sugar_block_thickness*norm2;
        tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);

      } else if(strncmp(res->GetResName(),"BMA",3)==0||strncmp(res->GetResName(),"MAN",3)==0||strncmp(res->GetResName(),"GAL",3)==0||strncmp(res->GetResName(),"GLC",3)==0||strncmp(res->GetResName(),"BGC",3)==0){

        Cartesian cp = centre - cat1;
        cp.normalize();
        Cartesian cpx = Cartesian::CrossProduct(cp,normal);
        cpx.normalize();
        Cartesian cpp = Cartesian::CrossProduct(normal,cpx);
        cpp.normalize(right.length());
        Cartesian pp = centre - cpp;
        res->PutUDData(udd_C1X,pp.get_x()); res->PutUDData(udd_C1Y,pp.get_y()); res->PutUDData(udd_C1Z,pp.get_z());

        cp = centre - cat3;
        cp.normalize();
        cpx = Cartesian::CrossProduct(cp,normal);
        cpx.normalize();
        cpp = Cartesian::CrossProduct(normal,cpx);
        cpp.normalize(right.length());
        pp = centre - cpp;
        res->PutUDData(udd_C3X,pp.get_x()); res->PutUDData(udd_C3Y,pp.get_y()); res->PutUDData(udd_C3Z,pp.get_z());

        cp = centre - cat4;
        cp.normalize();
        cpx = Cartesian::CrossProduct(cp,normal);
        cpx.normalize();
        cpp = Cartesian::CrossProduct(normal,cpx);
        cpp.normalize(right.length());
        pp = centre - cpp;
        res->PutUDData(udd_C4X,pp.get_x()); res->PutUDData(udd_C4Y,pp.get_y()); res->PutUDData(udd_C4Z,pp.get_z());

        cp = centre - cat5;
        cp.normalize();
        cpx = Cartesian::CrossProduct(cp,normal);
        cpx.normalize();
        cpp = Cartesian::CrossProduct(normal,cpx);
        cpp.normalize(right.length());
        pp = centre - cpp;
        res->PutUDData(udd_C5X,pp.get_x()); res->PutUDData(udd_C5Y,pp.get_y()); res->PutUDData(udd_C5Z,pp.get_z());

        cp = centre - cat2;
        cp.normalize();
        cpx = Cartesian::CrossProduct(cp,normal);
        cpx.normalize();
        cpp = Cartesian::CrossProduct(normal,cpx);
        cpp.normalize(right.length());
        pp = centre - cpp;
        res->PutUDData(udd_C2X,pp.get_x()); res->PutUDData(udd_C2Y,pp.get_y()); res->PutUDData(udd_C2Z,pp.get_z());

        for(int j=0;j<360;j=j+360/nsectors){
          double theta = (double)j/360.0 * PIBY2;
          double theta2 = (double)(j+360/nsectors)/360.0 * PIBY2;
          double x1 = cos(theta);
          double y1 = sin(theta);
          double x2 = cos(theta2);
          double y2 = sin(theta2);
          std::vector<Cartesian> cartsTri = std::vector<Cartesian>(3);
          cartsTri[0] = x1*up + y1*right + centre - sugar_block_thickness*normal;
          cartsTri[1] = x2*up + y2*right + centre - sugar_block_thickness*normal;
          cartsTri[2] = centre - sugar_block_thickness*normal;
          TriangleElement *tri = new TriangleElement(cartsTri,col1,centre,1.0);
          tris->add_primitive(tri);
          cartsTri[0] = x2*up + y2*right + centre + sugar_block_thickness*normal;
          cartsTri[1] = x1*up + y1*right + centre + sugar_block_thickness*normal;
          cartsTri[2] = centre + sugar_block_thickness*normal;
          tri = new TriangleElement(cartsTri,col1,centre,1.0);
          tris->add_primitive(tri);
        }
      } else if(strncmp(res->GetResName(),"FCA",3)==0
      ||strncmp(res->GetResName(),"FCB",3)==0
      ||strncmp(res->GetResName(),"FUC",3)==0
     ){
        /* Draw a triangle */
        std::cout << "Draw a triangle" << std::endl;
        Cartesian c1c4 = cat1-cat4;
        c1c4.normalize();
        Cartesian c1c4Perp = Cartesian::CrossProduct(normal,c1c4);
        c1c4Perp.normalize();
        Cartesian norm2 = Cartesian::CrossProduct(c1c4,c1c4Perp);
        c1c4.normalize(sugar_block_radius);
        c1c4Perp.normalize(sugar_block_radius);

        double ytri = -0.5;
        double xtri = -0.8660254037844387;
        Cartesian ptri = ytri*c1c4Perp + xtri*c1c4;
        Cartesian mtri = ytri*c1c4Perp - xtri*c1c4;

        Cartesian ls = centre+mtri;
        Cartesian le = centre+c1c4Perp;
        double t = DistanceBetweenPointAndLine(ls,le,cat1)[1];
        Cartesian pp = ls + t*(le-ls);
        res->PutUDData(udd_C1X,pp.get_x()); res->PutUDData(udd_C1Y,pp.get_y()); res->PutUDData(udd_C1Z,pp.get_z());

        t = DistanceBetweenPointAndLine(ls,le,cat2)[1];
        pp = ls + t*(le-ls);
        res->PutUDData(udd_C2X,pp.get_x()); res->PutUDData(udd_C2Y,pp.get_y()); res->PutUDData(udd_C2Z,pp.get_z());

        ls = centre+ptri;
        le = centre+c1c4Perp;
        t = DistanceBetweenPointAndLine(ls,le,cat3)[1];
        pp = ls + t*(le-ls);
        res->PutUDData(udd_C3X,pp.get_x()); res->PutUDData(udd_C3Y,pp.get_z()); res->PutUDData(udd_C3Z,pp.get_z());

        t = DistanceBetweenPointAndLine(ls,le,cat4)[1];
        pp = ls + t*(le-ls);
        res->PutUDData(udd_C4X,pp.get_x()); res->PutUDData(udd_C4Y,pp.get_y()); res->PutUDData(udd_C4Z,pp.get_z());

        ls = centre+ptri;
        le = centre+mtri;
        t = DistanceBetweenPointAndLine(ls,le,cat5)[1];
        pp = ls + t*(le-ls);
        res->PutUDData(udd_C5X,pp.get_x()); res->PutUDData(udd_C5Y,pp.get_y()); res->PutUDData(udd_C5Z,pp.get_z());

        std::vector<Cartesian> cartsTri = std::vector<Cartesian>(3);
        cartsTri[0] = centre + c1c4Perp - sugar_block_thickness*norm2;
        cartsTri[1] = centre + mtri - sugar_block_thickness*norm2;
        cartsTri[2] = centre + ptri - sugar_block_thickness*norm2;
        TriangleElement *tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);

        cartsTri[0] = centre + c1c4Perp + sugar_block_thickness*norm2;
        cartsTri[2] = centre + mtri + sugar_block_thickness*norm2;
        cartsTri[1] = centre + ptri + sugar_block_thickness*norm2;
        tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);

        cartsTri[0] = centre + mtri + sugar_block_thickness*norm2;
        cartsTri[2] = centre + mtri - sugar_block_thickness*norm2;
        cartsTri[1] = centre + ptri - sugar_block_thickness*norm2;
        tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);

        cartsTri[0] = centre + mtri + sugar_block_thickness*norm2;
        cartsTri[1] = centre + ptri + sugar_block_thickness*norm2;
        cartsTri[2] = centre + ptri - sugar_block_thickness*norm2;
        tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);

        cartsTri[0] = centre + c1c4Perp + sugar_block_thickness*norm2;
        cartsTri[2] = centre + ptri + sugar_block_thickness*norm2;
        cartsTri[1] = centre + ptri - sugar_block_thickness*norm2;
        tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);

        cartsTri[0] = centre + c1c4Perp + sugar_block_thickness*norm2;
        cartsTri[1] = centre + c1c4Perp - sugar_block_thickness*norm2;
        cartsTri[2] = centre + ptri - sugar_block_thickness*norm2;
        tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);

        cartsTri[0] = centre + c1c4Perp + sugar_block_thickness*norm2;
        cartsTri[1] = centre + mtri + sugar_block_thickness*norm2;
        cartsTri[2] = centre + mtri - sugar_block_thickness*norm2;
        tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);

        cartsTri[0] = centre + c1c4Perp + sugar_block_thickness*norm2;
        cartsTri[2] = centre + c1c4Perp - sugar_block_thickness*norm2;
        cartsTri[1] = centre + mtri - sugar_block_thickness*norm2;
        tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);

      } else if(strncmp(res->GetResName(),"BEM",3)==0
      ||strncmp(res->GetResName(),"GTR",3)==0
      ||strncmp(res->GetResName(),"ADA",3)==0
      ||strncmp(res->GetResName(),"DGU",3)==0
      ||strncmp(res->GetResName(),"KDN",3)==0
      ||strncmp(res->GetResName(),"SI3",3)==0
      ||strncmp(res->GetResName(),"NCC",3)==0

      ||strncmp(res->GetResName(),"IDR",3)==0
      ||strncmp(res->GetResName(),"GC4",3)==0
      ||strncmp(res->GetResName(),"GCD",3)==0
      ||strncmp(res->GetResName(),"GCU",3)==0
      ||strncmp(res->GetResName(),"GCV",3)==0
      ||strncmp(res->GetResName(),"GCW",3)==0
      ||strncmp(res->GetResName(),"IDS",3)==0
      ||strncmp(res->GetResName(),"REL",3)==0

     ){
        bool horizontal = false;
        bool vertical = false;
        bool invert_colour = false;
        double white[] = {1.0,1.0,1.0,1.0};
        white[3] = col1[3];
        if(strncmp(res->GetResName(),"IDR",3)==0
         ||strncmp(res->GetResName(),"GC4",3)==0
         ||strncmp(res->GetResName(),"GCD",3)==0
         ||strncmp(res->GetResName(),"GCU",3)==0
         ||strncmp(res->GetResName(),"GCV",3)==0
         ||strncmp(res->GetResName(),"GCW",3)==0
         ||strncmp(res->GetResName(),"IDS",3)==0
         ||strncmp(res->GetResName(),"REL",3)==0) {
         horizontal = true;
        }
        if(strncmp(res->GetResName(),"GTR",3)==0
         ||strncmp(res->GetResName(),"ADA",3)==0
         ||strncmp(res->GetResName(),"DGU",3)==0
         ||strncmp(res->GetResName(),"BEM",3)==0) {
         vertical = true;
        }
        if(strncmp(res->GetResName(),"IDR",3)==0
         ||strncmp(res->GetResName(),"BEM",3)==0) {
         invert_colour = true;
        }

        two_colour = two_colour && (horizontal||vertical);

        /* Draw a diamond */
        std::cout << "Draw a diamond" << std::endl;
        Cartesian c1c4 = cat1-cat4;
        c1c4.normalize();
        Cartesian c1c4Perp = Cartesian::CrossProduct(normal,c1c4);
        c1c4Perp.normalize();
        Cartesian norm2 = Cartesian::CrossProduct(c1c4,c1c4Perp);
        c1c4.normalize(sugar_block_radius);
        c1c4Perp.normalize(sugar_block_radius);
        // FIXME - TOTALLY UNTESTED
        Cartesian ls = centre+c1c4;
        Cartesian le = centre-c1c4Perp;
        double t = DistanceBetweenPointAndLine(ls,le,cat1)[1];
        Cartesian pp = ls + t*(le-ls);
        res->PutUDData(udd_C1X,pp.get_x()); res->PutUDData(udd_C1Y,pp.get_y()); res->PutUDData(udd_C1Z,pp.get_z());

        t = DistanceBetweenPointAndLine(ls,le,cat2)[1];
        pp = ls + t*(le-ls);
        res->PutUDData(udd_C2X,pp.get_x()); res->PutUDData(udd_C2Y,pp.get_y()); res->PutUDData(udd_C2Z,pp.get_z());

        ls = centre-c1c4;
        le = centre-c1c4Perp;
        t = DistanceBetweenPointAndLine(ls,le,cat3)[1];
        pp = ls + t*(le-ls);
        res->PutUDData(udd_C3X,pp.get_x()); res->PutUDData(udd_C3Y,pp.get_z()); res->PutUDData(udd_C3Z,pp.get_z());

        t = DistanceBetweenPointAndLine(ls,le,cat4)[1];
        pp = ls + t*(le-ls);
        res->PutUDData(udd_C4X,pp.get_x()); res->PutUDData(udd_C4Y,pp.get_y()); res->PutUDData(udd_C4Z,pp.get_z());

        ls = centre-c1c4;
        le = centre+c1c4Perp;
        t = DistanceBetweenPointAndLine(ls,le,cat5)[1];
        pp = ls + t*(le-ls);
        res->PutUDData(udd_C5X,pp.get_x()); res->PutUDData(udd_C5Y,pp.get_y()); res->PutUDData(udd_C5Z,pp.get_z());

        std::vector<Cartesian> cartsTri= std::vector<Cartesian>(3);
        TriangleElement *tri;

        if(horizontal){
          cartsTri[1] = centre - c1c4 - sugar_block_thickness*norm2;
          cartsTri[0] = centre + c1c4 - sugar_block_thickness*norm2;
          cartsTri[2] = centre + c1c4Perp - sugar_block_thickness*norm2;
          if(two_colour&&invert_colour)
            tri = new TriangleElement(cartsTri,white,centre,1.0);
          else
            tri = new TriangleElement(cartsTri,col1,centre,1.0);
          tris->add_primitive(tri);

          cartsTri[0] = centre + c1c4 - sugar_block_thickness*norm2;
          cartsTri[1] = centre - c1c4Perp - sugar_block_thickness*norm2;
          cartsTri[2] = centre - c1c4 - sugar_block_thickness*norm2;
          if(two_colour&&!invert_colour)
            tri = new TriangleElement(cartsTri,white,centre,1.0);
          else
            tri = new TriangleElement(cartsTri,col1,centre,1.0);
          tris->add_primitive(tri);

          cartsTri[1] = centre - c1c4 + sugar_block_thickness*norm2;
          cartsTri[2] = centre + c1c4 + sugar_block_thickness*norm2;
          cartsTri[0] = centre + c1c4Perp + sugar_block_thickness*norm2;
          if(two_colour&&invert_colour)
            tri = new TriangleElement(cartsTri,white,centre,1.0);
          else
            tri = new TriangleElement(cartsTri,col1,centre,1.0);
          tris->add_primitive(tri);

          cartsTri[0] = centre + c1c4 + sugar_block_thickness*norm2;
          cartsTri[2] = centre - c1c4Perp + sugar_block_thickness*norm2;
          cartsTri[1] = centre - c1c4 + sugar_block_thickness*norm2;
          if(two_colour&&!invert_colour)
            tri = new TriangleElement(cartsTri,white,centre,1.0);
          else
            tri = new TriangleElement(cartsTri,col1,centre,1.0);
          tris->add_primitive(tri);
        } else {
          cartsTri[0] = centre - c1c4 - sugar_block_thickness*norm2;
          cartsTri[1] = centre + c1c4Perp - sugar_block_thickness*norm2;
          cartsTri[2] = centre - c1c4Perp - sugar_block_thickness*norm2;
          if(two_colour&&!invert_colour)
            tri = new TriangleElement(cartsTri,white,centre,1.0);
          else
            tri = new TriangleElement(cartsTri,col1,centre,1.0);
          tris->add_primitive(tri);

          cartsTri[0] = centre + c1c4 - sugar_block_thickness*norm2;
          cartsTri[1] = centre - c1c4Perp - sugar_block_thickness*norm2;
          cartsTri[2] = centre + c1c4Perp - sugar_block_thickness*norm2;
          if(two_colour&&invert_colour)
            tri = new TriangleElement(cartsTri,white,centre,1.0);
          else
            tri = new TriangleElement(cartsTri,col1,centre,1.0);
          tris->add_primitive(tri);

          cartsTri[0] = centre - c1c4 + sugar_block_thickness*norm2;
          cartsTri[2] = centre + c1c4Perp + sugar_block_thickness*norm2;
          cartsTri[1] = centre - c1c4Perp + sugar_block_thickness*norm2;
          if(two_colour&&!invert_colour)
            tri = new TriangleElement(cartsTri,white,centre,1.0);
          else
            tri = new TriangleElement(cartsTri,col1,centre,1.0);
          tris->add_primitive(tri);

          cartsTri[0] = centre + c1c4 + sugar_block_thickness*norm2;
          cartsTri[2] = centre - c1c4Perp + sugar_block_thickness*norm2;
          cartsTri[1] = centre + c1c4Perp + sugar_block_thickness*norm2;
          if(two_colour&&invert_colour)
            tri = new TriangleElement(cartsTri,white,centre,1.0);
          else
            tri = new TriangleElement(cartsTri,col1,centre,1.0);
          tris->add_primitive(tri);
        }


        cartsTri[1] = centre + c1c4 + sugar_block_thickness*norm2;
        cartsTri[0] = centre + c1c4 - sugar_block_thickness*norm2;
        cartsTri[2] = centre - c1c4Perp - sugar_block_thickness*norm2;
        tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);

        cartsTri[0] = centre + c1c4 + sugar_block_thickness*norm2;
        cartsTri[1] = centre - c1c4Perp + sugar_block_thickness*norm2;
        cartsTri[2] = centre - c1c4Perp - sugar_block_thickness*norm2;
        tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);

        cartsTri[1] = centre - c1c4 + sugar_block_thickness*norm2;
        cartsTri[0] = centre - c1c4 - sugar_block_thickness*norm2;
        cartsTri[2] = centre + c1c4Perp - sugar_block_thickness*norm2;
        tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);

        cartsTri[0] = centre - c1c4 + sugar_block_thickness*norm2;
        cartsTri[1] = centre + c1c4Perp + sugar_block_thickness*norm2;
        cartsTri[2] = centre + c1c4Perp - sugar_block_thickness*norm2;
        tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);

        cartsTri[0] = centre - c1c4 + sugar_block_thickness*norm2;
        cartsTri[1] = centre - c1c4 - sugar_block_thickness*norm2;
        cartsTri[2] = centre - c1c4Perp - sugar_block_thickness*norm2;
        tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);

        cartsTri[1] = centre - c1c4 + sugar_block_thickness*norm2;
        cartsTri[0] = centre - c1c4Perp + sugar_block_thickness*norm2;
        cartsTri[2] = centre - c1c4Perp - sugar_block_thickness*norm2;
        tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);

        cartsTri[0] = centre + c1c4 + sugar_block_thickness*norm2;
        cartsTri[1] = centre + c1c4 - sugar_block_thickness*norm2;
        cartsTri[2] = centre + c1c4Perp - sugar_block_thickness*norm2;
        tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);

        cartsTri[1] = centre + c1c4 + sugar_block_thickness*norm2;
        cartsTri[0] = centre + c1c4Perp + sugar_block_thickness*norm2;
        cartsTri[2] = centre + c1c4Perp - sugar_block_thickness*norm2;
        tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);


      } else if(strncmp(res->GetResName(),"XLS",3)==0
      ||strncmp(res->GetResName(),"CXY",3)==0
      ||strncmp(res->GetResName(),"RBY",3)==0
      ||strncmp(res->GetResName(),"TDX",3)==0
      ||strncmp(res->GetResName(),"XYL",3)==0 /* Should all be orange */
      ||strncmp(res->GetResName(),"XYS",3)==0
     ){
        /* Draw a star! */
        std::cout << "Draw a star!" << std::endl;
        Cartesian c1c4 = cat1-cat4;
        c1c4.normalize();
        Cartesian c1c4Perp = Cartesian::CrossProduct(normal,c1c4);
        c1c4Perp.normalize();
        Cartesian norm2 = Cartesian::CrossProduct(c1c4,c1c4Perp);
        c1c4.normalize(sugar_block_radius);
        c1c4Perp.normalize(sugar_block_radius);

        // http://mathworld.wolfram.com/Pentagon.html
        double c1 = cos(2.0*M_PI/5.);
        double c2 = cos(M_PI/5.);
        double s1 = sin(2.0*M_PI/5.);
        double s2 = sin(4.0*M_PI/5.);
        //Cartesian p1(0,1,0);
        //Cartesian p2(s1,c1,0);
        //Cartesian p3(s2,-c2,0);
        //Cartesian p4(-s2,-c2,0);
        //Cartesian p5(-s1,c1,0);
        Cartesian p1 = c1c4Perp;
        Cartesian p2 = s1 * c1c4 + c1 * c1c4Perp;
        Cartesian p3 = s2 * c1c4 - c2 * c1c4Perp;
        Cartesian p4 = -s2 * c1c4 - c2 * c1c4Perp;
        Cartesian p5 = -s1 * c1c4 + c1 * c1c4Perp;

        double t6 = DistanceBetweenTwoLines(p1,p3,p5,p2)[1];
        double t7 = DistanceBetweenTwoLines(p1,p3,p2,p4)[1];
        double t8 = DistanceBetweenTwoLines(p2,p4,p3,p5)[1];
        double t9 = DistanceBetweenTwoLines(p4,p1,p3,p5)[1];
        double t10 = DistanceBetweenTwoLines(p5,p2,p4,p1)[1];

        Cartesian p6 = p1 + t6*(p3-p1);
        Cartesian p7 = p1 + t7*(p3-p1);
        Cartesian p8 = p2 + t8*(p4-p2);
        Cartesian p9 = p4 + t9*(p1-p4);
        Cartesian p10 = p5 + t6*(p2-p5);

        // FIXME - TOTALLY UNTESTED
        res->PutUDData(udd_C1X,(centre+p2).get_x());  res->PutUDData(udd_C1Y,(centre+p2).get_y());  res->PutUDData(udd_C1Z,(centre+p2).get_z());
        res->PutUDData(udd_C2X,(centre+p10).get_x()); res->PutUDData(udd_C2Y,(centre+p10).get_y()); res->PutUDData(udd_C2Z,(centre+p10).get_z());
        res->PutUDData(udd_C3X,(centre+p6).get_x());  res->PutUDData(udd_C3Y,(centre+p6).get_z());  res->PutUDData(udd_C3Z,(centre+p6).get_z());
        res->PutUDData(udd_C4X,(centre+p5).get_x());  res->PutUDData(udd_C4Y,(centre+p5).get_y());  res->PutUDData(udd_C4Z,(centre+p5).get_z());
        res->PutUDData(udd_C5X,(centre+p3).get_x());  res->PutUDData(udd_C5Y,(centre+p3).get_y());  res->PutUDData(udd_C5Z,(centre+p3).get_z());

        std::vector<Cartesian> cartsTri = std::vector<Cartesian>(3);
        cartsTri[1] = centre + p6 + sugar_block_thickness*norm2;
        cartsTri[0] = centre + p7 + sugar_block_thickness*norm2;
        cartsTri[2] = centre + sugar_block_thickness*norm2;
        TriangleElement *tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);
        
        cartsTri[1] = centre + p7 + sugar_block_thickness*norm2;
        cartsTri[0] = centre + p8 + sugar_block_thickness*norm2;
        cartsTri[2] = centre + sugar_block_thickness*norm2;
        tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);

        cartsTri[1] = centre + p8 + sugar_block_thickness*norm2;
        cartsTri[0] = centre + p9 + sugar_block_thickness*norm2;
        cartsTri[2] = centre + sugar_block_thickness*norm2;
        tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);

        cartsTri[1] = centre + p9 + sugar_block_thickness*norm2;
        cartsTri[0] = centre + p10 + sugar_block_thickness*norm2;
        cartsTri[2] = centre + sugar_block_thickness*norm2;
        tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);

        cartsTri[1] = centre + p10 + sugar_block_thickness*norm2;
        cartsTri[0] = centre + p6 + sugar_block_thickness*norm2;
        cartsTri[2] = centre + sugar_block_thickness*norm2;
        tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);

        cartsTri[0] = centre + p1 + sugar_block_thickness*norm2;
        cartsTri[1] = centre + p10 + sugar_block_thickness*norm2;
        cartsTri[2] = centre + p6 + sugar_block_thickness*norm2;
        tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);

        cartsTri[0] = centre + p2 + sugar_block_thickness*norm2;
        cartsTri[1] = centre + p6 + sugar_block_thickness*norm2;
        cartsTri[2] = centre + p7 + sugar_block_thickness*norm2;
        tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);

        cartsTri[0] = centre + p3 + sugar_block_thickness*norm2;
        cartsTri[1] = centre + p7 + sugar_block_thickness*norm2;
        cartsTri[2] = centre + p8 + sugar_block_thickness*norm2;
        tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);

        cartsTri[0] = centre + p4 + sugar_block_thickness*norm2;
        cartsTri[1] = centre + p8 + sugar_block_thickness*norm2;
        cartsTri[2] = centre + p9 + sugar_block_thickness*norm2;
        tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);

        cartsTri[0] = centre + p5 + sugar_block_thickness*norm2;
        cartsTri[1] = centre + p9 + sugar_block_thickness*norm2;
        cartsTri[2] = centre + p10 + sugar_block_thickness*norm2;
        tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);

        cartsTri[0] = centre + p6 - sugar_block_thickness*norm2;
        cartsTri[1] = centre + p7 - sugar_block_thickness*norm2;
        cartsTri[2] = centre - sugar_block_thickness*norm2;
        tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);
        
        cartsTri[0] = centre + p7 - sugar_block_thickness*norm2;
        cartsTri[1] = centre + p8 - sugar_block_thickness*norm2;
        cartsTri[2] = centre - sugar_block_thickness*norm2;
        tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);

        cartsTri[0] = centre + p8 - sugar_block_thickness*norm2;
        cartsTri[1] = centre + p9 - sugar_block_thickness*norm2;
        cartsTri[2] = centre - sugar_block_thickness*norm2;
        tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);

        cartsTri[0] = centre + p9 - sugar_block_thickness*norm2;
        cartsTri[1] = centre + p10 - sugar_block_thickness*norm2;
        cartsTri[2] = centre - sugar_block_thickness*norm2;
        tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);

        cartsTri[0] = centre + p10 - sugar_block_thickness*norm2;
        cartsTri[1] = centre + p6 - sugar_block_thickness*norm2;
        cartsTri[2] = centre - sugar_block_thickness*norm2;
        tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);

        cartsTri[1] = centre + p1 - sugar_block_thickness*norm2;
        cartsTri[0] = centre + p10 - sugar_block_thickness*norm2;
        cartsTri[2] = centre + p6 - sugar_block_thickness*norm2;
        tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);

        cartsTri[1] = centre + p2 - sugar_block_thickness*norm2;
        cartsTri[0] = centre + p6 - sugar_block_thickness*norm2;
        cartsTri[2] = centre + p7 - sugar_block_thickness*norm2;
        tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);

        cartsTri[1] = centre + p3 - sugar_block_thickness*norm2;
        cartsTri[0] = centre + p7 - sugar_block_thickness*norm2;
        cartsTri[2] = centre + p8 - sugar_block_thickness*norm2;
        tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);

        cartsTri[1] = centre + p4 - sugar_block_thickness*norm2;
        cartsTri[0] = centre + p8 - sugar_block_thickness*norm2;
        cartsTri[2] = centre + p9 - sugar_block_thickness*norm2;
        tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);

        cartsTri[1] = centre + p5 - sugar_block_thickness*norm2;
        cartsTri[0] = centre + p9 - sugar_block_thickness*norm2;
        cartsTri[2] = centre + p10 - sugar_block_thickness*norm2;
        tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);

        cartsTri[0] = centre + p1 - sugar_block_thickness*norm2;
        cartsTri[1] = centre + p1 + sugar_block_thickness*norm2;
        cartsTri[2] = centre + p6 + sugar_block_thickness*norm2;
        tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);

        cartsTri[0] = centre + p1 - sugar_block_thickness*norm2;
        cartsTri[1] = centre + p6 + sugar_block_thickness*norm2;
        cartsTri[2] = centre + p6 - sugar_block_thickness*norm2;
        tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);

        cartsTri[0] = centre + p6 - sugar_block_thickness*norm2;
        cartsTri[1] = centre + p6 + sugar_block_thickness*norm2;
        cartsTri[2] = centre + p2 + sugar_block_thickness*norm2;
        tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);

        cartsTri[0] = centre + p6 - sugar_block_thickness*norm2;
        cartsTri[1] = centre + p2 + sugar_block_thickness*norm2;
        cartsTri[2] = centre + p2 - sugar_block_thickness*norm2;
        tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);

        cartsTri[0] = centre + p2 - sugar_block_thickness*norm2;
        cartsTri[1] = centre + p2 + sugar_block_thickness*norm2;
        cartsTri[2] = centre + p7 + sugar_block_thickness*norm2;
        tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);

        cartsTri[0] = centre + p2 - sugar_block_thickness*norm2;
        cartsTri[1] = centre + p7 + sugar_block_thickness*norm2;
        cartsTri[2] = centre + p7 - sugar_block_thickness*norm2;
        tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);

        cartsTri[0] = centre + p7 - sugar_block_thickness*norm2;
        cartsTri[1] = centre + p7 + sugar_block_thickness*norm2;
        cartsTri[2] = centre + p3 + sugar_block_thickness*norm2;
        tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);

        cartsTri[0] = centre + p7 - sugar_block_thickness*norm2;
        cartsTri[1] = centre + p3 + sugar_block_thickness*norm2;
        cartsTri[2] = centre + p3 - sugar_block_thickness*norm2;
        tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);

        cartsTri[0] = centre + p3 - sugar_block_thickness*norm2;
        cartsTri[1] = centre + p3 + sugar_block_thickness*norm2;
        cartsTri[2] = centre + p8 + sugar_block_thickness*norm2;
        tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);

        cartsTri[0] = centre + p3 - sugar_block_thickness*norm2;
        cartsTri[1] = centre + p8 + sugar_block_thickness*norm2;
        cartsTri[2] = centre + p8 - sugar_block_thickness*norm2;
        tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);

        cartsTri[0] = centre + p8 - sugar_block_thickness*norm2;
        cartsTri[1] = centre + p8 + sugar_block_thickness*norm2;
        cartsTri[2] = centre + p4 + sugar_block_thickness*norm2;
        tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);

        cartsTri[0] = centre + p8 - sugar_block_thickness*norm2;
        cartsTri[1] = centre + p4 + sugar_block_thickness*norm2;
        cartsTri[2] = centre + p4 - sugar_block_thickness*norm2;
        tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);

        cartsTri[0] = centre + p4 - sugar_block_thickness*norm2;
        cartsTri[1] = centre + p4 + sugar_block_thickness*norm2;
        cartsTri[2] = centre + p9 + sugar_block_thickness*norm2;
        tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);

        cartsTri[0] = centre + p4 - sugar_block_thickness*norm2;
        cartsTri[1] = centre + p9 + sugar_block_thickness*norm2;
        cartsTri[2] = centre + p9 - sugar_block_thickness*norm2;
        tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);

        cartsTri[0] = centre + p9 - sugar_block_thickness*norm2;
        cartsTri[1] = centre + p9 + sugar_block_thickness*norm2;
        cartsTri[2] = centre + p5 + sugar_block_thickness*norm2;
        tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);

        cartsTri[0] = centre + p9 - sugar_block_thickness*norm2;
        cartsTri[1] = centre + p5 + sugar_block_thickness*norm2;
        cartsTri[2] = centre + p5 - sugar_block_thickness*norm2;
        tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);

        cartsTri[0] = centre + p5 - sugar_block_thickness*norm2;
        cartsTri[1] = centre + p5 + sugar_block_thickness*norm2;
        cartsTri[2] = centre + p10 + sugar_block_thickness*norm2;
        tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);

        cartsTri[0] = centre + p5 - sugar_block_thickness*norm2;
        cartsTri[1] = centre + p10 + sugar_block_thickness*norm2;
        cartsTri[2] = centre + p10 - sugar_block_thickness*norm2;
        tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);

        cartsTri[0] = centre + p10 - sugar_block_thickness*norm2;
        cartsTri[1] = centre + p10 + sugar_block_thickness*norm2;
        cartsTri[2] = centre + p1 + sugar_block_thickness*norm2;
        tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);

        cartsTri[0] = centre + p10 - sugar_block_thickness*norm2;
        cartsTri[1] = centre + p1 + sugar_block_thickness*norm2;
        cartsTri[2] = centre + p1 - sugar_block_thickness*norm2;
        tri = new TriangleElement(cartsTri,col1,centre,1.0);
        tris->add_primitive(tri);

      }
    }
  }
  return retval;
}

std::vector<Cartesian> DrawSugarBlock(TriangleCollection *tris, CylinderCollection *cyls, mmdb::Residue* res1, const double *col1, const CParamsManager &params, const CParamsManager &global_params, int selHnd, mmdb::Manager *molHnd, bool two_colour, float sugar_block_thickness, float sugar_block_scale){
  //std::cout << "DrawSugarBlock " << res1->GetResName() << "\n";

  // FIXME - We should almost certainly check that residue names match as well. This requires user being able
  // to provide custom saccharide residues....
  // We also need to be checking so that we know what shape to draw. Maybe Jon's code would be helpful here.

  //Complicated ?
  // RIB     C1',...,C4',O4' (5-membered)
  // F6P     C2,...,C5,O5 (5-membered)
  // GLO     O1,C1 ... C5(O5),C6,O6 (linear! Ignored?)
  // SUC     C1,...,C5,O5;C1',...,C4',O4' (6,5-membered)
  // DMU     C1,...,C5,O5;    C10,C5,C7,C8,C9,01 (6,6-membered)

  //Simpler
  // NAG     C1,...,C5,O5; (6-membered)
  // FFC     C1A,...,C5A,O5A; C1B,...,C5B,O5B; (6,6-membered)
  // LAT,LBT C1,...,C5,O5;    C1',...,C5',O5' (6,6-membered)
  // LMT     C1',...,C5',O5'; C1B,...,C5B,O5B; (6,6-membered)
  // TRE     C1P,...,C5P,O5P; C1B,...,C5B,O5B; (6,6-membered)

  std::vector<std::vector<std::vector<std::string> > > sugarAtomNames;
  std::vector<std::vector<std::string> > ribAtomNames;
  std::vector<std::vector<std::string> > nagAtomNames;
  std::vector<std::vector<std::string> > ffcAtomNames;
  std::vector<std::vector<std::string> > latAtomNames;
  std::vector<std::vector<std::string> > lmtAtomNames;
  std::vector<std::vector<std::string> > treAtomNames;
  std::vector<std::vector<std::string> > dmuAtomNames;
  std::vector<std::string> ribR1AtomNames;
  std::vector<std::string> nagR1AtomNames;
  std::vector<std::string> ffcR1AtomNames;
  std::vector<std::string> latR1AtomNames;
  std::vector<std::string> lmtR1AtomNames;
  std::vector<std::string> treR1AtomNames;
  std::vector<std::string> dmuR1AtomNames;
  std::vector<std::string> ffcR2AtomNames;
  std::vector<std::string> latR2AtomNames;
  std::vector<std::string> lmtR2AtomNames;
  std::vector<std::string> treR2AtomNames;
  std::vector<std::string> dmuR2AtomNames;

  ffcR1AtomNames.push_back(std::string("C1A"));
  ffcR1AtomNames.push_back(std::string("C2A"));
  ffcR1AtomNames.push_back(std::string("C3A"));
  ffcR1AtomNames.push_back(std::string("C4A"));
  ffcR1AtomNames.push_back(std::string("C5A"));
  ffcR1AtomNames.push_back(std::string("O5A"));
  ffcR2AtomNames.push_back(std::string("C1B"));
  ffcR2AtomNames.push_back(std::string("C2B"));
  ffcR2AtomNames.push_back(std::string("C3B"));
  ffcR2AtomNames.push_back(std::string("C4B"));
  ffcR2AtomNames.push_back(std::string("C5B"));
  ffcR2AtomNames.push_back(std::string("O5B"));
  ffcAtomNames.push_back(ffcR1AtomNames);
  ffcAtomNames.push_back(ffcR2AtomNames);
  sugarAtomNames.push_back(ffcAtomNames);

  latR1AtomNames.push_back(std::string("C1"));
  latR1AtomNames.push_back(std::string("C2"));
  latR1AtomNames.push_back(std::string("C3"));
  latR1AtomNames.push_back(std::string("C4"));
  latR1AtomNames.push_back(std::string("C5"));
  latR1AtomNames.push_back(std::string("O5"));
  latR2AtomNames.push_back(std::string("C1'"));
  latR2AtomNames.push_back(std::string("C2'"));
  latR2AtomNames.push_back(std::string("C3'"));
  latR2AtomNames.push_back(std::string("C4'"));
  latR2AtomNames.push_back(std::string("C5'"));
  latR2AtomNames.push_back(std::string("O5'"));
  latAtomNames.push_back(latR1AtomNames);
  latAtomNames.push_back(latR2AtomNames);
  sugarAtomNames.push_back(latAtomNames);

  lmtR1AtomNames.push_back(std::string("C1'"));
  lmtR1AtomNames.push_back(std::string("C2'"));
  lmtR1AtomNames.push_back(std::string("C3'"));
  lmtR1AtomNames.push_back(std::string("C4'"));
  lmtR1AtomNames.push_back(std::string("C5'"));
  lmtR1AtomNames.push_back(std::string("O5'"));
  lmtR2AtomNames.push_back(std::string("C1B"));
  lmtR2AtomNames.push_back(std::string("C2B"));
  lmtR2AtomNames.push_back(std::string("C3B"));
  lmtR2AtomNames.push_back(std::string("C4B"));
  lmtR2AtomNames.push_back(std::string("C5B"));
  lmtR2AtomNames.push_back(std::string("O5B"));
  lmtAtomNames.push_back(lmtR1AtomNames);
  lmtAtomNames.push_back(lmtR2AtomNames);
  sugarAtomNames.push_back(lmtAtomNames);

  treR1AtomNames.push_back(std::string("C1P"));
  treR1AtomNames.push_back(std::string("C2P"));
  treR1AtomNames.push_back(std::string("C3P"));
  treR1AtomNames.push_back(std::string("C4P"));
  treR1AtomNames.push_back(std::string("C5P"));
  treR1AtomNames.push_back(std::string("O5P"));
  treR2AtomNames.push_back(std::string("C1B"));
  treR2AtomNames.push_back(std::string("C2B"));
  treR2AtomNames.push_back(std::string("C3B"));
  treR2AtomNames.push_back(std::string("C4B"));
  treR2AtomNames.push_back(std::string("C5B"));
  treR2AtomNames.push_back(std::string("O5B"));
  treAtomNames.push_back(treR1AtomNames);
  treAtomNames.push_back(treR2AtomNames);
  sugarAtomNames.push_back(treAtomNames);

  ribR1AtomNames.push_back(std::string("C1'"));
  ribR1AtomNames.push_back(std::string("C2'"));
  ribR1AtomNames.push_back(std::string("C3'"));
  ribR1AtomNames.push_back(std::string("C4'"));
  ribR1AtomNames.push_back(std::string("O4'"));
  ribAtomNames.push_back(ribR1AtomNames);
  sugarAtomNames.push_back(ribAtomNames);

  nagR1AtomNames.push_back(std::string("C1"));
  nagR1AtomNames.push_back(std::string("C2"));
  nagR1AtomNames.push_back(std::string("C3"));
  nagR1AtomNames.push_back(std::string("C4"));
  nagR1AtomNames.push_back(std::string("C5"));
  nagR1AtomNames.push_back(std::string("O5"));
  nagAtomNames.push_back(nagR1AtomNames);
  sugarAtomNames.push_back(nagAtomNames);

  for(unsigned isugartype=0;isugartype<sugarAtomNames.size();isugartype++){
   for(unsigned iring=0;iring<sugarAtomNames[isugartype].size();iring++){
    if(sugarAtomNames[isugartype][iring].size()>4){
    const char* c1name = sugarAtomNames[isugartype][iring][0].c_str();
    const char* c2name = sugarAtomNames[isugartype][iring][1].c_str();
    const char* c3name = sugarAtomNames[isugartype][iring][2].c_str();
    const char* c4name = sugarAtomNames[isugartype][iring][3].c_str();
    const char* c5name = sugarAtomNames[isugartype][iring][4].c_str();
    const char* o5name = NULL;
    mmdb::Atom* o5 = NULL;
    if(sugarAtomNames[isugartype][iring].size()==6){
      o5name = sugarAtomNames[isugartype][iring][5].c_str();
      o5 = res1->GetAtom(o5name);
      //std::cout << "Looking for " << c1name << " " << c2name << " " << c3name << " " << c4name << " " << c5name << " " << o5name << "\n";
    } else {
      //std::cout << "Looking for " << c1name << " " << c2name << " " << c3name << " " << c4name << " " << c5name << "\n";
    }
    mmdb::Atom* c1 = res1->GetAtom(c1name);
    mmdb::Atom* c2 = res1->GetAtom(c2name);
    mmdb::Atom* c3 = res1->GetAtom(c3name);
    mmdb::Atom* c4 = res1->GetAtom(c4name);
    mmdb::Atom* c5 = res1->GetAtom(c5name);
    if(c1&&c2&&c3&&c4&&c5){
      //std::cout << "Have them all\n";
      std::vector<Cartesian> res = DrawSugarBlockInt(tris, cyls, res1, col1, params, global_params, selHnd, c1,c2,c3,c4,c5,o5,res1, molHnd, two_colour, sugar_block_thickness, sugar_block_scale );
      if(iring==sugarAtomNames[isugartype].size()-1&&res.size()==2) return res;
    }

    c1 = res1->GetAtom(c1name,NULL,"A");
    c2 = res1->GetAtom(c2name,NULL,"A");
    c3 = res1->GetAtom(c3name,NULL,"A");
    c4 = res1->GetAtom(c4name,NULL,"A");
    c5 = res1->GetAtom(c5name,NULL,"A");
    if(sugarAtomNames[isugartype][iring].size()==6){
      o5 = res1->GetAtom(o5name,NULL,"A");
    }
    if(c1&&c2&&c3&&c4&&c5){
      if(c1->isInSelection(selHnd)&&c2->isInSelection(selHnd)&&c3->isInSelection(selHnd)&&c4->isInSelection(selHnd)&&c5->isInSelection(selHnd)&&o5->isInSelection(selHnd)){
        std::vector<Cartesian> res = DrawSugarBlockInt(tris, cyls, res1, col1, params, global_params, selHnd, c1,c2,c3,c4,c5,o5,res1, molHnd, two_colour, sugar_block_thickness, sugar_block_scale );
        if(iring==sugarAtomNames[isugartype].size()-1&&res.size()==2) return res;
      }
    }
    c1 = res1->GetAtom(c1name,NULL,"B");
    c2 = res1->GetAtom(c2name,NULL,"B");
    c3 = res1->GetAtom(c3name,NULL,"B");
    c4 = res1->GetAtom(c4name,NULL,"B");
    c5 = res1->GetAtom(c5name,NULL,"B");
    if(sugarAtomNames[isugartype][iring].size()==6){
      o5 = res1->GetAtom(o5name,NULL,"B");
    }
    if(c1&&c2&&c3&&c4&&c5){
      if(c1->isInSelection(selHnd)&&c2->isInSelection(selHnd)&&c3->isInSelection(selHnd)&&c4->isInSelection(selHnd)&&c5->isInSelection(selHnd)&&o5->isInSelection(selHnd)){
        std::vector<Cartesian> res = DrawSugarBlockInt(tris, cyls, res1, col1, params, global_params, selHnd, c1,c2,c3,c4,c5,o5,res1, molHnd, two_colour, sugar_block_thickness, sugar_block_scale );
        if(iring==sugarAtomNames[isugartype].size()-1&&res.size()==2) return res;
      }
    }
    c1 = res1->GetAtom(c1name,NULL,"C");
    c2 = res1->GetAtom(c2name,NULL,"C");
    c3 = res1->GetAtom(c3name,NULL,"C");
    c4 = res1->GetAtom(c4name,NULL,"C");
    c5 = res1->GetAtom(c5name,NULL,"C");
    if(sugarAtomNames[isugartype][iring].size()==6){
      o5 = res1->GetAtom(o5name,NULL,"C");
    }
    if(c1&&c2&&c3&&c4&&c5){
      if(c1->isInSelection(selHnd)&&c2->isInSelection(selHnd)&&c3->isInSelection(selHnd)&&c4->isInSelection(selHnd)&&c5->isInSelection(selHnd)&&o5->isInSelection(selHnd)){
        std::vector<Cartesian> res = DrawSugarBlockInt(tris, cyls, res1, col1, params, global_params, selHnd, c1,c2,c3,c4,c5,o5,res1, molHnd, two_colour, sugar_block_thickness, sugar_block_scale );
        if(iring==sugarAtomNames[isugartype].size()-1&&res.size()==2) return res;
      }
     }
    }
   }
  }
  std::vector<Cartesian> retval;
  return retval;
}

void DrawSugarBlocks(Displayobject &obj, CMMANManager *molHnd, int selHnd, mmdb::Atom** selAtoms, int nSelAtoms, const AtomColourVector &atom_colour_vector,const CParamsManager &params, const CParamsManager &global_params, double *hb_params_array){

  int two_colour = params.GetInt("two_glycoblock_colour"); 
  std::cout << "DrawSugarBlocks, two colour: " << two_colour << "\n";

  double stick_col[4];
  bool have_stick_col = false;
  if(params.GetString("glycoblock_stick_colour")!=std::string("")&&params.GetString("glycoblock_stick_colour")!=std::string("default")){
    std::vector<double> stick_col_v = RGBReps::GetColour(params.GetString("glycoblock_stick_colour"));
    stick_col[0] = stick_col_v[0];
    stick_col[1] = stick_col_v[1];
    stick_col[2] = stick_col_v[2];
    stick_col[3] = stick_col_v[3];
    have_stick_col = true;
  }

  int first_nmr_model=-1;
  int n_nmr_models=molHnd->GetNumberOfModels();
  if(n_nmr_models>0){
    for(int i=1;i<=n_nmr_models;i++){
      if(molHnd->GetModel(i)){
        first_nmr_model = i;
        break;
      }
    }
  } else {
    // Fallback.
    first_nmr_model=1;
  }
  const char *amino_acid = "GLY,ALA,VAL,PRO,SER,THR,LEU,ILE,CYS,ASP,GLU,ASN,GLN,ARG,LYS,MET,MSE,HIS,PHE,TYR,TRP,HCS,ALO,PDD";
  int proteinSelHnd = molHnd->NewSelection();
  molHnd->SelectAminoNotHet(proteinSelHnd,mmdb::STYPE_ATOM,0,"*",mmdb::ANY_RES,"*",mmdb::ANY_RES,"*",amino_acid,"*","*","*",mmdb::SKEY_NEW);
  int nProteinAtoms;
  mmdb::Atom** ProteinAtoms=0;
  molHnd->GetSelIndex(proteinSelHnd,ProteinAtoms,nProteinAtoms);
  std::cout << "first_nmr_model " << first_nmr_model << std::endl;
  std::cout << "Selected " << nProteinAtoms << " atoms" << std::endl;

  bool drawProteinInteractions = false;
  bool drawCovalentProteinInteractions = false;
  bool labelProteinInteractions = false;
  double interaction_cylinder_radius = 0.1;
  int interaction_line_width = 1;
  int interaction_style = DASHCYLINDER;

  double cylinders_accu = 20;
  int quality = global_params.GetInt("solid_quality");
  if(quality==0){
    cylinders_accu = 20;
  } else if(quality==1){
    cylinders_accu = 40;
  } else if(quality==2){
    cylinders_accu = 60;
  }

  double interaction_colour[4] = {0.0,0.0,0.0,1.0};
  float interaction_dash_length = 0.36;
  float sugar_block_thickness = 0.1;
  float sugar_block_scale = 1.0;

  if(params.GetString("glycoblock_interaction_colour")!=std::string("")&&params.GetString("glycoblock_interaction_colour")!=std::string("default")){
    std::vector<double> stick_col_v = RGBReps::GetColour(params.GetString("glycoblock_interaction_colour"));
    interaction_colour[0] = stick_col_v[0];
    interaction_colour[1] = stick_col_v[1];
    interaction_colour[2] = stick_col_v[2];
    interaction_colour[3] = stick_col_v[3];
  }
  sugar_block_thickness = params.GetFloat("glycoblock_thickness");
  sugar_block_scale = params.GetFloat("glycoblock_scale");
  interaction_cylinder_radius = params.GetFloat("glycoblock_interaction_cylinder_radius");
  interaction_line_width = params.GetInt("glycoblock_interaction_line_width");
  std::string interaction_style_pref = params.GetString("glycoblock_interaction_style");
  if(interaction_style_pref=="Dashed line"){
    interaction_style = DASHLINE;
  } else {
    interaction_style = DASHCYLINDER;
  }
  interaction_dash_length = params.GetFloat("glycoblock_interaction_dash_length");
  std::cout << "glycoblock_draw_protein_interactions " << params.GetInt("glycoblock_draw_protein_interactions") << std::endl;
  std::cout << "glycoblock_draw_cov_protein_interactions " << params.GetInt("glycoblock_draw_cov_protein_interactions") << std::endl;
  drawProteinInteractions = params.GetInt("glycoblock_draw_protein_interactions");
  drawCovalentProteinInteractions = params.GetInt("glycoblock_draw_cov_protein_interactions");
  labelProteinInteractions = params.GetInt("glycoblock_label_protein_interactions");

  //FIXME - Not cylinders in general!
  SphereCollection *spheres = new SphereCollection();
  CylinderCollection *cyls = new CylinderCollection();
  TriangleCollection *tris = new TriangleCollection();
  DashCylinderCollection *dash_cylinders = new DashCylinderCollection();
  DashLinesCollection *dash_lines_collection = new DashLinesCollection();
  char ResidueID1[30];
  char ResidueID2[30];
  int C1sel = molHnd->NewSelection();
  // FFC     C1A,...,C5A,O5A;C1B,...,C5B,O5B; (6,6-membered)
  // RIB     C1',...,C4',O4' (5-membered)
  // F6P     C2,...,C5,O5 (5-membered)
  // SUC     C1,...,C5,O5;C1',...,C4',O4' (6,5-membered)
  // LAT,LBT C1,...,C5,O5;C1',...,C5',O5' (6,6-membered)
  // GLO     O1,C1 ... C5(O5),C6,O6 (linear!)
  // DMU     C1,...,C5,O5;C10,C5,C7,C8,C9,01 (6,6-membered)
  // LMT     C1',...,C5',O5';C1B,...,C5B,O5B; (6,6-membered)
  // TRE     C1P,...,C5P,O5P;C1B,...,C5B,O5B; (6,6-membered)

  int udd_C1X = molHnd->GetUDDHandle ( UDR_RESIDUE,"GLYCO_BLOCK_C1X" );
  if (udd_C1X<=0) udd_C1X = molHnd->RegisterUDReal(UDR_RESIDUE,"GLYCO_BLOCK_C1X");

  int udd_C1Y = molHnd->GetUDDHandle ( UDR_RESIDUE,"GLYCO_BLOCK_C1Y" );
  if (udd_C1Y<=0) udd_C1Y = molHnd->RegisterUDReal(UDR_RESIDUE,"GLYCO_BLOCK_C1Y");

  int udd_C1Z = molHnd->GetUDDHandle ( UDR_RESIDUE,"GLYCO_BLOCK_C1Z" );
  if (udd_C1Z<=0) udd_C1Z = molHnd->RegisterUDReal(UDR_RESIDUE,"GLYCO_BLOCK_C1Z");

  int udd_C3X = molHnd->GetUDDHandle ( UDR_RESIDUE,"GLYCO_BLOCK_C3X" );
  if (udd_C3X<=0) udd_C3X = molHnd->RegisterUDReal(UDR_RESIDUE,"GLYCO_BLOCK_C3X");

  int udd_C3Y = molHnd->GetUDDHandle ( UDR_RESIDUE,"GLYCO_BLOCK_C3Y" );
  if (udd_C3Y<=0) udd_C3Y = molHnd->RegisterUDReal(UDR_RESIDUE,"GLYCO_BLOCK_C3Y");

  int udd_C3Z = molHnd->GetUDDHandle ( UDR_RESIDUE,"GLYCO_BLOCK_C3Z" );
  if (udd_C3Z<=0) udd_C3Z = molHnd->RegisterUDReal(UDR_RESIDUE,"GLYCO_BLOCK_C3Z");

  int udd_C4X = molHnd->GetUDDHandle ( UDR_RESIDUE,"GLYCO_BLOCK_C4X" );
  if (udd_C4X<=0) udd_C4X = molHnd->RegisterUDReal(UDR_RESIDUE,"GLYCO_BLOCK_C4X");

  int udd_C4Y = molHnd->GetUDDHandle ( UDR_RESIDUE,"GLYCO_BLOCK_C4Y" );
  if (udd_C4Y<=0) udd_C4Y = molHnd->RegisterUDReal(UDR_RESIDUE,"GLYCO_BLOCK_C4Y");

  int udd_C4Z = molHnd->GetUDDHandle ( UDR_RESIDUE,"GLYCO_BLOCK_C4Z" );
  if (udd_C4Z<=0) udd_C4Z = molHnd->RegisterUDReal(UDR_RESIDUE,"GLYCO_BLOCK_C4Z");

  int udd_C5X = molHnd->GetUDDHandle ( UDR_RESIDUE,"GLYCO_BLOCK_C5X" );
  if (udd_C5X<=0) udd_C5X = molHnd->RegisterUDReal(UDR_RESIDUE,"GLYCO_BLOCK_C5X");

  int udd_C5Y = molHnd->GetUDDHandle ( UDR_RESIDUE,"GLYCO_BLOCK_C5Y" );
  if (udd_C5Y<=0) udd_C5Y = molHnd->RegisterUDReal(UDR_RESIDUE,"GLYCO_BLOCK_C5Y");

  int udd_C5Z = molHnd->GetUDDHandle ( UDR_RESIDUE,"GLYCO_BLOCK_C5Z" );
  if (udd_C5Z<=0) udd_C5Z = molHnd->RegisterUDReal(UDR_RESIDUE,"GLYCO_BLOCK_C5Z");

  int udd_C2X = molHnd->GetUDDHandle ( UDR_RESIDUE,"GLYCO_BLOCK_C2X" );
  if (udd_C2X<=0) udd_C2X = molHnd->RegisterUDReal(UDR_RESIDUE,"GLYCO_BLOCK_C2X");

  int udd_C2Y = molHnd->GetUDDHandle ( UDR_RESIDUE,"GLYCO_BLOCK_C2Y" );
  if (udd_C2Y<=0) udd_C2Y = molHnd->RegisterUDReal(UDR_RESIDUE,"GLYCO_BLOCK_C2Y");

  int udd_C2Z = molHnd->GetUDDHandle ( UDR_RESIDUE,"GLYCO_BLOCK_C2Z" );
  if (udd_C2Z<=0) udd_C2Z = molHnd->RegisterUDReal(UDR_RESIDUE,"GLYCO_BLOCK_C2Z");

  molHnd->Select(C1sel,STYPE_ATOM,0,"*",ANY_RES,"*",ANY_RES,"*","*","C1","C","*",SKEY_NEW);
  molHnd->Select(C1sel,STYPE_ATOM,0,"*",ANY_RES,"*",ANY_RES,"*","*","C2","C","*",SKEY_NEW);
  molHnd->Select(C1sel,STYPE_ATOM,0,"*",ANY_RES,"*",ANY_RES,"*","*","C1A","C","*",SKEY_OR);
  molHnd->Select(C1sel,STYPE_ATOM,0,"*",ANY_RES,"*",ANY_RES,"*","*","C1P","C","*",SKEY_OR);
  molHnd->Select(C1sel,STYPE_ATOM,0,"*",ANY_RES,"*",ANY_RES,"*","*","C1'","C","*",SKEY_OR);
  AtomBond *AtomBond;
  int nAtomBonds;

  CHBond proteinSugarHBonds(molHnd,selHnd,molHnd,proteinSelHnd);
  CHBond sugarSugarHBonds(molHnd,selHnd,molHnd,selHnd);

  if(hb_params_array) proteinSugarHBonds.SetParams(8,hb_params_array);
  proteinSugarHBonds.Calculate();

  if(hb_params_array) sugarSugarHBonds.SetParams(8,hb_params_array);
  sugarSugarHBonds.Calculate();

  //std::cout << proteinSugarHBonds.Print() << std::endl;

  std::map<std::string,std::vector<Cartesian> > idPlaneMap;
  for(int ii=0;ii<nSelAtoms;ii++){
    if(selAtoms[ii]->isInSelection(C1sel)){
      mmdb::Residue* res = selAtoms[ii]->GetResidue();
      if(res){
        int restype = molHnd->GetRestypeCode(res);
        //std::cout << res->name << " " << restype << " " << res->seqNum << std::endl;
        if(restype==RESTYPE_NONPOLY||restype==RESTYPE_SACH||restype==RESTYPE_DSACH||restype==RESTYPE_LSACH){
           const double *col = atom_colour_vector.GetRGB(ii);
           std::vector<Cartesian> conn_carts = DrawSugarBlock(tris, cyls,res,col, params, global_params, selHnd, molHnd, two_colour,sugar_block_thickness, sugar_block_scale);
           if(drawProteinInteractions&&conn_carts.size()>1){
             //std::cout << "drawProteinInteractions" << std::endl;
             std::vector<mmdb::Atom*>::iterator hbondIter = proteinSugarHBonds.hbonds.pAtom2.begin();
             std::vector<mmdb::Atom*>::iterator hbondIter1 = proteinSugarHBonds.hbonds.pAtom1.begin();
             while(hbondIter!=proteinSugarHBonds.hbonds.pAtom2.end()&&hbondIter1!=proteinSugarHBonds.hbonds.pAtom1.end()){
               if((*hbondIter)->GetResidue()&&(*hbondIter1)->GetResidue()&&(*hbondIter1)->GetResidue()==res){
                 mmdb::Atom* proteinCA = (*hbondIter)->GetResidue()->GetAtom("CA");
                 if(proteinCA){
                    Cartesian CAcart(proteinCA->x,proteinCA->y,proteinCA->z);
                    std::vector<Cartesian> interCarts;
                    interCarts.push_back(conn_carts[0]);
                    interCarts.push_back(CAcart);
                    if(interaction_style==DASHLINE){
                      dash_lines_collection->add(interCarts,interaction_colour,0);
                    } else {
                      DashCylinderElement *line = new DashCylinderElement(interCarts,interaction_colour,interCarts[0],interaction_cylinder_radius);
                      dash_cylinders->add_primitive(line);
                    }
                    if(labelProteinInteractions){
                      double length = (Cartesian((*hbondIter)->x,(*hbondIter)->y,(*hbondIter)->z)-Cartesian((*hbondIter1)->x,(*hbondIter1)->y,(*hbondIter1)->z)).length();
                      std::string length_text = FloatToString(length,"%.3f");
                      Text *text = new Text (Cartesian::MidPoint(interCarts), length_text,Cartesian::MidPoint(interCarts));
                      text->SetDefaultColour();
                      obj.add_text_primitive(text);
                    }
                 }
               }
               hbondIter++;
               hbondIter1++;
             }
             hbondIter = sugarSugarHBonds.hbonds.pAtom2.begin();
             hbondIter1 = sugarSugarHBonds.hbonds.pAtom1.begin();
             //std::cout << "Sugar-sugar hbonds: " << sugarSugarHBonds.hbonds.pAtom1.size() << " " << sugarSugarHBonds.hbonds.pAtom2.size() << std::endl;

             while(hbondIter!=sugarSugarHBonds.hbonds.pAtom2.end()&&hbondIter1!=proteinSugarHBonds.hbonds.pAtom1.end()){

               mmdb::Residue* res1 = (*hbondIter)->GetResidue();
               mmdb::Residue* res2 = (*hbondIter1)->GetResidue();

               bool ignoreHBond = false;
               for(int iat1=0;iat1<res1->GetNumberOfAtoms()&&(!ignoreHBond);iat1++){
                 res1->GetAtom(iat1)->GetBonds(AtomBond,nAtomBonds);
                 for(int iatB=0;iatB<nAtomBonds;iatB++){
                   if(AtomBond[iatB].atom->residue==res2){
                      ignoreHBond = true;
                      break;
                   }
                 }
               }

               if((!ignoreHBond)&&(*hbondIter)->GetResidue()&&(*hbondIter1)->GetResidue()&&(*hbondIter1)->GetResidue()==res){
                 mmdb::Atom* proteinC1 = (*hbondIter)->GetResidue()->GetAtom("C1");
                 mmdb::Atom* proteinC2 = (*hbondIter)->GetResidue()->GetAtom("C2");
                 mmdb::Atom* proteinC3 = (*hbondIter)->GetResidue()->GetAtom("C3");
                 mmdb::Atom* proteinC4 = (*hbondIter)->GetResidue()->GetAtom("C4");
                 mmdb::Atom* proteinC5 = (*hbondIter)->GetResidue()->GetAtom("C5");
                 mmdb::Atom* proteinO5 = (*hbondIter)->GetResidue()->GetAtom("O5");
                 mmdb::Atom* proteinO6 = (*hbondIter)->GetResidue()->GetAtom("O6");
                 std::vector<Cartesian> centreCarts;
                 if(proteinC1){
                   centreCarts.push_back(Cartesian(proteinC1->x,proteinC1->y,proteinC1->z));
                 }
                 if(proteinC2){
                   centreCarts.push_back(Cartesian(proteinC2->x,proteinC2->y,proteinC2->z));
                 }
                 if(proteinC3){
                   centreCarts.push_back(Cartesian(proteinC3->x,proteinC3->y,proteinC3->z));
                 }
                 if(proteinC4){
                   centreCarts.push_back(Cartesian(proteinC4->x,proteinC4->y,proteinC4->z));
                 }
                 if(proteinC5){
                   centreCarts.push_back(Cartesian(proteinC5->x,proteinC5->y,proteinC5->z));
                 }
                 if(proteinO5){
                   centreCarts.push_back(Cartesian(proteinO5->x,proteinO5->y,proteinO5->z));
                 }
                 if(proteinO6){
                   centreCarts.push_back(Cartesian(proteinO6->x,proteinO6->y,proteinO6->z));
                 }
                 if(centreCarts.size()>2){
                    std::vector<Cartesian> interCarts;
                    interCarts.push_back(conn_carts[0]);
                    interCarts.push_back(Cartesian::MidPoint(centreCarts));
                    if(interaction_style==DASHLINE){
                      dash_lines_collection->add(interCarts,interaction_colour,0);
                    } else {
                      DashCylinderElement *line = new DashCylinderElement(interCarts,interaction_colour,interCarts[0],interaction_cylinder_radius);
                      dash_cylinders->add_primitive(line);
                    }
                    if(labelProteinInteractions){
                      double length = (Cartesian((*hbondIter)->x,(*hbondIter)->y,(*hbondIter)->z)-Cartesian((*hbondIter1)->x,(*hbondIter1)->y,(*hbondIter1)->z)).length();
                      std::string length_text = FloatToString(length,"%.3f");
                      Text *text = new Text (Cartesian::MidPoint(interCarts), length_text,Cartesian::MidPoint(interCarts));
                      text->SetDefaultColour();
                      obj.add_text_primitive(text);
                    }
                 }
               }
               hbondIter++;
               hbondIter1++;
             }
           }
           delete [] col;
        }
      }
    }
  }
  for(int ii=0;ii<nSelAtoms;ii++){
    if(selAtoms[ii]->isInSelection(C1sel)){
      mmdb::Residue* res = selAtoms[ii]->GetResidue();
      if(res){
        int restype = molHnd->GetRestypeCode(res);
        //std::cout << res->name << " " << restype << " " << res->seqNum << std::endl;
        std::vector<Cartesian> conn_carts(2);
        for(int iat=0;iat<res->nAtoms;iat++){
          if(res->atom[iat]){
            if(strncmp(res->atom[iat]->name," C1 ",4)==0){
              // Not atom coords! UDDData for C1.
              mmdb::realtype UDDX, UDDY, UDDZ;
              res->GetUDData(udd_C1X,UDDX); res->GetUDData(udd_C1Y,UDDY); res->GetUDData(udd_C1Z,UDDZ);
              conn_carts[0] = Cartesian(UDDX,UDDY,UDDZ);
              res->atom[iat]->GetBonds(AtomBond,nAtomBonds);
              for (int i=0;i<nAtomBonds;i++){
                if(AtomBond[i].atom->residue&&AtomBond[i].atom->residue->seqNum!=res->seqNum){
                  int restype2 = molHnd->GetRestypeCode(AtomBond[i].atom->residue);
                  if(drawCovalentProteinInteractions&&(restype2==RESTYPE_PEPTIDE||restype2==RESTYPE_DPEPTIDE||restype2==RESTYPE_LPEPTIDE)){
                    if(AtomBond[i].atom->GetResidue()&&AtomBond[i].atom->GetResidue()->GetAtom("CA")){
                      mmdb::Atom* proteinCA = AtomBond[i].atom->GetResidue()->GetAtom("CA");
                      Cartesian CAcart(proteinCA->x,proteinCA->y,proteinCA->z);
                      std::vector<Cartesian> interCarts;
                      interCarts.push_back(conn_carts[0]);
                      interCarts.push_back(CAcart);
                      CylinderElement *line = new CylinderElement(interCarts,interaction_colour,interCarts[0],interaction_cylinder_radius);
                      cyls->add_primitive(line);
                      SphereElement *sphere = new SphereElement(conn_carts[0],interaction_colour,conn_carts[0],interaction_cylinder_radius,1.0,1);
                      spheres->add_primitive(sphere);
                    }
                  }
                  if(restype2==RESTYPE_SACH||restype2==RESTYPE_DSACH||restype2==RESTYPE_LSACH){
                    const double *col = atom_colour_vector.GetRGB(ii);
                    if(strncmp(AtomBond[i].atom->name," O4 ",4)==0){
                      AtomBond[i].atom->GetResidue()->GetUDData(udd_C4X,UDDX);
                      AtomBond[i].atom->GetResidue()->GetUDData(udd_C4Y,UDDY);
                      AtomBond[i].atom->GetResidue()->GetUDData(udd_C4Z,UDDZ);
                      conn_carts[1] = Cartesian(UDDX,UDDY,UDDZ);
                    } else if(strncmp(AtomBond[i].atom->name," O3 ",4)==0){
                      AtomBond[i].atom->GetResidue()->GetUDData(udd_C3X,UDDX);
                      AtomBond[i].atom->GetResidue()->GetUDData(udd_C3Y,UDDY);
                      AtomBond[i].atom->GetResidue()->GetUDData(udd_C3Z,UDDZ);
                      conn_carts[1] = Cartesian(UDDX,UDDY,UDDZ);
                    } else if(strncmp(AtomBond[i].atom->name," O6 ",4)==0){
                      AtomBond[i].atom->GetResidue()->GetUDData(udd_C5X,UDDX);
                      AtomBond[i].atom->GetResidue()->GetUDData(udd_C5Y,UDDY);
                      AtomBond[i].atom->GetResidue()->GetUDData(udd_C5Z,UDDZ);
                      conn_carts[1] = Cartesian(UDDX,UDDY,UDDZ);
                    } else if(strncmp(AtomBond[i].atom->name," O2 ",4)==0){
                      AtomBond[i].atom->GetResidue()->GetUDData(udd_C2X,UDDX);
                      AtomBond[i].atom->GetResidue()->GetUDData(udd_C2Y,UDDY);
                      AtomBond[i].atom->GetResidue()->GetUDData(udd_C2Z,UDDZ);
                      conn_carts[1] = Cartesian(UDDX,UDDY,UDDZ);
                    }
                    conn_carts[1] = Cartesian::MidPoint(conn_carts);
                    CylinderElement *cylinder;
                    SphereElement *sphere;
                    if(have_stick_col){
                      cylinder = new CylinderElement(conn_carts,stick_col,conn_carts[1],interaction_cylinder_radius,1.0);
                      sphere = new SphereElement(conn_carts[0],stick_col,conn_carts[0],interaction_cylinder_radius,1.0,1);
                    } else {
                      cylinder = new CylinderElement(conn_carts,col,conn_carts[1],interaction_cylinder_radius,1.0);
                      sphere = new SphereElement(conn_carts[0],col,conn_carts[0],interaction_cylinder_radius,1.0,1);
                    }
                    cyls->add_primitive(cylinder);
                    spheres->add_primitive(sphere);
                  }
                } 
              } 
            }
            if(strncmp(res->atom[iat]->element," O",2)==0){
              res->atom[iat]->GetBonds(AtomBond,nAtomBonds);
              mmdb::realtype UDDX, UDDY, UDDZ;
              if(strncmp(res->atom[iat]->name," O4 ",4)==0){
                res->GetUDData(udd_C4X,UDDX); res->GetUDData(udd_C4Y,UDDY); res->GetUDData(udd_C4Z,UDDZ);
                conn_carts[0] = Cartesian(UDDX,UDDY,UDDZ);
              } else if(strncmp(res->atom[iat]->name," O3 ",4)==0) {
                res->GetUDData(udd_C3X,UDDX); res->GetUDData(udd_C3Y,UDDY); res->GetUDData(udd_C3Z,UDDZ);
                conn_carts[0] = Cartesian(UDDX,UDDY,UDDZ);
              } else if(strncmp(res->atom[iat]->name," O6 ",4)==0) {
                res->GetUDData(udd_C5X,UDDX); res->GetUDData(udd_C5Y,UDDY); res->GetUDData(udd_C5Z,UDDZ);
                conn_carts[0] = Cartesian(UDDX,UDDY,UDDZ);
              } else if(strncmp(res->atom[iat]->name," O2 ",4)==0) {
                res->GetUDData(udd_C2X,UDDX); res->GetUDData(udd_C2Y,UDDY); res->GetUDData(udd_C2Z,UDDZ);
                conn_carts[0] = Cartesian(UDDX,UDDY,UDDZ);
              }
              for (int i=0;i<nAtomBonds;i++){
                if(AtomBond[i].atom->residue&&AtomBond[i].atom->residue->seqNum!=res->seqNum){
                  int restype2 = molHnd->GetRestypeCode(AtomBond[i].atom->residue);
                  if(restype2==RESTYPE_SACH||restype2==RESTYPE_DSACH||restype2==RESTYPE_LSACH){
                    const double *col = atom_colour_vector.GetRGB(ii);
                    // Not atom coords! UDDData for C1.
                    AtomBond[i].atom->residue->GetUDData(udd_C1X,UDDX);
                    AtomBond[i].atom->residue->GetUDData(udd_C1Y,UDDY);
                    AtomBond[i].atom->residue->GetUDData(udd_C1Z,UDDZ);
                    conn_carts[1] = Cartesian(UDDX,UDDY,UDDZ);
                    conn_carts[1] = Cartesian::MidPoint(conn_carts);
                    CylinderElement *cylinder;
                    SphereElement *sphere;
                    if(have_stick_col){
                      cylinder = new CylinderElement(conn_carts,stick_col,conn_carts[1],interaction_cylinder_radius,1.0);
                      sphere = new SphereElement(conn_carts[0],stick_col,conn_carts[0],interaction_cylinder_radius,1.0,1);
                    } else {
                      cylinder = new CylinderElement(conn_carts,col,conn_carts[1],interaction_cylinder_radius,1.0);
                      sphere = new SphereElement(conn_carts[0],col,conn_carts[0],interaction_cylinder_radius,1.0,1);
                    }
                    cyls->add_primitive(cylinder);
                    spheres->add_primitive(sphere);
                    delete [] col;
                  }
                } 
              } 
            }
          }
        }
      }
    }
  }

  cyls->SetAlpha(obj.GetAlpha());
  tris->SetAlpha(obj.GetAlpha());
  spheres->SetAlpha(obj.GetAlpha());
  cyls->SetAccu(cylinders_accu);
  obj.add_primitive(cyls);
  obj.add_primitive(spheres);
  obj.add_primitive(tris);
  dash_cylinders->SetDashLength(interaction_dash_length);
  dash_cylinders->SetAccu(cylinders_accu);
  dash_cylinders->SetCapped(true);
  dash_lines_collection->SetSize(interaction_line_width);
  dash_lines_collection->SetDashLength(interaction_dash_length);
  obj.add_primitive(dash_cylinders);
  obj.add_primitive(dash_lines_collection);
}

void DrawBaseBlocks(Displayobject &obj, CMMANManager *molHnd, int selHnd, mmdb::Atom** selAtoms, int nSelAtoms, const AtomColourVector &atom_colour_vector,const CParamsManager &params ){
  TriangleCollection *polys = new TriangleCollection();
  int C5sel = molHnd->NewSelection();
  molHnd->Select(C5sel,STYPE_ATOM,0,"*",ANY_RES,"*",ANY_RES,"*","*","C5*","C","*",SKEY_NEW);
  molHnd->Select(C5sel,STYPE_ATOM,0,"*",ANY_RES,"*",ANY_RES,"*","*","C5'","C","*",SKEY_OR);
  for(int ii=0;ii<nSelAtoms;ii++){
    if(selAtoms[ii]->isInSelection(C5sel)){
      mmdb::Residue* res = selAtoms[ii]->GetResidue();
      if(res){
        int restype = molHnd->GetRestypeCode(res);
	//std::cout << res->name << " " << restype << std::endl;
        if(restype==RESTYPE_NUCL||restype==RESTYPE_DNA||restype==RESTYPE_RNA){
           const double *col = atom_colour_vector.GetRGB(ii);
           DrawBaseBlock(polys,res,col, params, selHnd);
           delete [] col;
        }
      }
    }
  }
  polys->SetAlpha(obj.GetAlpha());
  obj.add_primitive(polys);
}

//-----------------------------------------------------------------
void DrawBlobEllipsoids(Displayobject &obj,int selHnd,
     const Connectivity &connectivity, 
     const AtomColourVector &atom_colour_vector,
     const CParamsManager &params,
     const CParamsManager &global_params){
//-----------------------------------------------------------------
  int spheres_accu = global_params.GetInt("solid_quality")+2; 
  mmdb::Atom** atoms = connectivity.GetAtoms();
  int natoms = connectivity.GetNumberOfAtoms();
  if(natoms<1) return;
  SpheroidCollection *polys = new SpheroidCollection();
  polys->SetAccu(spheres_accu);

  Manager *molHnd = atoms[0]->GetChain()->GetCoordHierarchy();
  int previousTotalAtoms = 0;
  for(int i=0;i<molHnd->GetNumberOfChains(1);i++){
    mmdb::Chain* ch = molHnd->GetChain(1,i);
    mmdb::Atom** chainAtoms;
    int nChainAtoms;
    int chainSelHnd = molHnd->NewSelection();
    molHnd->SelectAtoms(chainSelHnd, 0,ch->GetChainID(),ANY_RES,"*",ANY_RES,"*","*","*","*","*",SKEY_NEW );
    molHnd->Select(chainSelHnd, STYPE_ATOM,selHnd,SKEY_AND );
    molHnd->GetSelIndex(chainSelHnd,chainAtoms,nChainAtoms);
    if(nChainAtoms>0){
      std::vector <Cartesian> carts = CartesiansFromAtoms(chainAtoms,nChainAtoms);
      std::vector<Cartesian> pca = Cartesian::PrincipalComponentAnalysis(carts);
      Cartesian centre = Cartesian::MidPoint(carts);
      const double *col = atom_colour_vector.GetRGB(previousTotalAtoms);
      SpheroidElement *spheroid = new SpheroidElement(centre,col,centre,2.0*sqrt(pca[4].get_x())+.7,2.0*sqrt(pca[4].get_y())+.7,2.0*sqrt(pca[4].get_z())+.7,pca[0],pca[1],pca[2],1.0,spheres_accu);
      polys->add_primitive(spheroid);
    }
    previousTotalAtoms += nChainAtoms;
    molHnd->DeleteSelection(chainSelHnd);
  }

  obj.add_primitive(polys);
}


//-----------------------------------------------------------------
void DrawImposterSpheres(Displayobject &obj,
     const Connectivity &connectivity, 
     const AtomColourVector &atom_colour_vector,
     const CParamsManager &params,
     const CParamsManager &global_params,
     const std::vector<double> &atomRadii){
//-----------------------------------------------------------------

  std::vector<std::vector<std::vector<int> > > conn_order_lists = connectivity.GetConnectivityLists();
  std::vector<std::vector<int> > conn_lists = conn_order_lists[0];
  std::vector<std::vector<int> > ext_conn_lists = connectivity.GetExternalConnectivityLists();

  mmdb::Atom** selAtoms = connectivity.GetAtoms();
  int nSelAtoms = connectivity.GetNumberOfAtoms();

  std::vector <Cartesian> int_carts = CartesiansFromAtoms(selAtoms,nSelAtoms);

  ImposterSphereCollection *collection = new ImposterSphereCollection();
  for(unsigned i=0;i<conn_lists.size();i++){
    const double *col = atom_colour_vector.GetRGB(i);
    ImposterSphereElement *p = new ImposterSphereElement(int_carts[i],col,int_carts[i],atomRadii[i],1,0);
    collection->add_primitive(p);
  }
  obj.add_primitive(collection);
}

//-----------------------------------------------------------------
void DrawCircles(Displayobject &obj,
     const Connectivity &connectivity, 
     const AtomColourVector &atom_colour_vector,
     const CParamsManager &params,
     const CParamsManager &global_params,
     int stick_colour) {
//-----------------------------------------------------------------

  std::vector<std::vector<std::vector<int> > > conn_order_lists = connectivity.GetConnectivityLists();
  std::vector<std::vector<int> > conn_lists = conn_order_lists[0];
  std::vector<std::vector<int> > ext_conn_lists = connectivity.GetExternalConnectivityLists();

  mmdb::Atom** selAtoms = connectivity.GetAtoms();
  int nSelAtoms = connectivity.GetNumberOfAtoms();

  double *stick_colour_array=0;
  if ( stick_colour > 0 ) stick_colour_array = RGBReps::GetColourP(stick_colour);

  std::vector <Cartesian> int_carts = CartesiansFromAtoms(selAtoms,nSelAtoms);

  // These need to be got from prefs
  double size = 0.3;
  double cylinders_size = 0.1;

  std::vector<Cartesian> carts(2);
  for(unsigned i=0;i<conn_lists.size();i++){
    const double *col = atom_colour_vector.GetRGB(i);
    bool warning = false;
    for(unsigned j=0;j<conn_lists[i].size();j++){
      Cartesian midpoint = Cartesian::MidPoint(int_carts[i],int_carts[conn_lists[i][j]]);
      // 0.85 is a fudge factor which seems to get rid of gaps at ends of bonds.
      double frac = 1-size*0.85/(int_carts[i]-midpoint).length();
      carts[1] = midpoint;
      carts[0] = midpoint+frac*(int_carts[i]-midpoint);
      Cartesian tmp_v = carts[0] -  carts[1];
      if (tmp_v.length() < 0.001||tmp_v.length() > 6   ) {
        if (!warning) {
          warning = true;
          std::cout << "Error drawing FLAT CYLINDER \n";
        }
        continue;
      }
      // Now we have to implement BillboardCylinderElement.
      BillboardCylinderElement *line;
      if ( stick_colour > 0 )
        line = new BillboardCylinderElement(carts,stick_colour_array,Cartesian::MidPoint(carts),cylinders_size,1.0);
      else
        line = new BillboardCylinderElement(carts,col,Cartesian::MidPoint(carts),cylinders_size,1.0);
      line->SetTextureStartXCoord(0.003);
      line->SetTextureEndXCoord(0.125);
      line->SetTextureStartYCoord(0.0);
      line->SetTextureEndYCoord(0.125);
      obj.add_billboard_primitive(line);
    }
    BillboardLabelledCircleElement *p = new BillboardLabelledCircleElement(int_carts[i],col,int_carts[i],size,1,0,"Foo");
    std::vector<double> coords = TexCoords::GetCoords(std::string(selAtoms[i]->element));
    p->SetTextureStartXCoord(coords[0]);
    p->SetTextureEndXCoord(coords[1]);
    p->SetTextureStartYCoord(coords[2]);
    p->SetTextureEndYCoord(coords[3]);
    obj.add_billboard_primitive(p);
  }
  if ( stick_colour > 0 && stick_colour_array ) delete [] stick_colour_array;
}


//-----------------------------------------------------------------
void DrawAnisoU(Displayobject &obj,
     int selHndin, mmdb::Atom** selAtoms, int nSelAtoms, 
     const AtomColourVector &atom_colour_vector, const std::vector<double> &atomRadii,
     int style, double scale, 
     const CParamsManager &global_params) {
//-----------------------------------------------------------------

  std::string label;
  Spheroid *ellipse;
  double a,b,c;
  Cartesian a_axis,b_axis,c_axis;

  int accu = global_params.GetInt("solid_quality");
  int show_axes=0, show_solid=0;
  if (style == SPHEROID_AXES) show_axes=1;
  if (style == SPHEROID_SOLID) show_solid=1;
  if (style == SPHEROID_SOLIDAXES) {
    show_axes=1;
    show_solid=1;
  }

  SpheroidCollection *polys = new SpheroidCollection();
  polys->SetAccu(accu);
  for ( int i = 0; i < nSelAtoms; i++ ) {
    Cartesian vertex = Cartesian(selAtoms[i]->x, 
                             selAtoms[i]->y,selAtoms[i]->z);
    const double *col = atom_colour_vector.GetRGB(i);
    if (selAtoms[i]->WhatIsSet&ASET_Anis_tFac) {
      double array[9] = { 
	  selAtoms[i]->u11,selAtoms[i]->u12, selAtoms[i]->u13, 
	  selAtoms[i]->u12,selAtoms[i]->u22, selAtoms[i]->u23, 
	  selAtoms[i]->u13,selAtoms[i]->u23, selAtoms[i]->u33 };
      matrix U(3,3,array);
      matrix A = atomRadii[i]*scale*U.Cholesky();
      matrix B(4,4);
      B(0,0) = A(0,0);
      B(1,0) = A(1,0);
      B(1,1) = A(1,1);
      B(2,0) = A(2,0);
      B(2,1) = A(2,1);
      B(2,2) = A(2,2);
      B(3,3) = 1.0;
      ellipse = new Spheroid(vertex, col, vertex,
		    B,
         	    1.0, accu , show_axes, show_solid );
      polys->add_primitive(ellipse);
    } else {
      matrix B(4,4);
      double Bfac = atomRadii[i]*scale * selAtoms[i]->tempFactor/(8*M_PI*M_PI) ;
      if(fabs(Bfac)>1e-5){
        B(0,0) = B(1,1) = B(2,2) = Bfac;
        B(3,3) = 1.0;
        ellipse = new Spheroid(vertex, col, vertex,
		      B,
         	      1.0, accu , show_axes, show_solid );
        polys->add_primitive(ellipse);
      }
    }
    delete [] col;
  }
  obj.add_primitive(polys);

}


Cartesian GetClosestSplinePoint(const Cartesian &cart, const SplineInfo &splineinfo){

  //double min_t = -1.0;
  double min_dist = 1e8;
  int idx = 0;
  int chain = 0;

  for(unsigned j=0;j<splineinfo.nasplines.size();j++){
    for(unsigned i=0;i<splineinfo.nasplines[j].size();i++){
      double dist = (cart-splineinfo.nasplines[j][i]).length();
      if(dist<min_dist&&dist>0.0){
        min_dist = dist;
        idx = i;
        chain = j;
      }
    }
  }
  return splineinfo.nasplines[chain][idx];
}

Cartesian GetClosestSplinePoint(const std::vector<Cartesian> &carts, const SplineInfo &splineinfo){
  Cartesian closest;

  Cartesian ls = carts[0];
  Cartesian le = carts[1];

  int idx = 0;
  int chain = 0;
  int step = 4;

  double min_t = -1.0;
  double min_dist = 1e8;

  for(unsigned j=0;j<splineinfo.nasplines.size();j++){
    for(unsigned i=0;i<splineinfo.nasplines[j].size();i+=step){
      std::vector<double> dist = DistanceBetweenPointAndLine(ls,le,splineinfo.nasplines[j][i]);
      if(dist[0]<min_dist&&dist[1]>0.0){
        min_dist = dist[0];
        idx = i;
        chain = j;
	min_t = dist[1];
      }
    }
  }

  if(idx<step) idx = step;
  if(idx>(int)(splineinfo.nasplines[chain].size())-step-1) idx = splineinfo.nasplines[chain].size()-step-1;

  /*
  if(idx<step){
    std::vector<double> dist = DistanceBetweenPointAndLine(ls,le,splineinfo.nasplines[chain][0]);
    min_dist = dist[0];
    min_t = dist[1];
    idx = 0;
    std::cout << "Beginning\n";
  } else if(idx>splineinfo.nasplines[chain].size()-step-1) {
    std::vector<double> dist = DistanceBetweenPointAndLine(ls,le,splineinfo.nasplines[chain].back());
    min_dist = dist[0];
    min_t = dist[1];
    idx = splineinfo.nasplines[chain].size()-1;
    std::cout << "End (" << splineinfo.nasplines[chain].size() << "\n";
  } else {
  */
    int idx2 = idx;
    for(int i=idx-step;i<idx+step;i++){
       std::vector<double> dist = DistanceBetweenPointAndLine(ls,le,splineinfo.nasplines[chain][i]);
       if(dist[0]<min_dist){
          min_dist = dist[0];
	  idx2 = i;
	  min_t = dist[1];
       }
    }
  //}
  idx = idx2;

  Cartesian MP=le-ls;
  double length = MP.length();
  MP.normalize();
  Cartesian n1 = splineinfo.n1_nasplines[chain][idx];
  n1.normalize(); // Already normalized ?
  double l = min_dist*(Cartesian::DotProduct(n1,MP));
  l = l/length;

  Cartesian sp = splineinfo.nasplines[chain][idx]; 
  Cartesian new_le = ls+min_t*(le-ls);
  Cartesian l1 = new_le - sp;
  double ang_test = Cartesian::DotProduct(l1,n1);

  if(ang_test>0&&(idx<step||idx>(int)(splineinfo.nasplines[chain].size())-step-1)){
    min_t += l;
  }else{
    min_t -= l;
  }

  if(min_t>0.0)
    return ls+min_t*(le-ls);

  return le;
}

void DrawBasePairs(Displayobject &obj, const CNABasePairs &bp, const SplineInfo &splineinfo, const CParamsManager &params){
  std::vector<std::pair<mmdb::Residue*,mmdb::Residue*> > base_pairs = bp.GetPairs();
  std::vector<std::pair<const double*,const double*> > colours = bp.GetColours();
  CylinderCollection *polys = new CylinderCollection();
  double cylinders_size = params.GetFloat("nucleic_stick_width");
  if (cylinders_size<0.02)cylinders_size=0.1;

  int cylinders_accu = 8;

  //std::cout << "DrawBasePairs " << base_pairs.size() << "\n";

  for(unsigned ii=0;ii<base_pairs.size();ii++){
    mmdb::Residue* res1 = base_pairs[ii].first;
    mmdb::Residue* res2 = base_pairs[ii].second;
    const double *col1 = colours[ii].first;
    const double *col2 = colours[ii].second;

    std::vector<Cartesian> cartsend = GetBasePairEnds(res1,res2,splineinfo);
    //std::cout << res1->name << " " << res2->name << "\n";
    //std::cout << cartsend[0] << "\n" << cartsend[1] << " " << (cartsend[0]-cartsend[1]).length() << "\n\n";
    if(cartsend.size()>0){
      std::vector<Cartesian> carts(2);
      std::vector<Cartesian> carts1(2);
      std::vector<Cartesian> carts2(2);
      Cartesian Mnew = Cartesian::MidPoint(cartsend[0],cartsend[1]);
      carts[0] = Mnew;
      carts[1] = cartsend[0];
      //carts2[0] = Mnew;
      //carts2[1] = cartsend[1];

      CylinderElement *line = new CylinderElement(carts,col1,carts[1],cylinders_size,1.0,cylinders_accu);
      polys->add_primitive(line);

      //line = new CylinderElement(carts2,col2,carts2[1],cylinders_size,1.0,cylinders_accu);
      //polys->add_primitive(line);

    }
  }
  polys->SetAlpha(obj.GetAlpha());
  obj.add_primitive(polys);
}

void DrawLipids(Displayobject &obj, const std::vector<MMUTLipid> &lipids, CMMANManager *molHnd, const AtomColourVector &atom_colour_vector, const CParamsManager &params, const CParamsManager &global_params){
#ifdef _DO_TIMINGS_
  clock_t t1 = clock();
#endif
  /* Need colours, etc. */

  if(lipids.size()<1) return;

  SpheroidCollection *polys = new SpheroidCollection();
  int spheres_accu=1;
  spheres_accu = global_params.GetInt("solid_quality");
  int spline_accu = 1+global_params.GetInt("solid_quality");
  int ribbon_accu = 4+global_params.GetInt("solid_quality")*4;

  //std::cout << "In routine to DRAWLIPIDS\n";

  //std::cout << "ribbon_accu: " << ribbon_accu << "\n";
  //std::cout << "spline_accu: " << spline_accu << "\n";
  std::cout << "There are " << lipids.size() << "\n";
  if(lipids.size()<1) return;

  // This requires that the whole lipid vector is from one
  // selHnd. It *should* be, but we probably want a container
  // class rather than std::vector of lipids.
  int selHnd = lipids[0].GetMainSelectionHandle();
  mmdb::Atom** atomTable_main_selection=0;
  int nAtoms_main_selection;
  molHnd->GetSelIndex ( selHnd, atomTable_main_selection, nAtoms_main_selection);
  //std::cout << nAtoms_main_selection << " atoms in main selection\n";
  if(nAtoms_main_selection==0) return;

  std::vector<int> serNums;
  for(int jj=0;jj<nAtoms_main_selection;jj++){
    serNums.push_back(atomTable_main_selection[jj]->serNum);
  }
  std::sort(serNums.begin(),serNums.end());

  for(unsigned ilipid=0;ilipid<lipids.size();ilipid++){
  MMUTLipid lipid = lipids[ilipid];
  std::vector<std::vector<Cartesian> > sorted_tail_carts = lipid.GetTailCartesians();
  std::vector<std::vector<int> > sorted_tail_serNums = lipid.GetTailSerNums();
  std::vector<std::vector<Cartesian> > all_head_carts = lipid.GetHeadCartesians();
  std::vector<std::vector<int> > all_head_serNums = lipid.GetHeadSerNums();

  //std::cout << "Heads:\n";
  std::vector<Cartesian> head_centres;
  for(unsigned ii=0;ii<all_head_carts.size();ii++){
     std::vector<int> head_serNums = all_head_serNums[ii];
     std::vector<Cartesian> head_carts = all_head_carts[ii];
     if(head_carts.size()>0){ //Hopefully
       Cartesian centre = Cartesian::MidPoint(head_carts);
       head_centres.push_back(centre);

       int atom_map=0;
       std::vector<int>::iterator ub = std::lower_bound(serNums.begin(),serNums.end(),head_serNums[0]);
        if(ub!=serNums.end()){
          //atom_map = (*ub)-1;
         atom_map = std::distance(serNums.begin(),ub);
        }
       //else std::cout << "Bummer\n";
       //std::cout << atom_map << "\n";
       /*
       for(int j=0;j<nAtoms_main_selection;j++) {
         if(head_serNums[0]==serNums[j]){
           atom_map = j;
           break;
         }
       }
       */

       //std::cout << "atom_map " << atom_map << "\n";
       const double *col = atom_colour_vector.GetRGB(atom_map);
       if(head_carts.size()>3){
       std::vector<Cartesian> pca = Cartesian::PrincipalComponentAnalysis(head_carts);
         SpheroidElement *spheroid = new SpheroidElement(centre,col,centre,2.0*sqrt(pca[4].get_x())+.7,2.0*sqrt(pca[4].get_y())+.7,2.0*sqrt(pca[4].get_z())+.7,pca[0],pca[1],pca[2],1.0,spheres_accu);
         polys->add_primitive(spheroid);
       }else if(head_carts.size()>1){
         double max_dist = 0.0;
         for(unsigned ih=0;ih<head_carts.size();ih++){
           if(LineLength(head_carts[ih],centre)>max_dist) max_dist = LineLength(head_carts[ih],centre);
         }
         SphereElement *spheroid = new SphereElement(centre,col,centre,max_dist+.7,1.0,spheres_accu);
         polys->add_primitive(spheroid);
       }else{
         SphereElement *spheroid = new SphereElement(centre,col,centre,1.0,1.0,spheres_accu);
         polys->add_primitive(spheroid);
       }
       delete [] col;
     //delete [] atoms;
     }
  }

  //std::cout << "Tails:\n";
  for(unsigned ii=0;ii<sorted_tail_serNums.size();ii++){
     std::vector<Cartesian> tail_carts = sorted_tail_carts[ii];
     std::vector<int> tail_serNums = sorted_tail_serNums[ii];
     std::vector<int> atom_map;
     for(unsigned i=0;i<tail_serNums.size();i++) {
       int serNum = tail_serNums[i];
       std::vector<int>::iterator ub = std::lower_bound(serNums.begin(),serNums.end(),serNum);
       if(ub!=serNums.end()) {
         //atom_map.push_back((*ub)-1);
         atom_map.push_back(std::distance(serNums.begin(),ub));
         //std::cout << "Index (1.1) " << std::distance(serNums.begin(),ub) << "\n";
       }
       /*
       for(int j=0;j<nAtoms_main_selection;j++) {
	  if(serNum==serNums[j]){
            atom_map.push_back(j);
            std::cout << "Index (2) " << j << "\n";
            break;
	  }
       }
       */
     }

     if(tail_carts.size()>3) {
       // Drawing from hereonin about ??% time
       int head_offset = 0;
       if(head_centres.size()>0){
	  head_offset = -1;
	  Cartesian posn = tail_carts[0];
	  Cartesian extra_posn = head_centres[0];
	  int first_t = 1;
	  double min_dist = LineLength(extra_posn,tail_carts[0]);
	  for(unsigned ic=1;ic<head_centres.size();ic++){
            double dist = LineLength(head_centres[ic],tail_carts[0]);
            if(dist<min_dist){
              dist = min_dist;             
	      extra_posn = head_centres[ic];
            }
	  }
	  for(unsigned ic=0;ic<head_centres.size();ic++){
            double dist = LineLength(head_centres[ic],tail_carts.back());
            if(dist<min_dist){
              dist = min_dist;             
	      extra_posn = head_centres[ic];
	      posn = tail_carts.back();
	      first_t = 0;
            }
	  }
	  Cartesian dirn = extra_posn-posn;
	  dirn.normalize();
	  if(first_t){
	    double wanted_length = LineLength(tail_carts[1],posn);
            tail_carts.insert(tail_carts.begin(),tail_carts[0]+wanted_length*dirn);
	  }else{
	    head_offset = 1;
	    double wanted_length = LineLength(tail_carts[tail_carts.size()-2],posn);
            tail_carts.insert(tail_carts.end(),tail_carts.back()+wanted_length*dirn);
	  }
       }

       std::vector<Cartesian> tail_n1;
       std::vector<Cartesian>::iterator t_iter = tail_carts.begin()+1;
       while(t_iter!=tail_carts.end()-1){
         tail_n1.push_back(Cartesian::CrossProduct((*(t_iter)-*(t_iter-1)),(*(t_iter+1)-*(t_iter))));
	 tail_n1.back().normalize();
         t_iter++;
       }
       tail_n1.push_back(tail_n1.back());
       tail_n1.insert(tail_n1.begin(),tail_n1[0]);

       t_iter = tail_n1.begin()+1;
       while(t_iter!=tail_n1.end()){
         if(Cartesian::DotProduct(*t_iter,*(t_iter-1))<0.0){
           *(t_iter) = -(*t_iter);
         }
         t_iter++;
       }
       //std::cout << "Size of spline inputs: " << tail_carts.size() << " " << tail_n1.size() << "\n";
       //std::cout.flush();
       // And this is what we draw.
       std::vector<Cartesian> tail_spline = BezierCurve(tail_carts,spline_accu);
       std::vector<Cartesian> tail_n1_spline = BezierCurve(tail_n1,spline_accu);
       std::vector<Cartesian> tail_n2_spline;
       std::vector<Cartesian> colour_vector;
       const double *col = atom_colour_vector.GetRGB(serNums[0]-1);
       if(head_offset>0){
         for(int iat=0;iat<(int)tail_carts.size()-1;iat++){
           const double *atcol = atom_colour_vector.GetRGB(atom_map[iat]);
	   for(int iacc=0;iacc<spline_accu;iacc++)
             colour_vector.push_back(Cartesian(atcol));
           delete [] atcol;
         }
         const double *atcol = atom_colour_vector.GetRGB(atom_map[tail_serNums.size()-1]);
	 for(int iacc=0;iacc<spline_accu;iacc++)
           colour_vector.push_back(Cartesian(atcol));
         delete [] atcol;
       } else {
         if(head_offset<0){
           const double *atcol = atom_colour_vector.GetRGB(atom_map[0]);
	   for(int iacc=0;iacc<spline_accu;iacc++)
             colour_vector.push_back(Cartesian(atcol));
           delete [] atcol;
         }
         for(int iat=-head_offset;iat<(int)tail_carts.size();iat++){
           const double *atcol = atom_colour_vector.GetRGB(atom_map[iat+head_offset]);
	   for(int iacc=0;iacc<spline_accu;iacc++)
             colour_vector.push_back(Cartesian(atcol));
           delete [] atcol;
         }
       }
       delete [] col;
       for(int iacc=0;iacc<spline_accu;iacc++)
         colour_vector.push_back(colour_vector.back());
       for(unsigned ispl=0;ispl<tail_spline.size()-1;ispl++){
	 Cartesian ts = tail_spline[ispl]-tail_spline[ispl+1];
	 Cartesian tn1s = tail_n1_spline[ispl];
	 //std::cout << Cartesian::DotProduct(ts,tn1s) << "\n";
	 ts.normalize();
	 tn1s.normalize();
         tail_n2_spline.push_back(Cartesian::CrossProduct(ts,tn1s));
	 tail_n2_spline.back().normalize();
       }
       tail_n2_spline.push_back(tail_n2_spline.back());
       float worm_width = params.GetFloat("worm_width");
       Ribbon *ribbon = new Ribbon(tail_spline,tail_n2_spline,tail_n1_spline,col,tail_spline[0],colour_vector,worm_width,worm_width,worm_width,1.0,ribbon_accu,spline_accu);
       obj.add_primitive(ribbon);
     }
  }
  }
  polys->SetAccu(spheres_accu);
  obj.add_primitive(polys);
#ifdef _DO_TIMINGS_
  clock_t t2 = clock();
  std::cout << "Time for lipid draw " << ((t2-t1)*1000.0/CLOCKS_PER_SEC)/1000.0<< "\n";
#endif
}

