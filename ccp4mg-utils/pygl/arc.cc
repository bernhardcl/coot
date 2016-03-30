/*
     pygl/arc.cc: CCP4MG Molecular Graphics Program
     Copyright (C) 2001-2008 University of York, CCLRC
     Copyright (C) 2009-2010 University of York
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

#include <math.h>
#ifndef M_PI
#define M_PI 3.141592653589793238462643
#endif

#define PIBY2 (M_PI * 2)
#include "cprimitive.h"

MGArc::MGArc() : ArcElement(){}
ArcElement::ArcElement() : Primitive(){vertices.push_back(Cartesian());}

MGArc::MGArc(const Cartesian &vertex, const double *colour_in, const Cartesian &origin_in, double radius, double width_in, double alpha_in, int textured_in) : ArcElement(vertex, colour_in, origin_in, radius, width_in, alpha_in, textured_in){};

void MGArc::draw(const double *override_colour, int selective_override){
  glDisable(GL_LIGHTING);
  glLineWidth(getWidth());
  if (useVBO) {
    drawArrays(override_colour,selective_override);
    return;
  }
  glBegin(GL_LINES);
  ArcElement::draw(override_colour,selective_override);
  glEnd();
}

ArcElement::ArcElement(const Cartesian &vertex, const double *colour_in, const Cartesian &origin_in, double size_in, double width_in, double alpha_in, int textured_in) : Primitive(){
  width = width_in;
  size = size_in;
  alpha = alpha_in;
  vertices.clear();
  vertices.push_back(vertex);
  colour[0] = colour_in[0];
  colour[1] = colour_in[1];
  colour[2] = colour_in[2];
  origin = origin_in;
  textured = textured_in;
  startAngle = 0;
  sweepAngle = 360;
  dash_length = 0.2;
  isDashed = false;

  double arcLength = sweepAngle*M_PI/180. * size;

  arrowLength = 0.1 * arcLength;
  arrowWidth = 0.5 * arrowLength;

  arrow_head = 0;
}

ArcElement::~ArcElement(){
}

MGArc::~MGArc(){
}

void ArcElement::set_draw_colour(const GLfloat *col){
  if(col)
    set_draw_colour_line_override(col);
  else
    set_draw_colour_line();
}

void ArcElement::draw(const double *override_colour, int selective_override){

  if (useVBO) {
    drawArrays(override_colour,selective_override);
    return;
  }
  int nsectors = 36;

  if(getIsDashed()){
    nsectors = int(round(2*M_PI*size / GetDashLength()));
  }
  if(nsectors>360) nsectors = 360;

  bool useMatrix = false;
  matrix m = getMatrix();
  if(m.get_rows()==4&&m.get_rows()==4) useMatrix = true;

  if(GetArrowHead()==2||GetArrowHead()==3){
    double cosTheta = 1.0 - (arrowLength*arrowLength)/(2.0*size*size);
    double theta = acos(cosTheta)*180./M_PI;
    if(sweepAngle<0) theta = -theta;
    double coneEndX = size * cos(startAngle*M_PI/180.);
    double coneEndY = size * sin(startAngle*M_PI/180.);
    double coneStartX = size * cos((startAngle+theta)*M_PI/180.);
    double coneStartY = size * sin((startAngle+theta)*M_PI/180.);
    Cartesian coneEnd = Cartesian(coneEndX,coneEndY,0);
    Cartesian coneStart = Cartesian(coneStartX,coneStartY,0);
    if(useMatrix){
      coneEnd = m.Transpose() * coneEnd;
      coneStart = m.Transpose() * coneStart;
    }

    coneStart += vertices[0];
    coneEnd += vertices[0];
    Cartesian vec = (coneEnd-coneStart);
    vec.normalize();

    Cartesian atmp,btmp;
    Cartesian xaxis(1,0,0);
    Cartesian yaxis(0,1,0);
    Cartesian zaxis(0,0,1);
    atmp = vec.CrossProduct(vec,zaxis);
    if(atmp.length() < 0.000000001){
      atmp = vec.CrossProduct(vec,yaxis);
      if(atmp.length() < 0.000000001){
        atmp = vec.CrossProduct(vec,xaxis);
      }
    }
    btmp = vec.CrossProduct(vec,atmp);
    btmp.normalize(arrowWidth);
    atmp.normalize(arrowWidth);
    Cartesian p1 = coneStart + btmp;
    Cartesian p2 = coneStart - btmp;
    Cartesian p3 = coneStart + atmp;
    Cartesian p4 = coneStart - atmp;
    glVertex3f(p1.get_x(),p1.get_y(),p1.get_z());
    glVertex3f(coneEnd.get_x(),coneEnd.get_y(),coneEnd.get_z());
    glVertex3f(p2.get_x(),p2.get_y(),p2.get_z());
    glVertex3f(coneEnd.get_x(),coneEnd.get_y(),coneEnd.get_z());
    glVertex3f(p3.get_x(),p3.get_y(),p3.get_z());
    glVertex3f(coneEnd.get_x(),coneEnd.get_y(),coneEnd.get_z());
    glVertex3f(p4.get_x(),p4.get_y(),p4.get_z());
    glVertex3f(coneEnd.get_x(),coneEnd.get_y(),coneEnd.get_z());
  }

  if(GetArrowHead()==1||GetArrowHead()==3){
    double cosTheta = 1.0 - (arrowLength*arrowLength)/(2.0*size*size);
    double theta = acos(cosTheta)*180./M_PI;
    if(sweepAngle<0) theta = -theta;
    double coneEndX = size * cos((startAngle+sweepAngle)*M_PI/180.);
    double coneEndY = size * sin((startAngle+sweepAngle)*M_PI/180.);
    double coneStartX = size * cos((startAngle+sweepAngle-theta)*M_PI/180.);
    double coneStartY = size * sin((startAngle+sweepAngle-theta)*M_PI/180.);
    Cartesian coneEnd = Cartesian(coneEndX,coneEndY,0);
    Cartesian coneStart = Cartesian(coneStartX,coneStartY,0);
    if(useMatrix){
      coneEnd = m.Transpose() * coneEnd;
      coneStart = m.Transpose() * coneStart;
    }
    coneStart += vertices[0];
    coneEnd += vertices[0];

    Cartesian vec = (coneEnd-coneStart);
    vec.normalize();

    Cartesian atmp,btmp;
    Cartesian xaxis(1,0,0);
    Cartesian yaxis(0,1,0);
    Cartesian zaxis(0,0,1);
    atmp = vec.CrossProduct(vec,zaxis);
    if(atmp.length() < 0.000000001){
      atmp = vec.CrossProduct(vec,yaxis);
      if(atmp.length() < 0.000000001){
        atmp = vec.CrossProduct(vec,xaxis);
      }
    }
    btmp = vec.CrossProduct(vec,atmp);
    btmp.normalize(arrowWidth);
    atmp.normalize(arrowWidth);
    Cartesian p1 = coneStart + btmp;
    Cartesian p2 = coneStart - btmp;
    Cartesian p3 = coneStart + atmp;
    Cartesian p4 = coneStart - atmp;
    glVertex3f(p1.get_x(),p1.get_y(),p1.get_z());
    glVertex3f(coneEnd.get_x(),coneEnd.get_y(),coneEnd.get_z());
    glVertex3f(p2.get_x(),p2.get_y(),p2.get_z());
    glVertex3f(coneEnd.get_x(),coneEnd.get_y(),coneEnd.get_z());
    glVertex3f(p3.get_x(),p3.get_y(),p3.get_z());
    glVertex3f(coneEnd.get_x(),coneEnd.get_y(),coneEnd.get_z());
    glVertex3f(p4.get_x(),p4.get_y(),p4.get_z());
    glVertex3f(coneEnd.get_x(),coneEnd.get_y(),coneEnd.get_z());
  }

  int iloop = 0;

  float sa,ea;
  if(sweepAngle>0){
    sa = startAngle;
    ea = startAngle+sweepAngle;
  }else{
   sa = startAngle+sweepAngle;
   ea = startAngle;
  }
  for(float j=sa;j<ea;j=j+360/nsectors,iloop++){
    double phi = (double)j/360.0 * PIBY2;
    double phi2 = (double)(j+360.0/nsectors)/360.0 * PIBY2;
    if(sweepAngle>0&&j+360.0/nsectors>startAngle+sweepAngle)  phi2 = (double)(startAngle+sweepAngle)/360.0 * PIBY2;
    if(sweepAngle<0&&j+360.0/nsectors>startAngle) phi2 = (double)(startAngle)/360.0 * PIBY2;

    if(isDashed&&(iloop%2==0)) continue;
    double x = size * cos(phi);
    double y = size * sin(phi);
    double z = 0;
    if(useMatrix) {
      Cartesian cart(x,y,z);
      cart = m.Transpose()*cart+vertices[0];
      glVertex3f(cart.get_x(),cart.get_y(),cart.get_z());
    } else {
      glVertex3f(x+vertices[0].get_x(),y+vertices[0].get_y(),z+vertices[0].get_z());
    }
    x = size * cos(phi2);
    y = size * sin(phi2);
    if(useMatrix) {
      Cartesian cart(x,y,z);
      cart = m.Transpose()*cart+vertices[0];
      glVertex3f(cart.get_x(),cart.get_y(),cart.get_z());
    } else {
      glVertex3f(x+vertices[0].get_x(),y+vertices[0].get_y(),z+vertices[0].get_z());
    }
  }


}

std::vector<Primitive*> ArcElement::GetSimplePrimitivesArrays(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin, int start, int end){
  std::vector<Primitive*> a;

  int nsectors = 36;

  if(getIsDashed()){
    nsectors = int(round(2*M_PI*size / GetDashLength()));
  }
  if(nsectors>360) nsectors = 360;

  nVerts = nsectors*2;
  if(GetArrowHead()==1||GetArrowHead()==2)
    nVerts += 8;
  if(GetArrowHead()==3)
    nVerts += 16;
  nRealVerts = nVerts * 3;
  bufferVertices = new float[nRealVerts];
  bufferNormals = 0;
  bufferColors = new float[nRealVerts*4/3];
  bufferIndices = new GLuint[nVerts];

  bool useMatrix = false;
  matrix m = getMatrix();
  if(m.get_rows()==4&&m.get_rows()==4) useMatrix = true;

  int iloop = 0;

  float sa,ea;
  if(sweepAngle>0){
    sa = startAngle;
    ea = startAngle+sweepAngle;
  }else{
   sa = startAngle+sweepAngle;
   ea = startAngle;
  }
  int maxInd = 0;
  int il = 0;
  for(float j=sa;j<ea;j=j+360/nsectors,iloop++){
    double phi = (double)j/360.0 * PIBY2;
    double phi2 = (double)(j+360.0/nsectors)/360.0 * PIBY2;
    if(sweepAngle>0&&j+360.0/nsectors>startAngle+sweepAngle)  phi2 = (double)(startAngle+sweepAngle)/360.0 * PIBY2;
    if(sweepAngle<0&&j+360.0/nsectors>startAngle) phi2 = (double)(startAngle)/360.0 * PIBY2;

    if(isDashed&&(iloop%2==0)) continue;
    double x = size * cos(phi);
    double y = size * sin(phi);
    double z = 0;
    bufferIndices[il] = il;
    bufferIndices[il+1] = il+1;
    
    bufferColors[4*il]   = colour[0];
    bufferColors[4*il+1] = colour[1];
    bufferColors[4*il+2] = colour[2];
    bufferColors[4*il+3] = alpha;
    bufferColors[4*(il+1)]   = colour[0];
    bufferColors[4*(il+1)+1] = colour[1];
    bufferColors[4*(il+1)+2] = colour[2];
    bufferColors[4*(il+1)+3] = alpha;

    if(useMatrix) {
      Cartesian cart(x,y,z);
      cart = m.Transpose()*cart+vertices[0];
      //glVertex3f(cart.get_x(),cart.get_y(),cart.get_z());
      bufferVertices[3*il]   = cart.get_x();
      bufferVertices[3*il+1] = cart.get_y();
      bufferVertices[3*il+2] = cart.get_z();
    } else {
      //glVertex3f(x+vertices[0].get_x(),y+vertices[0].get_y(),z+vertices[0].get_z());
      bufferVertices[3*il]   = x+vertices[0].get_x();
      bufferVertices[3*il+1] = y+vertices[0].get_y();
      bufferVertices[3*il+2] = z+vertices[0].get_z();
    }
    x = size * cos(phi2);
    y = size * sin(phi2);
    if(useMatrix) {
      Cartesian cart(x,y,z);
      cart = m.Transpose()*cart+vertices[0];
      //glVertex3f(cart.get_x(),cart.get_y(),cart.get_z());
      bufferVertices[3*(il+1)]   = cart.get_x();
      bufferVertices[3*(il+1)+1] = cart.get_y();
      bufferVertices[3*(il+1)+2] = cart.get_z();
    } else {
      //glVertex3f(x+vertices[0].get_x(),y+vertices[0].get_y(),z+vertices[0].get_z());
      bufferVertices[3*(il+1)]   = x+vertices[0].get_x();
      bufferVertices[3*(il+1)+1] = y+vertices[0].get_y();
      bufferVertices[3*(il+1)+2] = z+vertices[0].get_z();
    }
    maxInd = il;
    il+=2;
  }

  if(GetArrowHead()==2||GetArrowHead()==3){
    double cosTheta = 1.0 - (arrowLength*arrowLength)/(2.0*size*size);
    double theta = acos(cosTheta)*180./M_PI;
    if(sweepAngle<0) theta = -theta;
    double coneEndX = size * cos(startAngle*M_PI/180.);
    double coneEndY = size * sin(startAngle*M_PI/180.);
    double coneStartX = size * cos((startAngle+theta)*M_PI/180.);
    double coneStartY = size * sin((startAngle+theta)*M_PI/180.);
    Cartesian coneEnd = Cartesian(coneEndX,coneEndY,0);
    Cartesian coneStart = Cartesian(coneStartX,coneStartY,0);
    if(useMatrix){
      coneEnd = m.Transpose() * coneEnd;
      coneStart = m.Transpose() * coneStart;
    }

    coneStart += vertices[0];
    coneEnd += vertices[0];
    Cartesian vec = (coneEnd-coneStart);
    vec.normalize();

    Cartesian atmp,btmp;
    Cartesian xaxis(1,0,0);
    Cartesian yaxis(0,1,0);
    Cartesian zaxis(0,0,1);
    atmp = vec.CrossProduct(vec,zaxis);
    if(atmp.length() < 0.000000001){
      atmp = vec.CrossProduct(vec,yaxis);
      if(atmp.length() < 0.000000001){
        atmp = vec.CrossProduct(vec,xaxis);
      }
    }
    btmp = vec.CrossProduct(vec,atmp);
    btmp.normalize(arrowWidth);
    atmp.normalize(arrowWidth);
    Cartesian p1 = coneStart + btmp;
    Cartesian p2 = coneStart - btmp;
    Cartesian p3 = coneStart + atmp;
    Cartesian p4 = coneStart - atmp;
    bufferIndices[il] = il;
      bufferVertices[3*(il)]   = p1.get_x();
      bufferVertices[3*(il)+1] = p1.get_y();
      bufferVertices[3*(il)+2] = p1.get_z();
      bufferColors[4*(il)]   = colour[0];
      bufferColors[4*(il)+1] = colour[1];
      bufferColors[4*(il)+2] = colour[2];
      bufferColors[4*(il)+3] = alpha;
    bufferIndices[il+1] = il+1;
      bufferVertices[3*(il+1)]   = coneEnd.get_x();
      bufferVertices[3*(il+1)+1] = coneEnd.get_y();
      bufferVertices[3*(il+1)+2] = coneEnd.get_z();
      bufferColors[4*(il+1)]   = colour[0];
      bufferColors[4*(il+1)+1] = colour[1];
      bufferColors[4*(il+1)+2] = colour[2];
      bufferColors[4*(il+1)+3] = alpha;
    bufferIndices[il+2] = il+2;
      bufferVertices[3*(il+2)]   = p2.get_x();
      bufferVertices[3*(il+2)+1] = p2.get_y();
      bufferVertices[3*(il+2)+2] = p2.get_z();
      bufferColors[4*(il+2)]   = colour[0];
      bufferColors[4*(il+2)+1] = colour[1];
      bufferColors[4*(il+2)+2] = colour[2];
      bufferColors[4*(il+2)+3] = alpha;
    bufferIndices[il+3] = il+3;
      bufferVertices[3*(il+3)]   = coneEnd.get_x();
      bufferVertices[3*(il+3)+1] = coneEnd.get_y();
      bufferVertices[3*(il+3)+2] = coneEnd.get_z();
      bufferColors[4*(il+3)]   = colour[0];
      bufferColors[4*(il+3)+1] = colour[1];
      bufferColors[4*(il+3)+2] = colour[2];
      bufferColors[4*(il+3)+3] = alpha;
    bufferIndices[il+4] = il+4;
      bufferVertices[3*(il+4)]   = p3.get_x();
      bufferVertices[3*(il+4)+1] = p3.get_y();
      bufferVertices[3*(il+4)+2] = p3.get_z();
      bufferColors[4*(il+4)]   = colour[0];
      bufferColors[4*(il+4)+1] = colour[1];
      bufferColors[4*(il+4)+2] = colour[2];
      bufferColors[4*(il+4)+3] = alpha;
    bufferIndices[il+5] = il+5;
      bufferVertices[3*(il+5)]   = coneEnd.get_x();
      bufferVertices[3*(il+5)+1] = coneEnd.get_y();
      bufferVertices[3*(il+5)+2] = coneEnd.get_z();
      bufferColors[4*(il+5)]   = colour[0];
      bufferColors[4*(il+5)+1] = colour[1];
      bufferColors[4*(il+5)+2] = colour[2];
      bufferColors[4*(il+5)+3] = alpha;
    bufferIndices[il+6] = il+6;
      bufferVertices[3*(il+6)]   = p4.get_x();
      bufferVertices[3*(il+6)+1] = p4.get_y();
      bufferVertices[3*(il+6)+2] = p4.get_z();
      bufferColors[4*(il+6)]   = colour[0];
      bufferColors[4*(il+6)+1] = colour[1];
      bufferColors[4*(il+6)+2] = colour[2];
      bufferColors[4*(il+6)+3] = alpha;
    bufferIndices[il+7] = il+7;
    maxInd = il+7;
      bufferVertices[3*(il+7)]   = coneEnd.get_x();
      bufferVertices[3*(il+7)+1] = coneEnd.get_y();
      bufferVertices[3*(il+7)+2] = coneEnd.get_z();
      bufferColors[4*(il+7)]   = colour[0];
      bufferColors[4*(il+7)+1] = colour[1];
      bufferColors[4*(il+7)+2] = colour[2];
      bufferColors[4*(il+7)+3] = alpha;
  }

  int arrowSkip = 0;
  if(GetArrowHead()==3) arrowSkip = 8;
  if(GetArrowHead()==1||GetArrowHead()==3){
    double cosTheta = 1.0 - (arrowLength*arrowLength)/(2.0*size*size);
    double theta = acos(cosTheta)*180./M_PI;
    if(sweepAngle<0) theta = -theta;
    double coneEndX = size * cos((startAngle+sweepAngle)*M_PI/180.);
    double coneEndY = size * sin((startAngle+sweepAngle)*M_PI/180.);
    double coneStartX = size * cos((startAngle+sweepAngle-theta)*M_PI/180.);
    double coneStartY = size * sin((startAngle+sweepAngle-theta)*M_PI/180.);
    Cartesian coneEnd = Cartesian(coneEndX,coneEndY,0);
    Cartesian coneStart = Cartesian(coneStartX,coneStartY,0);
    if(useMatrix){
      coneEnd = m.Transpose() * coneEnd;
      coneStart = m.Transpose() * coneStart;
    }
    coneStart += vertices[0];
    coneEnd += vertices[0];

    Cartesian vec = (coneEnd-coneStart);
    vec.normalize();

    Cartesian atmp,btmp;
    Cartesian xaxis(1,0,0);
    Cartesian yaxis(0,1,0);
    Cartesian zaxis(0,0,1);
    atmp = vec.CrossProduct(vec,zaxis);
    if(atmp.length() < 0.000000001){
      atmp = vec.CrossProduct(vec,yaxis);
      if(atmp.length() < 0.000000001){
        atmp = vec.CrossProduct(vec,xaxis);
      }
    }
    btmp = vec.CrossProduct(vec,atmp);
    btmp.normalize(arrowWidth);
    atmp.normalize(arrowWidth);
    Cartesian p1 = coneStart + btmp;
    Cartesian p2 = coneStart - btmp;
    Cartesian p3 = coneStart + atmp;
    Cartesian p4 = coneStart - atmp;
    bufferIndices[il+arrowSkip] = il+arrowSkip;
      bufferVertices[3*(il+arrowSkip)]   = p1.get_x();
      bufferVertices[3*(il+arrowSkip)+1] = p1.get_y();
      bufferVertices[3*(il+arrowSkip)+2] = p1.get_z();
      bufferColors[4*(il+arrowSkip)]   = colour[0];
      bufferColors[4*(il+arrowSkip)+1] = colour[1];
      bufferColors[4*(il+arrowSkip)+2] = colour[2];
      bufferColors[4*(il+arrowSkip)+3] = alpha;
    bufferIndices[il+arrowSkip+1] = il+arrowSkip+1;
      bufferVertices[3*(il+arrowSkip+1)]   = coneEnd.get_x();
      bufferVertices[3*(il+arrowSkip+1)+1] = coneEnd.get_y();
      bufferVertices[3*(il+arrowSkip+1)+2] = coneEnd.get_z();
      bufferColors[4*(il+arrowSkip+1)]   = colour[0];
      bufferColors[4*(il+arrowSkip+1)+1] = colour[1];
      bufferColors[4*(il+arrowSkip+1)+2] = colour[2];
      bufferColors[4*(il+arrowSkip+1)+3] = alpha;
    bufferIndices[il+arrowSkip+2] = il+arrowSkip+2;
      bufferVertices[3*(il+arrowSkip+2)]   = p2.get_x();
      bufferVertices[3*(il+arrowSkip+2)+1] = p2.get_y();
      bufferVertices[3*(il+arrowSkip+2)+2] = p2.get_z();
      bufferColors[4*(il+arrowSkip+2)]   = colour[0];
      bufferColors[4*(il+arrowSkip+2)+1] = colour[1];
      bufferColors[4*(il+arrowSkip+2)+2] = colour[2];
      bufferColors[4*(il+arrowSkip+2)+3] = alpha;
    bufferIndices[il+arrowSkip+3] = il+arrowSkip+3;
      bufferVertices[3*(il+arrowSkip+3)]   = coneEnd.get_x();
      bufferVertices[3*(il+arrowSkip+3)+1] = coneEnd.get_y();
      bufferVertices[3*(il+arrowSkip+3)+2] = coneEnd.get_z();
      bufferColors[4*(il+arrowSkip+3)]   = colour[0];
      bufferColors[4*(il+arrowSkip+3)+1] = colour[1];
      bufferColors[4*(il+arrowSkip+3)+2] = colour[2];
      bufferColors[4*(il+arrowSkip+3)+3] = alpha;
    bufferIndices[il+arrowSkip+4] = il+arrowSkip+4;
      bufferVertices[3*(il+arrowSkip+4)]   = p3.get_x();
      bufferVertices[3*(il+arrowSkip+4)+1] = p3.get_y();
      bufferVertices[3*(il+arrowSkip+4)+2] = p3.get_z();
      bufferColors[4*(il+arrowSkip+4)]   = colour[0];
      bufferColors[4*(il+arrowSkip+4)+1] = colour[1];
      bufferColors[4*(il+arrowSkip+4)+2] = colour[2];
      bufferColors[4*(il+arrowSkip+4)+3] = alpha;
    bufferIndices[il+arrowSkip+5] = il+arrowSkip+5;
      bufferVertices[3*(il+arrowSkip+5)]   = coneEnd.get_x();
      bufferVertices[3*(il+arrowSkip+5)+1] = coneEnd.get_y();
      bufferVertices[3*(il+arrowSkip+5)+2] = coneEnd.get_z();
      bufferColors[4*(il+arrowSkip+5)]   = colour[0];
      bufferColors[4*(il+arrowSkip+5)+1] = colour[1];
      bufferColors[4*(il+arrowSkip+5)+2] = colour[2];
      bufferColors[4*(il+arrowSkip+5)+3] = alpha;
    bufferIndices[il+arrowSkip+6] = il+arrowSkip+6;
      bufferVertices[3*(il+arrowSkip+6)]   = p4.get_x();
      bufferVertices[3*(il+arrowSkip+6)+1] = p4.get_y();
      bufferVertices[3*(il+arrowSkip+6)+2] = p4.get_z();
      bufferColors[4*(il+arrowSkip+6)]   = colour[0];
      bufferColors[4*(il+arrowSkip+6)+1] = colour[1];
      bufferColors[4*(il+arrowSkip+6)+2] = colour[2];
      bufferColors[4*(il+arrowSkip+6)+3] = alpha;
    bufferIndices[il+arrowSkip+7] = il+arrowSkip+7;
    maxInd = il+arrowSkip+7;
      bufferVertices[3*(il+arrowSkip+7)]   = coneEnd.get_x();
      bufferVertices[3*(il+arrowSkip+7)+1] = coneEnd.get_y();
      bufferVertices[3*(il+arrowSkip+7)+2] = coneEnd.get_z();
      bufferColors[4*(il+arrowSkip+7)]   = colour[0];
      bufferColors[4*(il+arrowSkip+7)+1] = colour[1];
      bufferColors[4*(il+arrowSkip+7)+2] = colour[2];
      bufferColors[4*(il+arrowSkip+7)+3] = alpha;
  }
  nVerts = maxInd + 2;
  nRealVerts = nVerts * 3;

  polygonType = GL_LINES;

  return a;
}

std::vector<Primitive*> ArcElement::GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin, int start, int end) const {
  std::cout << "ArcElement::GetSimplePrimitives not implemented\n";
  std::vector<Primitive*> a;
  return a;
}

void ArcElement::DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const Volume &v){}

void ArcElement::DrawPostScript(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, double xoff, double yoff, double xscale, double yscale, double xscaleps, const Volume &v){}
