/*
     pygl/torus.cc: CCP4MG Molecular Graphics Program
     Copyright (C) 2001-2008 University of York, CCLRC
     Copyright (C) 2009-2010 University of York

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

Torus::Torus() : TorusElement(){}
TorusElement::TorusElement() : Primitive(){vertices.push_back(Cartesian());}

Torus::Torus(const Cartesian &vertex, const double *colour_in, const Cartesian &origin_in, double majorRadius, double minorRadius,  double alpha_in, int textured_in) : TorusElement(vertex, colour_in, origin_in, majorRadius, minorRadius, alpha_in, textured_in){};

void Torus::draw(const double *override_colour, int selective_override){
  glEnable(GL_LIGHTING);
  if (useVBO) {
    drawArrays(override_colour,selective_override);
    return;
  }
  TorusElement::draw(override_colour,selective_override);
}

TorusElement::TorusElement(const Cartesian &vertex, const double *colour_in, const Cartesian &origin_in, double majorRadius_in, double minorRadius_in,  double alpha_in, int textured_in) : Primitive(){
  majorRadius = majorRadius_in;
  minorRadius = minorRadius_in;
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
}

TorusElement::~TorusElement(){
}

Torus::~Torus(){
}

void TorusElement::set_draw_colour(const GLfloat *col){
  if(col)
    set_draw_colour_poly_override(col);
  else
    set_draw_colour_poly();
}

void TorusElement::draw(const double *override_colour, int selective_override){

  int nsectors = 36;

  bool useMatrix = false;
  matrix m = getMatrix();
  if(m.get_rows()==4&&m.get_rows()==4) useMatrix = true;

  glPushMatrix();
  glTranslatef(vertices[0].get_x(),vertices[0].get_y(),vertices[0].get_z());

  if(useMatrix) {
    glMultMatrixd(getMatrix().to_dp());
  }

  if((360-sweepAngle)>1e-6){
      GLUquadric *qobj = gluNewQuadric();
      double j = startAngle;
      double phi = (double)j/360.0 * PIBY2;
      double x = majorRadius*cos(phi);
      double y = majorRadius*sin(phi);
      double z = 0.0;
      glPushMatrix();
      glTranslatef(x,y,z);
      Quat q1(0,1,0,1,90);
      Quat q2(0,0,1,1,startAngle-90);
      glMultMatrixd(q2.getInvMatrix().to_dp());
      glMultMatrixd(q1.getInvMatrix().to_dp());
      gluDisk(qobj,0.0,minorRadius,nsectors,16);
      glPopMatrix();
      j = startAngle+sweepAngle;
      phi = (double)j/360.0 * PIBY2;
      x = majorRadius*cos(phi);
      y = majorRadius*sin(phi);
      z = 0.0;
      glPushMatrix();
      glTranslatef(x,y,z);
      Quat q1end(0,1,0,1,-90);
      Quat q2end(0,0,1,1,(startAngle+sweepAngle)-90);
      glMultMatrixd(q2end.getInvMatrix().to_dp());
      glMultMatrixd(q1end.getInvMatrix().to_dp());
      gluDisk(qobj,0.0,minorRadius,nsectors,16);
      glPopMatrix();
      gluDeleteQuadric(qobj);
  }

  if(getIsDashed()){
    nsectors = int(round(2*M_PI*getMajorRadius() / GetDashLength()));
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
    glBegin(GL_QUAD_STRIP);
    for(int i=0;i<=360;i=i+360/nsectors){

      double theta = (double)i/360.0 * PIBY2;

      double x = (majorRadius +  minorRadius * cos(theta)) * cos(phi);
      double y = (majorRadius +  minorRadius * cos(theta)) * sin(phi);
      double z = minorRadius * sin(theta);
      double norm_x = cos(theta) * cos(phi);
      double norm_y = cos(theta) * sin(phi);
      double norm_z = sin(theta);
      glNormal3f(norm_x,norm_y,norm_z);
      glVertex3f(x,y,z);
      x = (majorRadius +  minorRadius * cos(theta)) * cos(phi2);
      y = (majorRadius +  minorRadius * cos(theta)) * sin(phi2);
      z = minorRadius * sin(theta);
      norm_x = cos(theta) * cos(phi2);
      norm_y = cos(theta) * sin(phi2);
      glNormal3f(norm_x,norm_y,norm_z);
      glVertex3f(x,y,z);

    }
    glEnd();
  }

  glPopMatrix();
   
}

std::vector<Primitive*> TorusElement::GetSimplePrimitivesArrays(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin, int start, int end){

  int nsectors = 36;

  float sa,ea;
  if(sweepAngle>0){
    sa = startAngle;
    ea = startAngle+sweepAngle;
  }else{
   sa = startAngle+sweepAngle;
   ea = startAngle;
  }

  int nquads = 0;
  for(float j=sa;j<ea;j=j+360/nsectors){
    for(int i=0;i<=360;i=i+360/nsectors){
       nquads++;
    }
  }
  nVerts = 4*nquads;
  nRealVerts = nVerts * 3;
  nVerts = 6*nquads;
  bufferVertices = new float[nRealVerts];
  bufferNormals = new float[nRealVerts];
  bufferColors = new float[nRealVerts*4/3];
  bufferIndices = new GLuint[nVerts];
  int ivert = 0;
  int indx = 0;

  bool useMatrix = false;
  matrix m = getMatrix();
  if(m.get_rows()==4&&m.get_rows()==4) useMatrix = true;

  for(float j=sa;j<ea;j=j+360/nsectors){
    double phi = (double)j/360.0 * PIBY2;
    double phi2 = (double)(j+360.0/nsectors)/360.0 * PIBY2;
    if(sweepAngle>0&&j+360.0/nsectors>startAngle+sweepAngle)  phi2 = (double)(startAngle+sweepAngle)/360.0 * PIBY2;
    if(sweepAngle<0&&j+360.0/nsectors>startAngle) phi2 = (double)(startAngle)/360.0 * PIBY2;
    for(int i=0;i<=360;i=i+360/nsectors){

      double theta = (double)i/360.0 * PIBY2;
      double theta2 = (double)(i+360.0/nsectors)/360.0 * PIBY2;

      double x,y,z;
      double norm_x,norm_y,norm_z;
      if(useMatrix){
        x = (majorRadius +  minorRadius * cos(theta)) * cos(phi);
        y = (majorRadius +  minorRadius * cos(theta)) * sin(phi);
        z = minorRadius * sin(theta);
        Cartesian cart = m.Transpose()*Cartesian(x,y,z) + vertices[0];
        x = cart.get_x();
        y = cart.get_y();
        z = cart.get_z();
        norm_x = cos(theta) * cos(phi);
        norm_y = cos(theta) * sin(phi);
        norm_z = sin(theta);
        Cartesian norm = m.Transpose()*Cartesian(norm_x,norm_y,norm_z);
        norm_x = norm.get_x();
        norm_y = norm.get_y();
        norm_z = norm.get_z();
      } else {
        x = (majorRadius +  minorRadius * cos(theta)) * cos(phi) + vertices[0].get_x();
        y = (majorRadius +  minorRadius * cos(theta)) * sin(phi) + vertices[0].get_y();
        z = minorRadius * sin(theta) + vertices[0].get_z();
        norm_x = cos(theta) * cos(phi);
        norm_y = cos(theta) * sin(phi);
        norm_z = sin(theta);
      }
      bufferVertices[3*ivert]   = x;
      bufferVertices[3*ivert+1] = y;
      bufferVertices[3*ivert+2] = z;
      bufferNormals[3*ivert]    = norm_x;
      bufferNormals[3*ivert+1]  = norm_y;
      bufferNormals[3*ivert+2]  = norm_z;
      bufferColors[4*ivert]     = colour[0];
      bufferColors[4*ivert+1]   = colour[1];
      bufferColors[4*ivert+2]   = colour[2];
      bufferColors[4*ivert+3]   = alpha;
      ivert++;

      if(useMatrix){
        x = (majorRadius +  minorRadius * cos(theta)) * cos(phi2);
        y = (majorRadius +  minorRadius * cos(theta)) * sin(phi2);
        z = minorRadius * sin(theta);
        Cartesian cart = m.Transpose()*Cartesian(x,y,z) + vertices[0];
        x = cart.get_x();
        y = cart.get_y();
        z = cart.get_z();
        norm_x = cos(theta) * cos(phi2);
        norm_y = cos(theta) * sin(phi2);
        norm_z = sin(theta);
        Cartesian norm = m.Transpose()*Cartesian(norm_x,norm_y,norm_z);
        norm_x = norm.get_x();
        norm_y = norm.get_y();
        norm_z = norm.get_z();
      } else {
        x = (majorRadius +  minorRadius * cos(theta)) * cos(phi2) + vertices[0].get_x();
        y = (majorRadius +  minorRadius * cos(theta)) * sin(phi2) + vertices[0].get_y();
        z = minorRadius * sin(theta) + vertices[0].get_z();
        norm_x = cos(theta) * cos(phi2);
        norm_y = cos(theta) * sin(phi2);
      }
      bufferVertices[3*ivert]   = x;
      bufferVertices[3*ivert+1] = y;
      bufferVertices[3*ivert+2] = z;
      bufferNormals[3*ivert]    = norm_x;
      bufferNormals[3*ivert+1]  = norm_y;
      bufferNormals[3*ivert+2]  = norm_z;
      bufferColors[4*ivert]     = colour[0];
      bufferColors[4*ivert+1]   = colour[1];
      bufferColors[4*ivert+2]   = colour[2];
      bufferColors[4*ivert+3]   = alpha;
      ivert++;

      if(useMatrix){
        x = (majorRadius +  minorRadius * cos(theta2)) * cos(phi2);
        y = (majorRadius +  minorRadius * cos(theta2)) * sin(phi2);
        z = minorRadius * sin(theta2);
        Cartesian cart = m.Transpose()*Cartesian(x,y,z) + vertices[0];
        x = cart.get_x();
        y = cart.get_y();
        z = cart.get_z();
        norm_x = cos(theta2) * cos(phi2);
        norm_y = cos(theta2) * sin(phi2);
        norm_z = sin(theta2);
        Cartesian norm = m.Transpose()*Cartesian(norm_x,norm_y,norm_z);
        norm_x = norm.get_x();
        norm_y = norm.get_y();
        norm_z = norm.get_z();
      } else {
        x = (majorRadius +  minorRadius * cos(theta2)) * cos(phi2) + vertices[0].get_x();
        y = (majorRadius +  minorRadius * cos(theta2)) * sin(phi2) + vertices[0].get_y();
        z = minorRadius * sin(theta2) + vertices[0].get_z();
        norm_x = cos(theta2) * cos(phi2);
        norm_y = cos(theta2) * sin(phi2);
        norm_z = sin(theta2);
      }
      bufferVertices[3*ivert]   = x;
      bufferVertices[3*ivert+1] = y;
      bufferVertices[3*ivert+2] = z;
      bufferNormals[3*ivert]    = norm_x;
      bufferNormals[3*ivert+1]  = norm_y;
      bufferNormals[3*ivert+2]  = norm_z;
      bufferColors[4*ivert]     = colour[0];
      bufferColors[4*ivert+1]   = colour[1];
      bufferColors[4*ivert+2]   = colour[2];
      bufferColors[4*ivert+3]   = alpha;
      ivert++;

      if(useMatrix){
        x = (majorRadius +  minorRadius * cos(theta2)) * cos(phi);
        y = (majorRadius +  minorRadius * cos(theta2)) * sin(phi);
        z = minorRadius * sin(theta2);
        Cartesian cart = m.Transpose()*Cartesian(x,y,z) + vertices[0];
        x = cart.get_x();
        y = cart.get_y();
        z = cart.get_z();
        norm_x = cos(theta2) * cos(phi);
        norm_y = cos(theta2) * sin(phi);
        norm_z = sin(theta2);
        Cartesian norm = m.Transpose()*Cartesian(norm_x,norm_y,norm_z);
        norm_x = norm.get_x();
        norm_y = norm.get_y();
        norm_z = norm.get_z();
      } else {
        x = (majorRadius +  minorRadius * cos(theta2)) * cos(phi) + vertices[0].get_x();
        y = (majorRadius +  minorRadius * cos(theta2)) * sin(phi) + vertices[0].get_y();
        z = minorRadius * sin(theta2) + vertices[0].get_z();
        norm_x = cos(theta2) * cos(phi);
        norm_y = cos(theta2) * sin(phi);
      }
      bufferVertices[3*ivert]   = x;
      bufferVertices[3*ivert+1] = y;
      bufferVertices[3*ivert+2] = z;
      bufferNormals[3*ivert]    = norm_x;
      bufferNormals[3*ivert+1]  = norm_y;
      bufferNormals[3*ivert+2]  = norm_z;
      bufferColors[4*ivert]     = colour[0];
      bufferColors[4*ivert+1]   = colour[1];
      bufferColors[4*ivert+2]   = colour[2];
      bufferColors[4*ivert+3]   = alpha;
      ivert++;

      bufferIndices[indx++] = ivert-4;
      bufferIndices[indx++] = ivert-3;
      bufferIndices[indx++] = ivert-2;

      bufferIndices[indx++] = ivert-4;
      bufferIndices[indx++] = ivert-2;
      bufferIndices[indx++] = ivert-1;

    }
  }

  polygonType = GL_TRIANGLES;
  std::vector<Primitive*> a;
  return a;
}

std::vector<Primitive*> TorusElement::GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin, int start, int end) const {
	std::vector<Primitive*> a;
	return a;
}

void TorusElement::DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const Volume &v){}

void TorusElement::DrawPostScript(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, double xoff, double yoff, double xscale, double yscale, double xscaleps, const Volume &v){}
