#ifndef _MG_HELP_GLOBALS_H_
#define _MG_HELP_GLOBALS_H_
// These are  helper functions which don't depend on an OPENGL context.
#if defined (_WIN32) && not defined (WINDOWS_MINGW)
#define HELP_DLL __declspec(dllimport)
int __stdcall HELP_DLL findprimc(const std::vector<Cartesian> &xyzbox, const std::vector<Cartesian> &primorigin, const Cartesian &origin, const matrix &objrotmat);
int __stdcall HELP_DLL findprimc(double *xyzmpc, double *xyzmmc, double *xyzpmc, double *xyzppc, std::vector<Cartesian> primorigin, Cartesian origin, matrix objrotmat);
int __stdcall HELP_DLL findprimc_main(const Cartesian &xyzmpf, const Cartesian &xyzmpb, const Cartesian &xyzmmf, const Cartesian &xyzmmb, const Cartesian &xyzpmf, const Cartesian &xyzpmb, const Cartesian &xyzppf, const Cartesian &xyzppb, const std::vector<Cartesian> &primorigin,const Cartesian &origin,const matrix &objrotmat);
std::vector<Cartesian> __stdcall HELP_DLL getxyzc(double x,double y);
Volume __stdcall HELP_DLL GetClippingPlanes();
Volume __stdcall HELP_DLL GetFrontAndBackClippingPlanes();
#else
int findprimc(const std::vector<Cartesian> &xyzbox, const std::vector<Cartesian> &primorigin, const Cartesian &origin, const matrix &objrotmat);
int findprimc(double *xyzmpc, double *xyzmmc, double *xyzpmc, double *xyzppc, std::vector<Cartesian> primorigin, Cartesian origin, matrix objrotmat);
int findprimc_main(const Cartesian &xyzmpf, const Cartesian &xyzmpb, const Cartesian &xyzmmf, const Cartesian &xyzmmb, const Cartesian &xyzpmf, const Cartesian &xyzpmb, const Cartesian &xyzppf, const Cartesian &xyzppb, const std::vector<Cartesian> &primorigin,const Cartesian &origin,const matrix &objrotmat);
std::vector<Cartesian> getxyzc(double x,double y);
Volume GetClippingPlanes();
Volume GetFrontAndBackClippingPlanes();
#endif
void setTextureMatrix(void);
void SetMyTextureUnit(GLuint texUnit);
image_info get_pixdata(int trans=0);
void write_pixdata(const char *filename, int width=-1, int height=-1, int trans=0);
const double *GLf2f(const GLfloat *in, int size);
const GLfloat *f2GLf(double *in, int size);
const GLfloat *buildrotmatrix_from_c(matrix a);
const GLfloat *buildrotmatrix(
GLfloat a0, GLfloat a1, GLfloat a2, GLfloat a3, 
GLfloat a4, GLfloat a5, GLfloat a6, GLfloat a7, 
GLfloat a8, GLfloat a9, GLfloat a10, GLfloat a11, 
GLfloat a12, GLfloat a13, GLfloat a14, GLfloat a15);
int CheckIfStereoAvailable(void);
int CheckIfAlphaAvailable(void);
image_info_yuv_t get_yuvdata(int trans=0);
bool isPointInClippingVolume(const Cartesian &p, const Volume &v);
void SetupFBOBlending();
void DrawAxes(float y, int w, int h, const Quat &quat, float viewsize, const std::vector<double> &axes_position,
 unsigned tex_id_x_0, unsigned tex_id_x_1, unsigned tex_id_y_0, unsigned tex_id_y_1, unsigned tex_id_z_0, unsigned tex_id_z_1, 
 unsigned tex_image_x_0_w, unsigned tex_image_x_0_h, unsigned tex_image_y_0_w, unsigned tex_image_y_0_h, unsigned tex_image_z_0_w, unsigned tex_image_z_0_h, int descent
);
void DrawAxesText(float y, int w, int h, const Quat &quat, float viewsize, const std::vector<double> &axes_position,
 unsigned tex_id_x_0, unsigned tex_id_x_1, unsigned tex_id_y_0, unsigned tex_id_y_1, unsigned tex_id_z_0, unsigned tex_id_z_1, 
 unsigned tex_image_x_0_w, unsigned tex_image_x_0_h, unsigned tex_image_y_0_w, unsigned tex_image_y_0_h, unsigned tex_image_z_0_w, unsigned tex_image_z_0_h, int descent
);
void DrawScale(float y, int w, int h, float viewsize, unsigned tex_id_0, unsigned tex_id_1, unsigned tex_image_w, unsigned tex_image_h, int descent);
void DrawScaleText(float y, int w, int h, float viewsize, unsigned tex_id_0, unsigned tex_id_1, unsigned tex_image_w, unsigned tex_image_h, int descent);

#endif // _MG_HELP_GLOBALS_H_
