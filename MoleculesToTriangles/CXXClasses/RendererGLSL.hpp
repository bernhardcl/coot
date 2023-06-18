//
//  RendererGLSL.hpp
//  MoleculesToTriangles
//
//  Created by Martin Noble on 2/1/17.
//  Copyright © 2017 MartinNobleSoftware. All rights reserved.
//

#ifndef RendererGLSL_hpp
#define RendererGLSL_hpp

#ifdef WINDOWS_MINGW
typedef unsigned int uint;
#endif

#include <stdio.h>
#include "RendererGL.h"

class RendererGLSL : public RendererGL
{
private:
    int program;
    virtual void vboRenderVCN(VertexColorNormalPrimitive *prim);
    void myglMaterialfv(int faceEnum, int propertyEnum, float *values);
    void myglMaterialf(int faceEnum, int propertyEnum, float value);
    void myglEnableClientState(int clientStateEnum);
    void myglColorPointer(int size, int type, int stride, const void *pointer);
    void myglVertexPointer(int size, int type, int stride, const void *pointer);
    void myglNormalPointer(int type, int stride, const void *pointer);
    std::map<std::string, int>uniforms;
public:
    virtual void renderVertexColorNormalPrimitive(VertexColorNormalPrimitive *prim);
    virtual void init();
    void setProgram(int _program);
    void Use() const;
    void test(){};
    static std::string PointLightFragmentShaderText;
    static std::string PointLightVertexShaderText;
    static std::shared_ptr<Renderer> create();
    std::pair<bool, std::string> readShaderFile(std::string strFilename) const;
    bool loadShaderFile(std::string strFilename, uint iHandle);
    void loadShaders();
};



#endif /* RendererGLSL_hpp */
