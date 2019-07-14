
#include "monteCarloPathTracing.h"

void GLInit()
{
    glEnable(GL_DEPTH_TEST);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    // "-near" and "-far" define clipping planes
  //  glOrtho(-256, 256, -256, 256, -10000, 10000);

    gluPerspective(40, 1.0, 1.0, 10000.0);

    glMatrixMode(GL_MODELVIEW);


    glEnable(GL_NORMALIZE);
    glFrontFace(GL_CCW);    // CW: clock wise if we render teapot. Otherwise, use "GL_CCW"
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);
}

