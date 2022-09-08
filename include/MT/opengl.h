/*  Copyright (C) 2000, 2006  Marc Toussaint (mtoussai@inf.ed.ac.uk)
    under the terms of the GNU LGPL (http://www.gnu.org/copyleft/lesser.html)
    see the `std.h' file for a full copyright statement  */

/*! \file opengl.h
    \brief core file: defines the OpenGL interface to freeglut or Qt */

#ifndef MT_opengl_h
#define MT_opengl_h

#include "std.h"
#include "array.h"
#include "geo3d.h"
#ifdef MT_GL
#  include <GL/gl.h>
#  include <GL/glu.h>
#endif

#ifdef MT_QT
#  if !defined MT_Cygwin 
#    include <GL/glut.h>
#  endif
#  include <qapplication.h>
#  include <qpixmap.h>
#  include <qimage.h>
#  include <qgl.h>
#  include <qtimer.h>
#  if defined MT_Cygwin //|| defined MT_Linux
#    define GLformat DirectRendering | DepthBuffer | Rgba
#    define GLosformat DirectRendering | DepthBuffer | Rgba
#  else
#    define GLformat DirectRendering | DepthBuffer | Rgba
#    define GLosformat DirectRendering | DepthBuffer | Rgba
#  endif
#endif

#ifdef MT_FREEGLUT
#  define FREEGLUT_STATIC
#  include <GL/freeglut.h>
extern "C"{
  void FGAPIENTRY glutMainLoopMT( void );
  void fgDeinitialize( void );
}
#endif

#ifdef MT_GL2PS
#  include<gl2ps.h>
#endif


//===========================================================================
//
// utility functions
//

void glStandardLight();
void glColor(float r, float g, float b,float a=1.f);
void glColor(int col);
void glDrawText(const char* txt,float x,float y,float z);
//void glShadowTransform();
void glTransform (const double pos[3], const double R[12]);
void glDrawRect(float x1,float y1,float z1,float x2,float y2,float z2,
		float x3,float y3,float z3,float x4,float y4,float z4,
		float r,float g,float b);
void glDrawRect(float x1,float y1,float z1,float x2,float y2,float z2,
		float x3,float y3,float z3,float x4,float y4,float z4);
void glDrawRect(float x,float y,float z,float rad);
void glDrawFloor(float x,float r,float g,float b);
void glDrawBox(float x,float y,float z);
void glDrawDimond(float x,float y,float z);
void glDrawSphere(float radius);
void glDrawCylinder(float radius, float length,bool closed=true);
void glDrawCappedCylinder(float radius, float length);
void glDrawAxes(double scale);
void glDrawGridBox(float x);
void glDrawGridBox(float x1,float y1,float z1,float x2,float y2,float z2);
void glDrawKhepera();
void glMakeSquare(int num);
void glMakeStdSimplex(int num);
void glMakeTorus(int num);
void glDrawRobotArm(float a,float b,float c,float d,float e,float f);
uint glImageTexture(const byteA &img);
void glDrawTexQuad(uint texture,
		   float x1,float y1,float z1,float x2,float y2,float z2,
	           float x3,float y3,float z3,float x4,float y4,float z4,
		   float mulX=1.,float mulY=1.);
/*!\brief return the RGBA-image of scenery drawn just before; the image
  buffer has to have either 2 dimensions [width,height] for a
  gray-scale luminance image or 3 dimensions [width,height,4] for an
  RGBA-image. */
void glGrabImage(byteA& img);
/*!\brief return the depth map of the scenery drawn just before; the depth
    buffer has to be a 2-dimensional [width,height] and is filled with
    depth values between 0 and 1. */
void glGrabDepth(byteA& depth);
/*!\brief return the depth map of the scenery drawn just before; the depth
    buffer has to be a 2-dimensional [width,height] and is filled with
    depth values between 0 and 1. */
void glGrabDepth(floatA& depth);
void glRasterImage(int x,int y,byteA &img,float zoom=1.);
void glWatchImage(byteA &img,bool wait=true,float zoom=1.);
void glDisplayGrey(const arr &x,uint d0=0,uint d1=0,bool wait=false,uint win=0);
void glDisplayRedBlue(const arr &x,uint d0=0,uint d1=0,bool wait=false,uint win=0);
void glColor(float r, float g, float b,float a);


//===========================================================================
//
// basic gui routines - wrapped for compatibility
//

void MTprocessEvents();
void MTenterLoop();
void MTexitLoop();


//===========================================================================
//
// OpenGL class
//

#ifndef MT_QT
class QGLWidget{};
#endif


/*!\brief A class to display and control 3D scenes using OpenGL and Qt.

    Minimal use: call \ref add to add routines or objects to be drawn
    and \ref update or \ref watch to start the display. */
class OpenGL:private MT::InitQt,public QGLWidget{
  Q_OBJECT
public:
  // little structs to store objects to be drawn by one OpenGL instance, etc
  struct GLDrawer   { void *classP; void (*call)(void*); };
  struct GLInitCall { void *classP; bool (*call)(void*,OpenGL*); };
  struct GLHoverCall{ void *classP; bool (*call)(void*,OpenGL*); };
  struct GLClickCall{ void *classP; bool (*call)(void*,OpenGL*); };
  struct GLEvent{ int button,key,x,y; float dx,dy; void set(int b,int k,int _x,int _y,float _dx,float _dy){ button=b; key=k; x=_x; y=_y; dx=_dx; dy=_dy; } };
  struct GLSelect{ int name; double dmin,dmax; };

  //global window list
  static MT::Array<OpenGL*> glwins;

  int windowID; //!< id of this window in the global glwins list

  MT::Array<GLDrawer> drawers;         //!< list of draw routines
  MT::Array<GLInitCall> initCalls;     //!< list of initialization routines
  MT::Array<GLHoverCall> hoverCalls;   //!< list of hover callbacks
  MT::Array<GLClickCall> clickCalls;   //!< list of click callbacks

  String text; //!< the text to be drawn as title within the opengl frame
  geo3d::Camera camera; //!< the camera used for projection
  float clearR,clearG,clearB,clearA;  //!< colors of the beackground (called in glClearColor(...))
  bool reportEvents,reportSelects; //!< flags for verbosity
  int pressedkey;          //!< stores the key pressed
  int mouse_button;        //!< stores which button was pressed
  int mouseposx,mouseposy; //!< current x- and y-position of mouse
  MT::Array<GLSelect> selection; //!< list of all selected objects
  GLSelect *topSelection;        //!< top selected object


private:
  GLEvent lastEvent;
  geo3d::Vector downVec,downPos;
  geo3d::Rotation downRot;
  static GLuint selectionBuffer[1000];
#ifdef MT_QT
  bool quitLoopOnTimer;
  QPixmap *osPixmap;      // the paint device for off-screen rendering
  QGLContext *osContext;  //the GL context for off-screen rendering
#endif


private:
  void init(); //initializes camera etc


public://!@name constructors & destructors (seperately for freeglut and Qt)
  //! constructor
  OpenGL(const char* title="MT::OpenGL",int w=400,int h=400,int posx=100,int posy=100);
#ifdef MT_QT
  //! Qt constructor when window is parent of another window
  OpenGL(QWidget *parent,const char* title,int width=400,int height=400,int posx=0,int posy=0);
#endif
  //! destructor
  ~OpenGL();


public://!@name adding routines and callbacks
  //! add a draw routine
  void add(void (*call)(void*),const void* classP=0);
  //! add a class or struct with a staticDraw routine
  template<class T> void add(const T& x){ add(x.staticDraw,&x); }
  //! clear the list of all draw routines
  void clear();
  //! add a hover callback
  void addHoverCall(bool (*call)(void*,OpenGL*),const void* classP=0);
  //! add a click callback
  void addClickCall(bool (*call)(void*,OpenGL*),const void* classP=0);


public://!@name the core draw routines (actually only for internal use)
  /*!\brief basic OpenGL drawing routine -- first arranges the perspective
      with respect to the geo3d::Camera camera - then calls the external
      drawing rountines in the \c drawers list */
  void Draw(geo3d::Camera& c,int w,int h);

  /*!\brief similarly to the Draw() routine, this is called in the
    GL_SELECT mode; however the camera is a narrow window around the
    mouseposx and mouseposy; the selection and topSelection buffers
    contain the outcome of selection */
  void Select();  


private: //public://!@name internal callbacks provided for freeglut or Qt
  //general callbacks
  void Key(unsigned char key, int x, int y);
  void Mouse(int button, int updown, int x, int y);
  void Motion(int x, int y);
  void PassiveMotion(int x, int y);
  void Close(){ }
  void Reshape(int w,int h);

#ifdef MT_FREEGLUT
  //hooks for FREEGLUT (static callbacks)
public:
  static void _Void(){ }
  static void _Draw(){ OpenGL *gl=glwins(glutGetWindow()); gl->Draw(gl->camera,gl->width(),gl->height()); glutSwapBuffers(); }
private:
  static void _Key(unsigned char key, int x, int y){ glwins(glutGetWindow())->Key(key,x,y); }
  static void _Mouse(int button, int updown, int x, int y){ glwins(glutGetWindow())->Mouse(button,updown,x,y); }
  static void _Motion(int x, int y){ glwins(glutGetWindow())->Motion(x,y); }
  static void _PassiveMotion(int x, int y){ glwins(glutGetWindow())->PassiveMotion(x,y); }
  static void _Close(){ glwins(glutGetWindow())->Close(); }
  static void _Reshape(int w,int h){ glwins(glutGetWindow())->Reshape(w,h); }
#endif
#ifdef MT_QT
  //hooks for Qt (overloading virtuals)
  void paintGL(){ Draw(camera,width(),height()); }
  void initializeGL(){ }
  void resizeGL(int w,int h){ Reshape(w,h); }
  void keyPressEvent(QKeyEvent *e){ pressedkey=e->ascii(); Key(pressedkey,mouseposx,mouseposy); }
  void timerEvent(QTimerEvent *e){ if(quitLoopOnTimer) MTexitLoop(); }
  void mouseMoveEvent(QMouseEvent* e){ if(mouse_button<0) PassiveMotion(e->x(),e->y()); else Motion(e->x(),e->y()); }
  void mousePressEvent(QMouseEvent* e){
    if(e->button()==LeftButton) { Mouse(0,0,e->x(),e->y()); }
    if(e->button()==MidButton)  { Mouse(1,0,e->x(),e->y()); }
    if(e->button()==RightButton){ Mouse(2,0,e->x(),e->y()); }
  }
  void mouseReleaseEvent(QMouseEvent* e){
    if(e->button()==LeftButton) { Mouse(0,1,e->x(),e->y()); }
    if(e->button()==MidButton)  { Mouse(1,1,e->x(),e->y()); }
    if(e->button()==RightButton){ Mouse(2,1,e->x(),e->y()); }
  }
#endif


public://!@name showing, updating, and watching
  //! update the view (in Qt: also starts displaying the window)
  bool update();
  /*!\brief watch in interactive mode and wait for an exiting event
    (key pressed or right mouse) */
  int watch();
  //! waits some msecons before updating
  int timedupdate(uint msec);
  //! resize the window
  void resize(int w,int h);
  //! set the four clear colors
  void setClearColors(float r,float g,float b,float a);
  /*!\brief inverse projection: given a 2D+depth coordinates in the
    camera view (e.g. as a result of selection) computes the world 3D
    coordinates */
  void unproject(double &x,double &y,double &z);


#ifdef MT_QT
public://!@name off-screen rendering (only in Qt yet!!)
  /*!\brief creates a off-screen rendering context for future backround
      rendering routines -- the off-screen context cannot be
      resized... */
  void createOffscreen(int width,int height);
  /*!\brief return the RGBA-image of the given perspective; rendering is done
      off-screen (on an internal QPixmap) */
  void offscreenGrab(byteA& image);
  /*!\brief return the RGBA-image of the given perspective; rendering
      is done off-screen (on an internal QPixmap) */
  void offscreenGrab(byteA& image,byteA& depth);
  /*!\brief return only the depth gray-scale map of given perspective;
      rendering is done off-screen (on an internal QPixmap) */
  void offscreenGrabDepth(byteA& depth);
  /*!\brief return only the depth gray-scale map of given perspective;
      rendering is done off-screen (on an internal QPixmap) */
  void offscreenGrabDepth(floatA& depth);


private:
  void setOffscreen(int width,int height);
#endif

//#if !defined MT_QT && defined MT_MSVC
//  HWND__* winId(){ NIY; return 0; }
//#endif


public://!@name info & I/O
  //! print some info on the selection buffer
  void reportSelection();
#if defined MT_FREEGLUT || !defined MT_GL
  //! get width
  int width();
  //! get height
  int height();
#endif
  /*!\brief generates a ps from the current OpenGL display, using gl2ps */
  void saveEPS(const char *filename);
  /*!\brief report on the OpenGL capabilities (the QGLFormat) */
  void about(std::ostream& os=std::cout);
};


//===========================================================================
//
// simple UI
//

#include <MT/list.h>
//#include <MT/png.h>

class glUI{
public:
  int top;
  struct Button{ byteA img1,img2; bool hover; uint x,y,w,h; const char* name; };
  MT::List<Button> buttons;

  glUI(){ top=-1; }

  void addButton(uint x,uint y,const char *name,const char *img1=0,const char *img2=0);
  void glDraw();
  bool checkMouse(int _x,int _y);
};

void glDrawUI(void *p);
bool glHoverUI(void *p,OpenGL *gl);
bool glClickUI(void *p,OpenGL *gl);
//--------- implementation


#ifdef MT_IMPLEMENTATION
#  include "opengl.cpp"
#endif

#endif
