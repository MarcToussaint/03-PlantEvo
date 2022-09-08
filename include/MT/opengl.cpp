/*  Copyright (C) 2000, 2006  Marc Toussaint (mtoussai@inf.ed.ac.uk)
    under the terms of the GNU LGPL (http://www.gnu.org/copyleft/lesser.html)
    see the `std.h' file for a full copyright statement  */

#include "opengl.h"
//#include "target.h"

//===========================================================================
//
// static objects
//

OpenGL *staticgl [10]; //ten pointers to be potentially used as display windows
uint nrWins=0;
GLuint OpenGL::selectionBuffer[1000];
MT::Array<OpenGL*> OpenGL::glwins;



//===========================================================================
//
// basic gui routines - wrapped for compatibility
//

#ifdef MT_FREEGLUT
void MTprocessEvents(){ glutMainLoopEvent(); }
void MTenterLoop(){     glutMainLoopMT(); }
void MTexitLoop(){      glutLeaveMainLoop(); }
#endif
#ifdef MT_QT
void MTprocessEvents(){ qApp->processEvents(); }
void MTenterLoop(){     qApp->enter_loop(); }
void MTexitLoop(){      qApp->exit_loop(); }
#endif
#ifndef MT_GL
void MTprocessEvents(){ }
void MTenterLoop(){     }
void MTexitLoop(){      }
#endif

//===========================================================================
//
// utility implementations
//

void glStandardLight(){
  glEnable(GL_LIGHTING);

  static GLfloat diffuse[]   = { 1.0, 1.0, 1.0, 1.0 };
  static GLfloat specular[]  = { 1.0, 1.0, 1.0, 1.0 };
  static GLfloat position[]  = { 1000.0, -800.0, 1000.0, 1.0 };
  static GLfloat direction[] = { -1.0, .8, -1.0 };
  glLightfv(GL_LIGHT0, GL_POSITION, position);
  glLightfv(GL_LIGHT0, GL_SPOT_DIRECTION, direction);
  glLighti (GL_LIGHT0, GL_SPOT_CUTOFF,   90);
  glLighti (GL_LIGHT0, GL_SPOT_EXPONENT, 10);
  glLightfv(GL_LIGHT0, GL_DIFFUSE,  diffuse);
  glLightfv(GL_LIGHT0, GL_SPECULAR, specular);
  glEnable(GL_LIGHT0);

  static GLfloat diffuse1[]   = { 0.5, 0.5, 0.5, 1.0 };
  static GLfloat specular1[]  = { 0.1, 0.1, 0.1, 1.0 };
  static GLfloat position1[]  = { -100.0, 20.0, -100.0, 1.0 };
  static GLfloat direction1[] = { 1.0, -.2, 1.0 };
  glLightfv(GL_LIGHT1, GL_POSITION, position1);
  glLightfv(GL_LIGHT1, GL_SPOT_DIRECTION, direction1);
  glLighti (GL_LIGHT1, GL_SPOT_CUTOFF,   90);
  glLighti (GL_LIGHT1, GL_SPOT_EXPONENT, 10);
  glLightfv(GL_LIGHT1, GL_DIFFUSE,  diffuse1);
  glLightfv(GL_LIGHT1, GL_SPECULAR, specular1);
  glEnable(GL_LIGHT1);
}

// Henryk's openGL color tab, 4. Mar 06 (hh)
void glColor(int col){
  static const GLfloat colorsTab[6][4] = {{0.2, 0.2, 1.0, 1.0}, // blue
					{1.0, 0.8, 0.0, 1.0}, // gold
					{1.0, 0.0, 0.0, 1.0}, // red
					{0.7, 0.7, 0.7, 1.0}, // gray
					{1.0, 1.0, 1.0, 1.0}, // white
					{0.2, 1.0, 0.2, 1.0}}; // green

  if (col<0) col=0; if (col>5) col=5;
  glColor(colorsTab[col][0],colorsTab[col][1],colorsTab[col][2],colorsTab[col][3]);
}

void glColor(float r, float g, float b, float alpha){
  float amb=1.f,diff=1.f,spec=1.f;
  GLfloat ambient[4],diffuse[4],specular[4];
  ambient[0] = r*amb;
  ambient[1] = g*amb;
  ambient[2] = b*amb;
  ambient[3] = alpha;
  diffuse[0] = r*diff;
  diffuse[1] = g*diff;
  diffuse[2] = b*diff;
  diffuse[3] = alpha;
  specular[0] = spec*.5*(1.+r);
  specular[1] = spec*.5*(1.+g);
  specular[2] = spec*.5*(1.+b);
  specular[3] = alpha;
#if 0
  glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, diffuse);
#else
  glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, ambient);
  glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse);
  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular);
  glMaterialf (GL_FRONT_AND_BACK, GL_SHININESS, 50.0f);
#endif
  glColor4f(r,g,b,alpha);
}

void glColor(float *rgb){ glColor(rgb[0],rgb[1],rgb[2],1.); }

/* // shadows do not work with a light source;
   // thus, we need to leave this out. 4. Mar 06 (hh)
void glShadowTransform()
{
  GLfloat matrix[16];
  for (int i=0; i<16; i++) matrix[i] = 0;
  matrix[0]=1;
  matrix[5]=1;
  matrix[8]=-1;  //light_x
  matrix[9]=-1;  //light_y
  matrix[14]=.02; //ground offset
  matrix[15]=1;
  glPushMatrix();
  glMultMatrixf (matrix);
}
*/

void glTransform(const double pos[3], const double R[12]){
    GLfloat matrix[16];
    matrix[0]=R[0];
    matrix[1]=R[4];
    matrix[2]=R[8];
    matrix[3]=0;
    matrix[4]=R[1];
    matrix[5]=R[5];
    matrix[6]=R[9];
    matrix[7]=0;
    matrix[8]=R[2];
    matrix[9]=R[6];
    matrix[10]=R[10];
    matrix[11]=0;
    matrix[12]=pos[0];
    matrix[13]=pos[1];
    matrix[14]=pos[2];
    matrix[15]=1;
    glPushMatrix();
    glMultMatrixf(matrix);
}

static GLboolean glLightIsOn = false;
void glPushLightOff(){ glGetBooleanv(GL_LIGHTING,&glLightIsOn); glDisable(GL_LIGHTING); }
void glPopLight(){ if(glLightIsOn) glEnable(GL_LIGHTING); }

void glDrawText(const char* txt,float x,float y,float z){
#if defined __glut_h__ || (defined __FREEGLUT_H__ && !defined MT_Cygwin)
  glPushLightOff();
  glRasterPos3f(x,y,z);
  void *font=GLUT_BITMAP_HELVETICA_12;
  while(*txt){ 
    switch(*txt){
    case '\n':
      y+=15;
      glRasterPos3f(x,y,z);
      break;
    case '\b':
      if(font==GLUT_BITMAP_HELVETICA_12) font=GLUT_BITMAP_HELVETICA_18;
      else font=GLUT_BITMAP_HELVETICA_12;
      break;
    default:
      glutBitmapCharacter(font,*txt); 
    }
    txt++;
  }
  glPopLight();
#endif
}

void glDrawRect(float x1,float y1,float z1,float x2,float y2,float z2,
	        float x3,float y3,float z3,float x4,float y4,float z4){
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glBegin(GL_POLYGON);
  glVertex3f(x1,y1,z1);
  glVertex3f(x2,y2,z2);
  glVertex3f(x3,y3,z3);
  glVertex3f(x4,y4,z4);
  glVertex3f(x1,y1,z1);
  glEnd();
}

void glDrawRect(float x1,float y1,float z1,float x2,float y2,float z2,
	        float x3,float y3,float z3,float x4,float y4,float z4,
		float r,float g,float b){
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glBegin(GL_POLYGON);
  glColor(r,g,b);
  glVertex3f(x1,y1,z1);
  glVertex3f(x2,y2,z2);
  glVertex3f(x3,y3,z3);
  glVertex3f(x4,y4,z4);
  glVertex3f(x1,y1,z1);
  glEnd();
}

void glDrawRect(float x,float y,float z,float r){
  glDrawRect(x-r,y-r,z,x-r,y+r,z,x+r,y+r,z,x+r,y-r,z);
}

void glDrawFloor(float x,float r,float g,float b){
  x/=2.;
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glColor(r,g,b);
  glBegin(GL_POLYGON);
  glNormal3f(0,0,1);
  glVertex3f(-x,-x,0.);
  glVertex3f( x,-x,0.);
  glVertex3f( x, x,0.);
  glVertex3f(-x, x,0.);
  glVertex3f(-x,-x,0.);
  glEnd();

  float grey=0.75;
  glColor(grey,grey,grey);
  for(float i=0.0; i<= 1.1;i= i+0.1){
    for(float j=0.0; j <= 1.1;j=j+0.1){
      glBegin(GL_LINES);
      glVertex3f(i,-1,0.001);
      glVertex3f(i,1,0.001);
      glEnd();
      glBegin(GL_LINES);
      glVertex3f(-i,-1,0.001);
      glVertex3f(-i,1,0.001);
      glEnd();
      glBegin(GL_LINES);
      glVertex3f(-1,i,0.001);
      glVertex3f(1,i,0.001);
      glEnd();
      glBegin(GL_LINES);
      glVertex3f(-1,-i,0.001);
      glVertex3f(1,-i,0.001);
      glEnd();
    }
  }

  grey=0.25;
  glColor(grey,grey,grey);
  glBegin(GL_LINE_STRIP);
  glVertex3f(-1,-1,0.002);
  glVertex3f(-1,1,0.002);
  glVertex3f(1,1,0.002);
  glVertex3f(1,-1,0.002);
  glVertex3f(-1,-1,0.002);
  glEnd();
}

// used Henryk's drawBox routine to get the 
// correct shading,  4. Mar 06 (hh)
void glDrawBox(float x,float y,float z)
{
  static GLfloat n[6][3] =
    {
      {-1.0, 0.0, 0.0},
      {0.0, 1.0, 0.0},
      {1.0, 0.0, 0.0},
      {0.0, -1.0, 0.0},
      {0.0, 0.0, 1.0},
      {0.0, 0.0, -1.0}
    };

  static GLint faces[6][4] =
    {
      {0, 1, 2,3},
      {3, 2, 6, 7},
      {7, 6, 5, 4},
      {4, 5, 1, 0},
      {5, 6, 2, 1},
      {7, 4, 0, 3}
    };
  GLfloat v[8][3];
  GLint i;

  v[0][0] = v[1][0] = v[2][0] = v[3][0] = -x / 2;
  v[4][0] = v[5][0] = v[6][0] = v[7][0] =  x / 2;
  v[0][1] = v[1][1] = v[4][1] = v[5][1] =  -y / 2;
  v[2][1] = v[3][1] = v[6][1] = v[7][1] =  y / 2;
  v[0][2] = v[3][2] = v[4][2] = v[7][2] =  -z / 2;
  v[1][2] = v[2][2] = v[5][2] = v[6][2] =  z / 2;

  for (i = 5; i >= 0; i--) {
    glBegin(GL_QUADS);
    glNormal3fv(&n[i][0]);
    glVertex3fv(&v[faces[i][0]][0]);
    glVertex3fv(&v[faces[i][1]][0]);
    glVertex3fv(&v[faces[i][2]][0]);
    glVertex3fv(&v[faces[i][3]][0]);
    glEnd();
  }
}

void glDrawDimond(float x,float y,float z){
  glBegin(GL_TRIANGLE_FAN);
  glVertex3f(.0,.0, z);
  glVertex3f( x,.0,.0);
  glVertex3f(.0, y,.0);
  glVertex3f(-x,.0,.0);
  glVertex3f(.0,-y,.0);
  glVertex3f( x,.0,.0);
  glEnd();
  glBegin(GL_TRIANGLE_FAN);
  glVertex3f(.0,.0,-z);
  glVertex3f( x,.0,.0);
  glVertex3f(.0, y,.0);
  glVertex3f(-x,.0,.0);
  glVertex3f(.0,-y,.0);
  glVertex3f( x,.0,.0);
  glEnd();
}

void glDrawAxes(double scale){
  GLUquadric *style=gluNewQuadric();

  for(uint i=0;i<3;i++){
    glPushMatrix();
    glScalef(scale,scale,scale);
    switch(i){
    case 0:  glColor(1,0,0);  break;
    case 1:  glColor(0,1,0);  glRotatef(90,0,0,1);  break;
    case 2:  glColor(0,0,1);  glRotatef(90,0,-1,0);  break;
    }
    glBegin(GL_LINES);
    glVertex3f(0,0,0);
    glVertex3f(.95,0,0);
    glEnd();
    glTranslatef(.9,0,0);
    glRotatef(90,0,1,0);
    gluCylinder(style,.1,0,.3,20,1); 
    glPopMatrix();
  }
  
  gluDeleteQuadric(style);
}

void drawCoordinateFrame(void* motionObj)
{
	//motion *m = (motion *)motionObj;

	// x-axis = green
	glBegin(GL_LINE_STRIP);
		glColor(4.0,0.0,0.0);
		glVertex3f(0.0,0.0,0.001);
		glVertex3f(1.0,0.0,0.001);
		glVertex3f(0.0,0.0,0.002);
		glVertex3f(1.0,0.0,0.002);
		glVertex3f(0.0,0.0,0.003);
		glVertex3f(1.0,0.0,0.003);
		glVertex3f(0.0,0.0,0.004);
		glVertex3f(1.0,0.0,0.004);
    glEnd();

	// y-axis = green
	glBegin(GL_LINE_STRIP);
		glColor(0.0,4.0,0.0);
		glVertex3f(0.0,0.0,0.001);
		glVertex3f(0.0,1.0,0.001);
		glVertex3f(0.0,0.0,0.002);
		glVertex3f(0.0,1.0,0.002);
		glVertex3f(0.0,0.0,0.003);
		glVertex3f(0.0,1.0,0.003);
		glVertex3f(0.0,0.0,0.004);
		glVertex3f(0.0,1.0,0.004);
    glEnd();

	// z-axis = blue
	glBegin(GL_LINE_STRIP);
		glColor(0.0,0.0,4.0);
		glVertex3f(0.001,0.0,0.0);
		glVertex3f(0.001,0.0,1.0);
		glVertex3f(0.002,0.0,0.0);
		glVertex3f(0.002,0.0,1.0);
		glVertex3f(0.003,0.0,0.0);
		glVertex3f(0.003,0.0,1.0);
		glVertex3f(0.004,0.0,0.0);
		glVertex3f(0.004,0.0,1.0);
    glEnd();

}

#if 0
void drawTarget(void* tar){
	geo3d::Vector color;
	geo3d::Vector position;
	float radius;
	int  stacks, slices;
	
	target *_target = (target*)tar;

	color = _target->getRGB();
	position = _target->getPosition();
	radius = _target->getRadius();
	stacks = _target->getStacks();
	slices = _target->getSlices();

	glPushMatrix();
	glColor(color[0],color[1],color[2]);
	
	GLUquadric *style=gluNewQuadric();
	glTranslatef(position[0],position[1],position[2]);
	gluSphere(style,radius,slices,stacks); // last two value for detail
	gluDeleteQuadric(style);
	glPopMatrix();


	
	glDrawGridBox(0.0,0.0,0.0,position[0],position[1],position[2]);
}
#endif

void glDrawSphere(float radius){
  GLUquadric *style=gluNewQuadric();
  gluSphere(style,radius,10,10); // last two value for detail
  gluDeleteQuadric(style);
}

// add cylinder, 6. Mar 06 (hh)
void glDrawCylinder(float radius, float length,bool closed){
  GLUquadric *style=gluNewQuadric();
  glTranslatef(0,0,-length/2);
  gluCylinder(style,radius,radius,length,20,1); 
  if(closed){
    gluDisk(style,0,radius,20,1); 
    glTranslatef(0,0,length);
    gluDisk(style,0,radius,20,1); 
    glTranslatef(0,0,-length/2);
  }
  gluDeleteQuadric(style);
}


// add capped cylinder, 6. Mar 06 (hh)
void glDrawCappedCylinder(float radius, float length)
{
  GLUquadric *style1=gluNewQuadric();
  GLUquadric *style2=gluNewQuadric();
  GLUquadric *style3=gluNewQuadric();

  glTranslatef(0,0,-length/2);
  gluCylinder(style1,radius,radius,length,20,1);
  // add two spheres to the cylinder:
  glTranslatef(0,0,length);
  gluSphere(style2,radius,10,10);
  glTranslatef(0,0,-length);
  gluSphere(style3,radius,10,10);

  gluDeleteQuadric(style1);
  gluDeleteQuadric(style2);
  gluDeleteQuadric(style3);
}

void glDrawGridBox(float x=10.){
  x/=2.;
  glBegin(GL_LINE_LOOP);
  glVertex3f(-x,-x,-x);
  glVertex3f(-x,-x, x);
  glVertex3f(-x, x, x);
  glVertex3f(-x, x,-x);
  glEnd();
  glBegin(GL_LINE_LOOP);
  glVertex3f( x,-x,-x);
  glVertex3f( x,-x, x);
  glVertex3f( x, x, x);
  glVertex3f( x, x,-x);
  glEnd();
  glBegin(GL_LINES);
  glVertex3f( x, x, x);
  glVertex3f(-x, x, x);
  glVertex3f( x,-x, x);
  glVertex3f(-x,-x, x);
  glVertex3f( x, x,-x);
  glVertex3f(-x, x,-x);
  glVertex3f( x,-x,-x);
  glVertex3f(-x,-x,-x);
  glEnd();
}

void glDrawGridBox(float x1,float y1,float z1,float x2,float y2,float z2){
  //glDrawFloor(x);
  glBegin(GL_LINE_STRIP);
  glVertex3f(x1,y1,z1+0.001);
  glVertex3f(x2,y1,z1+0.001);
  glVertex3f(x2,y2,z1+0.001);
  glVertex3f(x1,y2,z1+0.001);
  glVertex3f(x1,y1,z1+0.001);
  glEnd();
 
  glBegin(GL_LINES);
   glVertex3f( x2, y2,z1 +0.001);
  glVertex3f( x2, y2,z2 +0.001);
  glEnd();
}

void glDrawKhepera(){
  GLUquadric *style=gluNewQuadric();
  glPushMatrix();
  glRotatef(-90,1,0,0);
  glColor3f (.3,.3,.3);
  gluCylinder(style,1.5,1.5,2,20,1);
  glPopMatrix();

  glColor3f (1,0,0);
  glBegin(GL_LINES);
  glVertex3f(0,2,0);
  glVertex3f(0,2,-2.5);
  glEnd();
  gluDeleteQuadric(style);
}

void glMakeSquare(int num){
  glNewList(num,GL_COMPILE);
  glColor3f(1.,0.,0.);
  glBegin( GL_LINE_LOOP );
  glVertex3f( -1 , -1 , 0. );
  glVertex3f( -1 , +1 , 0. );
  glVertex3f( +1 , +1 , 0. );
  glVertex3f( +1 , -1 , 0. );
  glEnd();
  glEndList();
}

void glMakeStdSimplex(int num){
  glNewList(num,GL_COMPILE);
  //glPolygonMode(GL_BACK,GL_FILL);
  glShadeModel(GL_SMOOTH);
  glBegin(GL_TRIANGLE_FAN);
  glColor3f (1.,1.,1.);
  glVertex3f(0.,0.,0.);
  glColor3f (1.,0.,0.);
  glVertex3f(1.,0.,0.);
  glColor3f (0.,1.,0.);
  glVertex3f(0.,1.,0.);
  glColor3f (0.,0.,1.);
  glVertex3f(0.,0.,1.);
  glColor3f (1.,0.,0.);
  glVertex3f(1.,0.,0.);
  glEnd();
  /*
    glColor4f(.5,.5,.5,.9);
    glBegin(GL_POLYGON);
    glVertex3f( 1. , 0. , 0. );
    glVertex3f( 0. , 1. , 0. );
    glVertex3f( 0. , 0. , 1. );
    glEnd();
  */
  glEndList();
}

void glMakeTorus(int num){
  glNewList(num,GL_COMPILE);

  GLint i, j, rings, sides;
  float theta1, phi1, theta2, phi2;
  float v0[03], v1[3], v2[3], v3[3];
  float t0[03], t1[3], t2[3], t3[3];
  float n0[3], n1[3], n2[3], n3[3];
  float innerRadius=0.4;
  float outerRadius=0.8;
  float scalFac;

  rings = 8;
  sides = 10;
  scalFac=1/(outerRadius*2);

  for(i=0;i<rings;i++){
    theta1 = (float)i * 2.0 * MT_PI / rings;
    theta2 = (float)(i + 1) * 2.0 * MT_PI / rings;
    for(j=0;j<sides;j++){
      phi1 = (float)j * 2.0 * MT_PI / sides;
      phi2 = (float)(j + 1) * 2.0 * MT_PI / sides;
	
      v0[0] = cos(theta1) * (outerRadius + innerRadius * cos(phi1));
      v0[1] =-sin(theta1) * (outerRadius + innerRadius * cos(phi1));
      v0[2] = innerRadius * sin(phi1);

      v1[0] = cos(theta2) * (outerRadius + innerRadius * cos(phi1));
      v1[1] =-sin(theta2) * (outerRadius + innerRadius * cos(phi1));
      v1[2] = innerRadius * sin(phi1);

      v2[0] = cos(theta2) * (outerRadius + innerRadius * cos(phi2));
      v2[1] =-sin(theta2) * (outerRadius + innerRadius * cos(phi2));
      v2[2] = innerRadius * sin(phi2);

      v3[0] = cos(theta1) * (outerRadius + innerRadius * cos(phi2));
      v3[1] =-sin(theta1) * (outerRadius + innerRadius * cos(phi2));
      v3[2] = innerRadius * sin(phi2);

      n0[0] = cos(theta1) * (cos(phi1));
      n0[1] =-sin(theta1) * (cos(phi1));
      n0[2] = sin(phi1);

      n1[0] = cos(theta2) * (cos(phi1));
      n1[1] =-sin(theta2) * (cos(phi1));
      n1[2] = sin(phi1);

      n2[0] = cos(theta2) * (cos(phi2));
      n2[1] =-sin(theta2) * (cos(phi2));
      n2[2] = sin(phi2);

      n3[0] = cos(theta1) * (cos(phi2));
      n3[1] =-sin(theta1) * (cos(phi2));
      n3[2] = sin(phi2);

      t0[0] = v0[0]*scalFac + 0.5;
      t0[1] = v0[1]*scalFac + 0.5;
      t0[2] = v0[2]*scalFac + 0.5;

      t1[0] = v1[0]*scalFac + 0.5;
      t1[1] = v1[1]*scalFac + 0.5;
      t1[2] = v1[2]*scalFac + 0.5;

      t2[0] = v2[0]*scalFac + 0.5;
      t2[1] = v2[1]*scalFac + 0.5;
      t2[2] = v2[2]*scalFac + 0.5;

      t3[0] = v3[0]*scalFac + 0.5;
      t3[1] = v3[1]*scalFac + 0.5;
      t3[2] = v3[2]*scalFac + 0.5;

      if((i+j)%2) glColor3f(0.,1.,0.);
      else glColor3f(0.,0.,1.);

      glBegin(GL_POLYGON);
      glNormal3fv(n3); glTexCoord3fv(t3); glVertex3fv(v3);
      glNormal3fv(n2); glTexCoord3fv(t2); glVertex3fv(v2);
      glNormal3fv(n1); glTexCoord3fv(t1); glVertex3fv(v1);
      glNormal3fv(n0); glTexCoord3fv(t0); glVertex3fv(v0);
      glEnd();
    }
  }
  glEndList();
}

void glDrawRobotArm(MT::Array<double>& x){
  CHECK(x.N==6,"wrong DOFs for robot arm draw");

  MT::constrain(x(0),0,90);
  MT::constrain(x(1),-30,90);
  MT::constrain(x(2),-30,90);
  MT::constrain(x(3),0,160);
  MT::constrain(x(4),0,160);
  MT::constrain(x(5),-90,90);

  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glColor3f(.3,.3,.3);
  glPushMatrix();
  glTranslatef(0,0,10);
  glDrawBox(10,10,20);

  glColor3f(.8,.0,.0);
  glTranslatef(5,0,5);
  glRotatef(x(0),0,1,0);
  glRotatef(x(1),0,0,1);
  glRotatef(x(2),1,0,0);
  glTranslatef(5,0,0);
  glDrawBox(10,5,5);

  glColor3f(.0,.0,.8);
  glTranslatef(5,0,0);
  glRotatef(x(3),0,0,1);
  glTranslatef(5,0,0);
  glDrawBox(10,5,5);

  glColor3f(.8,.8,.0);
  glTranslatef(5,0,0);
  glRotatef(x(4),1,0,0);
  glRotatef(x(5),0,1,0);
  glTranslatef(2,0,0);
  glDrawBox(4,5,2);

  glPopMatrix();
}

#if 0
// dm 08.06.2006
void glDrawTrimesh(geo3d::TriMesh mesh) 
{

		glBegin(GL_TRIANGLES);
	glShadeModel(GL_SMOOTH);							
	for(int i=1;i<mesh.T.d0;i++){
			
		  if(mesh.trinormals) glNormal3dv(mesh.Tn(i));
		  if(!mesh.trinormals) glNormal3dv(mesh.Tn(i));//(mesh.Vn(mesh.T(i).nindices[0]));
		////  if(mesh.colored) glColor3fv(mesh.C(T(i,0)));
		  glVertex3dv(mesh.V(mesh.T(i).vindices[0]));

		  if(mesh.trinormals) glNormal3dv(mesh.Tn(i));
		  if(!mesh.trinormals) glNormal3dv(mesh.Tn(i));//glNormal3dv(mesh.Vn(mesh.T(i).nindices[1]));
		////  if(mesh.colored) glColor3fv(mesh.C(T(i,1)));
		  glVertex3dv(mesh.V(mesh.T(i).vindices[1]));

		  if(mesh.trinormals) glNormal3dv(mesh.Tn(i));
		  if(!mesh.trinormals) glNormal3dv(mesh.Tn(i));//glNormal3dv(mesh.Vn(mesh.T(i).nindices[2]));
		////  if(mesh.colored) glColor3fv(mesh.C(T(i,0)));
		  glVertex3dv(mesh.V(mesh.T(i).vindices[2]));
		  
		
	}
	glEnd();

}
#endif

uint glImageTexture(const byteA &img){
  GLuint texName;

  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
  glGenTextures(1, &texName);

  glBindTexture(GL_TEXTURE_2D, texName);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  switch(img.d2){
  case 0:
  case 1:
    glTexImage2D(GL_TEXTURE_2D, 0, 4, img.d1, img.d0, 0, GL_LUMINANCE, GL_UNSIGNED_BYTE, img.p);
    break;
  case 2:
    glTexImage2D(GL_TEXTURE_2D, 0, 4, img.d1, img.d0, 0, GL_LUMINANCE_ALPHA, GL_UNSIGNED_BYTE, img.p);
    break;
  case 3:
    glTexImage2D(GL_TEXTURE_2D, 0, 4, img.d1, img.d0, 0, GL_RGB, GL_UNSIGNED_BYTE, img.p);
    break;
  case 4:
    glTexImage2D(GL_TEXTURE_2D, 0, 4, img.d1, img.d0, 0, GL_RGBA, GL_UNSIGNED_BYTE, img.p);
    break;
  default:
    HALT("no image fomat");
  }
  return texName;
}

void glDrawTexQuad(uint texture,
		   float x1,float y1,float z1,float x2,float y2,float z2,
	           float x3,float y3,float z3,float x4,float y4,float z4,
		   float mulX,float mulY){

  glEnable(GL_TEXTURE_2D);
  glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
  glBindTexture(GL_TEXTURE_2D, texture);
  glBegin(GL_QUADS);
  glTexCoord2f(0.0,  0.0);  glVertex3f(x1,y1,z1);
  glTexCoord2f(0.0,  mulY); glVertex3f(x2,y2,z2);
  glTexCoord2f(mulX, mulY); glVertex3f(x3,y3,z3);
  glTexCoord2f(mulX, 0.0);  glVertex3f(x4,y4,z4);
  glEnd();
  glFlush();
  glDisable(GL_TEXTURE_2D);
}

void glGrabImage(byteA& image){
  CHECK(image.nd==2 ||image.nd==3,"not an image format");
  GLint w=image.d1,h=image.d0;

  //glPixelStorei(GL_PACK_SWAP_BYTES,0); 
  switch(image.d2){
  case 0:
  case 1:
    glPixelTransferf(GL_RED_SCALE,.3333);
    glPixelTransferf(GL_GREEN_SCALE,.3333);
    glPixelTransferf(GL_BLUE_SCALE,.3333);
    glReadPixels(0,0,w,h,GL_LUMINANCE,GL_UNSIGNED_BYTE,image.p);
    glPixelTransferf(GL_RED_SCALE,1.);
    glPixelTransferf(GL_GREEN_SCALE,1.);
    glPixelTransferf(GL_BLUE_SCALE,1.);
    break;
  case 2:
    //glReadPixels(0,0,w,h,GL_GA,GL_UNSIGNED_BYTE,image.p);
    break;
  case 3:
    glReadPixels(0,0,w,h,GL_RGB,GL_UNSIGNED_BYTE,image.p);
    break;
  case 4:
#if defined MT_SunOS
    glReadPixels(0,0,w,h,GL_ABGR_EXT,GL_UNSIGNED_BYTE,image.p);
#else
#if defined MT_Cygwin
    glReadPixels(0,0,w,h,GL_RGBA,GL_UNSIGNED_BYTE,image.p);
#else
    glReadPixels(0,0,w,h,GL_BGRA_EXT,GL_UNSIGNED_BYTE,image.p);
#endif
#endif
    break;
  default: HALT("wrong image format");
  }
}

void glGrabDepth(byteA& depth){
  CHECK(depth.nd==2,"depth buffer has to be either 2-dimensional");
  GLint w=depth.d1,h=depth.d0;
  glReadPixels(0,0,w,h,GL_DEPTH_COMPONENT,GL_UNSIGNED_BYTE,depth.p);
}

void glGrabDepth(floatA& depth){
  CHECK(depth.nd==2,"depth buffer has to be either 2-dimensional");
  GLint w=depth.d1,h=depth.d0;
  glReadPixels(0,0,w,h,GL_DEPTH_COMPONENT,GL_FLOAT,depth.p);
}

void glRasterImage(int x,int y,byteA &img,float zoom){
  glRasterPos2i(x,y); //(int)(y+zoom*img.d0)); (the latter was necessary for other pixel/raster coordinates)
  glPixelZoom(zoom,-zoom);
  if(img.d2<2 && (img.d1%4)){
    img.insColumns(img.d1,4-img.d1%4);
  }
  switch(img.d2){
  case 0:
  case 1:  glDrawPixels(img.d1,img.d0,GL_LUMINANCE,GL_UNSIGNED_BYTE,img.p);        break;
  case 2:  glDrawPixels(img.d1,img.d0,GL_LUMINANCE_ALPHA,GL_UNSIGNED_BYTE,img.p);  break;
  case 3:  glDrawPixels(img.d1,img.d0,GL_RGB,GL_UNSIGNED_BYTE,img.p);              break;
  case 4:  glDrawPixels(img.d1,img.d0,GL_RGBA,GL_UNSIGNED_BYTE,img.p);             break;
  default: HALT("no image format");
  };
}

void glWatchImage(byteA &img,bool wait,float zoom){
#ifdef MT_FREEGLUT
  int w=(int)(zoom*img.d1+10),h=(int)(zoom*img.d0+10);
  glutReshapeWindow(w,h);
  glutDisplayFunc(OpenGL::_Void);
  glutPostRedisplay();
  MTprocessEvents();
  glViewport(0,0,w,h);
  glDisable(GL_DEPTH_TEST);
  glClearColor(.5,.9,.9,0.);
  glClear(GL_COLOR_BUFFER_BIT);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glOrtho(-5,5+zoom*img.d1,5+zoom*img.d0,-5,-1.,1.);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glRasterImage(0,0,img,zoom);
  glutSwapBuffers();
  if(wait) MTenterLoop();
  glutDisplayFunc(OpenGL::_Draw);
#else
  NIY
#endif
}

void glDrawUI(void *p){
  glPushName(0x10);
  ((glUI*)p)->glDraw();
  glPopName();
}

bool glHoverUI(void *p,OpenGL *gl){
  //bool b=
  ((glUI*)p)->checkMouse(gl->mouseposx,gl->mouseposy);
  //if(b) glutPostRedisplay();
  return true;
}

bool glClickUI(void *p,OpenGL *gl){
  bool b=((glUI*)p)->checkMouse(gl->mouseposx,gl->mouseposy);
  if(b) gl->update();
  int t=((glUI*)p)->top;
  if(t!=-1){
    cout <<"CLICK! on button #" <<t <<endl;
    MTexitLoop(); //glutLeaveMainLoop();
    return false;
  }
  return true;
}

void glSelectWin(uint win){
  if(!staticgl[win]) staticgl[win]=new OpenGL;
  glutSetWindow(staticgl[win]->windowID);
}

void glWatchImage(const floatA &x,bool wait,float zoom){
  double ma=x.max();
  double mi=x.min();
  if(wait) cout <<"watched image min/max = " <<mi <<' ' <<ma <<endl;
  byteA img;
  img.resize(x.d0*x.d1);
  img.setZero();
  for(uint i=0;i<x.N;i++){
    img(i)=(byte)(255.*(x.elem(i)-mi)/(ma-mi));
  }
  img.reshape(x.d0,x.d1);
  glWatchImage(img,wait,20);
}

void glDisplayGrey(const arr &x,uint d0,uint d1,bool wait,uint win){
  if(!staticgl[win]) staticgl[win]=new OpenGL;
  if(!d0) d0=x.d0;
  if(!d1) d1=x.d1;
  glutSetWindow(staticgl[win]->windowID);
  double ma=x.max();
  staticgl[win]->text.clr() <<"display"<<win<<" max="<<ma<<endl;
  byteA img;
  img.resize(d0*d1);
  img.setZero();
  for(uint i=0;i<x.N;i++){
    if(x.elem(i)>0.) img(i)=(byte)(255.*x.elem(i)/ma);
    if(x.elem(i)<0.) MT_MSG("warning: negative entry");
  }
  img.reshape(d0,d1);
  glWatchImage(img,wait,20);
}

void glDisplayRedBlue(const arr &x,uint d0,uint d1,bool wait,uint win){
  if(!staticgl[win]) staticgl[win]=new OpenGL;
  if(!d0) d0=x.d0;
  if(!d1) d1=x.d1;
  glutSetWindow(staticgl[win]->windowID);
  double mi=x.min(),ma=x.max();
  staticgl[win]->text.clr() <<"display"<<win<<" max="<<ma<<"min="<<mi<<endl;
  byteA img;
  img.resize(d0*d1,4);
  img.setZero();
  for(uint i=0;i<x.N;i++){
    if(x.elem(i)>0.) img(i,0)=(byte)(255.*x.elem(i)/ma);
    if(x.elem(i)<0.) img(i,2)=(byte)(255.*x.elem(i)/mi);
  }
  img.reshape(d0,d1,4);
  glWatchImage(img,wait,20);
}

//===========================================================================
//
// OpenGL implementations
//

#ifdef MT_FREEGLUT
OpenGL::OpenGL(const char* title,int w,int h,int posx,int posy){
  init();

  if(!nrWins){
    int argc=1;
    char *argv[1]={"x"};
    glutInit(&argc, argv);
    glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE,GLUT_ACTION_CONTINUE_EXECUTION);
  }
  nrWins++;

  glutInitWindowSize(w,h);
  glutInitWindowPosition(posx,posy);
  glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);

  windowID = glutCreateWindow(title);

  //OpenGL initialization
  //two optional thins:
  glEnable(GL_DEPTH_TEST);
  //glEnable(GL_CULL_FACE); glFrontFace(GL_CW);
  glDepthFunc(GL_LESS);
  glShadeModel(GL_SMOOTH);

  if(glwins.N<(uint)windowID+1) glwins.resizeCopy(windowID+1);
  glwins(windowID) = this;

  glutDisplayFunc( _Draw );
  glutKeyboardFunc( _Key );
  glutMouseFunc ( _Mouse ) ;
  glutMotionFunc ( _Motion ) ;
  glutPassiveMotionFunc ( _PassiveMotion ) ;
  glutCloseFunc ( _Close ) ;
  glutReshapeFunc( _Reshape );
  //glutMouseFunc (MouseButton);
  //glutMotionFunc (MouseMotion);

  //  glutSpecialFunc( Special );
  //  glutVisibilityFunc( Visibility );
  //  glutKeyboardUpFunc( KeyUp );
  //  glutSpecialUpFunc( SpecialUp );
  //  glutJoystickFunc( Joystick, 100 );
  //  glutMouseWheelFunc ( MouseWheel ) ;
  //  glutEntryFunc ( Entry ) ;
}
// freeglut destructor
OpenGL::~OpenGL(){
  nrWins--;
  if(!nrWins){
    fgDeinitialize();
  }
}
#endif
#ifdef MT_QT
OpenGL::OpenGL(const char* title,int width,int height,int posx,int posy)
  :QGLWidget(QGLFormat(GLformat),0,title){
  QGLWidget::move(posx,posy);
  QGLWidget::resize(width,height);
  QWidget::setMouseTracking(true);
  QWidget::setCaption(title);
  init();
  windowID=(int)winId();
}
OpenGL::OpenGL(QWidget *parent,const char* title,int width,int height,int posx,int posy)
  :QGLWidget(QGLFormat(GLformat),parent,title){
  QGLWidget::move(posx,posy);
  QGLWidget::resize(width,height);
  QWidget::setMouseTracking(true);
  QWidget::setCaption(title);
  init();
  windowID=(int)winId();
}
OpenGL::~OpenGL(){
  if(osContext) delete osContext;
  if(osPixmap) delete osPixmap;
};
#endif
#ifndef MT_GL
OpenGL::OpenGL(const char* title,int w,int h,int posx,int posy){
  MT_MSG("Warning: creating dummy OpenGL without GL support");
}
#endif

void OpenGL::init(){
  camera.setPosition(0.,0.,10.);
  camera.focus(0,0,0);
  camera.setZRange(1.,1000.);
  camera.setHeightAngle(12.);

  clearR=clearG=clearB=.8; clearA=0.;
  drawers.memMove=true;
  mouseposx=mouseposy=0;
  mouse_button=-1;

  reportEvents=false;
  reportSelects=false;

#ifdef MT_QT
  osPixmap=0;
  osContext=0;
  quitLoopOnTimer=reportEvents=false;
#endif
};


void OpenGL::add(void (*call)(void*),const void* classP){
  CHECK(call!=0,"OpenGL: NULL pointer to drawing routine");
  GLDrawer d; d.classP=(void*)classP; d.call=call;
  drawers.append(d);
}
void OpenGL::clear(){ drawers.resize(0); }
void OpenGL::addHoverCall(bool (*call)(void*,OpenGL*),const void* classP){
  CHECK(call!=0,"OpenGL: NULL pointer to drawing routine");
  GLHoverCall c; c.classP=(void*)classP; c.call=call;
  hoverCalls.append(c);
}
void OpenGL::addClickCall(bool (*call)(void*,OpenGL*),const void* classP){
  CHECK(call!=0,"OpenGL: NULL pointer to drawing routine");
  GLClickCall c; c.classP=(void*)classP; c.call=call;
  clickCalls.append(c);
}


#ifdef MT_GL
void OpenGL::Draw(geo3d::Camera& c,int w,int h){
  //clear bufferer
  glViewport(0,0,w,h);
  glClearColor(clearR,clearG,clearB,clearA);
  if(drawers.N)
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  //OpenGL initialization
  //two optional thins:
  glEnable(GL_DEPTH_TEST);
  //glEnable(GL_CULL_FACE); glFrontFace(GL_CW);
  glDepthFunc(GL_LESS);
  glShadeModel(GL_SMOOTH);
  //glEnable(GL_DEPTH_TEST);
  //glEnable(GL_CULL_FACE);

  //projection
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  c.glSetProjectionMatrix();

  //draw objects
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  if(!drawers.N){ MT_MSG("nothing to be drawn"); if(!text.N()) text <<"<nothing to draw>"; }
  for(uint i=0;i<drawers.N;i++) (*drawers(i).call)(drawers(i).classP);

  //check matrix stack
  GLint s;
  glGetIntegerv(GL_MODELVIEW_STACK_DEPTH,&s);
  CHECK(s==1,"OpenGL matrix stack has not depth 1 (pushs>pops)");

  //draw text
  glGetIntegerv(GL_RENDER_MODE,&s);
  if(text.N() && s!=GL_SELECT){
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glOrtho(0.,(double)w,(double)h,.0,-1.,1.);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    if(clearR+clearG+clearB>1.){
      glColor(0.0,0.0,0.0,1.0); // clear text  5. Mar 06 (hh)
    }else{
      glColor(1.0,1.0,1.0,1.0); // white colored text,  5. Mar 06 (hh)
    }
    glDrawText(text,10,20,0);
  }
}
void OpenGL::Select(){
  uint i,k;
  int j;

  glSelectBuffer(1000,selectionBuffer);
  glRenderMode(GL_SELECT);

  //projection
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  GLint viewport[4];
  viewport[0]=0; viewport[1]=0; viewport[2]=width(); viewport[3]=height();
  gluPickMatrix((GLdouble) mouseposx, (GLdouble) (height()-mouseposy), 2., 2., viewport);
  camera.glSetProjectionMatrix();

  //draw objects
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glInitNames();
  glPushName(0);
  for(i=0;i<drawers.N;i++){
    glLoadName(i);
    (*drawers(i).call)(drawers(i).classP);
    glGetIntegerv(GL_NAME_STACK_DEPTH,&j);
    CHECK(j==1,"OpenGL name stack has not depth 1 (pushs>pops)");
  }

  GLint n;
  n=glRenderMode(GL_RENDER);
  selection.resize(n);

  GLuint *obj,maxD=(GLuint)(-1);
  topSelection=0;
  for(j=0,i=0;i<(uint)n;i++){
    obj=selectionBuffer+j;
    j+=3+obj[0];

    //get name as superposition of all names
    selection(i).name = 0;
    for(k=0;k<obj[0];k++) selection(i).name |= obj[3+k];

    //get dim and dmax
    selection(i).dmin=(double)obj[1]/maxD;  camera.glConvertToTrueDepth(selection(i).dmin);
    selection(i).dmax=(double)obj[2]/maxD;  camera.glConvertToTrueDepth(selection(i).dmax);

    //get top-most selection
    if(!topSelection || selection(i).dmin < topSelection->dmin) topSelection = &selection(i);
  }
  if(reportSelects) reportSelection();
}

#else //MT_GL
void OpenGL::Draw(geo3d::Camera& c,int w,int h){}
void OpenGL::Select(){}
#endif

bool OpenGL::update(){
  pressedkey=0;
#ifdef MT_FREEGLUT
  glutSetWindow(windowID);
  glutPostRedisplay () ;
#endif
#ifdef MT_QT
  show();
  QGLWidget::update();
#endif
  MTprocessEvents();
  return !pressedkey;
}
int OpenGL::watch(){
  pressedkey=0;
  update();
  MTenterLoop();
  MTprocessEvents();
  return pressedkey;
}
int OpenGL::timedupdate(uint msec){
#ifdef MT_FREEGLUT
  int r=update();
  MT::wait(.001*msec);
  return r;
#endif
#ifdef MT_QT
  int i;
  quitLoopOnTimer=true;
  i=startTimer(msec);
  MTenterLoop();
  killTimer(i);
  return update();
#endif
}
void OpenGL::resize(int w,int h){
#ifdef MT_FREEGLUT
  glutSetWindow(windowID);
  glutReshapeWindow(w,h);
#elif defined MT_QT
  QGLWidget::resize(w,h);
#endif
  MTprocessEvents();
}
void OpenGL::setClearColors(float r,float g,float b,float a){
  clearR=r; clearG=g; clearB=b; clearA=a;
}
void OpenGL::unproject(double &x,double &y,double &z){
  double _x,_y,_z;
  GLdouble modelMatrix[16],projMatrix[16];
  GLint viewPort[4];
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  camera.glSetProjectionMatrix();
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glGetDoublev(GL_MODELVIEW_MATRIX,modelMatrix);
  glGetDoublev(GL_PROJECTION_MATRIX,projMatrix);
  glGetIntegerv(GL_VIEWPORT,viewPort);
  gluUnProject(x,y,z, modelMatrix, projMatrix, viewPort, &_x, &_y, &_z);
  x=_x; y=_y; z=_z;
}

#ifdef MT_QT
void OpenGL::createOffscreen(int width,int height){
  if(osContext && (width>osPixmap->width() || height>osPixmap->height())){
    delete osContext;
    delete osPixmap;
    osContext=0;
  }
  if(!osContext){
    osPixmap=new QPixmap(width,height);
    if(!osPixmap) MT_MSG("can't create off-screen Pixmap");
    osContext=new QGLContext(QGLFormat(GLosformat),osPixmap);
    if(!osContext->create()) MT_MSG("can't create off-screen OpenGL context");
  }
}
void OpenGL::offscreenGrab(byteA& image){
  if(image.nd==3){ CHECK(image.d2==4,"3rd dim of image has to be 4 for RGBA");}
  else{ CHECK(image.nd==2,"image has to be either 2- or 3(for RGBA)-dimensional");}
  setOffscreen(image.d1,image.d0);
  Draw(camera,image.d1,image.d0);
  glGrabImage(image);
}
void OpenGL::offscreenGrab(byteA& image,byteA& depth){
  if(image.nd==3){ CHECK(image.d2==4,"3rd dim of image has to be 4 for RGBA");}
  else{ CHECK(image.nd==2,"image has to be either 2- or 3(for RGBA)-dimensional");}
  CHECK(depth.nd==2,"depth buffer has to be either 2-dimensional");
  setOffscreen(image.d1,image.d0);
  Draw(camera,image.d1,image.d0);
  glGrabImage(image);
  glGrabDepth(depth);
}
void OpenGL::offscreenGrabDepth(byteA& depth){
  CHECK(depth.nd==2,"depth buffer has to be either 2-dimensional");
  setOffscreen(depth.d1,depth.d0);
  Draw(camera,depth.d1,depth.d0);
  glGrabDepth(depth);
}
void OpenGL::offscreenGrabDepth(floatA& depth){
  CHECK(depth.nd==2,"depth buffer has to be either 2-dimensional");
  setOffscreen(depth.d1,depth.d0);
  Draw(camera,depth.d1,depth.d0);
  glGrabDepth(depth);
}
void OpenGL::setOffscreen(int width,int height){
  createOffscreen(width,height);
  CHECK(width<=osPixmap->width() && height<=osPixmap->height(),
	"width ("<<width<<") or height ("<<height
	<<") too large for the created pixmap - create and set size earlier!");
  osContext->makeCurrent();
  //if(initRoutine) (*initRoutine)();
}
#endif


void OpenGL::reportSelection(){
  uint i;
  std::cout <<"selection report: mouse=" <<mouseposx <<" " <<mouseposy <<" -> #selections=" <<selection.N <<std::endl;
  for(i=0;i<selection.N;i++){
    if(topSelection == &selection(i)) std::cout <<"  TOP: "; else std::cout <<"       ";
    std::cout
      <<"name = 0x" <<std::hex <<selection(i).name <<std::dec
      <<" min-depth:" <<selection(i).dmin <<" max-depth:" <<selection(i).dmax
      <<endl;
  }
}

#ifdef MT_FREEGLUT
int OpenGL::width(){  glutSetWindow(windowID); return glutGet(GLUT_WINDOW_WIDTH); }
int OpenGL::height(){ glutSetWindow(windowID); return glutGet(GLUT_WINDOW_HEIGHT); }
#endif
#ifndef MT_GL
int OpenGL::width(){ return 0; }
int OpenGL::height(){ return 0; }
#endif

#ifdef MT_GL2PS
void OpenGL::saveEPS(const char *filename){
  FILE *fp = fopen(filename, "wb");
  GLint buffsize = 0, state = GL2PS_OVERFLOW;
  GLint viewport[4];
  glGetIntegerv(GL_VIEWPORT, viewport);
  while(state==GL2PS_OVERFLOW){
    buffsize+=1024*1024;
    gl2psBeginPage( "Marc Toussaint", "MT", viewport,
		    GL2PS_EPS, GL2PS_BSP_SORT, GL2PS_SILENT |
		    GL2PS_SIMPLE_LINE_OFFSET | GL2PS_NO_BLENDING |
		    GL2PS_OCCLUSION_CULL | GL2PS_BEST_ROOT,
		    GL_RGBA, 0, NULL, 0, 0, 0, buffsize,
		    fp, filename );
    Draw(camera,width(),height());
    state = gl2psEndPage();
  }
  fclose(fp);
}
#else
void OpenGL::saveEPS(const char *filename){
  MT_MSG("WARNING: OpenGL::saveEPS was called without MT_GL2PS configured!")
    }
#endif

#ifdef MT_QT
void OpenGL::about(std::ostream& os){
  os <<"Widget's OpenGL capabilities:\n";
  QGLFormat f=format();
  os <<"direct rendering: " <<f.directRendering() <<"\n"
     <<"double buffering: " <<f.doubleBuffer()  <<"\n"
     <<"depth:            " <<f.depth() <<"\n"
     <<"rgba:             " <<f.rgba() <<"\n"
     <<"alpha:            " <<f.alpha() <<"\n"
     <<"accum:            " <<f.accum() <<"\n"
     <<"stencil:          " <<f.stencil() <<"\n"
     <<"stereo:           " <<f.stereo() <<"\n"
     <<"overlay:          " <<f.hasOverlay() <<"\n"
     <<"plane:            " <<f.plane() <<std::endl;

  if(!osContext){
    os <<"no off-screen context created yet" <<std::endl;
  }else{
    os <<"Off-screen pixmaps's OpenGL capabilities:\n";
    f=osContext->format();
    os <<"direct rendering: " <<f.directRendering() <<"\n"
       <<"double buffering: " <<f.doubleBuffer()  <<"\n"
       <<"depth:            " <<f.depth() <<"\n"
       <<"rgba:             " <<f.rgba() <<"\n"
       <<"alpha:            " <<f.alpha() <<"\n"
       <<"accum:            " <<f.accum() <<"\n"
       <<"stencil:          " <<f.stencil() <<"\n"
       <<"stereo:           " <<f.stereo() <<"\n"
       <<"overlay:          " <<f.hasOverlay() <<"\n"
       <<"plane:            " <<f.plane() <<std::endl;
  }
}
#else
void OpenGL::about(std::ostream& os){ MT_MSG("NIY"); }
#endif


//===========================================================================
//
// modified main loop routine
// only change to free glut: don't deinitialize the glut engine!
//

#if defined MT_FREEGLUT //&& !defined MT_Cygwin
extern "C"{
#include "opengl_freeglut_modified.c"
}
#endif


//===========================================================================
//
// callbacks
//

#if 1
#  define CALLBACK_DEBUG(x) if(reportEvents) x
#else
#  define CALLBACK_DEBUG(x)
#endif

void getSphereVector(geo3d::Vector& vec,int _x,int _y,int width,int height){
  int minwh = width<height?width:height;
  double x,y;
  x=(double)_x;  x=x-.5*width;  x*= 2./minwh;
  y=(double)_y;  y=y-.5*height; y*=-2./minwh;
  vec(0)=x;
  vec(1)=y;
  vec(2)=.5-(x*x+y*y);
  if(vec(2)<0.) vec(2)=0.;
}

void OpenGL::Reshape(int width, int height){
  CALLBACK_DEBUG(printf("Window %d Reshape Callback:  %d %d\n", windowID, width, height ));
  camera.setWHRatio((double)width/height);
  //update();
}

void OpenGL::Key(unsigned char key, int x, int y){
  CALLBACK_DEBUG(printf("Window %d Keyboard Callback:  %d (`%c') %d %d\n", windowID, key, (char)key, x, y ));
  pressedkey=key;
  if(key==13 || key==32 || key==27) MTexitLoop();
}

void OpenGL::Mouse(int button, int updown, int _x, int _y){
  CALLBACK_DEBUG(printf("Window %d Mouse Click Callback:  %d %d %d %d\n", windowID, button, updown, _x, _y ));
  mouse_button=1+button;
  if(updown) mouse_button=-1-mouse_button;
  mouseposx=_x; mouseposy=_y;
  geo3d::Vector vec;
  getSphereVector(vec,_x,_y,width(),height());
  lastEvent.set(mouse_button,-1,_x,_y,0.,0.);
  bool cont=true;
  if(updown==0){
    for(uint i=0;i<clickCalls.N;i++) cont=cont && (*clickCalls(i).call)(clickCalls(i).classP,this);
  }
  if(cont){
    if(mouse_button>0){
      downVec=vec;
      downRot=camera.f->r;
      downPos=camera.f->p - camera.foc;
      if(reportSelects) Select();
    }
    if(mouse_button==3){ //right down
      MTexitLoop();
    }
  }
}

void OpenGL::Motion(int _x, int _y){
  CALLBACK_DEBUG(printf("Window %d Mouse Motion Callback:  %d %d\n", windowID, _x, _y ));
  geo3d::Vector vec;
  getSphereVector(vec,_x,_y,width(),height());
  lastEvent.set(mouse_button,-1,_x,_y,vec(0)-downVec(0),vec(1)-downVec(1));
  //mouseposx=vec(0); mouseposy=vec(1);
  mouseposx=_x; mouseposy=_y;
  if(mouse_button==1){
    geo3d::Rotation rot;
    if(downVec(2)<.1){
      rot.setDiff(vec,downVec);  //consider imaged shere rotation of mouse-move
    }else{
      rot.setVec((vec-downVec) ^ geo3d::Vector(0,0,1)); //consider only xy-mouse-move
    }
    rot = downRot * rot / downRot; //interpret rotation relative to current viewing
    camera.f->p = camera.foc + rot * downPos;   //rotate camera's position
    camera.f->r = rot * downRot;   //rotate camera's direction
    update();
  }
  if(mouse_button==2){
    double dy = vec(1) - downVec(1);
    if(dy<-.99) dy = -.99;
    camera.f->p = camera.foc + downPos / (dy+1.);
    update();
  }
}

void OpenGL::PassiveMotion(int x, int y){
  CALLBACK_DEBUG(printf("Window %d Mouse Passive Motion Callback:  %d %d\n", windowID, x, y ));
  mouseposx=x; mouseposy=y;
  if(reportSelects) Select();
  bool ud=false;
  for(uint i=0;i<hoverCalls.N;i++) ud=ud || (*hoverCalls(i).call)(hoverCalls(i).classP,this);
  if(ud) update();
}


//===========================================================================
//
// GUI implementation
//

void glUI::addButton(uint x,uint y,const char *name,const char *img1,const char *img2){
  MT::List<Button>::node b=buttons.new_node();
  byteA img;
  b->hover=false;
  b->x=x; b->y=y; b->name=name;
  //read_png(img,tex1);
  if(img1){
    read_ppm(img,img1,true);
  }else{
    img.resize(18,strlen(name)*9+10,3);
    img=255;
  }
  b->w=img.d1; b->h=img.d0;
  b->img1=img;    add_alpha_channel(b->img1,100);
  if(img2){
    read_ppm(img,img1,true);
    CHECK(img.d1==b->w && img.d0==b->h,"mismatched size");
  }
  b->img2=img;    add_alpha_channel(b->img2,200);
}

void glUI::glDraw(){
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  GLint viewPort[4];
  glGetIntegerv(GL_VIEWPORT,viewPort);
  glOrtho(0.,viewPort[2],viewPort[3],.0,-1.,1.);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
    
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  MT::List<Button>::node b;
  float x,y,w,h;
  forNodes(b,buttons){
    x=b->x-b->w/2.;
    y=b->y-b->h/2.;
    w=b->w;
    h=b->h;
    glColor(0,0,0,1);
    glDrawText(b->name,x+5,y+h-5);
    if((int)b->index==top) glRasterImage((int)x,(int)y,b->img2);
    else	             glRasterImage((int)x,(int)y,b->img1);
  }
}

bool glUI::checkMouse(int _x,int _y){
  float x,y,w,h;
  MT::List<Button>::node b;
  int otop=top;
  top=-1;
  forNodes(b,buttons){
    x=b->x-b->w/2.;
    y=b->y-b->h/2.;
    w=b->w;
    h=b->h;
    if(_x>=x && _x <=x+w && _y>=y && _y<=y+h) top = b->index;
  }
  if(otop==top) return false;
  //glutPostRedisplay();
  return true;
}


#ifdef MT_QT

#if   defined MT_MSVC
#  include"opengl_MSVC.moccpp"
#elif defined MT_SunOS
#  include"opengl_SunOS.moccpp"
#elif defined MT_Linux
#  include"opengl_Linux.moccpp"
#elif defined MT_Cygwin
#  include"opengl_Cygwin.moccpp"
#endif

#endif
