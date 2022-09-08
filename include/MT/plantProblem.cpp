#include"plantProblem.h"

byteA col;
floatA dep;
OpenGL *gl1=0,*gl2=0;

#define MT_maxPlantSymbol 9

//===========================================================================
//
// some parameters read from MT.cfg
//

struct PlantProblemParameters{
  uint xyResolution;
  uint boxRange;
  double weightFactor;
  geo3d::Vector cameraPos;

  PlantProblemParameters(){
    xyResolution=Parameter<uint>("xyResolution");
    boxRange=Parameter<uint>("boxRange");
    weightFactor=Parameter<double>("weightFactor");
    cameraPos=geo3d::Vector(-.0, -.6, .5); //cameraPos: [1:3] -1. -1.3 .8
  }
}PPP;


//===========================================================================
//
// the Plant class
//

/*! A plant in the sense of Prusinkiewicz: A string with turtle
    commands and a method to draw it in OpenGL. The \ref ProPlant
    class handles such plant object to display or evaluate them.*/
class Plant{
public:
  
  Parameter<int> degree;
  Parameter<float> leaf;

  /*! the string of turtle commands (including a final zero) */
  MT::Array<char> str;
  /*! the weight associated to each element of the string (except the last zero) */
  MT::Array<uint> weights;

  /*! the green value of its leaves (to distinguish from other plants) */
  byte color;
  /*! the number of elements (is set when glDraw() was called) */
  uint elements;
  /*! plant's weight (is set when weighing() was called) */
  uint weight;
  /*! plant's fitness (can be set externally) */
  double fitness;

  Plant():
    degree("plantDegree",20),
    leaf("plantLeaf",.3)
  {
    color=0xff;
    elements=0;
  }


public://{ Assignment
  /*! set the plants string */
  void set(const charA& a){ str=a; translate(); }
  /*! set the plants string */
  void set(const byteA& x){ str.copy(x); translate(); };
  /*! set the plants string */
  void set(const char* s,int N){ str.resize(N); for(int i=0;i<N;i++) str(i)=s[i]; translate(); }
  /*! set the plants string */
  void set(String& s){ set(s,s.N()); translate(); }

//private:
  void translate(){
    str.resizeCopy(str.N+1);
    uint i;
    for(i=0;i<str.N-1;i++){
      switch(str(i)){
      case 0: str(i)='F'; break; //A
      case 1: str(i)='.'; break; //B
      case 2: str(i)='+'; break; //C
      case 3: str(i)='-'; break; //D
      case 4: str(i)='&'; break; //E
      case 5: str(i)='^'; break; //F
      case 6: str(i)='\\'; break;//G
      case 7: str(i)='/'; break; //H
      case 8: str(i)='['; break; //I
      case 9: str(i)=']'; break; //J
      default: str(i)='.'; break;
      }        
    }
    str(i)=0;
  }


public://{ Evaluation
  /*! the OpenGL drawing routine; calculates also the number of
      elements (=stems+leaves) */
  void glDraw(){
    elements=0;
    glPushMatrix();
    glRotatef(180,1,0,0);
    int i=0;
    while(i!=(int)str.N) glDraw(i);
    glPopMatrix();
  }

  /*! `weighs' the each stem and leave of the plant; leaves have
    weight 1, stems have weight 1 plus the weight `resting on them from
    above' */
  uint weighing(){
    weight=0;
    int i=str.N-1;
    weights.resize(i);
    weights=0;
    weighing(i);
    for(i=0;i<(int)weights.N;i++) weight+=weights(i);
    return weight;
  }
private:
    uint weighing(int& i,uint topweight=0){
    static uint calls=0;
    CHECK(calls==0 || topweight!=0,"Plant::weighing::calls is not zero!");
    if(calls++>10000){ i--; calls--; return 1+topweight; }

    uint weight=topweight;
    for(;i--;){
      switch(str(i)){
      case 'F':	weight++; weights(i)=weight; break;
      case '[':	if(topweight){ calls--; return weight; } break;
      case ']':	weights(i)=1; weight+=weighing(i,1); break;
      default:	break;
      }
    }
    i=0;
    calls--;
    return weight;
  }


public:
  static void staticDraw(void *p){ ((Plant*)p)->glDraw(); }
  void glDraw(int& i){
    static uint calls=0;
    CHECK(calls==0 || i!=0,"Plant::draw::calls is not zero!");
    if(calls++>10000){ calls--; return; }

    glPushMatrix();	
    for(;i<(int)str.N;i++){
      switch(str(i)){
      case 'F':	drawStem(); elements++; break;
      case 'f':	glTranslatef(0,0,-1); elements++; break;
      case '+':	glRotatef( degree(),0,1,0);	break;
      case '-':	glRotatef(-degree(),0,1,0);	break;
      case '&':	glRotatef(-degree(),1,0,0);	break;
      case '^':	glRotatef( degree(),1,0,0);	break;
      case '\\':glRotatef(-degree(),0,0,1);	break;
      case '/':	glRotatef(-degree(),0,0,1);	break;
      case '[':	i++; glDraw(i); i--; break;
      case ']':	drawLeaf(); elements++; i++; glPopMatrix(); calls--; return;
      case '.':	break;
      case 0:	break;
      //default:	HALT("wrong symbol '" <<char(str(i)) <<"' in plant string");
      default:	HALT("wrong symbol in plant string");
      }
    }
    glPopMatrix();
    //drawLeaf();
    calls--;
  }
  void drawStem(){
    glColor3ub(0xff,0,0);
    glBegin(GL_LINES);
    glVertex3f(0,0,0);
    glVertex3f(0,0,-1.);
    glEnd();
    glTranslatef(0,0,-1);
  }
  void drawLeaf(){
    glColor3ub(0,color,0);
    glBegin(GL_POLYGON);
    glVertex3f(0,0,0);
    glVertex3f( leaf(),0,-1); //.5!!
    glVertex3f(-leaf(),0,-1); //.5!!
    glVertex3f(0,0,0);
    glEnd();
  }


public://{ I/O
  void write(std::ostream& os) const{ os <<str.p; }
};
stdOutPipe(Plant);


Plant plant;

//===========================================================================
//
// evaluate and display
//

void evaluatePlant(double info[3],byte* x,uint x_size){
  byteA   X(x,x_size);
  doubleA Info(info,3);
  evaluatePlant(Info,X);
}

void evaluatePlant(arr& info,const byteA& x){
  col.resize(PPP.xyResolution,PPP.xyResolution,4);
  dep.resize(PPP.xyResolution,PPP.xyResolution);

  plant.set(x);

  if(!gl1){
    gl1=new OpenGL("plant evaluation",PPP.xyResolution,PPP.xyResolution);
    gl1->clear();
    gl1->add(plant);
  }
  gl1->camera.reset();
  gl1->setClearColors(0.,0.,0.,0.);
  gl1->camera.setPosition(0,0,PPP.boxRange);
  gl1->camera.watchDirection(0,0,-1);
  gl1->camera.setHeightAngle(0);
  gl1->camera.setHeightAbs(PPP.boxRange);
  gl1->camera.setZRange(0,PPP.boxRange);

  gl1->text.clr();
  gl1->update();
  glGrabImage(col);
  glGrabDepth(dep);
  
  uint i,j;
  double height,integral=0;
  for(i=0;i<PPP.xyResolution;i++) for(j=0;j<PPP.xyResolution;j++){
    if(col(i,j,1)==plant.color){ //check for plant's green! color
      height=((float)255-dep(i,j))/255.;
      integral+=height;
    }
  }
  double weight=plant.weighing();
  
  plant.fitness=integral/(PPP.xyResolution*PPP.xyResolution) - PPP.weightFactor*weight;

  info.resize(3);
  info(0)=plant.fitness;
  info(1)=integral;
  info(2)=weight;
}

void drawenv(void*){
  glDrawFloor(PPP.boxRange,.6,.6,.6);
}

void displayPlant(byte* x,uint x_size,bool watch,char* filename){
  byteA  X(x,x_size);
  displayPlant(X,watch,filename);
}

void displayPlant(byteA& x,bool watch,char* filename){
  arr info;
  plant.set(x);
  evaluatePlant(info,x);

  if(!gl2){
    gl2=new OpenGL("plant visualization",300,300);
    gl2->clear();
    gl2->add(plant);
    gl2->add(drawenv,0);
    gl2->camera.reset();
    gl2->setClearColors(1.,1.,1.,0.);
    gl2->camera.setPosition(10.*PPP.cameraPos(0)*PPP.boxRange,
      10.*PPP.cameraPos(1)*PPP.boxRange,
      10.*PPP.cameraPos(2)*PPP.boxRange);
    gl2->camera.focus(0,0,.5*PPP.boxRange);
    gl2->camera.setZRange(1.,1000.);
    gl2->camera.setHeightAngle(12.);
  }
  gl2->text.clr() <<"value=" <<info(0) <<" green=" <<info(1) <<" weight=" <<info(2);
  gl2->update();
  if(watch) gl2->watch();

  //glWatchImage(col);
  //glRasterImage(0,100,dep);

  if(filename){
    byteA img;
    img.resize(gl2->height(),gl2->width(),3);
    glGrabImage(img);
    write_ppm(img,filename);
  }
}



//===========================================================================
//
// rules
//

void readGenoFile(std::istream& is, MT::Array<byteRule>& rules, arr& reals){
  uint n;
  char c;
  is >>"<genotype [" >>n >>"]>";
  rules.resize(n+1);
  rules(0).code=0;
  rules(0).str.resize(0);
  is >>"(A):<";
  for(;;){ is.get(c); if(c=='>') break; rules(0).str.append(c-'A'); }
  for(uint i=1;i<rules.N;i++){
    is >>"(" >>n >>"):<";
    CHECK(i==n+1,"rule numbering is wrong");
    is.get(c);
    rules(i).code=c-'A';
    rules(i).str.resize(0);
    for(;;){ is.get(c); if(c=='>') break; rules(i).str.append(c-'A'); }
  }
  reals.resize(3);
  is >>"(reals): [1:3]" >>reals(0) >>reals(1) >>reals(2);
}

void writeByteRule(std::ostream& os, const MT::Array<byteRule>& rules){
  uint i,j;
  for(i=0;i<rules.N;i++){
    os <<"(" <<i <<"):<" <<(char)(rules(i).code+'A') <<"|";
    for(j=0;j<rules(i).str.N;j++) os <<(char)(rules(i).str(j)+'A');
    os <<">" <<endl;
  }
}

void decodeRules(byteA& x,const MT::Array<byteRule>& rules,uint maxT,uint maxN){
  uint t,i,j;
  bool change;
  for(i=0;i<rules.N;i++){ rules(i).use=0; rules(i).hierarchy=0; }
  x=rules(0).str;
  for(t=0;t<maxT;t++){
    change=false;
    for(i=1;i<rules.N;i++){
      if(x.N>maxN) return;
      for(j=0;j<x.N;j++) if(x(j)==rules(i).code){
	//cout <<t <<' ' <<i <<' ' <<j <<' ' <<x+'A' <<rules(i).str+'A';
	if(x.N-1+rules(i).str.N<maxN){
	  x.replace(j,1,rules(i).str);
	  j+=rules(i).str.N-1;
	  rules(i).use++;
	  if(!rules(i).hierarchy) rules(i).hierarchy=t;
	  change=true;
	}
	//cout <<x+'A' <<endl;
      }
    }
    if(!change) return;
  }
}

