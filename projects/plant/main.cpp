#define MT_IMPLEMENTATION

#include<MT/opengl.h>
#include<MT/plantProblem.h>

using namespace std;


//scan a genotype file and display/store phenotype pictures
void movie(){
  Parameter<String> genofile("genofile");
  Parameter<uint> maxN("maxResource");
  Parameter<uint> maxT("maxBuildTime");
  //Parameter<String> moviepath ("moviepath");
  //Parameter<uint> movierate ("movierate");
  //Parameter<String> repfile("repfile");


  //ProPlant pro;
  MT::Array<byteRule> rules;
  byteA str;
  arr reals;

  //setup input and output file
  ifstream gen;
  MT::open(gen,genofile());
  //ofstream rep;
  //MT::open(rep,repfile());

  //setup window size and image class
  /*
  floatA img; //img.mirrorDisplay=true;
  img.resize(pro.opengl->width(),pro.opengl->height(),4);
  String outfile;
  QPixmap pix(600,600);
  */

  //read buffers
  uint g;
  double fit=0.;
  ARRAYOLDREAD=1;

  //loop
  while(gen.good()){
    gen >>"generation:" >>g;
    readGenoFile(gen,rules,reals);
    gen >>"fitness:" >>fit;

    cout <<"generation: " <<g <<endl;
    writeByteRule(cout,rules);
    cout <<"reals: " <<reals <<endl;
    cout <<"fitness:" <<fit <<endl;

    //if(!(g%movierate())){
      //axiom.referToSubRange(gene.axiom.str,1,-1);
      decodeRules(str,rules,maxT,maxN);
      //cout <<"decoded: "; for(uint i=0;i<str.N;i++) cout <<(char)(str(i)+'A'); cout <<endl;
      displayPlant(str,false);

      /*
      if(moviepath()!="none"){
#if 0
	pro.opengl->grepImage(img);
	outfile=moviepath()+"/pic"+g+".jpg";
	cout <<"saving image `" <<outfile <<"'" <<endl;
	img.save(outfile,"JPEG");
#else
	pro.opengl->grabImage(img);
	pix=img.getQPixmap();
	pro.paintOn(&pix,String()()<<g,1.5);
	img.setQPixmap(pix);
	img.display();
	outfile=moviepath()+"/pic"+g+".bmp";
	cout <<"saving image `" <<outfile <<"'" <<endl;
	img.save(outfile().str(),"BMP");
#endif
	MT::wait();
      }
      rep <<g <<" " <<fit <<" "
	  <<pro.plant.fitness <<" " <<pro.plant.elements <<" " <<pro.plant.weight<<" "
	  <<gene.reals(0) <<" " <<gene.totUsage() <<" "
	  <<gene.totSize() <<" " <<gene.rule.N() <<" "
	  <<endl;
      */
  }
}


void single(){
  arr reals;
  MT::Array<byteRule> rules;
  byteA str;
  std::ifstream is("prun");
  readGenoFile(is,rules,reals);
  decodeRules(str,rules,4,1000000);

  displayPlant(str);
}

int main(int argc,char **argv){
  MT::init(argc,argv);
  cout.precision(3);
  
  //  single();
  movie();

  return 0;
}



