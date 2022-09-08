#include "geo3d.h"

/*void geo3d::Trimesh::getOBJ(char* filename){
  if (!glm){
  glm = glmReadOBJ(filename);
  glmReverseWinding(glm);
  }
  
  ////glmUnitize(glm);
  glmFacetNormals(glm);
  glmVertexNormals(glm, 90.0);
  
  // creates a display list for the OBJ
  ////	g._pmodel_displaylist = glmList(glm, GLM_SMOOTH | GLM_MATERIAL);
  }*/

// dm 07.06.2006
/*!\ initialises the ascii-obj file "filename"*/
void geo3d::Trimesh::readOBJ(char* filename){
  FILE* file;
  
  // open the file 
  file = fopen(filename, "r");
  if (!file) {
    fprintf(stderr, "readOBJ() failed: can't open data file \"%s\".\n",
	    filename);
    exit(1);
  }
  else {
    printf("%s opened.\n ", filename);
  }
  // allocate a new model
  //	pathname = strdup(filename);
  
  //uint nV=0,nT=0,nVN=0,nFN=0,nTex=0;
  
  // make a first pass through the file to get a count of the number
  // of vertices, normals, texcoords & triangles 
  firstPass(file);
  
  
  // allocate memory 
  
  // rewind to beginning of file and read in the data this pass 
  rewind(file);
  
  secondPass(file);
  
  // close the file 
  fclose(file);
 
}


/*!\First pass at a Wavefront OBJ file that gets all the
  statistics of the model (such as #vertices, #normals, etc)

  file  - (fopen'd) file descriptor 
*/
void geo3d::Trimesh::firstPass(FILE* file){
  uint nV;        
  uint nN;        
  uint nTex;       
  uint nT;       
  ////    GLMgroup* group;           // current group 
  int v, n, t;
  char buf[128];
  
  /* make a default group */
  /*   group = glmAddGroup(model, "default");*/
  
  nV = nN = nTex = nT = 0;
  while(fscanf(file, "%s", buf) != EOF) {
    switch(buf[0]) {
    case '#':               // comment
      // eat up rest of line 
      fgets(buf, sizeof(buf), file);
      break;
    case 'v':               // v, vn, vt 
      switch(buf[1]) {
      case '\0':          // vertex 
	// eat up rest of line 
	fgets(buf, sizeof(buf), file);
	nV++;
	break;
      case 'n':           // normal 
	// eat up rest of line 
	fgets(buf, sizeof(buf), file);
	nN++;
	break;
      case 't':           // texcoord 
	// eat up rest of line 
	fgets(buf, sizeof(buf), file);
	nTex++;
	break;
      default:
	printf("firstPass(): Unknown token \"%s\".\n", buf);
	exit(1);
	break;
      }
      break;
    case 'm':
      // TODO: implement mat
      /*fgets(buf, sizeof(buf), file);
	sscanf(buf, "%s %s", buf, buf);
	mtllibname = strdup(buf);
	glmReadMTL(model, buf);*/
      break;
    case 'u':
      // eat up rest of line 
      fgets(buf, sizeof(buf), file);
      break;
    case 'g':               // group 
      // eat up rest of line 
      fgets(buf, sizeof(buf), file);
#if SINGLE_STRING_GROUP_NAMES
      sscanf(buf, "%s", buf);
#else
      buf[strlen(buf)-1] = '\0';  //nuke '\n' 
#endif
      ////group = glmAddGroup(model, buf);
      break;
    case 'f':               // face 
      v = n = t = 0;
      fscanf(file, "%s", buf);
      // can be one of %d, %d//%d, %d/%d, %d/%d/%d %d//%d 
      if (strstr(buf, "//")) {
	// v//n 
	sscanf(buf, "%d//%d", &v, &n);
	fscanf(file, "%d//%d", &v, &n);
	fscanf(file, "%d//%d", &v, &n);
	nT++;
	//// group->numtriangles++;
	while(fscanf(file, "%d//%d", &v, &n) > 0) nT++;
      } else if (sscanf(buf, "%d/%d/%d", &v, &t, &n) == 3) {
	// v/t/n 
	fscanf(file, "%d/%d/%d", &v, &t, &n);
	fscanf(file, "%d/%d/%d", &v, &t, &n);
	nT++;
	//// group->numtriangles++;
	while(fscanf(file, "%d/%d/%d", &v, &t, &n) > 0) nT++;
      } else if (sscanf(buf, "%d/%d", &v, &t) == 2) {
	// v/t 
	fscanf(file, "%d/%d", &v, &t);
	fscanf(file, "%d/%d", &v, &t);
	nT++;
	////group->numtriangles++;
	while(fscanf(file, "%d/%d", &v, &t) > 0) nT++;
      } else {
	// v 
	fscanf(file, "%d", &v);
	fscanf(file, "%d", &v);
	nT++;
	while(fscanf(file, "%d", &v) > 0) nT++;
      }
      break;
	
    default:
      // eat up rest of line 
      fgets(buf, sizeof(buf), file);
      break;
    }
  }

  V.resize(nV,3);
  T.resize(nT,3);
  N.resize(nN,3);
  //if(nVN) N.resize(nVN,3);
  //if(nTex) Tex = (double*)malloc(sizeof(double) * 2 * (numtexcoords));
}
  
  
/*!\glmSecondPass: second pass at a Wavefront OBJ file that gets all
  the data.
  
  file  - (fopen'd) file descriptor 
*/
void geo3d::Trimesh::secondPass( FILE* file) 
{
  uint nV;        /* number of vertices in model */
  uint nN;         /* number of normals in model */
  uint nTex;       /* number of texcoords in model */
  uint nT;       /* number of triangles in model */
  int v, n, t;
  char buf[128];
      
  /* on the second pass through the file, read all the data into the
     allocated arrays */
  nV = nN = nTex = nT = 0; //changed from 1 to 0, Sep 1 06 (mt)
  ////_material = 0;
  while(fscanf(file, "%s", buf) != EOF) {
	
    switch(buf[0]) {
    case '#':               // comment 
      // eat up rest of line /
      fgets(buf, sizeof(buf), file);
      break;
    case 'v':               // v, vn, vt 
      switch(buf[1]) {
      case '\0':          // vertex 
	fscanf(file, "%lf %lf %lf", &V(nV,0),&V(nV,1),&V(nV,2));
	nV++;
	break;
      case 'n':           // normal 
	fscanf(file, "%lf %lf %lf", &N(nN,0),&N(nN,1),&N(nN,2));
	nN++;
	break;
      case 't':           // texcoord 
	/*    fscanf(file, "%f %f", 
	      &texcoords[2 * numtexcoords + 0],
	      &texcoords[2 * numtexcoords + 1]);
	      numtexcoords++;*/
	break;
      }
      break;
    case 'u':
      /* fgets(buf, sizeof(buf), file);
	 sscanf(buf, "%s %s", buf, buf);
	 group->material = material = glmFindMaterial(model, buf);*/
      break;
    case 'g':               // group 
      // eat up rest of line 
      fgets(buf, sizeof(buf), file);
#if SINGLE_STRING_GROUP_NAMES
      sscanf(buf, "%s", buf);
#else
      buf[strlen(buf)-1] = '\0';  // nuke '\n' 
#endif
      ////  group = glmFindGroup(model, buf);
      ////  group->material = material;
      break;
    case 'f':               // face 
      v = n = t = 0;
      fscanf(file, "%s", buf);
      //can be one of %d, %d//%d, %d/%d, %d/%d/%d %d//%d 
      if (strstr(buf, "//")) {
	// v//n 
	sscanf(buf, "%d//%d", &v, &n);

	T(nT,0) = v < 0 ? v + nV : v;
	Tni(nT,0) = n < 0 ? n + nN : n;
	fscanf(file, "%d//%d", &v, &n);
	T(nT,1) = v < 0 ? v + nV : v;
	Tni(nT,1) = n < 0 ? n + nN : n;
	fscanf(file, "%d//%d", &v, &n);
	T(nT,2) = v < 0 ? v + nV : v;
	Tni(nT,2) = n < 0 ? n + nN : n;
	//// group->triangles[group->nT++] = nT;
	nT++;
	while(fscanf(file, "%d//%d", &v, &n) > 0) {
	  T(nT,0) = T(nT-1,0);
	  Tni(nT,0) = Tni(nT-1,0);
	  T(nT,1) = T(nT-1,2);
	  Tni(nT,1) = Tni(nT-1,2);
	  T(nT,2) = v < 0 ? v + nV : v;
	  Tni(nT,2) = n < 0 ? n + nN : n;
	  //// group->triangles[group->numtriangles++] = numtriangles;
	  nT++;
	}
      } else if (sscanf(buf, "%d/%d/%d", &v, &t, &n) == 3) {
	// v/t/n 
	T(nT,0) = v < 0 ? v + nV : v;
	Tti(nT,0) = t < 0 ? t + nTex : t;
	Tni(nT,0) = n < 0 ? n + nN : n;
	fscanf(file, "%d/%d/%d", &v, &t, &n);
	T(nT,1) = v < 0 ? v + nV : v;
	Tti(nT,1) = t < 0 ? t + nTex : t;
	Tni(nT,1) = n < 0 ? n + nN : n;
	fscanf(file, "%d/%d/%d", &v, &t, &n);
	T(nT,2) = v < 0 ? v + nV : v;
	Tti(nT,2) = t < 0 ? t + nTex : t;
	Tni(nT,2) = n < 0 ? n + nN : n;
	//// group->triangles[group->numtriangles++] = numtriangles;
	nT++;
	while(fscanf(file, "%d/%d/%d", &v, &t, &n) > 0) {
	  T(nT,0) = T(nT-1,0);
	  Tti(nT,0) = Tti(nT-1,0);
	  Tni(nT,0) = Tni(nT-1,0);
	  T(nT,1) = T(nT-1,2);
	  Tti(nT,1) = Tti(nT-1,2);
	  Tni(nT,1) = Tni(nT-1,2);
	  T(nT,2) = v < 0 ? v + nV : v;
	  Tti(nT,2) = t < 0 ? t + nTex : t;
	  Tni(nT,2) = n < 0 ? n + nN : n;
	  //// group->triangles[group->numtriangles++] = numtriangles;
	  nT++;
	}
      } else if (sscanf(buf, "%d/%d", &v, &t) == 2) {
	// v/t 
	      
	T(nT,0) = v < 0 ? v + nV : v;
	Tti(nT,0) = t < 0 ? t + nTex : t;
	fscanf(file, "%d/%d", &v, &t);
	T(nT,1) = v < 0 ? v + nV : v;
	Tti(nT,1) = t < 0 ? t + nTex : t;
	fscanf(file, "%d/%d", &v, &t);
	T(nT,2) = v < 0 ? v + nV : v;
	Tti(nT,2) = t < 0 ? t + nTex : t;
	//// group->triangles[group->numtriangles++] = numtriangles;
	nT++;
	while(fscanf(file, "%d/%d", &v, &t) > 0) {
	  T(nT,0) = T(nT-1,0);
	  Tti(nT,0) = Tti(nT-1,0);
	  T(nT,1) = T(nT-1,2);
	  Tti(nT,1) = Tti(nT-1,2);
	  T(nT,2) = v < 0 ? v + nV : v;
	  Tti(nT,2) = t < 0 ? t + nTex : t;
	  //// group->triangles[group->numtriangles++] = numtriangles;
	  nT++;
	}
      } else {
	// v 
	sscanf(buf, "%d", &v);
	T(nT,0) = v < 0 ? v + nV : v;
	fscanf(file, "%d", &v);
	T(nT,1) = v < 0 ? v + nV : v;
	fscanf(file, "%d", &v);
	T(nT,2) = v < 0 ? v + nV : v;
	//// group->triangles[group->numtriangles++] = nT;
	nT++; // fixed bug, 21. Jun 06 (hh)
	while(fscanf(file, "%d", &v) > 0) {
	  T(nT,0) = T(nT-1,0);
	  T(nT,1) = T(nT-1,2);
	  T(nT,2) = v < 0 ? v + nV : v;
	  //// group->triangles[group->numtriangles++] = numtriangles;
	  nT++;
	}
      }
      break;
	    
    default:
      // eat up rest of line 
      fgets(buf, sizeof(buf), file);
      break;
    }
  }

  //Sep 1, 06 (mt): start counting indices from 0
  T -= 1;
}
  
void geo3d::Trimesh::clear()
{
  V.resize(0); N.resize(0); T.resize(0); C.resize(0); strips.resize(0);
}
  
/*!\brief calculate the normals of all triangles (Tn) and the average
  normals of the vertices (N); average normals are averaged over
  all adjacent triangles that are in the triangle list or member of
  a strip */
void geo3d::Trimesh::calcNormals(){
  uint i;
  Vector a,b,c;
  N.resize(V.d0,3);
  N.setZero();
  //triangle normals and contributions    
  for(i=0;i<T.d0;i++){      
    a=&V(T(i,0),0); b=&V(T(i,1),0); c=&V(T(i,2),0);
    b-=a; c-=a; a=c^b; a.normalize();
    //Tn(i,0)=a[0];  Tn(i,1)=a[1];  Tn(i,2)=a[2];
    N(T(i,0),0)+=a[0];  N(T(i,0),1)+=a[1];  N(T(i,0),2)+=a[2];
    N(T(i,1),0)+=a[0];  N(T(i,1),1)+=a[1];  N(T(i,1),2)+=a[2];
    N(T(i,2),0)+=a[0];  N(T(i,2),1)+=a[1];  N(T(i,2),2)+=a[2];
  }
  // dm 21.06.2006: Removed because not used and caused crash: 
  Vector *d;
  for(i=0;i<N.d0;i++){ d=(Vector*)&N(i,0); d->normalize(); }
}
  
/*!\brief add triangles according to the given grid; grid has to be a 2D
  Array, the elements of which are indices referring to vertices in
  the vertex list (V) */
void geo3d::Trimesh::gridToTriangles(const uintA &grid){
  uint i,j,k=T.d0;
  T.resizeCopy(T.d0+2*(grid.d0-1)*(grid.d1-1),3);
  for(i=0;i<grid.d0-1;i++) for(j=0;j<grid.d1-1;j++){
    if((i+j)&1){
      T(k,0)=grid(i+1,j  );
      T(k,1)=grid(i  ,j  );
      T(k,2)=grid(i  ,j+1);
      k++;
      T(k,0)=grid(i+1,j  );
      T(k,1)=grid(i  ,j+1);
      T(k,2)=grid(i+1,j+1);
      k++;
    }else{
      T(k,0)=grid(i+1,j  );
      T(k,1)=grid(i  ,j  );
      T(k,2)=grid(i+1,j+1);
      k++;
      T(k,0)=grid(i+1,j+1);
      T(k,1)=grid(i  ,j  );
      T(k,2)=grid(i  ,j+1);
      k++;
    }
  }
}
  
/*!\brief add strips according to the given grid (sliced in strips along
  the x-axis (the first index)); grid has to be a 2D Array, the
  elements of which are indices referring to vertices in the vertex
  list (V) */
void geo3d::Trimesh::gridToStrips(const uintA& grid){
  CHECK(grid.d0>1 && grid.d1>1,"grid has to be at least 2x2");
  uint i,j,k=strips.N,l;
  strips.resizeCopy(strips.N+grid.d0-1);
  for(i=0;i<grid.d0-1;i++){
    strips(k).resize(2*grid.d1);
    l=0;
    for(j=0;j<grid.d1;j++){
      strips(k)(l)=grid(i+1,j); l++;
      strips(k)(l)=grid(i  ,j); l++;
    }
#if 0 //code to make it less symmetric
      //}else{
    strips(k)(l)=grid(i,0); l++;
    for(j=0;j<grid.d1;j++){
      strips(k)(l)=grid(i  ,j); l++;
      strips(k)(l)=grid(i+1,j); l++;
    }
#endif
    k++;
  }
}
  
/*!\brief add strips according to the given grid (sliced in strips along
  the x-axis (the first index)); it is assumed that the vertices in
  the list V linearly correspond to points in the XxY grid */
void geo3d::Trimesh::gridToStrips(uint X,uint Y){
  CHECK(X>1 && Y>1,"grid has to be at least 2x2");
  uint i,j,k=strips.N,l;
  strips.resizeCopy(strips.N+Y-1);
  for(j=0;j<Y-1;j++){
    strips(k).resize(2*X);
    l=0;
    for(i=0;i<X;i++){
      strips(k)(l)=(j+1)*X+i;
      l++;
      strips(k)(l)=    j*X+i;
      l++;
    }
    k++;
  }
}
  
/*!\brief add triangles according to the given grid; grid has to be a 2D
  Array, the elements of which are indices referring to vertices in
  the vertex list (V) */
void geo3d::Trimesh::gridToTriangles(uint X,uint Y){
  CHECK(X>1 && Y>1,"grid has to be at least 2x2");
  uint i,j,k=T.d0;
  T.resizeCopy(k+(Y-1)*2*(X-1),3);
  for(j=0;j<Y-1;j++){
    for(i=0;i<X-1;i++){
      T(k,0)=j*X+i; T(k,1)=(j+1)*X+i; T(k,2)=(j+1)*X+(i+1);
      k++;
      T(k,0)=j*X+i; T(k,1)=(j+1)*X+(i+1); T(k,2)=j*X+(i+1);
      k++;
    }
  }
}
  
/*!\brief delete all void triangles (with vertex indices (0,0,0)) and void
  vertices (not used for triangles or strips) */
void geo3d::Trimesh::deleteUnused(){
  if(!V.N) return;
  MT::Permutation p;
  boolA u;
  uintA newT;
  arr newV;
  uint i,j,Nused;
  //find proper permutation of triangles
  p.init(T.d0);
  j=p.N-1;
  for(i=0;i<=j;i++) if(!(T(i,0)|T(i,1)|T(i,2))){ p.permute(i,j); j--; break; }
  Nused=j+1;
  //permute triangle list
  newT.resize(Nused,3);
  for(i=0;i<Nused;i++){ newT(i,0)=T(p(i),0); newT(i,1)=T(p(i),1); newT(i,2)=T(p(i),2); }
  T.resize(Nused,3);
  for(i=0;i<Nused;i++){ T(i,0)=newT(i,0); T(i,1)=newT(i,1); T(i,2)=newT(i,2); }
  //mark used vertices
  u.resize(V.d0);
  u=false;
  for(i=0;i<T.d0;i++){ u(T(i,0))=true; u(T(i,1))=true; u(T(i,2))=true; }
  for(i=0;i<strips.N;i++) for(j=0;j<strips(i).N;j++) u(strips(i)(j))=true;
  //find proper permutation
  p.init(V.d0);
  j=p.N-1;
  for(i=0;i<j;i++) if(!u(i)){ p.permute(i,j); j--; }
  Nused=j+1;
  //permute vertices and triangles/strips indices
  newV.resize(Nused,3);
  for(i=0;i<V.d0;i++) if(p(i)<Nused) newV[p(i)]()=V[i];
  for(i=0;i<T.d0;i++){ T(i,0)=p(T(i,0)); T(i,1)=p(T(i,1)); T(i,2)=p(T(i,2)); }
  for(i=0;i<strips.N;i++) for(j=0;j<strips(i).N;j++) strips(i)(j)=p(strips(i)(j));
  //reassign vertex array
  V=newV;
}
  
#ifdef MT_GL
void glColor(float *rgb);//{ glColor(rgb[0],rgb[1],rgb[2],1.); }
//! internal draw routine for OpenGL
void geo3d::Trimesh::glDraw(){
#if 0 //use OpenGL's Arrays for fast drawing...
  uint j;
  glEnableClientState(GL_VERTEX_ARRAY);
  glEnableClientState(GL_COLOR_ARRAY);
  glEnableClientState(GL_NORMAL_ARRAY);
  glVertexPointer(3,GL_DOUBLE,0,V.p );
  glColorPointer (3,GL_FLOAT ,0,C.p );
  glNormalPointer(  GL_DOUBLE,0,N.p);
    
  //triangles
  if(T.N) glDrawElements(GL_TRIANGLES,T.N,GL_UNSIGNED_INT,T.p);
  //strips
  for(j=0;j<strips.N;j++)
    glDrawElements(GL_TRIANGLE_STRIP,strips(j).N,GL_UNSIGNED_INT,strips(j).p);
#elif 1
  int g;
  uint v,t,i,j;
  Vector w;
  glShadeModel(GL_SMOOTH);							
  if(!GT.N){
    glBegin(GL_TRIANGLES);
    for(t=0;t<T.d0;t++){
      for(j=0;j<3;j++){
	v=T(t,j);
	if(G.N) g=G(v); else g=-1;
	w=&N(v,0);  if(g!=-1) w=GF(g)->r*w;  glNormal3dv(w.v);
	if(colored) glColor3fv(C(v));
        w=&V(v,0);  if(g!=-1) w=GF(g)->p+GF(g)->r*w;  glVertex3dv(w.v);
      }
    }
    glEnd();
  }else{
    //faces that belong to one group only
    for(g=0;g<(int)GT.N-1;g++){
      glPushMatrix();
      GF(g)->glLoadMatrix();
      glBegin(GL_TRIANGLES);
      for(i=0;i<GT(g).N;i++){
	t=GT(g)(i);
	for(j=0;j<3;j++){
	  v=T(t,j);
          glNormal3dv(&N(v,0));
          if(colored) glColor3fv(C(v));
          glVertex3dv(&V(v,0));
	}
      }
      glEnd();
      glPopMatrix();
    }
    //faces with vertices from different groups (transform each vertex individually)
    glBegin(GL_TRIANGLES);
    for(i=0;i<GT(GT.N-1).N;i++){
      t=GT(GT.N-1)(i);
      for(j=0;j<3;j++){
	v=T(t,j);
	g=G(v);
	w=&N(v,0);  if(g!=-1) w=GF(g)->r*w;  glNormal3dv(w.v);
	if(colored) glColor3fv(C(v));
        w=&V(v,0);  if(g!=-1) w=GF(g)->p+GF(g)->r*w;  glVertex3dv(w.v);
      }
    }
    glEnd();
  }
  for(j=0;j<strips.N;j++){
    glBegin(GL_TRIANGLE_STRIP);
    for(i=0;i<strips(j).N;i++){
      glNormal3dv(&N(strips(j)(i),0));
      if(colored) glColor3fv(C(strips(j)(i)));
      glVertex3dv(&V(strips(j)(i),0));
    }
    glEnd();
  }
#else
  glBegin(GL_TRIANGLES);
  glShadeModel(GL_SMOOTH);							
  for(int i=0;i<T.d0;i++){
      
    if(trinormals) glNormal3dv(Tn(i));
    if(!trinormals) glNormal3dv(Tn(i));//(N(Tni(i,0)));
    ////  if(colored) glColor3fv(C(T(i,0)));
    glVertex3dv(V(T(i,0)));
      
    if(trinormals) glNormal3dv(Tn(i));
    if(!trinormals) glNormal3dv(Tn(i));//glNormal3dv(N(Tni(i,1)));
    ////  if(colored) glColor3fv(C(T(i,1)));
    glVertex3dv(V(T(i,1)));
      
    if(trinormals) glNormal3dv(Tn(i));
    if(!trinormals) glNormal3dv(Tn(i));//glNormal3dv(N(Tni(i,2)));
    ////  if(colored) glColor3fv(C(T(i,0)));
    glVertex3dv(V(T(i,2)));
      
      
  }
  glEnd();
#endif
}
#endif
