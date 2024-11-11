
#include "gft_marker3.h"


namespace gft{


  namespace Scene32{


    int *GetMarkers(sScene32 *label, sAdjRel3 *A){
      int *S = NULL;
      int p, q, i, k, nseeds;
      Voxel u,v;
      nseeds = 0;
      for(p=0; p<label->n; p++){
	if(label->data[p]!=NIL) nseeds++;
      }
      S = (int *)malloc((nseeds+1)*sizeof(int));
      if(S == NULL) return NULL;
      
      if(A == NULL){
	S[0] = nseeds;
	i = 1;
	for(p=0; p<label->n; p++){
	  if(label->data[p]!=NIL){
	    S[i] = p;
	    i++;
	  }
	}
      }
      else{
	i = 1;
	for(p=0; p<label->n; p++){
	  if(label->data[p]==NIL) continue;
	  u.c.x = Scene32::GetAddressX(label, p);
	  u.c.y = Scene32::GetAddressY(label, p);
	  u.c.z = Scene32::GetAddressZ(label, p);	  
	  for (k=1; k < A->n; k++){
	    v.v = u.v + A->d[k].v;
	    if (Scene32::IsValidVoxel(label, v)){
	      q = Scene32::GetVoxelAddress(label, v);
	      if(label->data[q] != label->data[p]){
		S[i] = p;
		i++;
		break;
	      }
	    }
	  }
	}
	S = (int *)realloc(S, i*sizeof(int));
	if(S == NULL) return NULL;
	S[0] = i-1;
      }
      return S;
    }



    void SaveMarkers(char *filename,
		     int *S,
		     sScene32 *label){
      int i,p,lb,x,y,z,id;
      FILE *fp;
      fp = fopen(filename, "w");
      if (fp == NULL) {
	printf("Error writing seeds.\n");
	exit(1);
      }
      fprintf(fp, " %d\n", S[0]);
      for(i = 1; i <= S[0]; i++){
	p = S[i];
	lb = label->data[p];
	x = Scene32::GetAddressX(label, p);
	y = Scene32::GetAddressY(label, p);
	z = Scene32::GetAddressZ(label, p);
	id = 1;
	if(Scene32::IsValidVoxel(label, x, y, z))
	  fprintf(fp, " %d %d %d %d %d\n", x, y, z, id, lb);
      }
      fclose(fp);
    }


    
    int *GetVolumeBorder(sScene32 *scn){
      return GetVolumeBorder(scn, 0);
    }


    
    int *GetVolumeBorder(sScene32 *scn,
			 uchar exception_flags){
      int *B = NULL;
      int p,x,y,z,n;
      n = 2*scn->xsize*scn->ysize + 2*scn->xsize*scn->zsize + 2*scn->ysize*scn->zsize;
      B = (int *)calloc((n + 1), sizeof(int));
      n = 0;
      //------YZ Planes------------------
      if( !(exception_flags & EXCEPT_BORDER_MIN_X) ){
	x = 0;
	for(y = 0; y < scn->ysize; y++){
	  for(z = 0; z < scn->zsize; z++){
	    p = Scene32::GetVoxelAddress(scn, x, y, z);
	    n++;
	    B[n] = p;
	  }
	}
      }
      
      if( !(exception_flags & EXCEPT_BORDER_MAX_X) ){
	x = scn->xsize-1;
	for(y = 0; y < scn->ysize; y++){
	  for(z = 0; z < scn->zsize; z++){
	    p = Scene32::GetVoxelAddress(scn, x, y, z);
	    n++;
	    B[n] = p;
	  }
	}
      }
      //------XZ Planes------------------
      if( !(exception_flags & EXCEPT_BORDER_MIN_Y) ){
	y = 0;
	for(x = 1; x < scn->xsize-1; x++){
	  for(z = 0; z < scn->zsize; z++){
	    p = Scene32::GetVoxelAddress(scn, x, y, z);
	    n++;
	    B[n] = p;
	  }  
	}
      }
      
      if( !(exception_flags & EXCEPT_BORDER_MAX_Y) ){
	y = scn->ysize-1;
	for(x = 1; x < scn->xsize-1; x++){
	  for(z = 0; z < scn->zsize; z++){
	    p = Scene32::GetVoxelAddress(scn, x, y, z);
	    n++;
	    B[n] = p;
	  }  
	}
      }
      //------XY Planes------------------  
      if( !(exception_flags & EXCEPT_BORDER_MIN_Z) ){
	z = 0;
	for(x = 1; x < scn->xsize-1; x++){
	  for(y = 1; y < scn->ysize-1; y++){
	    p = Scene32::GetVoxelAddress(scn, x, y, z);
	    n++;
	    B[n] = p;
	  }
	}
      }
      
      if( !(exception_flags & EXCEPT_BORDER_MAX_Z) ){
	z = scn->zsize-1;
	for(x = 1; x < scn->xsize-1; x++){
	  for(y = 1; y < scn->ysize-1; y++){
	    p = Scene32::GetVoxelAddress(scn, x, y, z);
	    n++;
	    B[n] = p;
	  }  
	}
      }
  
      B[0] = n;
      B = (int *)realloc(B, (n + 1)*sizeof(int));
      return B;
    }
    

    
  } /*end Scene32 namespace*/



  namespace Scene8{


    int *Convert2Array(sScene8 *mask){
      int *A;
      int p, n;
      n = 0;
      for(p = 0; p < mask->n; p++)
	if(mask->data[p] > 0)
	  n++;
      A = (int *)calloc((n + 1), sizeof(int));
      A[0] = n;
      n = 0;
      for(p = 0; p < mask->n; p++){
	if(mask->data[p] > 0){
	  n++;
	  A[n] = p;
	}
      }
      return A;
    }
    

    int *GetMarkers(sScene8 *label, sAdjRel3 *A){
      int *S = NULL;
      int p, q, i, k, nseeds;
      Voxel u,v;
      nseeds = 0;
      for(p=0; p<label->n; p++){
	if(label->data[p]!=NIL) nseeds++;
      }
      S = (int *)malloc((nseeds+1)*sizeof(int));
      if(S == NULL) return NULL;
      
      if(A == NULL){
	S[0] = nseeds;
	i = 1;
	for(p=0; p<label->n; p++){
	  if(label->data[p]!=NIL){
	    S[i] = p;
	    i++;
	  }
	}
      }
      else{
	i = 1;
	for(p=0; p<label->n; p++){
	  if(label->data[p]==NIL) continue;
	  u.c.x = Scene8::GetAddressX(label, p);
	  u.c.y = Scene8::GetAddressY(label, p);
	  u.c.z = Scene8::GetAddressZ(label, p);	  
	  for (k=1; k < A->n; k++){
	    v.v = u.v + A->d[k].v;
	    if (Scene8::IsValidVoxel(label, v)){
	      q = Scene8::GetVoxelAddress(label, v);
	      if(label->data[q] != label->data[p]){
		S[i] = p;
		i++;
		break;
	      }
	    }
	  }
	}
	S = (int *)realloc(S, i*sizeof(int));
	if(S == NULL) return NULL;
	S[0] = i-1;
      }
      return S;
    }

    

    void SaveMarkers(char *filename,
		     int *S,
		     sScene8 *label){
      int i,p,lb,x,y,z,id;
      FILE *fp;
      fp = fopen(filename, "w");
      if (fp == NULL) {
	printf("Error writing seeds.\n");
	exit(1);
      }
      fprintf(fp, " %d\n", S[0]);
      for(i = 1; i <= S[0]; i++){
	p = S[i];
	lb = label->data[p];
	x = Scene8::GetAddressX(label, p);
	y = Scene8::GetAddressY(label, p);
	z = Scene8::GetAddressZ(label, p);
	id = 1;
	if(Scene8::IsValidVoxel(label, x, y, z))
	  fprintf(fp, " %d %d %d %d %d\n", x, y, z, id, lb);
      }
      fclose(fp);
    }


    
    int *GetVolumeBorder(sScene8 *scn){
      return GetVolumeBorder(scn, 0);
    }


    
    int *GetVolumeBorder(sScene8 *scn,
			 uchar exception_flags){
      int *B = NULL;
      int p,x,y,z,n;
      n = 2*scn->xsize*scn->ysize + 2*scn->xsize*scn->zsize + 2*scn->ysize*scn->zsize;
      B = (int *)calloc((n + 1), sizeof(int));
      n = 0;
      //------YZ Planes------------------
      if( !(exception_flags & EXCEPT_BORDER_MIN_X) ){
	x = 0;
	for(y = 0; y < scn->ysize; y++){
	  for(z = 0; z < scn->zsize; z++){
	    p = Scene8::GetVoxelAddress(scn, x, y, z);
	    n++;
	    B[n] = p;
	  }
	}
      }
      
      if( !(exception_flags & EXCEPT_BORDER_MAX_X) ){
	x = scn->xsize-1;
	for(y = 0; y < scn->ysize; y++){
	  for(z = 0; z < scn->zsize; z++){
	    p = Scene8::GetVoxelAddress(scn, x, y, z);
	    n++;
	    B[n] = p;
	  }
	}
      }
      //------XZ Planes------------------
      if( !(exception_flags & EXCEPT_BORDER_MIN_Y) ){
	y = 0;
	for(x = 1; x < scn->xsize-1; x++){
	  for(z = 0; z < scn->zsize; z++){
	    p = Scene8::GetVoxelAddress(scn, x, y, z);
	    n++;
	    B[n] = p;
	  }  
	}
      }
      
      if( !(exception_flags & EXCEPT_BORDER_MAX_Y) ){
	y = scn->ysize-1;
	for(x = 1; x < scn->xsize-1; x++){
	  for(z = 0; z < scn->zsize; z++){
	    p = Scene8::GetVoxelAddress(scn, x, y, z);
	    n++;
	    B[n] = p;
	  }  
	}
      }
      //------XY Planes------------------  
      if( !(exception_flags & EXCEPT_BORDER_MIN_Z) ){
	z = 0;
	for(x = 1; x < scn->xsize-1; x++){
	  for(y = 1; y < scn->ysize-1; y++){
	    p = Scene8::GetVoxelAddress(scn, x, y, z);
	    n++;
	    B[n] = p;
	  }
	}
      }
      
      if( !(exception_flags & EXCEPT_BORDER_MAX_Z) ){
	z = scn->zsize-1;
	for(x = 1; x < scn->xsize-1; x++){
	  for(y = 1; y < scn->ysize-1; y++){
	    p = Scene8::GetVoxelAddress(scn, x, y, z);
	    n++;
	    B[n] = p;
	  }  
	}
      }
  
      B[0] = n;
      B = (int *)realloc(B, (n + 1)*sizeof(int));
      return B;
    }
    
    
  } //end Scene8 namespace

  

} /*end gft namespace*/
