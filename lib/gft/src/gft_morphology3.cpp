
#include "gft_morphology3.h"

namespace gft{
  namespace Scene8{

    sScene8 *Dilate(sScene8 *scn, sAdjRel3 *A){
      sScene8 *dil=NULL;
      int p,q,n,max,i;
      Voxel u,v;
      n = scn->n;
      dil = gft::Scene8::Create(scn);
      for(p = 0; p < n; p++){
	u.c.x = gft::Scene8::GetAddressX(scn, p);
	u.c.y = gft::Scene8::GetAddressY(scn, p);
	u.c.z = gft::Scene8::GetAddressZ(scn, p);
	max = INT_MIN;
	for(i=0; i < A->n; i++){
	  v.v = u.v + A->d[i].v;
	  if(IsValidVoxel(scn, v)){
	    q = GetVoxelAddress(scn, v);
	    if(scn->data[q] > max)
	      max = scn->data[q];
	  }
	}
	dil->data[p] = max;
      }
      return(dil);
    }

    
    sScene8 *Erode(sScene8 *scn, sAdjRel3 *A){
      sScene8 *ero=NULL;
      int p,q,n,min,i;
      Voxel u,v;
      n = scn->n;
      ero = gft::Scene8::Create(scn);
      for(p = 0; p < n; p++){
	u.c.x = gft::Scene8::GetAddressX(scn, p);
	u.c.y = gft::Scene8::GetAddressY(scn, p);
	u.c.z = gft::Scene8::GetAddressZ(scn, p);
	min = INT_MAX;
	for(i=0; i < A->n; i++){
	  v.v = u.v + A->d[i].v;
	  if(IsValidVoxel(scn, v)){
	    q = GetVoxelAddress(scn, v);
	    if(scn->data[q] < min)
	      min = scn->data[q];
	  }
	}
	ero->data[p] = min;
      }
      return(ero);
    }


    /* It assumes that the next operation is a erosion */
    sScene8 *DilateBin(sScene8 *bin, sSet **seed, float radius){
      sScene8 *dil=NULL,*boundr=NULL,*ero=NULL;
      Voxel u,v,w;
      sScene32 *cost=NULL,*root;
      sPQueue32 *Q=NULL;
      int i,p,q,n,sz;
      int *sq=NULL,tmp=INT_MAX,dx,dy,dz;
      sAdjRel3 *A=NULL;
      float dist;

      /* Compute seeds */
      
      if (*seed == NULL) {
	A      = gft::AdjRel3::Spheric(1.0);
	ero    = Erode(bin, A);
	boundr = gft::Scene8::Create(bin);
	n = boundr->n;
	for(p = 0; p < n; p++){
	  boundr->data[p] = bin->data[p] -  ero->data[p];
	}
	gft::Scene8::Destroy(&ero);
	gft::AdjRel3::Destroy(&A);

	for (p=0; p < n; p++)
	  if (boundr->data[p]==1)
	    gft::Set::Insert(seed, p);
	gft::Scene8::Destroy(&boundr);
      }
      
      /* Dilate image */
      
      dil  = gft::Scene8::Clone(bin);
      dist = (radius*radius);
      A  = gft::AdjRel3::Spheric(1.5);
      n  = MAX( MAX(dil->xsize, dil->ysize), dil->zsize);
      sq = gft::AllocIntArray(n);
      for (i=0; i < n; i++)
	sq[i]=i*i;
      
      cost = gft::Scene32::Create(dil->xsize, dil->ysize, dil->zsize);
      root = gft::Scene32::Create(dil->xsize, dil->ysize, dil->zsize);
      gft::Scene32::Fill(cost, INT_MAX);
      n    = dil->n;
      sz   = gft::AdjRel3::GetFrameSize(A);
      Q    = gft::PQueue32::Create(2*sz*(sz+dil->xsize+dil->ysize+dil->zsize),
				   n, cost->data);
      
      while (*seed != NULL){
	p = gft::Set::Remove(seed);
	cost->data[p]=0;
	root->data[p]=p;
	gft::PQueue32::InsertElem(&Q, p);
      }
      
      while(!gft::PQueue32::IsEmpty(Q)) {
	p = gft::PQueue32::RemoveMinFIFO(Q);
	if (cost->data[p] <= dist){
	  
	  dil->data[p] = 1;

	  u.c.x = gft::Scene8::GetAddressX(dil, p);
	  u.c.y = gft::Scene8::GetAddressY(dil, p);
	  u.c.z = gft::Scene8::GetAddressZ(dil, p);

	  w.c.x = gft::Scene8::GetAddressX(dil, root->data[p]);
	  w.c.y = gft::Scene8::GetAddressY(dil, root->data[p]);
	  w.c.z = gft::Scene8::GetAddressZ(dil, root->data[p]);
	  
	  for (i=1; i < A->n; i++){
	    v.v = u.v + A->d[i].v;
	    if (IsValidVoxel(dil, v)){
	      q = GetVoxelAddress(dil, v);
	      if ((cost->data[p] < cost->data[q])&&(dil->data[q]==0)){
		dx  = abs(v.c.x-w.c.x);
		dy  = abs(v.c.y-w.c.y);
		dz  = abs(v.c.z-w.c.z);
		tmp = sq[dx] + sq[dy] + sq[dz];
		if (tmp < cost->data[q]){
		  if (cost->data[q] == INT_MAX){
		    cost->data[q] = tmp;
		    gft::PQueue32::InsertElem(&Q, q);
		  }
		  else
		    gft::PQueue32::UpdateElem(&Q, q, tmp);

		  root->data[q] = root->data[p];
		}
	      }
	    }
	  }
	} else {  /* Seeds for erosion */
	  gft::Set::Insert(seed, p);
	}
      }
      
      free(sq);
      gft::PQueue32::Destroy(&Q);
      gft::Scene32::Destroy(&root);
      gft::Scene32::Destroy(&cost);
      gft::AdjRel3::Destroy(&A);
      
      return(dil);
    }
    
    
    /* It assumes that the next operation is a dilation, but it may
       be an erosion if you remove comments below. */

    sScene8 *ErodeBin(sScene8 *bin, sSet **seed, float radius){
      sScene8 *ero=NULL,*boundr=NULL,*dil=NULL;
      Voxel u,v,w;
      sScene32 *cost=NULL,*root;
      sPQueue32 *Q=NULL;
      int i,p,q,n,sz;
      int *sq=NULL,tmp=INT_MAX,dx,dy,dz;
      sAdjRel3 *A=NULL;
      float dist;

      /* Compute seeds */
      
      if (*seed == NULL) {
	A      = gft::AdjRel3::Spheric(1.0);
	dil    = Dilate(bin, A);
	boundr = gft::Scene8::Create(bin);
	n = boundr->n;
	for(p = 0; p < n; p++){
	  boundr->data[p] = dil->data[p] -  bin->data[p];
	}
	gft::Scene8::Destroy(&dil);
	gft::AdjRel3::Destroy(&A);

	for (p=0; p < n; p++)
	  if (boundr->data[p]==1)
	    gft::Set::Insert(seed, p);
	gft::Scene8::Destroy(&boundr);
      }
      
      /* Erode image */
      
      ero  = gft::Scene8::Clone(bin);
      dist = (radius*radius);
      A  = gft::AdjRel3::Spheric(1.5);
      n  = MAX( MAX(ero->xsize, ero->ysize), ero->zsize);
      sq = gft::AllocIntArray(n);
      for (i=0; i < n; i++)
	sq[i]=i*i;
      
      cost = gft::Scene32::Create(ero->xsize, ero->ysize, ero->zsize);
      root = gft::Scene32::Create(ero->xsize, ero->ysize, ero->zsize);
      gft::Scene32::Fill(cost, INT_MAX);
      n    = ero->n;
      sz   = gft::AdjRel3::GetFrameSize(A);
      Q    = gft::PQueue32::Create(2*sz*(sz+ero->xsize+ero->ysize+ero->zsize),
				   n, cost->data);
      
      while (*seed != NULL){
	p = gft::Set::Remove(seed);
	cost->data[p]=0;
	root->data[p]=p;
	gft::PQueue32::InsertElem(&Q, p);
      }
      
      while(!gft::PQueue32::IsEmpty(Q)) {
	p = gft::PQueue32::RemoveMinFIFO(Q);
	if (cost->data[p] <= dist){
	  
	  ero->data[p] = 0;

	  /* Seeds for erosion if we wanted to compute sequences of erosions
	     
	     if (((sq[Dx->data[p]+1]+sq[Dy->data[p]]) > dist)||
	     ((sq[Dx->data[p]]+sq[Dy->data[p]+1]) > dist)){
	     InsertSet(seed,p);
	     }
	     
	  */
	  
	  u.c.x = gft::Scene8::GetAddressX(ero, p);
	  u.c.y = gft::Scene8::GetAddressY(ero, p);
	  u.c.z = gft::Scene8::GetAddressZ(ero, p);

	  w.c.x = gft::Scene8::GetAddressX(ero, root->data[p]);
	  w.c.y = gft::Scene8::GetAddressY(ero, root->data[p]);
	  w.c.z = gft::Scene8::GetAddressZ(ero, root->data[p]);
	  
	  for (i=1; i < A->n; i++){
	    v.v = u.v + A->d[i].v;
	    if (IsValidVoxel(ero, v)){
	      q = GetVoxelAddress(ero, v);
	      if ((cost->data[p] < cost->data[q])&&(ero->data[q]==1)){
		dx  = abs(v.c.x-w.c.x);
		dy  = abs(v.c.y-w.c.y);
		dz  = abs(v.c.z-w.c.z);
		tmp = sq[dx] + sq[dy] + sq[dz];
		if (tmp < cost->data[q]){
		  if (cost->data[q] == INT_MAX){
		    cost->data[q] = tmp;
		    gft::PQueue32::InsertElem(&Q, q);
		  }
		  else
		    gft::PQueue32::UpdateElem(&Q, q, tmp);

		  root->data[q] = root->data[p];
		}
	      }
	    }
	  }
	} else {  /* Seeds for dilation */
	  gft::Set::Insert(seed, p);
	}
      }
      
      free(sq);
      gft::PQueue32::Destroy(&Q);
      gft::Scene32::Destroy(&root);
      gft::Scene32::Destroy(&cost);
      gft::AdjRel3::Destroy(&A);
      
      return(ero);
    }



    /* It assumes that the next operation is a erosion */
    sScene8 *DilateBin(sScene8 *bin, sSet **seed, float radius, sAdjRel3 *A){
      sScene8 *dil=NULL,*boundr=NULL,*ero=NULL;
      Voxel u,v,w;
      sScene32 *cost=NULL,*root;
      sPQueue32 *Q=NULL;
      int i,p,q,n,sz;
      int *sq=NULL,tmp=INT_MAX,dx,dy,dz;
      sAdjRel3 *A6 = NULL;
      float dist;

      /* Compute seeds */
      
      if (*seed == NULL) {
	A6      = gft::AdjRel3::Spheric(1.0);
	ero    = Erode(bin, A6);
	boundr = gft::Scene8::Create(bin);
	n = boundr->n;
	for(p = 0; p < n; p++){
	  boundr->data[p] = bin->data[p] -  ero->data[p];
	}
	gft::Scene8::Destroy(&ero);
	gft::AdjRel3::Destroy(&A6);

	for (p=0; p < n; p++)
	  if (boundr->data[p]==1)
	    gft::Set::Insert(seed, p);
	gft::Scene8::Destroy(&boundr);
      }
      
      /* Dilate image */
      
      dil  = gft::Scene8::Clone(bin);
      dist = (radius*radius);
      n  = MAX( MAX(dil->xsize, dil->ysize), dil->zsize);
      sq = gft::AllocIntArray(n);
      for (i=0; i < n; i++)
	sq[i]=i*i;
      
      cost = gft::Scene32::Create(dil->xsize, dil->ysize, dil->zsize);
      root = gft::Scene32::Create(dil->xsize, dil->ysize, dil->zsize);
      gft::Scene32::Fill(cost, INT_MAX);
      n    = dil->n;
      sz   = gft::AdjRel3::GetFrameSize(A);
      Q    = gft::PQueue32::Create(2*sz*(sz+dil->xsize+dil->ysize+dil->zsize),
				   n, cost->data);
      
      while (*seed != NULL){
	p = gft::Set::Remove(seed);
	cost->data[p]=0;
	root->data[p]=p;
	gft::PQueue32::InsertElem(&Q, p);
      }
      
      while(!gft::PQueue32::IsEmpty(Q)) {
	p = gft::PQueue32::RemoveMinFIFO(Q);
	if (cost->data[p] <= dist){
	  
	  dil->data[p] = 1;

	  u.c.x = gft::Scene8::GetAddressX(dil, p);
	  u.c.y = gft::Scene8::GetAddressY(dil, p);
	  u.c.z = gft::Scene8::GetAddressZ(dil, p);

	  w.c.x = gft::Scene8::GetAddressX(dil, root->data[p]);
	  w.c.y = gft::Scene8::GetAddressY(dil, root->data[p]);
	  w.c.z = gft::Scene8::GetAddressZ(dil, root->data[p]);
	  
	  for (i=1; i < A->n; i++){
	    v.v = u.v + A->d[i].v;
	    if (IsValidVoxel(dil, v)){
	      q = GetVoxelAddress(dil, v);
	      if ((cost->data[p] < cost->data[q])&&(dil->data[q]==0)){
		dx  = abs(v.c.x-w.c.x);
		dy  = abs(v.c.y-w.c.y);
		dz  = abs(v.c.z-w.c.z);
		tmp = sq[dx] + sq[dy] + sq[dz];
		if (tmp < cost->data[q]){
		  if (cost->data[q] == INT_MAX){
		    cost->data[q] = tmp;
		    gft::PQueue32::InsertElem(&Q, q);
		  }
		  else
		    gft::PQueue32::UpdateElem(&Q, q, tmp);

		  root->data[q] = root->data[p];
		}
	      }
	    }
	  }
	} else {  /* Seeds for erosion */
	  gft::Set::Insert(seed, p);
	}
      }
      
      free(sq);
      gft::PQueue32::Destroy(&Q);
      gft::Scene32::Destroy(&root);
      gft::Scene32::Destroy(&cost);
     
      return(dil);
    }



    /* It assumes that the next operation is a dilation, but it may
       be an erosion if you remove comments below. */

    sScene8 *ErodeBin(sScene8 *bin, sSet **seed, float radius, sAdjRel3 *A){
      sScene8 *ero=NULL,*boundr=NULL,*dil=NULL;
      Voxel u,v,w;
      sScene32 *cost=NULL,*root;
      sPQueue32 *Q=NULL;
      int i,p,q,n,sz;
      int *sq=NULL,tmp=INT_MAX,dx,dy,dz;
      sAdjRel3 *A6=NULL;
      float dist;

      /* Compute seeds */
      
      if (*seed == NULL) {
	A6     = gft::AdjRel3::Spheric(1.0);
	dil    = Dilate(bin, A6);
	boundr = gft::Scene8::Create(bin);
	n = boundr->n;
	for(p = 0; p < n; p++){
	  boundr->data[p] = dil->data[p] -  bin->data[p];
	}
	gft::Scene8::Destroy(&dil);
	gft::AdjRel3::Destroy(&A6);

	for (p=0; p < n; p++)
	  if (boundr->data[p]==1)
	    gft::Set::Insert(seed, p);
	gft::Scene8::Destroy(&boundr);
      }
      
      /* Erode image */
      
      ero  = gft::Scene8::Clone(bin);
      dist = (radius*radius);
      n  = MAX( MAX(ero->xsize, ero->ysize), ero->zsize);
      sq = gft::AllocIntArray(n);
      for (i=0; i < n; i++)
	sq[i]=i*i;
      
      cost = gft::Scene32::Create(ero->xsize, ero->ysize, ero->zsize);
      root = gft::Scene32::Create(ero->xsize, ero->ysize, ero->zsize);
      gft::Scene32::Fill(cost, INT_MAX);
      n    = ero->n;
      sz   = gft::AdjRel3::GetFrameSize(A);
      Q    = gft::PQueue32::Create(2*sz*(sz+ero->xsize+ero->ysize+ero->zsize),
				   n, cost->data);
      
      while (*seed != NULL){
	p = gft::Set::Remove(seed);
	cost->data[p]=0;
	root->data[p]=p;
	gft::PQueue32::InsertElem(&Q, p);
      }
      
      while(!gft::PQueue32::IsEmpty(Q)) {
	p = gft::PQueue32::RemoveMinFIFO(Q);
	if (cost->data[p] <= dist){
	  
	  ero->data[p] = 0;

	  /* Seeds for erosion if we wanted to compute sequences of erosions
	     
	     if (((sq[Dx->data[p]+1]+sq[Dy->data[p]]) > dist)||
	     ((sq[Dx->data[p]]+sq[Dy->data[p]+1]) > dist)){
	     InsertSet(seed,p);
	     }
	     
	  */
	  
	  u.c.x = gft::Scene8::GetAddressX(ero, p);
	  u.c.y = gft::Scene8::GetAddressY(ero, p);
	  u.c.z = gft::Scene8::GetAddressZ(ero, p);

	  w.c.x = gft::Scene8::GetAddressX(ero, root->data[p]);
	  w.c.y = gft::Scene8::GetAddressY(ero, root->data[p]);
	  w.c.z = gft::Scene8::GetAddressZ(ero, root->data[p]);
	  
	  for (i=1; i < A->n; i++){
	    v.v = u.v + A->d[i].v;
	    if (IsValidVoxel(ero, v)){
	      q = GetVoxelAddress(ero, v);
	      if ((cost->data[p] < cost->data[q])&&(ero->data[q]==1)){
		dx  = abs(v.c.x-w.c.x);
		dy  = abs(v.c.y-w.c.y);
		dz  = abs(v.c.z-w.c.z);
		tmp = sq[dx] + sq[dy] + sq[dz];
		if (tmp < cost->data[q]){
		  if (cost->data[q] == INT_MAX){
		    cost->data[q] = tmp;
		    gft::PQueue32::InsertElem(&Q, q);
		  }
		  else
		    gft::PQueue32::UpdateElem(&Q, q, tmp);

		  root->data[q] = root->data[p];
		}
	      }
	    }
	  }
	} else {  /* Seeds for dilation */
	  gft::Set::Insert(seed, p);
	}
      }
      
      free(sq);
      gft::PQueue32::Destroy(&Q);
      gft::Scene32::Destroy(&root);
      gft::Scene32::Destroy(&cost);
      
      return(ero);
    }

    
    
    sScene8 *Open(sScene8 *scn, sAdjRel3 *A){
      sScene8 *open, *ero;
      ero  = Erode(scn, A);
      open = Dilate(ero, A);
      gft::Scene8::Destroy(&ero);
      return (open);
    }

    
    sScene8 *Close3(sScene8 *scn, sAdjRel3 *A){
      sScene8 *close, *dil;
      dil   = Dilate(scn, A);
      close = Erode(dil, A);
      gft::Scene8::Destroy(&dil);
      return (close);
    }
    

    sScene8 *CloseBin(sScene8 *bin, float radius){
      sScene8 *close=NULL,*dil=NULL;
      gft::sSet *seed=NULL;
      dil   = DilateBin(bin,&seed,radius);
      close = ErodeBin(dil,&seed,radius);
      gft::Scene8::Destroy(&dil);
      gft::Set::Destroy(&seed);
      return(close);
    }

    sScene8 *OpenBin(sScene8 *bin, float radius){
      sScene8 *open=NULL,*ero=NULL;
      gft::sSet *seed=NULL;
      ero   = ErodeBin(bin,&seed,radius);
      open  = DilateBin(ero,&seed,radius);
      gft::Scene8::Destroy(&ero);
      gft::Set::Destroy(&seed);
      return(open);
    }


    sScene8 *CloseBin(sScene8 *bin, float radius, sAdjRel3 *A){
      sScene8 *close=NULL,*dil=NULL;
      gft::sSet *seed=NULL;
      dil   = DilateBin(bin,&seed,radius,A);
      close = ErodeBin(dil,&seed,radius,A);
      gft::Scene8::Destroy(&dil);
      gft::Set::Destroy(&seed);
      return(close);
    }


    sScene8 *CloseHoles(sScene8 *scn){
      sScene8 *value;
      sAdjRel3 *A;
      int *B = NULL;
      A = AdjRel3::Spheric(1.0);
      B = Scene8::GetVolumeBorder(scn);
      value = CloseHolesFromSeeds(scn, A, B);
      free(B);
      AdjRel3::Destroy(&A);
      return value;
    }


    sScene8 *CloseHolesFromSeeds(sScene8 *scn,
				 sAdjRel3 *A,
				 int *S){
      sPQueue8 *Q=NULL;
      int i,j,p,q,n;
      uchar Imax,cst;
      sScene8 *value;
      Voxel u,v;
      
      if(S == NULL)
	return NULL;
      
      value = Scene8::Create(scn);
      Scene8::Fill(value, UCHAR_MAX);
      n = scn->n;
      Imax = Scene8::GetMaxVal(scn);
      Q = PQueue8::Create(Imax+2, n, value->data);
      
      for(i=1; i<=S[0]; i++){
	p = S[i];
	value->data[p] = scn->data[p];
	PQueue8::FastInsertElem(Q, p);
      }
      
      while(!PQueue8::IsEmpty(Q)) {
	p = PQueue8::FastRemoveMinFIFO(Q);
	u.c.x = Scene8::GetAddressX(scn, p);
	u.c.y = Scene8::GetAddressY(scn, p);
	u.c.z = Scene8::GetAddressZ(scn, p);	
    
	for(i=1; i<A->n; i++){
	  v.v = u.v + A->d[i].v;
	  if(Scene8::IsValidVoxel(scn, v)){
	    q = Scene8::GetVoxelAddress(scn,v);
	    if(Q->L.elem[q].color != BLACK){
	      cst = MAX(value->data[p], scn->data[q]);
	      
	      if(cst < value->data[q]){
		if(Q->L.elem[q].color == GRAY)
		  PQueue8::FastRemoveElem(Q, q);
		value->data[q] = cst;
		PQueue8::FastInsertElem(Q, q);
	      }
	    }
	  }
	}
      }
      PQueue8::Destroy(&Q);
      return value;
    }

    
    
    
  } /*end Scene8 namespace*/



  
  namespace Scene32{
    
    sScene32 *Dilate(sScene32 *scn, sAdjRel3 *A){
      sScene32 *dil=NULL;
      int p,q,n,max,i;
      Voxel u,v;
      n = scn->n;
      dil = gft::Scene32::Create(scn);
      for(p = 0; p < n; p++){
	u.c.x = gft::Scene32::GetAddressX(scn, p);
	u.c.y = gft::Scene32::GetAddressY(scn, p);
	u.c.z = gft::Scene32::GetAddressZ(scn, p);
	max = INT_MIN;
	for(i=0; i < A->n; i++){
	  v.v = u.v + A->d[i].v;
	  if(IsValidVoxel(scn, v)){
	    q = GetVoxelAddress(scn, v);
	    if(scn->data[q] > max)
	      max = scn->data[q];
	  }
	}
	dil->data[p] = max;
      }
      return(dil);
    }
    
    
    sScene32 *Erode(sScene32 *scn, sAdjRel3 *A){
      sScene32 *ero=NULL;
      int p,q,n,min,i;
      Voxel u,v;
      n = scn->n;
      ero = gft::Scene32::Create(scn);
      for(p = 0; p < n; p++){
	u.c.x = gft::Scene32::GetAddressX(scn, p);
	u.c.y = gft::Scene32::GetAddressY(scn, p);
	u.c.z = gft::Scene32::GetAddressZ(scn, p);
	min = INT_MAX;
	for(i=0; i < A->n; i++){
	  v.v = u.v + A->d[i].v;
	  if(IsValidVoxel(scn, v)){
	    q = GetVoxelAddress(scn, v);
	    if(scn->data[q] < min)
	      min = scn->data[q];
	  }
	}
	ero->data[p] = min;
      }
      return(ero);
    }



    sScene32 *Open(sScene32 *scn, sAdjRel3 *A){
      sScene32 *open, *ero;
      ero  = Erode(scn, A);
      open = Dilate(ero, A);
      gft::Scene32::Destroy(&ero);
      return (open);
    }

    
    sScene32 *Close3(sScene32 *scn, sAdjRel3 *A){
      sScene32 *close, *dil;
      dil   = Dilate(scn, A);
      close = Erode(dil, A);
      gft::Scene32::Destroy(&dil);
      return (close);
    }
   


    sScene32 *CloseHoles(sScene32 *scn){
      sScene32 *value;
      sAdjRel3 *A;
      int *B = NULL;
      A = AdjRel3::Spheric(1.0);
      B = Scene32::GetVolumeBorder(scn);
      value = CloseHolesFromSeeds(scn, A, B);
      free(B);
      AdjRel3::Destroy(&A);
      return value;
    }
    

    sScene32 *CloseHolesSliceBySlice_z(sScene32 *scn){
      sScene32 *fscn;
      sImage32 *img, *fimg;
      int z;
      fscn = Scene32::Create(scn);
      for(z = 0; z < scn->zsize; z++){
	img  = Scene32::GetSliceZ(scn, z);
	fimg = Image32::CloseHoles(img);
	Scene32::PutSliceZ(fscn, fimg, z);
	Image32::Destroy(&fimg);
	Image32::Destroy(&img);
      }
      return fscn;
    }

    

    sScene32 *CloseHolesFromSeeds(sScene32 *scn,
				  sAdjRel3 *A,
				  int *S){
      sPQueue32 *Q=NULL;
      int i,j,p,q,n,Imax,cst;
      sScene32 *value;
      Voxel u,v;
      
      if(S == NULL)
	return NULL;

      value = Scene32::Create(scn);
      Scene32::Fill(value, INT_MAX);
      n = scn->n;
      Imax = Scene32::GetMaxVal(scn);
      Q = PQueue32::Create(Imax+2, n, value->data);
      
      for(i=1; i<=S[0]; i++){
	p = S[i];
	value->data[p] = scn->data[p];
	PQueue32::FastInsertElem(Q, p);
      }
      
      while(!PQueue32::IsEmpty(Q)) {
	p = PQueue32::FastRemoveMinFIFO(Q);
	u.c.x = Scene32::GetAddressX(scn, p);
	u.c.y = Scene32::GetAddressY(scn, p);
	u.c.z = Scene32::GetAddressZ(scn, p);	
	
	for(i=1; i<A->n; i++){
	  v.v = u.v + A->d[i].v;
	  if(Scene32::IsValidVoxel(scn, v)){
	    q = Scene32::GetVoxelAddress(scn,v);
	    if(Q->L.elem[q].color != BLACK){
	      cst = MAX(value->data[p], scn->data[q]);
	      
	      if(cst < value->data[q]){
		if(Q->L.elem[q].color == GRAY)
		  PQueue32::FastRemoveElem(Q, q);
		value->data[q] = cst;
		PQueue32::FastInsertElem(Q, q);
	      }
	    }
	  }
	}
      }
      PQueue32::Destroy(&Q);
      return value;
    }


    
  } //end Scene32 namespace
  

} /*end gft namespace*/

