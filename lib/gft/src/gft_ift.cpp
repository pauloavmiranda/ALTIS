
#include "gft_ift.h"
#include <queue>

namespace gft{
  namespace ift{

    /* P_sum = Predecessor map obtained by the IFT fsum.*/
    void SC_conquer_path(int p,
			 sImageGraph *sg,
			 sImage32 *P_sum, 
			 sImage32 *V,
			 sPQueue32 *Q,
			 sImage32 *label);
    
    /* P_sum = Predecessor map obtained by the IFT fsum.*/
    void SC_prune_tree(int p,
		       sImageGraph *sg,
		       sImage32 *P_sum, 
		       sImage32 *V,
		       sPQueue32 *Q,
		       sQueue *Qfifo,
		       sImage32 *label);
    
    sSet *COIFT_new_seeds(sImageGraph *sg,
			  sSet *Si1,
			  sSet *Si2,
			  sImage32 *E,
			  sImage32 *label_oift);

    //View Boundary Band
    /*
    Image32::Image32 *BB_OIFT_View(ImageGraph::ImageGraph *G,
				   Image32::Image32 *L,
				   Image32::Image32 *C,
				   float delta,
				   int band_id);
    */ 
   
    void BB_OIFT_Propagate(int p,
			   sImageGraph *sg,
			   sPQueue32 *Q,
			   sImage32 *V,
			   sImage32 *L,
			   int *i_inv);
    
    void BB_OIFT_Propagate_bkg(int p,
			       sImageGraph *sg,
			       sPQueue32 *Q,
			       sPQueue32 *Qe,
			       sImage32 *V,
			       sImage32 *L,
			       int *i_inv);
    
    void BB_OIFT_Propagate_obj(int p,
			       sImageGraph *sg,
			       sPQueue32 *Q,
			       sPQueue32 *Qi,
			       sImage32 *V,
			       sImage32 *L,
			       int *i_inv);

    sImage32 *BB_OIFT_GetLeafNodes(sImage32 *pred,
				   sAdjRel *A);
    
    //----------------------------------

    void IFT_fmax(sImageGraph *g,
		  int *S,
		  sImage32 *label,
		  sImage32 *cost,
		  sImage32 *pred){
      sPQueue32 *Q=NULL;
      int i,p,q,n,cst;
      Pixel u,v;
      sAdjRel *A = g->A;

      n = g->ncols*g->nrows;
      Q = PQueue32::Create(2+g->Wmax, n, cost->data);

      Image32::Set(pred, NIL);
      for(p = 0; p < n; p++){
	if(label->data[p]==NIL) cost->data[p] = INT_MAX;
	else                    cost->data[p] = 0;
      }
      
      if(S != NULL){
	for(i = 1; i <= S[0]; i++)
	  PQueue32::FastInsertElem(Q, S[i]);
      }
      else{
	for(p=0; p<n; p++)
	  if(label->data[p]!=NIL)
	    PQueue32::FastInsertElem(Q, p);	    
      }
      
      while(!PQueue32::IsEmpty(Q)) {
	p = PQueue32::FastRemoveMinFIFO(Q);
	
	u.x = p%label->ncols;
	u.y = p/label->ncols;
	for (i=1; i < A->n; i++){
	  v.x = u.x + A->dx[i];
	  v.y = u.y + A->dy[i];
	  //if (gft::Image32::IsValidPixel(label, v.x, v.y)){
	  if(v.x >= 0 && v.x < label->ncols &&
	     v.y >= 0 && v.y < label->nrows){	  
	    q = v.x + v.y*label->ncols;

	    cst = MAX( (g->n_link[p])[i], cost->data[p] );
	    
	    if(cst < cost->data[q]){	    
	      if(Q->L.elem[q].color == GRAY)
		gft::PQueue32::FastRemoveElem(Q, q);
	      cost->data[q]  = cst;
	      pred->data[q]  = p;
	      label->data[q] = label->data[p];
	      gft::PQueue32::FastInsertElem(Q, q);
	    }
	  }
	}
      }
      gft::PQueue32::Destroy(&Q);
    }

    

    void IFT_feuc(sImageGraph *g,
		  int *S,
		  sImage32 *label,
		  sImage32 *cost,
		  sImage32 *pred){
      sImage32 *Dx=NULL,*Dy=NULL;
      sPQueue32 *Q=NULL;
      int i,p,q,n,cst,dx,dy;
      Pixel u,v;
      sAdjRel *A = g->A;

      Dx = Image32::Create(cost->ncols, cost->nrows);
      Dy = Image32::Create(cost->ncols, cost->nrows);
      
      n = g->ncols*g->nrows;
      Q = PQueue32::Create(2*(label->ncols+label->nrows), n, cost->data);

      Image32::Set(pred, NIL);
      for(p = 0; p < n; p++){
	if(label->data[p]==NIL) cost->data[p] = INT_MAX;
	else                    cost->data[p] = 0;
      }
      
      if(S != NULL){
	for(i = 1; i <= S[0]; i++)
	  PQueue32::FastInsertElem(Q, S[i]);
      }
      else{
	for(p=0; p<n; p++)
	  if(label->data[p]!=NIL)
	    PQueue32::FastInsertElem(Q, p);	    
      }
      
      while(!PQueue32::IsEmpty(Q)) {
	p = PQueue32::FastRemoveMinFIFO(Q);
	
	u.x = p%label->ncols;
	u.y = p/label->ncols;
	for (i=1; i < A->n; i++){
	  v.x = u.x + A->dx[i];
	  v.y = u.y + A->dy[i];
	  //if (gft::Image32::IsValidPixel(label, v.x, v.y)){
	  if(v.x >= 0 && v.x < label->ncols &&
	     v.y >= 0 && v.y < label->nrows){	  
	    q = v.x + v.y*label->ncols;

	    dx  = Dx->data[p] + abs(A->dx[i]);
	    dy  = Dy->data[p] + abs(A->dy[i]);
	    cst = dx*dx + dy*dy; //(g->n_link[p])[i];
	    
	    if(cst < cost->data[q]){	    
	      if(Q->L.elem[q].color == GRAY)
		gft::PQueue32::FastRemoveElem(Q, q);
	      Dx->data[q] = dx;
	      Dy->data[q] = dy;
	      cost->data[q]  = cst;
	      pred->data[q]  = p;
	      label->data[q] = label->data[p];
	      gft::PQueue32::FastInsertElem(Q, q);
	    }
	  }
	}
      }
      gft::PQueue32::Destroy(&Q);
      Image32::Destroy(&Dx);
      Image32::Destroy(&Dy);      
    }


    
    void IFT_fw(sImageGraph *g,
		int *S,
		sImage32 *label,
		sImage32 *cost,
		sImage32 *pred){
      sPQueue32 *Q=NULL;
      int i,p,q,n,cst;
      Pixel u,v;
      sAdjRel *A = g->A;
      
      n = g->ncols*g->nrows;
      Q = PQueue32::Create(g->Wmax+2, n, cost->data);

      Image32::Set(pred, NIL);
      for(p = 0; p < n; p++){
	if(label->data[p]==NIL) cost->data[p] = INT_MAX;
	else                    cost->data[p] = 0;
      }
      
      if(S != NULL){
	for(i = 1; i <= S[0]; i++)
	  PQueue32::FastInsertElem(Q, S[i]);
      }
      else{
	for(p=0; p<n; p++)
	  if(label->data[p]!=NIL)
	    PQueue32::FastInsertElem(Q, p);	    
      }
      
      while(!PQueue32::IsEmpty(Q)) {
	p = PQueue32::FastRemoveMinFIFO(Q);
	
	u.x = p%label->ncols;
	u.y = p/label->ncols;
	for (i=1; i < A->n; i++){
	  v.x = u.x + A->dx[i];
	  v.y = u.y + A->dy[i];
	  //if (gft::Image32::IsValidPixel(label, v.x, v.y)){
	  if(v.x >= 0 && v.x < label->ncols &&
	     v.y >= 0 && v.y < label->nrows){	  
	    q = v.x + v.y*label->ncols;
	    
	    cst = (g->n_link[p])[i];
	    
	    if(cst < cost->data[q]){	    
	      if(Q->L.elem[q].color == GRAY)
		gft::PQueue32::FastRemoveElem(Q, q);
	      cost->data[q]  = cst;
	      pred->data[q]  = p;
	      label->data[q] = label->data[p];
	      gft::PQueue32::FastInsertElem(Q, q);
	    }
	  }
	}
      }
      gft::PQueue32::Destroy(&Q);
    }
    
    //----------------------------------
    
    int GetEnergy_Min(sImageGraph *sg,
		      sImage32 *label,
		      int lb){
      sAdjRel *A;
      int u_x,u_y,v_x,v_y,p,q,n,i;
      int energy,w;

      energy = INT_MAX;
      A = sg->A;
      n = sg->ncols*sg->nrows;
      for(p = 0; p < n; p++){
	u_x = p%sg->ncols; 
	u_y = p/sg->ncols; 
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  if(v_x >= 0 && v_x < sg->ncols &&
	     v_y >= 0 && v_y < sg->nrows){
	    q = v_x + sg->ncols*v_y;
	    if(label->data[p] == lb && label->data[q] != lb){
	      w = (sg->n_link[p])[i];
	      energy = MIN(energy, w);
	    }
	  }
	}
      }
      return energy;
    }


    int GetEnergy_Max(sImageGraph *sg,
		      sImage32 *label,
		      int lb){
      sAdjRel *A;
      int u_x,u_y,v_x,v_y,p,q,n,i;
      int energy,w;

      energy = INT_MIN;
      A = sg->A;
      n = sg->ncols*sg->nrows;
      for(p = 0; p < n; p++){
	u_x = p%sg->ncols; 
	u_y = p/sg->ncols; 
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  if(v_x >= 0 && v_x < sg->ncols &&
	     v_y >= 0 && v_y < sg->nrows){
	    q = v_x + sg->ncols*v_y;
	    if(label->data[p] == lb && label->data[q] != lb){
	      w = (sg->n_link[p])[i];
	      energy = MAX(energy, w);
	    }
	  }
	}
      }
      return energy;
    }
    

    long long GetEnergy_Sum(sImageGraph *sg,
			    sImage32 *label,
			    int lb){
      sAdjRel *A;
      int u_x,u_y,v_x,v_y,p,q,n,i;
      long long energy,w;

      energy = 0;
      A = sg->A;
      n = sg->ncols*sg->nrows;
      for(p = 0; p < n; p++){
	u_x = p%sg->ncols; 
	u_y = p/sg->ncols; 
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  if(v_x >= 0 && v_x < sg->ncols &&
	     v_y >= 0 && v_y < sg->nrows){
	    q = v_x + sg->ncols*v_y;
	    if(label->data[p] == lb && label->data[q] != lb){
	      w = (sg->n_link[p])[i];
	      energy += w;
	    }
	  }
	}
      }
      return energy;
    }



    int GetEnergy_Min(sGraph *graph,
		      int *label,
		      int lb){
      int p,q,i;
      int energy,w;

      energy = INT_MAX;
      for(p = 0; p < graph->nnodes; p++){
	for(i = 0; i < graph->nodes[p].outdegree; i++){
	  q = graph->nodes[p].adjList[i];
	  w = graph->nodes[p].Warcs[i];
	  if(label[p] == lb && label[q] != lb){
	    energy = MIN(energy, w);
	  }
	}
      }
      return energy;
    }



    float GetEnergy_Mean(sGraph *graph,
			 int *label,
			 int lb){
      int p,q,i;
      int w, n = 0;
      float sum = 0.0;
      for(p = 0; p < graph->nnodes; p++){
	for(i = 0; i < graph->nodes[p].outdegree; i++){
	  q = graph->nodes[p].adjList[i];
	  w = graph->nodes[p].Warcs[i];
	  if(label[p] == lb && label[q] != lb){
	    sum += w;
	    n++;
	  }
	}
      }
      return (sum/n);
    }

    
    
    
    /*Weighted Distance Transform.*/
    sImage32 *SC_Pred_fsum(sImageGraph *sg,
			   int *S,
			   float power){
      sHeap *Q=NULL;
      int i,p,q,n;
      float edge,tmp;
      float *cost=NULL;
      int u_x,u_y,v_x,v_y;
      sImage32 *pred;
      sAdjRel *A;
      float *Dpq;
      
      n    = sg->ncols*sg->nrows;
      pred = Image32::Create(sg->ncols, sg->nrows);
      cost = gft::AllocFloatArray(n);
      Q = Heap::Create(n, cost);
      A = sg->A;

      //--------------------
      Dpq = (float *)malloc(A->n*sizeof(float));
      for(i=1; i<A->n; i++){
	Dpq[i] = sqrtf(A->dx[i]*A->dx[i] + A->dy[i]*A->dy[i]);
      }
      //--------------------
      
      Image32::Set(pred, NIL);
      for(p = 0; p < n; p++)
	cost[p] = FLT_MAX;
      
      for(i=1; i<=S[0]; i++){
	cost[S[i]] = 0.0;
	Heap::Insert_MinPolicy(Q, S[i]);
      }
	
      while(!Heap::IsEmpty(Q)){
	Heap::Remove_MinPolicy(Q, &p);
	u_x = p%sg->ncols; 
	u_y = p/sg->ncols; 
	
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  if(v_x >= 0 && v_x < sg->ncols &&
	     v_y >= 0 && v_y < sg->nrows){
	    q = v_x + sg->ncols*v_y;
	    if(Q->color[q] != BLACK){
	      
	      edge = (sg->n_link[p])[i];
	      tmp  = cost[p] + powf(MAX(edge,1.0), power) - 1.0 + Dpq[i];
	      
	      if(tmp < cost[q]){
		Heap::Update_MinPolicy(Q, q, tmp);
		pred->data[q] = p;
	      }
	    }
	  }
	}
      }
      free(Dpq);
      gft::FreeFloatArray(&cost);
      Heap::Destroy(&Q);
      return pred;
    }


    
    sScene32 *SC_Pred_fsum(sGraph *graph,
			   sScene32 *scn,
			   int *S,
			   float power){
      sHeap *Q=NULL;
      int i,p,q,n;
      float edge,tmp;
      float *cost=NULL;
      int u_x,u_y,u_z,v_x,v_y,v_z;
      gft::sScene32 *pred;
      float Dpq,dx2,dy2,dz2;
      dx2 = scn->dx*scn->dx;
      dy2 = scn->dy*scn->dy;
      dz2 = scn->dz*scn->dz;
      n    = graph->nnodes;
      pred = gft::Scene32::Create(scn);
      cost = gft::AllocFloatArray(n);
      Q = Heap::Create(n, cost);
      
      Scene32::Fill(pred, NIL);
      for(p = 0; p < n; p++)
	cost[p] = FLT_MAX;
      
      for(i=1; i<=S[0]; i++){
	cost[S[i]] = 0.0;
	Heap::Insert_MinPolicy(Q, S[i]);
      }
	
      while(!Heap::IsEmpty(Q)){
	Heap::Remove_MinPolicy(Q, &p);
	u_x = gft::Scene32::GetAddressX(scn, p);
	u_y = gft::Scene32::GetAddressY(scn, p);
	u_z = gft::Scene32::GetAddressZ(scn, p);	
	for(i = 0; i < graph->nodes[p].outdegree; i++){
	  q = graph->nodes[p].adjList[i];
	  if(Q->color[q] != BLACK){
	    v_x = gft::Scene32::GetAddressX(scn, q);
	    v_y = gft::Scene32::GetAddressY(scn, q);
	    v_z = gft::Scene32::GetAddressZ(scn, q);
	    Dpq = sqrtf((u_x-v_x)*(u_x-v_x)*dx2 + 
			(u_y-v_y)*(u_y-v_y)*dy2 + 
			(u_z-v_z)*(u_z-v_z)*dz2);
	    edge = graph->nodes[p].Warcs[i];
	    tmp  = cost[p] + powf(MAX(edge,1.0), power) - 1.0 + Dpq;
	    
	    if(tmp < cost[q]){
	      Heap::Update_MinPolicy(Q, q, tmp);
	      pred->data[q] = p;
	    }
	  }
	}
      }
      gft::FreeFloatArray(&cost);
      Heap::Destroy(&Q);
      return pred;
    }
    


    sScene32 *SC_Pred_fsum(sScene32 *W,
			   sAdjRel3 *A,
			   int *S, float power){
      sHeap *Q=NULL;
      int i,p,q,n;
      float edge,tmp;
      float *cost=NULL;
      gft::Voxel u,v;
      sScene32 *pred;
      float *Dpq;
      
      n    = W->n;
      pred = Scene32::Create(W);
      cost = gft::AllocFloatArray(n);
      Q = Heap::Create(n, cost);

      //--------------------
      Dpq = (float *)malloc(A->n*sizeof(float));
      for(i=1; i<A->n; i++){
	Dpq[i] = sqrtf(A->d[i].axis.x*A->d[i].axis.x*W->dx*W->dx + 
		       A->d[i].axis.y*A->d[i].axis.y*W->dy*W->dy + 
		       A->d[i].axis.z*A->d[i].axis.z*W->dz*W->dz);
      }
      //--------------------
      
      Scene32::Fill(pred, NIL);
      for(p = 0; p < n; p++)
	cost[p] = FLT_MAX;
      
      for(i=1; i<=S[0]; i++){
	cost[S[i]] = 0.0;
	Heap::Insert_MinPolicy(Q, S[i]);
      }
	
      while(!Heap::IsEmpty(Q)){
	Heap::Remove_MinPolicy(Q, &p);
	u.c.x = gft::Scene32::GetAddressX(W, p);
	u.c.y = gft::Scene32::GetAddressY(W, p);
	u.c.z = gft::Scene32::GetAddressZ(W, p);
	for(i=1; i<A->n; i++){
	  v.v = u.v + A->d[i].v;
	  if(gft::Scene32::IsValidVoxel(W,v)){
	    q = gft::Scene32::GetVoxelAddress(W, v);
	    if(Q->color[q] != BLACK){
	      edge = W->data[p] + W->data[q];
	      tmp  = cost[p] + powf(MAX(edge,1.0), power) - 1.0 + Dpq[i];
	      
	      if(tmp < cost[q]){
		Heap::Update_MinPolicy(Q, q, tmp);
		pred->data[q] = p;
	      }
	    }
	  }
	}
      }
      free(Dpq);
      gft::FreeFloatArray(&cost);
      Heap::Destroy(&Q);
      return pred;
    }
    

    //---------------------------------------

    void ORFC(sImageGraph *sg,
	      int *S,
	      sImage32 *label){
      gft::sImage32 *tmp,*value_e;
      gft::sPQueue32 *QS=NULL;
      gft::sQueue *Q = gft::Queue::Create(sg->ncols*sg->nrows);
      gft::sAdjRel *A = sg->A;
      int p,q,j,i,k,energy,nsi;
      gft::Pixel u,v;
      int *s_energy,*s_pixel;

      nsi = 0;
      for(i = 1; i <= S[0]; i++){
	p = S[i];
	if(label->data[p] > 0)
	  nsi++;
      }
      s_energy = gft::AllocIntArray(nsi);
      s_pixel  = gft::AllocIntArray(nsi);
      QS = gft::PQueue32::Create(sg->Wmax+2, nsi, s_energy);
  
      gft::ImageGraph::Transpose(sg);
      
      tmp = gft::Image32::Clone(label);
      value_e = gft::ift::Cost_fmin(sg, S, 0, tmp);
      gft::Image32::Destroy(&tmp);

      gft::ImageGraph::Transpose(sg);
      
      k = 0;
      for(i = 1; i <= S[0]; i++){
	p = S[i];
	if(label->data[p] > 0){
	  s_energy[k] = value_e->data[p];
	  s_pixel[k] = p;
	  gft::PQueue32::InsertElem(&QS, k);
	  k++;
	  label->data[p] = NIL;
	}
      }

      while(!gft::PQueue32::IsEmpty(QS)){
	j = gft::PQueue32::RemoveMinFIFO(QS);
	p = s_pixel[j];
	energy = s_energy[j];
	
	if(label->data[p] != NIL) continue;
	
	gft::Queue::Reset(Q);
	gft::Queue::Push(Q, p);
	label->data[p] = 1;
	
	while(!gft::Queue::IsEmpty(Q)){
	  p = gft::Queue::Pop(Q);
	  u.x = p%sg->ncols;
	  u.y = p/sg->ncols;      
	  
	  for(i = 1; i < A->n; i++){
	    v.x = u.x + A->dx[i];
	    v.y = u.y + A->dy[i];
	    //if(gft::Image32::IsValidPixel(label, v.x, v.y)){
	    if(v.x >= 0 && v.x < label->ncols &&
	       v.y >= 0 && v.y < label->nrows){	    
	      q = v.x + v.y*sg->ncols;
	      if(energy < (sg->n_link[p])[i] &&
		 label->data[q] == NIL){
		label->data[q] = 1;
		gft::Queue::Push(Q, q);
	      }
	    }
	  }
	}
      }
      for(p=0; p<sg->ncols*sg->nrows; p++){
	if(label->data[p] == NIL)
	  label->data[p] = 0;
      }

      gft::Queue::Destroy(&Q);
      gft::Image32::Destroy(&value_e);
      gft::PQueue32::Destroy(&QS);
      gft::FreeIntArray(&s_energy);
      gft::FreeIntArray(&s_pixel);      
    }


    //---------------------------------------

    //Outer Cut:
    void OIFT(sImage32 *W,
	      sAdjRel *A,
	      sImage32 *img,
	      float per,
	      int *S,
	      sImage32 *label){
      sPQueue32 *Q=NULL;
      int i,j,p,q,n;
      int w,Wmax;
      sImage32 *value;
      int u_x,u_y,v_x,v_y;
      float per_pq;
      
      value = Image32::Create(W->ncols,
			      W->nrows);
      n = label->n;
      Wmax = gft::Image32::GetMaxVal(W)*2;
      Wmax *= (1.0 + fabsf(per)/100.0);     
      Q = PQueue32::Create(Wmax+2,n,value->data);

      for(p=0; p<n; p++){
	if(label->data[p]==NIL) value->data[p] = INT_MAX;
	else                    value->data[p] = 0;
      }
      
      if(S != NULL){
	for(i=1; i<=S[0]; i++)
	  PQueue32::FastInsertElem(Q, S[i]);
      }
      else{
	for(p=0; p<n; p++)
	  if(label->data[p]!=NIL)
	    PQueue32::FastInsertElem(Q, p);	    
      }

      while(!PQueue32::IsEmpty(Q)) {
	p = PQueue32::FastRemoveMinFIFO(Q);
	u_x = p%label->ncols; //PixelX(label, p);
	u_y = p/label->ncols; //PixelY(label, p);
	
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  //if(Image32::IsValidPixel(label,v_x,v_y)){
	  if(v_x >= 0 && v_x < label->ncols &&
	     v_y >= 0 && v_y < label->nrows){
	    q = v_x + label->ncols*v_y;
	    if(Q->L.elem[q].color != BLACK){

	      w = W->data[p] + W->data[q];
	      if(label->data[p] > 0)
		per_pq = per;
	      else
		per_pq = -per;

	      if(img->data[p] > img->data[q])
		w *= (1.0 + per_pq/100.0);
	      else if(img->data[p] < img->data[q])
		w *= (1.0 - per_pq/100.0);
	      
	      if(w < value->data[q]){
		if(Q->L.elem[q].color == GRAY)
		  PQueue32::FastRemoveElem(Q, q);
		value->data[q] = w;
		label->data[q] = label->data[p];
		PQueue32::FastInsertElem(Q, q);
	      }
	    }
	  }
	}
      }
      Image32::Destroy(&value);
      PQueue32::Destroy(&Q);
    }




    //Outer Cut:
    void OIFT(sAdjRel3 *A,
	      sScene32 *scn,
	      float per,
	      int *S,
	      sScene32 *label){
      sPQueue32 *Q=NULL;
      int i,j,p,q,n;
      int w,Wmax;
      sScene32 *value;
      gft::Voxel u,v;
      float per_pq;
      
      value = Scene32::Create(scn);
      n = label->n;
      Wmax = gft::Scene32::GetMaximumValue(scn);
      Wmax *= (1.0 + fabsf(per)/100.0);     
      Q = PQueue32::Create(Wmax+2,n,value->data);

      for(p=0; p<n; p++){
	if(label->data[p]==NIL) value->data[p] = INT_MAX;
	else                    value->data[p] = 0;
      }
      
      if(S != NULL){
	for(i=1; i<=S[0]; i++)
	  PQueue32::FastInsertElem(Q, S[i]);
      }
      else{
	for(p=0; p<n; p++)
	  if(label->data[p]!=NIL)
	    PQueue32::FastInsertElem(Q, p);	    
      }

      while(!PQueue32::IsEmpty(Q)) {
	p = PQueue32::FastRemoveMinFIFO(Q);
	u.c.x = gft::Scene32::GetAddressX(label, p);
	u.c.y = gft::Scene32::GetAddressY(label, p);
	u.c.z = gft::Scene32::GetAddressZ(label, p);	
	
	for(i=1; i<A->n; i++){
	  v.v = u.v + A->d[i].v;
	  if(gft::Scene32::IsValidVoxel(label, v)){
	    q = gft::Scene32::GetVoxelAddress(label,v);
	    if(Q->L.elem[q].color != BLACK){

	      w = abs(scn->data[p] - scn->data[q]);
	      if(label->data[p] > 0)
		per_pq = per;
	      else
		per_pq = -per;

	      if(scn->data[p] > scn->data[q])
		w *= (1.0 + per_pq/100.0);
	      else if(scn->data[p] < scn->data[q])
		w *= (1.0 - per_pq/100.0);
	      
	      if(w < value->data[q]){
		if(Q->L.elem[q].color == GRAY)
		  PQueue32::FastRemoveElem(Q, q);
		value->data[q] = w;
		label->data[q] = label->data[p];
		PQueue32::FastInsertElem(Q, q);
	      }
	    }
	  }
	}
      }
      Scene32::Destroy(&value);
      PQueue32::Destroy(&Q);
    }

    
    
    
    //Outer Cut:
    void OIFT(sScene32 *W,
	      sAdjRel3 *A,
	      sScene32 *scn,
	      float per,
	      int *S,
	      sScene32 *label){
      sPQueue32 *Q=NULL;
      int i,j,p,q,n;
      int w,Wmax;
      sScene32 *value;
      gft::Voxel u,v;
      float per_pq;
      
      value = Scene32::Create(W);
      n = label->n;
      Wmax = gft::Scene32::GetMaximumValue(W)*2;
      Wmax *= (1.0 + fabsf(per)/100.0);     
      Q = PQueue32::Create(Wmax+2,n,value->data);

      for(p=0; p<n; p++){
	if(label->data[p]==NIL) value->data[p] = INT_MAX;
	else                    value->data[p] = 0;
      }
      
      if(S != NULL){
	for(i=1; i<=S[0]; i++)
	  PQueue32::FastInsertElem(Q, S[i]);
      }
      else{
	for(p=0; p<n; p++)
	  if(label->data[p]!=NIL)
	    PQueue32::FastInsertElem(Q, p);	    
      }

      while(!PQueue32::IsEmpty(Q)) {
	p = PQueue32::FastRemoveMinFIFO(Q);
	u.c.x = gft::Scene32::GetAddressX(label, p);
	u.c.y = gft::Scene32::GetAddressY(label, p);
	u.c.z = gft::Scene32::GetAddressZ(label, p);	
	
	for(i=1; i<A->n; i++){
	  v.v = u.v + A->d[i].v;
	  if(gft::Scene32::IsValidVoxel(label, v)){
	    q = gft::Scene32::GetVoxelAddress(label,v);
	    if(Q->L.elem[q].color != BLACK){

	      w = W->data[p] + W->data[q];
	      if(label->data[p] > 0)
		per_pq = per;
	      else
		per_pq = -per;

	      if(scn->data[p] > scn->data[q])
		w *= (1.0 + per_pq/100.0);
	      else if(scn->data[p] < scn->data[q])
		w *= (1.0 - per_pq/100.0);
	      
	      if(w < value->data[q]){
		if(Q->L.elem[q].color == GRAY)
		  PQueue32::FastRemoveElem(Q, q);
		value->data[q] = w;
		label->data[q] = label->data[p];
		PQueue32::FastInsertElem(Q, q);
	      }
	    }
	  }
	}
      }
      Scene32::Destroy(&value);
      PQueue32::Destroy(&Q);
    }



    //Outer Cut:
    void OIFT(sScene32 *Wx,
	      sScene32 *Wy,
	      sScene32 *Wz,
	      sScene32 *scn,
	      float per,
	      int *S,
	      sScene32 *label){
      sPQueue32 *Q=NULL;
      sAdjRel3 *A;
      int i,j,p,q,n,t;
      int w,Wmax;
      sScene32 *value;
      sScene32 *W[3];
      gft::Voxel u,v;
      float per_pq;
      int T[7];
      W[0] = Wx;
      W[1] = Wy;
      W[2] = Wz;
      A = gft::AdjRel3::Spheric(1.0);
      for(i=1; i<A->n; i++){
	if(A->d[i].axis.x != 0) T[i] = 0;
	if(A->d[i].axis.y != 0) T[i] = 1;
	if(A->d[i].axis.z != 0) T[i] = 2;
      }
      value = Scene32::Create(scn);
      n = label->n;
      Wmax = MAX(gft::Scene32::GetMaximumValue(Wx),
		 MAX(gft::Scene32::GetMaximumValue(Wy),
		     gft::Scene32::GetMaximumValue(Wz)))*2;
      Wmax *= (1.0 + fabsf(per)/100.0);     
      Q = PQueue32::Create(Wmax+2,n,value->data);

      for(p=0; p<n; p++){
	if(label->data[p]==NIL) value->data[p] = INT_MAX;
	else                    value->data[p] = 0;
      }
      
      if(S != NULL){
	for(i=1; i<=S[0]; i++)
	  PQueue32::FastInsertElem(Q, S[i]);
      }
      else{
	for(p=0; p<n; p++)
	  if(label->data[p]!=NIL)
	    PQueue32::FastInsertElem(Q, p);	    
      }

      while(!PQueue32::IsEmpty(Q)) {
	p = PQueue32::FastRemoveMinFIFO(Q);
	u.c.x = gft::Scene32::GetAddressX(label, p);
	u.c.y = gft::Scene32::GetAddressY(label, p);
	u.c.z = gft::Scene32::GetAddressZ(label, p);	
	
	for(i=1; i<A->n; i++){
	  v.v = u.v + A->d[i].v;
	  if(gft::Scene32::IsValidVoxel(label, v)){
	    q = gft::Scene32::GetVoxelAddress(label,v);
	    if(Q->L.elem[q].color != BLACK){
	      t = T[i];
	      w = W[t]->data[p] + W[t]->data[q];
	      if(label->data[p] > 0)
		per_pq = per;
	      else
		per_pq = -per;

	      if(scn->data[p] > scn->data[q])
		w *= (1.0 + per_pq/100.0);
	      else if(scn->data[p] < scn->data[q])
		w *= (1.0 - per_pq/100.0);
	      
	      if(w < value->data[q]){
		if(Q->L.elem[q].color == GRAY)
		  PQueue32::FastRemoveElem(Q, q);
		value->data[q] = w;
		label->data[q] = label->data[p];
		PQueue32::FastInsertElem(Q, q);
	      }
	    }
	  }
	}
      }
      Scene32::Destroy(&value);
      PQueue32::Destroy(&Q);
      AdjRel3::Destroy(&A);
    }

    
    //---------------------------------------
    
    void OIFT_in(sImageGraph *sg,
		 int *S,
		 sImage32 *label){
      sPQueue32 *Q=NULL;
      int i,j,p,q,n;
      int w;
      sImage32 *value;
      int u_x,u_y,v_x,v_y;
      sAdjRel *A;
      int *i_inv;
      
      value = Image32::Create(sg->ncols,
			      sg->nrows);
      n = label->n;
      Q = PQueue32::Create(sg->Wmax+2,n,value->data);
      A = sg->A;

      i_inv = gft::AdjRel::InverseIndexes(A);
      
      for(p=0; p<n; p++){
	if(label->data[p]==NIL) value->data[p] = INT_MAX;
	else                    value->data[p] = 0;
      }

      if(S != NULL){
	for(i=1; i<=S[0]; i++)
	  PQueue32::FastInsertElem(Q, S[i]);
      }
      else{
	for(p=0; p<n; p++)
	  if(label->data[p]!=NIL)
	    PQueue32::FastInsertElem(Q, p);	    
      }
      
      while(!PQueue32::IsEmpty(Q)) {
	p = PQueue32::FastRemoveMinFIFO(Q);
	u_x = p%label->ncols; //PixelX(label, p);
	u_y = p/label->ncols; //PixelY(label, p);
	
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  //if(Image32::IsValidPixel(label,v_x,v_y)){
	  if(v_x >= 0 && v_x < label->ncols &&
	     v_y >= 0 && v_y < label->nrows){	  
	    q = v_x + label->ncols*v_y;
	    if(Q->L.elem[q].color != BLACK){
	      
	      if(label->data[p] != 0){
		j = i_inv[i]; //j = ImageGraph::get_edge_index(q, p, sg);
		w = (sg->n_link[q])[j];
	      }
	      else
		w = (sg->n_link[p])[i];
	      
	      if(w < value->data[q]){
		if(Q->L.elem[q].color == GRAY)
		  PQueue32::FastRemoveElem(Q, q);
		value->data[q] = w;
		label->data[q] = label->data[p];
		PQueue32::FastInsertElem(Q, q);
	      }
	    }
	  }
	}
      }
      Image32::Destroy(&value);
      PQueue32::Destroy(&Q);
      free(i_inv);
    }
    
    
    void OIFT(sImageGraph *sg,
	      int *S,
	      sImage32 *label){
      sPQueue32 *Q=NULL;
      int i,j,p,q,n;
      int w;
      sImage32 *value;
      int u_x,u_y,v_x,v_y;
      sAdjRel *A;
      int *i_inv;

      value = Image32::Create(sg->ncols,
			      sg->nrows);
      n = label->ncols*label->nrows;
      Q = PQueue32::Create(sg->Wmax+2,n,value->data);
      A = sg->A;

      i_inv = gft::AdjRel::InverseIndexes(A);
      
      for(p=0; p<n; p++){
	if(label->data[p]==NIL) value->data[p] = INT_MAX;
	else                    value->data[p] = 0;
      }
      
      if(S != NULL){
	for(i=1; i<=S[0]; i++)
	  PQueue32::FastInsertElem(Q, S[i]);
      }
      else{
	for(p=0; p<n; p++)
	  if(label->data[p]!=NIL)
	    PQueue32::FastInsertElem(Q, p);	    
      }

      while(!PQueue32::IsEmpty(Q)) {
	p = PQueue32::FastRemoveMinFIFO(Q);
	u_x = p%label->ncols; //PixelX(label, p);
	u_y = p/label->ncols; //PixelY(label, p);
	
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  //if(Image32::IsValidPixel(label,v_x,v_y)){
	  if(v_x >= 0 && v_x < label->ncols &&
	     v_y >= 0 && v_y < label->nrows){
	    q = v_x + label->ncols*v_y;
	    if(Q->L.elem[q].color != BLACK){
	      
	      if(label->data[p]==0){
		j = i_inv[i]; //j = ImageGraph::get_edge_index(q, p, sg);
		w = (sg->n_link[q])[j];
	      }
	      else
		w = (sg->n_link[p])[i];
	      
	      if(w < value->data[q]){
		if(Q->L.elem[q].color == GRAY)
		  PQueue32::FastRemoveElem(Q, q);
		value->data[q] = w;
		label->data[q] = label->data[p];
		PQueue32::FastInsertElem(Q, q);
	      }
	    }
	  }
	}
      }
      Image32::Destroy(&value);
      PQueue32::Destroy(&Q);
      free(i_inv);
    }


    void OIFT_MaxMin(sImageGraph *sg,
		     int *S,
		     sImage32 *label,
		     sImage32 *pred,
		     sImage32 *value,
		     int niter){
      sPQueue32 *Q=NULL;
      int i,j,p,q,n,it = 0;
      int w;
      int u_x,u_y,v_x,v_y;
      sAdjRel *A;
      int *i_inv;

      Image32::Set(pred, NIL);
      n = label->ncols*label->nrows;
      Q = PQueue32::Create(sg->Wmax*2+2+1,n,value->data);
      A = sg->A;

      i_inv = gft::AdjRel::InverseIndexes(A);
      
      for(p=0; p<n; p++){
	if(label->data[p]==NIL) value->data[p] = sg->Wmax*2+2; //sg->Wmax+1; //INT_MAX;
	else                    value->data[p] = 0; 
      }
      
      if(S != NULL){
	for(i=1; i<=S[0]; i++){
	  PQueue32::FastInsertElem(Q, S[i]);
	  //label->data[S[i]] += 2;
	}
      }
      else{
	for(p=0; p<n; p++){
	  if(label->data[p]!=NIL){
	    PQueue32::FastInsertElem(Q, p);
	    //label->data[p] += 2;
	  }
	}
      }

      while(!PQueue32::IsEmpty(Q)) {
	if(it == niter)
	  break;

	p = PQueue32::FastRemoveMinFIFO(Q);
	u_x = p%label->ncols; //PixelX(label, p);
	u_y = p/label->ncols; //PixelY(label, p);

	//label->data[p] -= 2;
	
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  //if(Image32::IsValidPixel(label,v_x,v_y)){
	  if(v_x >= 0 && v_x < label->ncols &&
	     v_y >= 0 && v_y < label->nrows){
	    q = v_x + label->ncols*v_y;
	    if(Q->L.elem[q].color != BLACK){
	      
	      if(label->data[p]==0){
		j = i_inv[i]; //j = ImageGraph::get_edge_index(q, p, sg);
		w = (sg->n_link[q])[j];
		w = MAX(value->data[p], w*2);
	      }
	      else{
		w = (sg->n_link[p])[i];
		w = MAX(value->data[p], w*2+1);
	      }
	      if(w < value->data[q]){
		if(Q->L.elem[q].color == GRAY)
		  PQueue32::FastRemoveElem(Q, q);
		value->data[q] = w;
		label->data[q] = label->data[p];
		pred->data[q] = p;
		PQueue32::FastInsertElem(Q, q);
		//label->data[q] += 2;
	      }
	    }
	  }
	}

	it++;
      }
      /*
      for(u_y = 0; u_y < label->nrows; u_y++){
	for(u_x = 0; u_x < label->ncols; u_x++)
	  printf("%d ", label->array[u_y][u_x]);
	printf("\n");
      }
      */      
      PQueue32::Destroy(&Q);
      free(i_inv);
    }

    

    void OIFT_MinMax(sImageGraph *sg,
		     int *S,
		     sImage32 *label,
		     sImage32 *pred,
		     sImage32 *value,
		     int niter){
      sPQueue32 *Q=NULL;
      int i,j,p,q,n,it = 0;
      int w;
      int u_x,u_y,v_x,v_y;
      sAdjRel *A;
      int *i_inv;

      Image32::Set(pred, NIL);
      n = label->ncols*label->nrows;
      Q = PQueue32::Create(sg->Wmax+2,n,value->data);
      A = sg->A;

      i_inv = gft::AdjRel::InverseIndexes(A);
      
      for(p=0; p<n; p++){
	if(label->data[p]==NIL) value->data[p] = 0;
	else                    value->data[p] = sg->Wmax+1; //INT_MAX;
      }
      
      if(S != NULL){
	for(i=1; i<=S[0]; i++){
	  PQueue32::FastInsertElem(Q, S[i]);
	  label->data[S[i]] += 2;
	}
      }
      else{
	for(p=0; p<n; p++){
	  if(label->data[p]!=NIL){
	    PQueue32::FastInsertElem(Q, p);
	    label->data[p] += 2;
	  }
	}
      }

      while(!PQueue32::IsEmpty(Q)) {
	p = PQueue32::FastRemoveMaxFIFO(Q);
	u_x = p%label->ncols; //PixelX(label, p);
	u_y = p/label->ncols; //PixelY(label, p);

	label->data[p] -= 2;
	
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  //if(Image32::IsValidPixel(label,v_x,v_y)){
	  if(v_x >= 0 && v_x < label->ncols &&
	     v_y >= 0 && v_y < label->nrows){
	    q = v_x + label->ncols*v_y;
	    if(Q->L.elem[q].color != BLACK){
	      
	      if(label->data[p]==0){
		j = i_inv[i]; //j = ImageGraph::get_edge_index(q, p, sg);
		w = (sg->n_link[q])[j];
	      }
	      else
		w = (sg->n_link[p])[i];
	      
	      if(w > value->data[q]){
		if(Q->L.elem[q].color == GRAY)
		  PQueue32::FastRemoveElem(Q, q);
		value->data[q] = w;
		label->data[q] = label->data[p];
		pred->data[q] = p;
		PQueue32::FastInsertElem(Q, q);
		label->data[q] += 2;
	      }
	    }
	  }
	}

	it++;
	if(it == niter)
	  break;
      }
      PQueue32::Destroy(&Q);
      free(i_inv);
    }



    void OIFT_TZ2Bkg(sImageGraph *sg,
		     int *S,
		     sImage32 *label){
      sPQueue32 *Q=NULL;
      int i,j,p,q,n;
      int w;
      sImage32 *value;
      int u_x,u_y,v_x,v_y;
      sAdjRel *A;
      int *i_inv;

      value = Image32::Create(sg->ncols,
			      sg->nrows);
      n = label->ncols*label->nrows;
      Q = PQueue32::Create(sg->Wmax*2+3,n,value->data);
      A = sg->A;

      i_inv = gft::AdjRel::InverseIndexes(A);
      
      for(p=0; p<n; p++){
	if(label->data[p]==NIL) value->data[p] = INT_MAX;
	else                    value->data[p] = 0;
      }
      
      if(S != NULL){
	for(i=1; i<=S[0]; i++)
	  PQueue32::FastInsertElem(Q, S[i]);
      }
      else{
	for(p=0; p<n; p++)
	  if(label->data[p]!=NIL)
	    PQueue32::FastInsertElem(Q, p);	    
      }

      while(!PQueue32::IsEmpty(Q)) {
	p = PQueue32::FastRemoveMinFIFO(Q);
	u_x = p%label->ncols; //PixelX(label, p);
	u_y = p/label->ncols; //PixelY(label, p);
	
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  //if(Image32::IsValidPixel(label,v_x,v_y)){
	  if(v_x >= 0 && v_x < label->ncols &&
	     v_y >= 0 && v_y < label->nrows){
	    q = v_x + label->ncols*v_y;
	    if(Q->L.elem[q].color != BLACK){
	      
	      if(label->data[p]==0){
		j = i_inv[i]; //j = ImageGraph::get_edge_index(q, p, sg);
		w = (sg->n_link[q])[j] * 2;
	      }
	      else
		w = (sg->n_link[p])[i] * 2 + 1;

	      w = MAX(w, value->data[p]);
	      if(w < value->data[q]){
		if(Q->L.elem[q].color == GRAY)
		  PQueue32::FastRemoveElem(Q, q);
		value->data[q] = w;
		label->data[q] = label->data[p];
		PQueue32::FastInsertElem(Q, q);
	      }
	    }
	  }
	}
      }
      Image32::Destroy(&value);
      PQueue32::Destroy(&Q);
      free(i_inv);
    }



    void OIFT_TZ2Obj(sImageGraph *sg,
		     int *S,
		     sImage32 *label){
      sPQueue32 *Q=NULL;
      int i,j,p,q,n;
      int w;
      sImage32 *value;
      int u_x,u_y,v_x,v_y;
      sAdjRel *A;
      int *i_inv;

      value = Image32::Create(sg->ncols,
			      sg->nrows);
      n = label->ncols*label->nrows;
      Q = PQueue32::Create(sg->Wmax*2+3,n,value->data);
      A = sg->A;

      i_inv = gft::AdjRel::InverseIndexes(A);
      
      for(p=0; p<n; p++){
	if(label->data[p]==NIL) value->data[p] = INT_MAX;
	else                    value->data[p] = 0;
      }
      
      if(S != NULL){
	for(i=1; i<=S[0]; i++)
	  PQueue32::FastInsertElem(Q, S[i]);
      }
      else{
	for(p=0; p<n; p++)
	  if(label->data[p]!=NIL)
	    PQueue32::FastInsertElem(Q, p);	    
      }

      while(!PQueue32::IsEmpty(Q)) {
	p = PQueue32::FastRemoveMinFIFO(Q);
	u_x = p%label->ncols; //PixelX(label, p);
	u_y = p/label->ncols; //PixelY(label, p);
	
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  //if(Image32::IsValidPixel(label,v_x,v_y)){
	  if(v_x >= 0 && v_x < label->ncols &&
	     v_y >= 0 && v_y < label->nrows){
	    q = v_x + label->ncols*v_y;
	    if(Q->L.elem[q].color != BLACK){
	      
	      if(label->data[p]==0){
		j = i_inv[i]; //j = ImageGraph::get_edge_index(q, p, sg);
		w = (sg->n_link[q])[j] * 2 + 1;
	      }
	      else
		w = (sg->n_link[p])[i] * 2;

	      w = MAX(w, value->data[p]);
	      if(w < value->data[q]){
		if(Q->L.elem[q].color == GRAY)
		  PQueue32::FastRemoveElem(Q, q);
		value->data[q] = w;
		label->data[q] = label->data[p];
		PQueue32::FastInsertElem(Q, q);
	      }
	    }
	  }
	}
      }
      Image32::Destroy(&value);
      PQueue32::Destroy(&Q);
      free(i_inv);
    }



    

    void OIFT_TZ(sImageGraph *sg,
		 int *S,
		 sImage32 *label){
      sImage32 *lb_tzbkg, *lb_tzobj;
      int p;
      lb_tzbkg = gft::Image32::Clone(label);
      OIFT_TZ2Bkg(sg, S, lb_tzbkg);
      gft::ImageGraph::Transpose(sg);
      lb_tzobj = gft::Image32::Create(label);
      for(p = 0; p < label->n; p++){
	if(label->data[p] == NIL)
	  lb_tzobj->data[p] = NIL;
	else if(label->data[p] == 0)
	  lb_tzobj->data[p] = 1;
	else
	  lb_tzobj->data[p] = 0;
      }
      OIFT_TZ2Bkg(sg, S, lb_tzobj);
      for(p = 0; p < label->n; p++){
	if(lb_tzbkg->data[p] == 1)
	  label->data[p] = 1;
	else if(lb_tzobj->data[p] == 1)
	  label->data[p] = 0;
	else
	  label->data[p] = 2; /*Tie-zone*/
      }
      gft::ImageGraph::Transpose(sg);
      gft::Image32::Destroy(&lb_tzbkg);
      gft::Image32::Destroy(&lb_tzobj);
    }



    bool isOIFT(sImageGraph *sg,
		int *S,
		sImage32 *Slabel,
		sImage32 *label,
		sImage32 *pred,
		sImage32 *ord,
		bool complete_check){
      bool seg, forest;
      seg = isOIFT_Segmentation(sg,
				S,
				Slabel,
				label);
      if(complete_check)
	forest = isOIFT_Forest(sg,
			       S,
			       Slabel,
			       pred,
			       ord);
      else
	forest = isForest(sg,
			  pred);
      return (seg && forest);
    }

    
    
    bool isOIFT_Segmentation(sImageGraph *sg,
			     int *S,
			     sImage32 *Slabel,
			     sImage32 *label){
      //*************
      //static int i = 0;
      char filename[512];
      //*************
      gft::sImage32 *tz;
      bool flag = true;
      bool energy_test;
      int p;
      tz = gft::Image32::Clone(Slabel);
      OIFT_TZ(sg, S, tz);
      //*************
      //sprintf(filename, "tiezone.pgm");
      //Image32::Write(tz, filename);
      //*************
      for(p = 0; p < label->n; p++){
	if(label->data[p] != tz->data[p] &&
	   tz->data[p] != 2){
	  flag = false;
	  break;
	}
      }
      energy_test = (GetEnergy_Min(sg, tz, 1) == GetEnergy_Min(sg, label, 1));
      printf("within tie-zone: %d, energy test: %d\n", flag, energy_test);
      gft::Image32::Destroy(&tz);
      return (flag && energy_test);
    }


    

    bool isOIFT_Forest_tmp(sImageGraph *sg,
			   sImage32  *label,
			   sImage32  *Tpred,
			   sImage32  *value,
			   sPQueue32 *Q,
			   int *i_inv,
			   sImage32  *pred,
			   sImage32 *ord){
      struct node_oift_info { int label; int pred; int value; int color; };
      struct node_oift_info *backup;
      int u_x,u_y,v_x,v_y;
      int i,j,k,p,q,n,np,x;
      int w,bucket;
      bool flag = true;
      sAdjRel *A;
      int *F = NULL;
      A = sg->A;
      n = label->n;
      if(PQueue32::IsEmpty(Q))
	return true;
      
      bucket = Q->C.minvalue;
      while(Q->C.first[bucket] == NIL)
	bucket++;
      
      np = 0;
      p = Q->C.first[bucket];
      do{
	if(Tpred->data[p] == pred->data[p])
	  np++;
	p = Q->L.elem[p].next;
      }while(p != NIL);

      //printf("np: %d\n", np);
      
      if(np == 0)
	return false;
      
      F = (int *)calloc(np, sizeof(int));	
      np = 0;
      p = Q->C.first[bucket];
      w = value->data[p];
      do{
	if(Tpred->data[p] == pred->data[p]){
	  F[np] = p;
	  np++;
	}
	p = Q->L.elem[p].next;
      }while(p != NIL);	

      //ordenacao por insercao:
      for(i = 0; i < np-1; i++){
	// Insere F[i+1] em F[0],...,F[i].
	x = F[i+1];
	j = i;
	while( j >= 0 && ord->data[F[j]] > ord->data[x]){
	  F[j+1] = F[j];
	  j -= 1;
	}
	F[j+1] = x;
      }
      
      backup = (struct node_oift_info *)calloc(A->n,
					       sizeof(struct node_oift_info));
      
      for(k = 0; k < np; k++){
	p = F[k];
	u_x = p%label->ncols; //PixelX(label, p);
	u_y = p/label->ncols; //PixelY(label, p);

	//Backup:
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  //if(Image32::IsValidPixel(label,v_x,v_y)){
	  if(v_x >= 0 && v_x < label->ncols &&
	     v_y >= 0 && v_y < label->nrows){
	    q = v_x + label->ncols*v_y;
	    backup[i].value = value->data[q];
	    backup[i].label = label->data[q];
	    backup[i].pred  = Tpred->data[q];
	    backup[i].color = Q->L.elem[q].color;
	  }
	}

	PQueue32::FastRemoveElem(Q, p);
	flag = true;
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  //if(Image32::IsValidPixel(label,v_x,v_y)){
	  if(v_x >= 0 && v_x < label->ncols &&
	     v_y >= 0 && v_y < label->nrows){
	    q = v_x + label->ncols*v_y;
	    if(Q->L.elem[q].color != BLACK){
	      
	      if(label->data[p]==0){
		j = i_inv[i]; //j = ImageGraph::get_edge_index(q, p, sg);
		w = (sg->n_link[q])[j];
	      }
	      else
		w = (sg->n_link[p])[i];
	      
	      if(w < value->data[q]){
		if(Q->L.elem[q].color == GRAY)
		  PQueue32::FastRemoveElem(Q, q);
		value->data[q] = w;
		label->data[q] = label->data[p];
		Tpred->data[q] = p;
		PQueue32::FastInsertElem(Q, q);
	      }
	      
	      if(Tpred->data[q] != p && p == pred->data[q]){
		flag = false;
		break;
	      }
	      
	    }
	  }
	}

	if(flag)
	  flag = isOIFT_Forest_tmp(sg, label, Tpred, value,
				   Q, i_inv, pred, ord);
	
	if(flag)
	  break;

	//Restore:
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  //if(Image32::IsValidPixel(label,v_x,v_y)){
	  if(v_x >= 0 && v_x < label->ncols &&
	     v_y >= 0 && v_y < label->nrows){
	    q = v_x + label->ncols*v_y;
	    label->data[q] = backup[i].label;
	    Tpred->data[q] = backup[i].pred;
	    if(backup[i].color == WHITE){
	      if(Q->L.elem[q].color == GRAY)
		PQueue32::FastRemoveElem(Q, q);
	      value->data[q] = backup[i].value;
	    }
	    else if(backup[i].color == GRAY){
	      if(Q->L.elem[q].color == GRAY)
		PQueue32::FastRemoveElem(Q, q);
	      value->data[q] = backup[i].value;
	      PQueue32::FastInsertElem(Q, q);
	    }
	    Q->L.elem[q].color = backup[i].color;
	  }
	}
	PQueue32::FastInsertElem(Q, p);
      }
      free(F);
      free(backup);
      
      return flag;
    }

    

    
    bool isOIFT_Forest(sImageGraph *sg,
		       int *S,
		       sImage32 *Slabel,
		       sImage32 *pred,
		       sImage32 *ord){
      sPQueue32 *Q=NULL;
      int i,p,n;
      sImage32 *value, *label, *Tpred;
      int *i_inv;
      bool flag = true;

      if(!isForest(sg, pred))
	return false;
      
      label = Image32::Clone(Slabel);
      Tpred = Image32::Create(sg->ncols,
			      sg->nrows);
      Image32::Set(Tpred, NIL);
      value = Image32::Create(sg->ncols,
			      sg->nrows);
      n = label->n;
      Q = PQueue32::Create(sg->Wmax+2,n,value->data);

      i_inv = gft::AdjRel::InverseIndexes(sg->A);
      
      for(p=0; p<n; p++){
	if(label->data[p]==NIL) value->data[p] = INT_MAX;
	else                    value->data[p] = 0;
      }
      
      if(S != NULL){
	for(i=1; i<=S[0]; i++)
	  PQueue32::FastInsertElem(Q, S[i]);
      }
      else{
	for(p=0; p<n; p++)
	  if(label->data[p]!=NIL)
	    PQueue32::FastInsertElem(Q, p);	    
      }

      flag = isOIFT_Forest_tmp(sg,
			       label,
			       Tpred,
			       value,
			       Q,
			       i_inv,
			       pred,
			       ord);

      Image32::Destroy(&value);
      Image32::Destroy(&label);
      Image32::Destroy(&Tpred);
      PQueue32::Destroy(&Q);
      free(i_inv);      
      return flag;
    }



    bool isForest(sImageGraph *sg,
		  sImage32 *pred){
      gft::sBMap *path,*acyclic;
      bool forest = true;
      int p,q;
      acyclic = gft::BMap::Create(pred->n);
      gft::BMap::Fill(acyclic, 0);
      path = gft::BMap::Create(pred->n);
      gft::BMap::Fill(path, 0);      
      for(p = 0; p < pred->n; p++){
	q = p;
	do{
	  if(gft::BMap::Get(path, q) == 1){
	    forest = false;
	    break;
	  }
	  else if(gft::BMap::Get(acyclic, q) == 1)
	    break;

	  gft::BMap::Set1(path, q);
	  
	  q = pred->data[q];
	}while(q != NIL);

	if(!forest) break;

	q = p;
	do{
	  if(gft::BMap::Get(acyclic, q) == 1)
	    break;
	  gft::BMap::Set0(path, q);
	  gft::BMap::Set1(acyclic, q);	  
	  q = pred->data[q];
	}while(q != NIL);
      }      
      gft::BMap::Destroy(&acyclic);
      gft::BMap::Destroy(&path);

      if(!forest)
	printf("Cycle detected.\n");
      
      return forest;
    }




    void OIFT_guided(sImageGraph *sg,
		     int *S,
		     sImage32 *label,
		     sImage32 *pred,
		     sImage32 *ord){
      sHeap32fi_lex *Q=NULL;
      int i,j,p,q,n;
      int w;
      float *value;
      int u_x,u_y,v_x,v_y;
      sAdjRel *A;
      int *i_inv;

      value = (float *)calloc(label->n, sizeof(float));
      n = label->ncols*label->nrows;
      Q = gft::Heap32fi_lex::Create(label->n, value, ord->data);
      A = sg->A;

      i_inv = gft::AdjRel::InverseIndexes(A);

      Image32::Set(pred, NIL);
      for(p=0; p<n; p++){
	if(label->data[p]==NIL) value[p] = FLT_MAX;
	else                    value[p] = 0.0;
      }
      
      if(S != NULL){
	for(i=1; i<=S[0]; i++)
	  gft::Heap32fi_lex::Insert_MinPolicy(Q, S[i]);
      }
      else{
	for(p=0; p<n; p++)
	  if(label->data[p]!=NIL)
	    gft::Heap32fi_lex::Insert_MinPolicy(Q, p);
      }

      while(!Heap32fi_lex::IsEmpty(Q)) {
	Heap32fi_lex::Remove_MinPolicy(Q, &p);
	u_x = p%label->ncols; //PixelX(label, p);
	u_y = p/label->ncols; //PixelY(label, p);
	
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  //if(Image32::IsValidPixel(label,v_x,v_y)){
	  if(v_x >= 0 && v_x < label->ncols &&
	     v_y >= 0 && v_y < label->nrows){
	    q = v_x + label->ncols*v_y;
	    if(Q->color[q] != BLACK){
	      
	      if(label->data[p]==0){
		j = i_inv[i];
		w = (sg->n_link[q])[j];
	      }
	      else
		w = (sg->n_link[p])[i];
	      
	      if(w < value[q]){
		//value[q] = w;
		label->data[q] = label->data[p];
		pred->data[q] = p;
		gft::Heap32fi_lex::Update_MinPolicy(Q, q, w, ord->data[q]);
	      }
	    }
	  }
	}
      }
      free(value);
      Heap32fi_lex::Destroy(&Q);
      free(i_inv);
    }


    
    
    
    void OIFT(sGraph *graph,
	      sGraph *transpose,
	      int *S,
	      int *label){
      sPQueue32 *Q=NULL;
      sGraph *g;
      int i,j,p,q,n;
      int w;
      int *value;
      int Wmax;
      Wmax = Graph::GetMaximumArc(graph);
      n = graph->nnodes;
      value = (int *)malloc(n*sizeof(int));
      Q = PQueue32::Create(Wmax+2, n, value);

      for(p = 0; p < n; p++){
	if(label[p]==NIL) value[p] = INT_MAX;
	else              value[p] = 0;
      }
      
      if(S != NULL){
	for(i = 1; i <= S[0]; i++)
	  PQueue32::FastInsertElem(Q, S[i]);
      }
      else{
	for(p = 0; p < n; p++)
	  if(label[p]!=NIL)
	    PQueue32::FastInsertElem(Q, p);
      }

      while(!PQueue32::IsEmpty(Q)) {
	p = PQueue32::FastRemoveMinFIFO(Q);

	if(label[p]==0) g = transpose;
	else   	        g = graph;

	for(i = 0; i < g->nodes[p].outdegree; i++){
	  q = g->nodes[p].adjList[i];

	  if(Q->L.elem[q].color != BLACK){

	    /*
	    if(label[p]==0)
	      w = Graph::GetArcWeight(graph, q, p);
	    else
            */
	    w = g->nodes[p].Warcs[i];
	    
	    if(w < value[q]){
	      if(Q->L.elem[q].color == GRAY)
		PQueue32::FastRemoveElem(Q, q);
	      value[q] = w;
	      label[q] = label[p];
	      PQueue32::FastInsertElem(Q, q);
	    }
	  }
	}
      }
      free(value);
      PQueue32::Destroy(&Q);
    }
    



    void OIFT_Heap(sImageGraph *sg,
		   int *S,
		   sImage32 *label){
      sHeap *Q=NULL;
      int i,j,p,q,n;
      int wi;
      float w;
      float *value;
      int u_x,u_y,v_x,v_y;
      sAdjRel *A;
      int *i_inv;

      n = label->ncols*label->nrows;
      value = gft::AllocFloatArray(n);
      Q = Heap::Create(n, value);
      A = sg->A;

      i_inv = gft::AdjRel::InverseIndexes(A);
      
      for(p=0; p<n; p++){
	if(label->data[p]==NIL) value[p] = FLT_MAX;
	else                    value[p] = 0.0;
      }
      
      if(S != NULL){
	for(i=1; i<=S[0]; i++)
	  Heap::Insert_MinPolicy(Q, S[i]);
      }
      else{
	for(p=0; p<n; p++)
	  if(label->data[p]!=NIL)
	    Heap::Insert_MinPolicy(Q, p);
      }

      while(!Heap::IsEmpty(Q)) {
	Heap::Remove_MinPolicy(Q, &p);
	u_x = p%label->ncols; //PixelX(label, p);
	u_y = p/label->ncols; //PixelY(label, p);
	
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  //if(Image32::IsValidPixel(label,v_x,v_y)){
	  if(v_x >= 0 && v_x < label->ncols &&
	     v_y >= 0 && v_y < label->nrows){
	    q = v_x + label->ncols*v_y;
	    if(Q->color[q] != BLACK){
      
	      if(label->data[p]==0){
		j = i_inv[i]; //j = ImageGraph::get_edge_index(q, p, sg);
		wi = (sg->n_link[q])[j];
	      }
	      else
		wi = (sg->n_link[p])[i];

	      if(wi == INT_MAX)
		w = FLT_MAX;
	      else
		w = (float)wi;
	      
	      if(w < value[q]){
		label->data[q] = label->data[p];
		Heap::Update_MinPolicy(Q, q, w);
	      }
	    }
	  }
	}
      }
      gft::FreeFloatArray(&value);
      Heap::Destroy(&Q);
      free(i_inv);
    }



    void OIFT_Heap(sGraph *graph,
		   sGraph *transpose,
		   int *S,
		   int *label){
      sHeap *Q=NULL;
      sGraph *g;
      int i,j,p,q,n;
      int wi;
      float w;
      float *value;
      int u_x,u_y,v_x,v_y;

      n = graph->nnodes;
      value = gft::AllocFloatArray(n);
      Q = Heap::Create(n, value);

      for(p = 0; p < n; p++){
	if(label[p]==NIL) value[p] = FLT_MAX;
	else              value[p] = 0.0;
      }
      
      if(S != NULL){
	for(i = 1; i <= S[0]; i++)
	  Heap::Insert_MinPolicy(Q, S[i]);
      }
      else{
	for(p = 0; p < n; p++)
	  if(label[p]!=NIL)
	    Heap::Insert_MinPolicy(Q, p);
      }

      while(!Heap::IsEmpty(Q)) {
	Heap::Remove_MinPolicy(Q, &p);

	if(label[p]==0) g = transpose;
	else   	        g = graph;
	
	for(i = 0; i < g->nodes[p].outdegree; i++){
	  q = g->nodes[p].adjList[i];

	  if(Q->color[q] != BLACK){

	    /*
	    if(label[p]==0)
	      wi = Graph::GetArcWeight(graph, q, p);
	    else
	    */
	    wi = g->nodes[p].Warcs[i];
	    
	    if(wi == INT_MAX)
	      w = FLT_MAX;
	    else
	      w = (float)wi;
	    
	    if(w < value[q]){
	      label[q] = label[p];
	      Heap::Update_MinPolicy(Q, q, w);
	    }
	  }
	}
      }
      gft::FreeFloatArray(&value);
      Heap::Destroy(&Q);
    }
    
    

    /*
    void EOIFT(sImageGraph *sg,
	       int *S,
	       sImage32 *label){
      sPQueue32 *Qobj=NULL, *Qbkg=NULL;
      int i,j,p,p_obj,p_bkg,q,n;
      int l_ant,e_obj,e_bkg,e_max;
      int w;
      sImage32 *value;
      int u_x,u_y,v_x,v_y;
      sAdjRel *A;
      int *i_inv;
      int *Q=NULL;
      int Qtop = -1;
      value = Image32::Create(sg->ncols,
			      sg->nrows);
      n = label->ncols*label->nrows;
      Qobj = PQueue32::Create(sg->Wmax+2,n,value->data);
      Qbkg = PQueue32::Create(sg->Wmax+2,n,value->data);
      Q	= gft::AllocIntArray(n);
      A = sg->A;

      i_inv = gft::AdjRel::InverseIndexes(A);

      if(S != NULL){
	for(i=1; i<=S[0]; i++){
	  value->data[S[i]] = 0;
	  if(label->data[S[i]] == 0)
	    PQueue32::FastInsertElem(Qbkg, S[i]);
	  else if(label->data[S[i]] != NIL)
	    PQueue32::FastInsertElem(Qobj, S[i]);
	}
	for(p=0; p<n; p++){
	  if(label->data[p]==NIL){
	    value->data[p] = INT_MAX;
	    label->data[p] = 0;
	  }
	  else
	    value->data[p] = 0;
	}
      }
      else{
	for(p=0; p<n; p++){
	  if(label->data[p] == 0){
	    value->data[p] = 0;
	    PQueue32::FastInsertElem(Qbkg, p);
	  }
	  else if(label->data[p] != NIL){
	    value->data[p] = 0;
	    PQueue32::FastInsertElem(Qobj, p);
	  }
	  else{
	    value->data[p] = INT_MAX;
	    label->data[p] = 0;
	  }
	}
      }

      l_ant = 0;
      while(!PQueue32::IsEmpty(Qobj) && !PQueue32::IsEmpty(Qbkg)) {
	p_obj = PQueue32::FastGetMinFIFO(Qobj);
	p_bkg = PQueue32::FastGetMinFIFO(Qbkg);

	e_obj = value->data[p_obj];
	e_bkg = value->data[p_bkg];
	if(e_obj < e_bkg){
	  e_max = e_bkg;
	  p = p_obj;
	  PQueue32::FastRemoveElem(Qobj, p);
	}
	else if(e_obj > e_bkg){
	  e_max = e_obj;
	  p = p_bkg;
	  PQueue32::FastRemoveElem(Qbkg, p);
	}
	else{
	  e_max = e_obj;
          if(l_ant == 0){
	    p = p_obj;
	    PQueue32::FastRemoveElem(Qobj, p);
	  }
          else{
	    p = p_bkg;
	    PQueue32::FastRemoveElem(Qbkg, p);
	  }
          l_ant = 1 - l_ant;
	}

	Qobj->L.elem[p].color = BLACK;
	Qbkg->L.elem[p].color = BLACK;

	Qtop++;
	Q[Qtop] = p;
	while(Qtop > -1){
	  p = Q[Qtop];
	  Qtop--;
	  u_x = p%label->ncols;
	  u_y = p/label->ncols;
	  for(i=1; i<A->n; i++){
	    v_x = u_x + A->dx[i];
	    v_y = u_y + A->dy[i];
	    if(Image32::IsValidPixel(label,v_x,v_y)){
	      q = v_x + label->ncols*v_y;
	      if(Qobj->L.elem[q].color != BLACK){
		
		if(label->data[p]==0){
		  j = i_inv[i];
		  w = (sg->n_link[q])[j];
		}
		else
		  w = (sg->n_link[p])[i];

                if(w < e_max){
		  if(Qobj->L.elem[q].color == GRAY)
		    PQueue32::FastRemoveElem(Qobj, q);
		  else if(Qbkg->L.elem[q].color == GRAY)
		    PQueue32::FastRemoveElem(Qbkg, q);
		  Qobj->L.elem[q].color = BLACK;
		  Qbkg->L.elem[q].color = BLACK;
		  label->data[q] = label->data[p];
		  Qtop++;
		  Q[Qtop] = q;
		}
		else if(w < value->data[q]){
		  if(Qobj->L.elem[q].color == GRAY)
		    PQueue32::FastRemoveElem(Qobj, q);
		  if(Qbkg->L.elem[q].color == GRAY)
		    PQueue32::FastRemoveElem(Qbkg, q);
		  Qobj->L.elem[q].color = WHITE;
		  Qbkg->L.elem[q].color = WHITE;
		  value->data[q] = w;
		  label->data[q] = label->data[p];
		  if(label->data[q] > 0)
		    PQueue32::FastInsertElem(Qobj, q);
		  else
		    PQueue32::FastInsertElem(Qbkg, q);
		}
	      }
	    }
	  }
	}
      }

      while(!PQueue32::IsEmpty(Qobj)){
	p = PQueue32::FastRemoveMinFIFO(Qobj);
	Qtop++;
	Q[Qtop] = p;
	while(Qtop > -1){
	  p = Q[Qtop];
	  Qtop--;
	  u_x = p%label->ncols;
	  u_y = p/label->ncols;
	  for(i=1; i<A->n; i++){
	    v_x = u_x + A->dx[i];
	    v_y = u_y + A->dy[i];
	    if(Image32::IsValidPixel(label,v_x,v_y)){
	      q = v_x + label->ncols*v_y;
	      if(Qobj->L.elem[q].color != BLACK){

		if(Qobj->L.elem[q].color == GRAY)
		  PQueue32::FastRemoveElem(Qobj, q);
		Qobj->L.elem[q].color = BLACK;
		label->data[q] = label->data[p];
		Qtop++;
		Q[Qtop] = q;
	      }
	    }
	  }
	}
      }
      Image32::Destroy(&value);
      PQueue32::Destroy(&Qobj);
      PQueue32::Destroy(&Qbkg);
      gft::FreeIntArray(&Q);
      free(i_inv);
    }
    */

    void EOIFT(sImageGraph *sg,
	       int *S,
	       sImage32 *label,
	       int e_max){
      sPQueue32 *Qobj=NULL, *Qbkg=NULL;
      int i,j,p,p_obj,p_bkg,q,n,lp;
      int l_ant,e_obj,e_bkg;
      int w;
      sImage32 *value;
      int u_x,u_y,v_x,v_y;
      sAdjRel *A;
      int *i_inv;
      int *Q=NULL;
      int Qtop = -1;
      value = Image32::Create(sg->ncols,
			      sg->nrows);
      n = label->ncols*label->nrows;
      Qobj = PQueue32::Create(sg->Wmax+2,n,value->data);
      Qbkg = PQueue32::Create(sg->Wmax+2,n,value->data);
      Q	= gft::AllocIntArray(n);
      A = sg->A;

      i_inv = gft::AdjRel::InverseIndexes(A);

      if(S != NULL){
	for(i=1; i<=S[0]; i++){
	  value->data[S[i]] = 0;
	  if(label->data[S[i]] == 0)
	    PQueue32::FastInsertElem(Qbkg, S[i]);
	  else if(label->data[S[i]] != NIL)
	    PQueue32::FastInsertElem(Qobj, S[i]);
	}
	for(p=0; p<n; p++){
	  if(label->data[p]==NIL){
	    value->data[p] = INT_MAX;
	    label->data[p] = 0;
	  }
	  else
	    value->data[p] = 0;
	}
      }
      else{
	for(p=0; p<n; p++){
	  if(label->data[p] == 0){
	    value->data[p] = 0;
	    PQueue32::FastInsertElem(Qbkg, p);
	  }
	  else if(label->data[p] != NIL){
	    value->data[p] = 0;
	    PQueue32::FastInsertElem(Qobj, p);
	  }
	  else{
	    value->data[p] = INT_MAX;
	    label->data[p] = 0;
	  }
	}
      }

      l_ant = 0;
      //while(!PQueue32::IsEmpty(Qobj) && !PQueue32::IsEmpty(Qbkg)){
      while(Qobj->nadded != 0 && Qbkg->nadded != 0){
	p_obj = PQueue32::FastGetMinFIFO(Qobj);
	p_bkg = PQueue32::FastGetMinFIFO(Qbkg);

	e_obj = value->data[p_obj];
	e_bkg = value->data[p_bkg];
	if(e_obj < e_bkg){
	  e_max = MAX(e_max, e_bkg);
	  p = p_obj;
	  PQueue32::FastRemoveElem(Qobj, p);
	}
	else if(e_obj > e_bkg){
	  e_max = MAX(e_max, e_obj);
	  p = p_bkg;
	  PQueue32::FastRemoveElem(Qbkg, p);
	}
	else{
	  e_max = MAX(e_max, e_obj);
          if(l_ant == 0){
	    p = p_obj;
	    PQueue32::FastRemoveElem(Qobj, p);
	  }
          else{
	    p = p_bkg;
	    PQueue32::FastRemoveElem(Qbkg, p);
	  }
          l_ant = 1 - l_ant;
	}

	Qobj->L.elem[p].color = BLACK;
	Qbkg->L.elem[p].color = BLACK;

	lp = label->data[p];

	goto label04;
	
	Qtop++;
	Q[Qtop] = p;
	while(Qtop > -1){
	  p = Q[Qtop];
	  Qtop--;

	label04:

	  u_x = p%label->ncols;
	  u_y = p/label->ncols;
	  for(i=1; i<A->n; i++){
	    v_x = u_x + A->dx[i];
	    v_y = u_y + A->dy[i];
	    if(v_x >= 0 && v_x < label->ncols &&
	       v_y >= 0 && v_y < label->nrows){
	      q = v_x + label->ncols*v_y;
	      if(Qobj->L.elem[q].color != BLACK){
		
		if(lp == 0){
		  j = i_inv[i];
		  w = (sg->n_link[q])[j];
		}
		else
		  w = (sg->n_link[p])[i];

                if(w < e_max){
		  if(Qobj->L.elem[q].color == GRAY)
		    PQueue32::FastRemoveElem(Qobj, q);
		  else if(Qbkg->L.elem[q].color == GRAY)
		    PQueue32::FastRemoveElem(Qbkg, q);
		  Qobj->L.elem[q].color = BLACK;
		  Qbkg->L.elem[q].color = BLACK;
		  label->data[q] = lp;
		  Qtop++;
		  Q[Qtop] = q;
		}
		else if(w < value->data[q]){
		  if(Qobj->L.elem[q].color == GRAY)
		    PQueue32::FastRemoveElem(Qobj, q);
		  if(Qbkg->L.elem[q].color == GRAY)
		    PQueue32::FastRemoveElem(Qbkg, q);
		  Qobj->L.elem[q].color = WHITE;
		  Qbkg->L.elem[q].color = WHITE;
		  value->data[q] = w;
		  label->data[q] = lp;
		  if(lp > 0)
		    PQueue32::FastInsertElem(Qobj, q);
		  else
		    PQueue32::FastInsertElem(Qbkg, q);
		}
	      }
	    }
	  }
	}
      }

      //while(!PQueue32::IsEmpty(Qobj)){
      while(Qobj->nadded != 0){
	p = PQueue32::FastRemoveMinFIFO(Qobj);
	Qtop++;
	Q[Qtop] = p;
	while(Qtop > -1){
	  p = Q[Qtop];
	  Qtop--;
	  u_x = p%label->ncols;
	  u_y = p/label->ncols;
	  for(i=1; i<A->n; i++){
	    v_x = u_x + A->dx[i];
	    v_y = u_y + A->dy[i];
	    if(v_x >= 0 && v_x < label->ncols &&
	       v_y >= 0 && v_y < label->nrows){
	      q = v_x + label->ncols*v_y;
	      if(Qobj->L.elem[q].color != BLACK){
		if(Qobj->L.elem[q].color == GRAY)
		  PQueue32::FastRemoveElem(Qobj, q);
		Qobj->L.elem[q].color = BLACK;
		label->data[q] = label->data[p];
		Qtop++;
		Q[Qtop] = q;
	      }
	    }
	  }
	}
      }
      Image32::Destroy(&value);
      PQueue32::Destroy(&Qobj);
      PQueue32::Destroy(&Qbkg);
      gft::FreeIntArray(&Q);
      free(i_inv);
    }
    


    //------------

    void EOIFT(sGraph *graph,
	       sGraph *transpose,
	       int *S,
	       int *label,
	       int e_max){
      sPQueue32 *Qobj=NULL, *Qbkg=NULL;
      sGraph *g;
      int i,j,p,p_obj,p_bkg,q,n,lp;
      int l_ant,e_obj,e_bkg;
      int w;
      int *value;
      int *Q=NULL;
      int Qtop = -1;
      int Wmax;
      Wmax = Graph::GetMaximumArc(graph);
      n = graph->nnodes;
      value = (int *)malloc(n*sizeof(int));

      Qobj = PQueue32::Create(Wmax+2, n, value);
      Qbkg = PQueue32::Create(Wmax+2, n, value);
      Q	= gft::AllocIntArray(n);

      if(S != NULL){
	for(i=1; i<=S[0]; i++){
	  value[S[i]] = 0;
	  if(label[S[i]] == 0)
	    PQueue32::FastInsertElem(Qbkg, S[i]);
	  else if(label[S[i]] != NIL)
	    PQueue32::FastInsertElem(Qobj, S[i]);
	}
	for(p=0; p<n; p++){
	  if(label[p]==NIL){
	    value[p] = INT_MAX;
	    label[p] = 0;
	  }
	  else
	    value[p] = 0;
	}
      }
      else{
	for(p=0; p<n; p++){
	  if(label[p] == 0){
	    value[p] = 0;
	    PQueue32::FastInsertElem(Qbkg, p);
	  }
	  else if(label[p] != NIL){
	    value[p] = 0;
	    PQueue32::FastInsertElem(Qobj, p);
	  }
	  else{
	    value[p] = INT_MAX;
	    label[p] = 0;
	  }
	}
      }

      l_ant = 0;
      //while(!PQueue32::IsEmpty(Qobj) && !PQueue32::IsEmpty(Qbkg)){
      while(Qobj->nadded != 0 && Qbkg->nadded != 0){
	p_obj = PQueue32::FastGetMinFIFO(Qobj);
	p_bkg = PQueue32::FastGetMinFIFO(Qbkg);

	e_obj = value[p_obj];
	e_bkg = value[p_bkg];
	if(e_obj < e_bkg){
	  e_max = MAX(e_max, e_bkg);
	  p = p_obj;
	  PQueue32::FastRemoveElem(Qobj, p);
	}
	else if(e_obj > e_bkg){
	  e_max = MAX(e_max, e_obj);
	  p = p_bkg;
	  PQueue32::FastRemoveElem(Qbkg, p);
	}
	else{
	  e_max = MAX(e_max, e_obj);
          if(l_ant == 0){
	    p = p_obj;
	    PQueue32::FastRemoveElem(Qobj, p);
	  }
          else{
	    p = p_bkg;
	    PQueue32::FastRemoveElem(Qbkg, p);
	  }
          l_ant = 1 - l_ant;
	}
	
	Qobj->L.elem[p].color = BLACK;
	Qbkg->L.elem[p].color = BLACK;

	lp = label[p];

	if(lp==0) g = transpose;
	else      g = graph;
	
	goto label05;

	Qtop++;
	Q[Qtop] = p;
	while(Qtop > -1){
	  p = Q[Qtop];
	  Qtop--;

	label05:
	  
	  for(i = 0; i < g->nodes[p].outdegree; i++){
	    q = g->nodes[p].adjList[i];
	    if(Qobj->L.elem[q].color != BLACK){
	      w = g->nodes[p].Warcs[i];
	      
	      if(w < e_max){
		if(Qobj->L.elem[q].color == GRAY)
		  PQueue32::FastRemoveElem(Qobj, q);
		if(Qbkg->L.elem[q].color == GRAY)
		    PQueue32::FastRemoveElem(Qbkg, q);
		Qobj->L.elem[q].color = BLACK;
		Qbkg->L.elem[q].color = BLACK;
		label[q] = lp;
		Qtop++;
		Q[Qtop] = q;
	      }
	      else if(w < value[q]){
		if(Qobj->L.elem[q].color == GRAY)
		  PQueue32::FastRemoveElem(Qobj, q);
		if(Qbkg->L.elem[q].color == GRAY)
		  PQueue32::FastRemoveElem(Qbkg, q);
		Qobj->L.elem[q].color = WHITE;
		Qbkg->L.elem[q].color = WHITE;
		value[q] = w;
		label[q] = lp;
		if(lp > 0)
		  PQueue32::FastInsertElem(Qobj, q);
		else
		  PQueue32::FastInsertElem(Qbkg, q);		  
	      }
	      
	    }
	  }
	}
      }
      
      //while(!PQueue32::IsEmpty(Qobj)){
      while(Qobj->nadded != 0){
	p = PQueue32::FastRemoveMinFIFO(Qobj);
	Qtop++;
	Q[Qtop] = p;
	while(Qtop > -1){
	  p = Q[Qtop];
	  Qtop--;
	  for(i = 0; i < graph->nodes[p].outdegree; i++){
	    q = graph->nodes[p].adjList[i];
	    
	    if(Qobj->L.elem[q].color != BLACK){
	      
	      if(Qobj->L.elem[q].color == GRAY)
		PQueue32::FastRemoveElem(Qobj, q);
	      Qobj->L.elem[q].color = BLACK;
	      label[q] = label[p];
	      Qtop++;
	      Q[Qtop] = q;
	    }
	  }
	}
      }
      free(value);
      PQueue32::Destroy(&Qobj);
      PQueue32::Destroy(&Qbkg);
      gft::FreeIntArray(&Q);
    }
    

    //----------------------------------------
    /*
    void EOIFT_Heap_2(sImageGraph *sg,
		      int *S,
		      sImage32 *label){
      sHeap *Qobj=NULL, *Qbkg=NULL;
      int i,j,p,p_obj,p_bkg,q,n;
      int l_ant;
      float e_obj,e_bkg,e_max;
      int wi;
      float w;
      float *value;
      int u_x,u_y,v_x,v_y;
      sAdjRel *A;
      int *i_inv;
      int *Q=NULL;
      int Qtop = -1;
      //----------
      n = label->n;
      value = gft::AllocFloatArray(n);
      Qobj = Heap::Create(n, value);
      Qbkg = Heap::Create(n, value);
      Q	= gft::AllocIntArray(n);
      A = sg->A;
      
      i_inv = gft::AdjRel::InverseIndexes(A);
      
      if(S != NULL){
	for(i=1; i<=S[0]; i++){
	  value[S[i]] = 0.0;
	  if(label->data[S[i]] == 0)
	    Heap::Insert_MinPolicy(Qbkg, S[i]);
	  else if(label->data[S[i]] != NIL)
	    Heap::Insert_MinPolicy(Qobj, S[i]);
	}
	for(p=0; p<n; p++){
	  if(label->data[p]==NIL){
	    value[p] = FLT_MAX;
	    label->data[p] = 0;
	  }
	  else
	    value[p] = 0.0;
	}
      }
      else{
	for(p=0; p<n; p++){
	  if(label->data[p] == 0){
	    value[p] = 0.0;
	    Heap::Insert_MinPolicy(Qbkg, p);
	  }
	  else if(label->data[p] != NIL){
	    value[p] = 0.0;
	    Heap::Insert_MinPolicy(Qobj, p);
	  }
	  else{
	    value[p] = FLT_MAX;
	    label->data[p] = 0;
	  }
	}
      }

      l_ant = 0;
      while(!Heap::IsEmpty(Qobj) && !Heap::IsEmpty(Qbkg)) {
	Heap::Get_MinPolicy(Qobj, &p_obj);
	Heap::Get_MinPolicy(Qbkg, &p_bkg);

	e_obj = value[p_obj];
	e_bkg = value[p_bkg];
	if(e_obj < e_bkg){
	  e_max = e_bkg;
	  p = p_obj;
	  Heap::Delete_MinPolicy(Qobj, p);
	}
	else if(e_obj > e_bkg){
	  e_max = e_obj;
	  p = p_bkg;
	  Heap::Delete_MinPolicy(Qbkg, p);
	}
	else{
	  e_max = e_obj;
          if(l_ant == 0){
	    p = p_obj;
	    Heap::Delete_MinPolicy(Qobj, p);
	  }
          else{
	    p = p_bkg;
	    Heap::Delete_MinPolicy(Qbkg, p);
	  }
          l_ant = 1 - l_ant;
	}

	//-----------------
	//printf("e_obj: %5d, e_bkg: %5d, ", ROUND(e_obj), ROUND(e_bkg));
	//printf("x: %4d, y: %4d\n", p%label->ncols, p/label->ncols);
	//-----------------
	
	Qobj->color[p] = BLACK;
	Qbkg->color[p] = BLACK;

	//Stack::Push(Q, p);
	Qtop++;
	Q[Qtop] = p;
	while(Qtop > -1){
	  //p = Stack::Pop(Q);
	  p = Q[Qtop];
	  Qtop--;
	  u_x = p%label->ncols; //PixelX(label, p);
	  u_y = p/label->ncols; //PixelY(label, p);
	  for(i=1; i<A->n; i++){
	    v_x = u_x + A->dx[i];
	    v_y = u_y + A->dy[i];
	    if(Image32::IsValidPixel(label,v_x,v_y)){
	      q = v_x + label->ncols*v_y;
	      if(Qobj->color[q] != BLACK){
		
		if(label->data[p]==0){
		  j = i_inv[i]; //j = ImageGraph::get_edge_index(q, p, sg);
		  wi = (sg->n_link[q])[j];
		}
		else
		  wi = (sg->n_link[p])[i];

		if(wi == INT_MAX)
		  w = FLT_MAX;
		else
		  w = (float)wi;
		
                if(w < e_max){
		  if(Qobj->color[q] == GRAY)
		    Heap::Delete_MinPolicy(Qobj, q);
		  else if(Qbkg->color[q] == GRAY)
		    Heap::Delete_MinPolicy(Qbkg, q);
		  Qobj->color[q] = BLACK;
		  Qbkg->color[q] = BLACK;
		  label->data[q] = label->data[p];
		  //Stack::Push(Q, q);
		  Qtop++;
		  Q[Qtop] = q;
		  //---------
		  //NbyQ++;
		}
		else if(w < value[q]){
		  label->data[q] = label->data[p];
		  if(label->data[q] > 0){
		    if(Qbkg->color[q] == GRAY)
		      Heap::Delete_MinPolicy(Qbkg, q);
		    Heap::Update_MinPolicy(Qobj, q, w);
		  }
		  else{
		    if(Qobj->color[q] == GRAY)
		      Heap::Delete_MinPolicy(Qobj, q);
		    Heap::Update_MinPolicy(Qbkg, q, w);
		  }
		}
	      }
	    }
	  }
	}
      }

      while(!Heap::IsEmpty(Qobj)){
	Heap::Remove_MinPolicy(Qobj, &p);
	//Stack::Push(Q, p);
	Qtop++;
	Q[Qtop] = p;
	while(Qtop > -1){
	  //p = Stack::Pop(Q);
	  p = Q[Qtop];
	  Qtop--;
	  u_x = p%label->ncols; //PixelX(label, p);
	  u_y = p/label->ncols; //PixelY(label, p);
	  for(i=1; i<A->n; i++){
	    v_x = u_x + A->dx[i];
	    v_y = u_y + A->dy[i];
	    if(Image32::IsValidPixel(label,v_x,v_y)){
	      q = v_x + label->ncols*v_y;
	      if(Qobj->color[q] != BLACK){

		if(Qobj->color[q] == GRAY)
		  Heap::Delete_MinPolicy(Qobj, q);
		Qobj->color[q] = BLACK;
		label->data[q] = label->data[p];
		//Stack::Push(Q, q);
		Qtop++;
		Q[Qtop] = q;
		//--------
		//NbyQ++;
	      }
	    }
	  }
	}
      }

      //printf("NbyQ: %d -> %f\n",NbyQ, (float)NbyQ/(float)n);
      
      gft::FreeFloatArray(&value);
      Heap::Destroy(&Qobj);
      Heap::Destroy(&Qbkg);
      gft::FreeIntArray(&Q);
      free(i_inv);
    }
    */
    //------------------------------------

    /*
    void EOIFT_Heap(sImageGraph *sg,
		    int *S,
		    sImage32 *label,
		    float e_max){
      sHeap *Qobj=NULL, *Qbkg=NULL;
      int i,j,p,p_obj,p_bkg,q,n,lp;
      int l_ant;
      float e_obj,e_bkg;
      int wi;
      float w;
      float *value;
      int u_x,u_y,v_x,v_y;
      sAdjRel *A;
      int *i_inv;
      int *Q=NULL;
      int Qtop = -1;
      //----------
      n = label->n;
      value = gft::AllocFloatArray(n);
      Qobj = Heap::Create(n, value);
      Qbkg = Heap::Create(n, value);
      Q	= gft::AllocIntArray(n);

      A = sg->A;
      
      i_inv = gft::AdjRel::InverseIndexes(A);
      
      if(S != NULL){
	for(i=1; i<=S[0]; i++){
	  value[S[i]] = 0.0;
	  if(label->data[S[i]] == 0)
	    Heap::Insert_MinPolicy(Qbkg, S[i]);
	  else if(label->data[S[i]] != NIL)
	    Heap::Insert_MinPolicy(Qobj, S[i]);
	}
	for(p=0; p<n; p++){
	  if(label->data[p]==NIL){
	    value[p] = FLT_MAX;
	    label->data[p] = 0;
	  }
	  else
	    value[p] = 0.0;
	}
      }
      else{
	for(p=0; p<n; p++){
	  if(label->data[p] == 0){
	    value[p] = 0.0;
	    Heap::Insert_MinPolicy(Qbkg, p);
	  }
	  else if(label->data[p] != NIL){
	    value[p] = 0.0;
	    Heap::Insert_MinPolicy(Qobj, p);
	  }
	  else{
	    value[p] = FLT_MAX;
	    label->data[p] = 0;
	  }
	}
      }

      l_ant = 0;
      while(Qobj->last > 0 && Qbkg->last > 0){
      //---------------------
	p_obj = Qobj->pixel[1];
	p_bkg = Qbkg->pixel[1];
	
	e_obj = value[p_obj];
	e_bkg = value[p_bkg];
	if(e_obj < e_bkg){
	  e_max = MAX(e_max, e_bkg);
	  Heap::Remove_MinPolicy(Qobj, &p);
	}
	else if(e_obj > e_bkg){
	  e_max = MAX(e_max, e_obj);
	  Heap::Remove_MinPolicy(Qbkg, &p);
	}
	else{
	  e_max = MAX(e_max, e_obj);
          if(l_ant == 0){
	    Heap::Remove_MinPolicy(Qobj, &p);
	  }
          else{
	    Heap::Remove_MinPolicy(Qbkg, &p);
	  }
          l_ant = 1 - l_ant;
	}

	Qobj->color[p] = BLACK;
	Qbkg->color[p] = BLACK;

	lp = label->data[p];
	
	goto label01;
	  
	Qtop++;
	Q[Qtop] = p;
	while(Qtop > -1){
	  p = Q[Qtop];
	  Qtop--;
	  
	label01:
	  
	  u_x = p%label->ncols;
	  u_y = p/label->ncols;
	  for(i=1; i<A->n; i++){
	    v_x = u_x + A->dx[i];
	    v_y = u_y + A->dy[i];
	    if(v_x >= 0 && v_x < label->ncols &&
	       v_y >= 0 && v_y < label->nrows){
	      q = v_x + label->ncols*v_y;
	      if(Qobj->color[q] != BLACK){
		
		if(lp == 0){
		  j = i_inv[i];
		  wi = (sg->n_link[q])[j];
		}
		else
		  wi = (sg->n_link[p])[i];

		if(wi == INT_MAX) continue;

		w = (float)wi;
		
                if(w < e_max){
		  if(Qobj->color[q] == GRAY)
		    Heap::Delete_MinPolicy(Qobj, q);
		  else if(Qbkg->color[q] == GRAY)
		    Heap::Delete_MinPolicy(Qbkg, q);
		  Qobj->color[q] = BLACK;
		  Qbkg->color[q] = BLACK;
		  label->data[q] = lp;
		  Qtop++;
		  Q[Qtop] = q;
		}
		else if(w < value[q]){
		  label->data[q] = lp;
		  if(lp > 0){
		    if(Qbkg->color[q] == GRAY)
		      Heap::Delete_MinPolicy(Qbkg, q);
		    Heap::Update_MinPolicy(Qobj, q, w);
		  }
		  else{
		    if(Qobj->color[q] == GRAY)
		      Heap::Delete_MinPolicy(Qobj, q);
		    Heap::Update_MinPolicy(Qbkg, q, w);
		  }
		}
	      }
	    }
	  }
	}
      }

      for(i = 1; i <= Qobj->last; i++){
	p = Qobj->pixel[i];
	Qtop++;
	Q[Qtop] = p;
	Qobj->color[p] = BLACK;
      }
      
      while(Qtop > -1){
	p = Q[Qtop];
	Qtop--;
	u_x = p%label->ncols;
	u_y = p/label->ncols;
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  if(v_x >= 0 && v_x < label->ncols &&
	     v_y >= 0 && v_y < label->nrows){
	    q = v_x + label->ncols*v_y;
	    if(Qobj->color[q] != BLACK){
	      Qobj->color[q] = BLACK;
	      label->data[q] = label->data[p];
	      Qtop++;
	      Q[Qtop] = q;
	    }
	  }
	}
      }
      
      gft::FreeFloatArray(&value);
      Heap::Destroy(&Qobj);
      Heap::Destroy(&Qbkg);
      gft::FreeIntArray(&Q);
      free(i_inv);
    }
    */
    //-----------------------------------------------------


    void EOIFT_Heap(sImageGraph *sg,
		    int *S,
		    sImage32 *label,
		    float e_max){
      sHeapPair *QP=NULL;
      int i,j,p,p_obj,p_bkg,q,n,lp;
      int l_ant;
      float e_obj,e_bkg;
      int wi;
      float w;
      float *value;
      int u_x,u_y,v_x,v_y;
      sAdjRel *A;
      int *i_inv;
      int *Q=NULL;
      int Qtop = -1;
      n = label->n;
      value = gft::AllocFloatArray(n);
      QP = HeapPair::Create(n, value);
      Q	= gft::AllocIntArray(n);
      A = sg->A;
      
      i_inv = gft::AdjRel::InverseIndexes(A);
      
      if(S != NULL){
	for(i=1; i<=S[0]; i++){
	  value[S[i]] = 0.0;
	  if(label->data[S[i]] == 0)
	    HeapPair::Insert_MinPolicy_0(QP, S[i]);
	  else if(label->data[S[i]] != NIL)
	    HeapPair::Insert_MinPolicy_1(QP, S[i]);
	}
	for(p=0; p<n; p++){
	  if(label->data[p]==NIL){
	    value[p] = FLT_MAX;
	    label->data[p] = 0;
	  }
	  else
	    value[p] = 0.0;
	}
      }
      else{
	for(p=0; p<n; p++){
	  if(label->data[p] == 0){
	    value[p] = 0.0;
	    HeapPair::Insert_MinPolicy_0(QP, p);
	  }
	  else if(label->data[p] != NIL){
	    value[p] = 0.0;
	    HeapPair::Insert_MinPolicy_1(QP, p);
	  }
	  else{
	    value[p] = FLT_MAX;
	    label->data[p] = 0;
	  }
	}
      }

      l_ant = 0;
      while(QP->last_0 > 0 && QP->last_1 <= QP->n){
	p_obj = QP->pixel[QP->n];
	p_bkg = QP->pixel[1];
	
	e_obj = value[p_obj];
	e_bkg = value[p_bkg];
	if(e_obj < e_bkg){
	  e_max = MAX(e_max, e_bkg);
	  HeapPair::Remove_MinPolicy_1(QP, &p);
	}
	else if(e_obj > e_bkg){
	  e_max = MAX(e_max, e_obj);
	  HeapPair::Remove_MinPolicy_0(QP, &p);
	}
	else{
	  e_max = MAX(e_max, e_obj);
          if(l_ant == 0)
	    HeapPair::Remove_MinPolicy_1(QP, &p);
          else
	    HeapPair::Remove_MinPolicy_0(QP, &p);
          l_ant = 1 - l_ant;
	}

	QP->color[p] = BLACK;

	lp = label->data[p];
	
	goto label03;
	  
	Qtop++;
	Q[Qtop] = p;
	while(Qtop > -1){
	  p = Q[Qtop];
	  Qtop--;
	  
	label03:
	  
	  u_x = p%label->ncols;
	  u_y = p/label->ncols;
	  for(i=1; i<A->n; i++){
	    v_x = u_x + A->dx[i];
	    v_y = u_y + A->dy[i];
	    if(v_x >= 0 && v_x < label->ncols &&
	       v_y >= 0 && v_y < label->nrows){
	      q = v_x + label->ncols*v_y;
	      if(QP->color[q] != BLACK){
		
		if(lp == 0){
		  j = i_inv[i];
		  wi = (sg->n_link[q])[j];
		}
		else
		  wi = (sg->n_link[p])[i];

		if(wi == INT_MAX) continue;

		w = (float)wi;
		
                if(w < e_max){
		  if(QP->color[q] == GRAY)
		    HeapPair::Delete_MinPolicy(QP, q);
		  QP->color[q] = BLACK;
		  label->data[q] = lp;
		  Qtop++;
		  Q[Qtop] = q;
		}
		else if(w < value[q]){
		  label->data[q] = lp;
		  if(lp > 0)
		    HeapPair::Update_MinPolicy_1(QP, q, w);
		  else
		    HeapPair::Update_MinPolicy_0(QP, q, w);
		}
	      }
	    }
	  }
	}
      }

      for(i = QP->last_1; i <= QP->n; i++){
	p = QP->pixel[i];
	Qtop++;
	Q[Qtop] = p;
	QP->color[p] = BLACK;
      }
      
      while(Qtop > -1){
	p = Q[Qtop];
	Qtop--;
	u_x = p%label->ncols;
	u_y = p/label->ncols;
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  if(v_x >= 0 && v_x < label->ncols &&
	     v_y >= 0 && v_y < label->nrows){
	    q = v_x + label->ncols*v_y;
	    if(QP->color[q] != BLACK){
	      QP->color[q] = BLACK;
	      label->data[q] = label->data[p];
	      Qtop++;
	      Q[Qtop] = q;
	    }
	  }
	}
      }
      
      gft::FreeFloatArray(&value);
      HeapPair::Destroy(&QP);
      gft::FreeIntArray(&Q);
      free(i_inv);
    }
    

    //-----------------------------------------------------
    
    void EOIFT_Heap_2(sImageGraph *sg,
		      int *S,
		      sImage32 *label,
		      float e_max){
      sHeapPair *QP=NULL;
      int i,j,p,p_obj,p_bkg,q,n,lp;
      int l_ant;
      float e_obj,e_bkg;
      int wi;
      float w;
      float *value;
      int u_x,u_y,v_x,v_y;
      sAdjRel *A;
      int *i_inv;
      int *Q=NULL;
      int Qtop = -1;
      int *T=NULL;
      float *Tv = NULL;
      int Ttop = -1;
      n = label->n;
      value = gft::AllocFloatArray(n);
      QP = HeapPair::Create(n, value);
      Q	= gft::AllocIntArray(n);
      T	= gft::AllocIntArray(n);
      Tv = gft::AllocFloatArray(n);

      A = sg->A;
      i_inv = gft::AdjRel::InverseIndexes(A);
      
      if(S != NULL){
	for(i=1; i<=S[0]; i++){
	  Tv[S[i]] = value[S[i]] = 0.0;
	  if(label->data[S[i]] == 0)
	    HeapPair::Insert_MinPolicy_0(QP, S[i]);
	  else if(label->data[S[i]] != NIL)
	    HeapPair::Insert_MinPolicy_1(QP, S[i]);
	}
	for(p=0; p<n; p++){
	  if(label->data[p]==NIL){
	    Tv[p] = value[p] = FLT_MAX;
	    label->data[p] = 0;
	  }
	  else
	    Tv[p] = value[p] = 0.0;
	}
      }
      else{
	for(p=0; p<n; p++){
	  if(label->data[p] == 0){
	    Tv[p] = value[p] = 0.0;
	    HeapPair::Insert_MinPolicy_0(QP, p);
	  }
	  else if(label->data[p] != NIL){
	    Tv[p] = value[p] = 0.0;
	    HeapPair::Insert_MinPolicy_1(QP, p);
	  }
	  else{
	    Tv[p] = value[p] = FLT_MAX;
	    label->data[p] = 0;
	  }
	}
      }

      l_ant = 0;
      while(QP->last_0 > 0 && QP->last_1 <= QP->n){
	p_obj = QP->pixel[QP->n];
	p_bkg = QP->pixel[1];
	
	e_obj = value[p_obj];
	e_bkg = value[p_bkg];
	if(e_obj < e_bkg){
	  e_max = MAX(e_max, e_bkg);
	  HeapPair::Remove_MinPolicy_1(QP, &p);
	  //------------------
	  //p_obj = QP->pixel[QP->n];
	  //while(QP->last_1 <= QP->n && value[p_obj] < e_max){
	  //  HeapPair::Remove_MinPolicy_1(QP, &p_obj);
	  //  Qtop++;
	  //  Q[Qtop] = p_obj;
	  //  p_obj = QP->pixel[QP->n];
	  //}
	  //------------------
	}
	else if(e_obj > e_bkg){
	  e_max = MAX(e_max, e_obj);
	  HeapPair::Remove_MinPolicy_0(QP, &p);
	  //------------------
	  //p_bkg = QP->pixel[1];
	  //while(QP->last_0 > 0 && value[p_bkg] < e_max){
	  //  HeapPair::Remove_MinPolicy_0(QP, &p_bkg);
	  //  Qtop++;
	  //  Q[Qtop] = p_bkg;
	  //  p_bkg = QP->pixel[1];
	  //}
	  //------------------
	}
	else{
	  e_max = MAX(e_max, e_obj);
          if(l_ant == 0) HeapPair::Remove_MinPolicy_1(QP, &p);
          else           HeapPair::Remove_MinPolicy_0(QP, &p);
          l_ant = 1 - l_ant;
	  //------------------
	  //if(e_obj < e_max){
	  //  p_obj = QP->pixel[QP->n];
	  //  while(QP->last_1 <= QP->n && value[p_obj] < e_max){
	  //    HeapPair::Remove_MinPolicy_1(QP, &p_obj);
	  //    Qtop++;
	  //    Q[Qtop] = p_obj;
	  //    p_obj = QP->pixel[QP->n];
	  //  }
	  //  p_bkg = QP->pixel[1];
	  //  while(QP->last_0 > 0 && value[p_bkg] < e_max){
	  //    HeapPair::Remove_MinPolicy_0(QP, &p_bkg);
	  //    Qtop++;
	  //    Q[Qtop] = p_bkg;
	  //    p_bkg = QP->pixel[1];
	  //  }
	  //}
	  //else e_max = e_obj;
	  //------------------
	}

	//QP->color[p] = BLACK;

	lp = label->data[p];
	
	goto label01;
	  
	Qtop++;
	Q[Qtop] = p;
	while(Qtop > -1){
	  p = Q[Qtop];
	  Qtop--;
	  
	label01:
	  
	  u_x = p%label->ncols;
	  u_y = p/label->ncols;
	  for(i=1; i<A->n; i++){
	    v_x = u_x + A->dx[i];
	    v_y = u_y + A->dy[i];
	    if(v_x >= 0 && v_x < label->ncols &&
	       v_y >= 0 && v_y < label->nrows){
	      q = v_x + label->ncols*v_y;
	      if(QP->color[q] != BLACK){
		
		if(lp == 0){
		  j = i_inv[i];
		  wi = (sg->n_link[q])[j];
		}
		else
		  wi = (sg->n_link[p])[i];

		if(wi == INT_MAX) continue;

		w = (float)wi;
		
                if(w < e_max){
		  if(QP->color[q] == GRAY)
		    HeapPair::Delete_MinPolicy(QP, q);
		  QP->color[q] = BLACK;
		  label->data[q] = lp;
		  Qtop++;
		  Q[Qtop] = q;
		}
		else if(w < Tv[q]){ //if(w < value[q]){
		  if(Tv[q] >= value[q]){
		    label->data[q] = lp;
		    Ttop++;
		    T[Ttop] = q;
		  }
		  Tv[q] = w;
		}
	      }
	    }
	  }
	}
	while(Ttop > -1){
	  p = T[Ttop];
	  w = Tv[p];
	  Ttop--;
	  if(QP->color[p] != BLACK){
	    if(label->data[p] > 0)
	      HeapPair::Update_MinPolicy_1(QP, p, w);
	    else
	      HeapPair::Update_MinPolicy_0(QP, p, w);
	  }
	}
      }

      for(i = QP->last_1; i <= QP->n; i++){
	p = QP->pixel[i];
	Qtop++;
	Q[Qtop] = p;
	QP->color[p] = BLACK;
      }
      
      while(Qtop > -1){
	p = Q[Qtop];
	Qtop--;
	u_x = p%label->ncols;
	u_y = p/label->ncols;
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  if(v_x >= 0 && v_x < label->ncols &&
	     v_y >= 0 && v_y < label->nrows){
	    q = v_x + label->ncols*v_y;
	    if(QP->color[q] != BLACK){
	      QP->color[q] = BLACK;
	      label->data[q] = label->data[p];
	      Qtop++;
	      Q[Qtop] = q;
	    }
	  }
	}
      }
      
      gft::FreeFloatArray(&value);
      HeapPair::Destroy(&QP);
      gft::FreeIntArray(&Q);
      gft::FreeIntArray(&T);
      gft::FreeFloatArray(&Tv);
      free(i_inv);
    }


    //------------------------------------
    /*
    void EOIFT_Heap(sGraph *graph,
		    sGraph *transpose,
		    int *S,
		    int *label){
      sHeap *Qobj=NULL, *Qbkg=NULL;
      sGraph *g;
      int i,j,p,p_obj,p_bkg,q,n;
      int l_ant;
      float e_obj,e_bkg,e_max;
      int wi;
      float w;
      float *value;
      int u_x,u_y,v_x,v_y;
      sStack *Q=NULL;
      
      n = graph->nnodes;
      value = gft::AllocFloatArray(n);
      Qobj = Heap::Create(n, value);
      Qbkg = Heap::Create(n, value);
      Q	= Stack::Create(n);

      for(p = 0; p < n; p++){
	if(label[p]==NIL) value[p] = FLT_MAX;
	else              value[p] = 0.0;
      }
      
      if(S != NULL){
	for(i = 1; i <= S[0]; i++)
	  if(label[S[i]] == 0)
	    Heap::Insert_MinPolicy(Qbkg, S[i]);
	  else if(label[S[i]] != NIL)
	    Heap::Insert_MinPolicy(Qobj, S[i]);
      }
      else{
	for(p = 0; p < n; p++)
	  if(label[p] == 0)
	    Heap::Insert_MinPolicy(Qbkg, p);
	  else if(label[p] != NIL)
	    Heap::Insert_MinPolicy(Qobj, p);
      }

      l_ant = 0;
      while(!Heap::IsEmpty(Qobj) && !Heap::IsEmpty(Qbkg)) {
	Heap::Get_MinPolicy(Qobj, &p_obj);
	Heap::Get_MinPolicy(Qbkg, &p_bkg);

	e_obj = value[p_obj];
	e_bkg = value[p_bkg];
	if(e_obj < e_bkg){
	  e_max = e_bkg;
	  p = p_obj;
	}
	else if(e_obj > e_bkg){
	  e_max = e_obj;
	  p = p_bkg;
	}
	else{
	  e_max = e_obj;
          if(l_ant == 0)
	    p = p_obj;
          else
	    p = p_bkg;
          l_ant = 1 - l_ant;
	}

	if(Qobj->color[p] == GRAY)
	  Heap::Delete_MinPolicy(Qobj, p);
	if(Qbkg->color[p] == GRAY)
	  Heap::Delete_MinPolicy(Qbkg, p);
	Qobj->color[p] = BLACK;
	Qbkg->color[p] = BLACK;

	Stack::Push(Q, p);
	while(!Stack::IsEmpty(Q)){
	  p = Stack::Pop(Q);
	  
	  if(label[p]==0) g = transpose;
	  else   	  g = graph;
	  
	  for(i = 0; i < g->nodes[p].outdegree; i++){
	    q = g->nodes[p].adjList[i];

	    if(Qobj->color[q] != BLACK){

	      
	      //if(label[p]==0)
		//wi = Graph::GetArcWeight(graph, q, p);
	      //else
	      wi = g->nodes[p].Warcs[i];

	      if(wi == INT_MAX)
		w = FLT_MAX;
	      else
		w = (float)wi;

	      if(w < e_max){
		if(Qobj->color[q] == GRAY)
		  Heap::Delete_MinPolicy(Qobj, q);
		if(Qbkg->color[q] == GRAY)
		  Heap::Delete_MinPolicy(Qbkg, q);
		Qobj->color[q] = BLACK;
		Qbkg->color[q] = BLACK;
		
		label[q] = label[p];
		Stack::Push(Q, q);
	      }
	      else if(w < value[q]){
		if(Qobj->color[q] == GRAY)
		  Heap::Delete_MinPolicy(Qobj, q);
		if(Qbkg->color[q] == GRAY)
		  Heap::Delete_MinPolicy(Qbkg, q);
		value[q] = w;
		label[q] = label[p];
		if(label[q] > 0)
		  Heap::Insert_MinPolicy(Qobj, q);
		else
		  Heap::Insert_MinPolicy(Qbkg, q);		  
	      }
	      
	    }
	  }

	}
	
      }

      while(!Heap::IsEmpty(Qobj)){
	Heap::Remove_MinPolicy(Qobj, &p);
	Stack::Push(Q, p);
	while(!Stack::IsEmpty(Q)){
	  p = Stack::Pop(Q);

	  for(i = 0; i < graph->nodes[p].outdegree; i++){
	    q = graph->nodes[p].adjList[i];
	    
	    if(Qobj->color[q] != BLACK){

	      if(Qobj->color[q] == GRAY)
		Heap::Delete_MinPolicy(Qobj, q);
	      Qobj->color[q] = BLACK;
	      label[q] = label[p];
	      Stack::Push(Q, q);
	    }
	  }
	}
      }

      while(!Heap::IsEmpty(Qbkg)){
	Heap::Remove_MinPolicy(Qbkg, &p);
	Stack::Push(Q, p);
	Qobj->color[p] = BLACK;
	while(!Stack::IsEmpty(Q)){
	  p = Stack::Pop(Q);

	  for(i = 0; i < transpose->nodes[p].outdegree; i++){
	    q = transpose->nodes[p].adjList[i];
	    
	    if(Qobj->color[q] != BLACK){
	      
	      if(Qbkg->color[q] == GRAY)
		Heap::Delete_MinPolicy(Qbkg, q);
	      Qobj->color[q] = BLACK;
	      label[q] = label[p];
	      Stack::Push(Q, q);
	    }
	  }
	}
      }

      gft::FreeFloatArray(&value);
      Heap::Destroy(&Qobj);
      Heap::Destroy(&Qbkg);
      Stack::Destroy(&Q);
    }
    */



    /*    
    void EOIFT_Heap(sGraph *graph,
		    sGraph *transpose,
		    int *S,
		    int *label,
		    float e_max){
      sHeap *Qobj=NULL, *Qbkg=NULL;
      sGraph *g;
      int i,j,p,p_obj,p_bkg,q,n,lp;
      int l_ant;
      float e_obj,e_bkg;
      int wi;
      float w;
      float *value;
      int u_x,u_y,v_x,v_y;
      int *Q=NULL;
      int Qtop = -1;
      
      n = graph->nnodes;
      value = gft::AllocFloatArray(n);
      Qobj = Heap::Create(n, value);
      Qbkg = Heap::Create(n, value);
      Q	= gft::AllocIntArray(n);
      
      if(S != NULL){
	for(i=1; i<=S[0]; i++){
	  value[S[i]] = 0.0;
	  if(label[S[i]] == 0)
	    Heap::Insert_MinPolicy(Qbkg, S[i]);
	  else if(label[S[i]] != NIL)
	    Heap::Insert_MinPolicy(Qobj, S[i]);
	}
	for(p=0; p<n; p++){
	  if(label[p]==NIL){
	    value[p] = FLT_MAX;
	    label[p] = 0;
	  }
	  else
	    value[p] = 0.0;
	}
      }
      else{
	for(p=0; p<n; p++){
	  if(label[p] == 0){
	    value[p] = 0.0;
	    Heap::Insert_MinPolicy(Qbkg, p);
	  }
	  else if(label[p] != NIL){
	    value[p] = 0.0;
	    Heap::Insert_MinPolicy(Qobj, p);
	  }
	  else{
	    value[p] = FLT_MAX;
	    label[p] = 0;
	  }
	}
      }

      l_ant = 0;
      while(Qobj->last > 0 && Qbkg->last > 0){
	p_obj = Qobj->pixel[1];
	p_bkg = Qbkg->pixel[1];
	
	e_obj = value[p_obj];
	e_bkg = value[p_bkg];
	if(e_obj < e_bkg){
	  e_max = MAX(e_max, e_bkg);
	  Heap::Remove_MinPolicy(Qobj, &p);
	}
	else if(e_obj > e_bkg){
	  e_max = MAX(e_max, e_obj);
	  Heap::Remove_MinPolicy(Qbkg, &p);
	}
	else{
	  e_max = MAX(e_max, e_obj);
          if(l_ant == 0)
	    Heap::Remove_MinPolicy(Qobj, &p);
          else
	    Heap::Remove_MinPolicy(Qbkg, &p);
          l_ant = 1 - l_ant;
	}

	Qobj->color[p] = BLACK;
	Qbkg->color[p] = BLACK;

	lp = label[p];
	
	if(lp==0) g = transpose;
	else   	  g = graph;
	
	goto label02;
	  
	Qtop++;
	Q[Qtop] = p;
	while(Qtop > -1){
	  p = Q[Qtop];
	  Qtop--;
	  
	
	  
	  for(i = 0; i < g->nodes[p].outdegree; i++){
	    q = g->nodes[p].adjList[i];

	    if(Qobj->color[q] != BLACK){
	      
	      wi = g->nodes[p].Warcs[i];

	      if(wi == INT_MAX) continue;

	      w = (float)wi;

	      if(w < e_max){
		if(Qobj->color[q] == GRAY)
		  Heap::Delete_MinPolicy(Qobj, q);
		else if(Qbkg->color[q] == GRAY)
		  Heap::Delete_MinPolicy(Qbkg, q);
		Qobj->color[q] = BLACK;
		Qbkg->color[q] = BLACK;
		label[q] = lp;
		Qtop++;
		Q[Qtop] = q;
	      }
	      else if(w < value[q]){
		label[q] = lp;
		if(lp > 0){
		  if(Qbkg->color[q] == GRAY)
		    Heap::Delete_MinPolicy(Qbkg, q);
		  Heap::Update_MinPolicy(Qobj, q, w);
		}
		else{
		  if(Qobj->color[q] == GRAY)
		    Heap::Delete_MinPolicy(Qobj, q);
		  Heap::Update_MinPolicy(Qbkg, q, w);
		}
	      }
	      
	    }
	  }

	}
	
      }

      for(i = 1; i <= Qobj->last; i++){
	p = Qobj->pixel[i];
	Qtop++;
	Q[Qtop] = p;
	Qobj->color[p] = BLACK;
      }

      while(Qtop > -1){
	p = Q[Qtop];
	Qtop--;      
	for(i = 0; i < graph->nodes[p].outdegree; i++){
	  q = graph->nodes[p].adjList[i];
	  if(Qobj->color[q] != BLACK){
	    Qobj->color[q] = BLACK;
	    label[q] = label[p];
	    Qtop++;
	    Q[Qtop] = q;
	  }
	}
      }
	
      gft::FreeFloatArray(&value);
      Heap::Destroy(&Qobj);
      Heap::Destroy(&Qbkg);
      gft::FreeIntArray(&Q);
    }
    */



    void EOIFT_Heap(sGraph *graph,
		    sGraph *transpose,
		    int *S,
		    int *label,
		    float e_max){
      sHeapPair *QP=NULL;
      sGraph *g;
      int i,j,p,p_obj,p_bkg,q,n,lp;
      int l_ant;
      float e_obj,e_bkg;
      int wi;
      float w;
      float *value;
      int *Q=NULL;
      int Qtop = -1;
      n = graph->nnodes;
      value = gft::AllocFloatArray(n);
      QP = HeapPair::Create(n, value);
      Q	= gft::AllocIntArray(n);
      
      if(S != NULL){
	for(i=1; i<=S[0]; i++){
	  value[S[i]] = 0.0;
	  if(label[S[i]] == 0)
	    HeapPair::Insert_MinPolicy_0(QP, S[i]);
	  else if(label[S[i]] != NIL)
	    HeapPair::Insert_MinPolicy_1(QP, S[i]);
	}
	for(p=0; p<n; p++){
	  if(label[p]==NIL){
	    value[p] = FLT_MAX;
	    label[p] = 0;
	  }
	  else
	    value[p] = 0.0;
	}
      }
      else{
	for(p=0; p<n; p++){
	  if(label[p] == 0){
	    value[p] = 0.0;
	    HeapPair::Insert_MinPolicy_0(QP, p);
	  }
	  else if(label[p] != NIL){
	    value[p] = 0.0;
	    HeapPair::Insert_MinPolicy_1(QP, p);
	  }
	  else{
	    value[p] = FLT_MAX;
	    label[p] = 0;
	  }
	}
      }

      l_ant = 0;
      while(QP->last_0 > 0 && QP->last_1 <= QP->n){
	p_obj = QP->pixel[QP->n];
	p_bkg = QP->pixel[1];
	
	e_obj = value[p_obj];
	e_bkg = value[p_bkg];
	if(e_obj < e_bkg){
	  e_max = MAX(e_max, e_bkg);
	  HeapPair::Remove_MinPolicy_1(QP, &p);
	}
	else if(e_obj > e_bkg){
	  e_max = MAX(e_max, e_obj);
	  HeapPair::Remove_MinPolicy_0(QP, &p);
	}
	else{
	  e_max = MAX(e_max, e_obj);
          if(l_ant == 0) HeapPair::Remove_MinPolicy_1(QP, &p);
          else   	 HeapPair::Remove_MinPolicy_0(QP, &p);
          l_ant = 1 - l_ant;
	}

	QP->color[p] = BLACK;

	lp = label[p];
	
	if(lp==0) g = transpose;
	else   	  g = graph;
	
	goto label06;
	  
	Qtop++;
	Q[Qtop] = p;
	while(Qtop > -1){
	  p = Q[Qtop];
	  Qtop--;
	  
	label06:
	  
	  for(i = 0; i < g->nodes[p].outdegree; i++){
	    q = g->nodes[p].adjList[i];
	    if(QP->color[q] != BLACK){
	      
	      wi = g->nodes[p].Warcs[i];

	      if(wi == INT_MAX) continue;

	      w = (float)wi;

	      if(w < e_max){
		if(QP->color[q] == GRAY)
		  HeapPair::Delete_MinPolicy(QP, q);
		QP->color[q] = BLACK;
		label[q] = lp;
		Qtop++;
		Q[Qtop] = q;
	      }
	      else if(w < value[q]){
		label[q] = lp;
		if(lp > 0)
		  HeapPair::Update_MinPolicy_1(QP, q, w);
		else
		  HeapPair::Update_MinPolicy_0(QP, q, w);
	      }
	    }
	  }
	}
      }

      for(i = QP->last_1; i <= QP->n; i++){
	p = QP->pixel[i];
	Qtop++;
	Q[Qtop] = p;
	QP->color[p] = BLACK;
      }

      while(Qtop > -1){
	p = Q[Qtop];
	Qtop--;      
	for(i = 0; i < graph->nodes[p].outdegree; i++){
	  q = graph->nodes[p].adjList[i];
	  if(QP->color[q] != BLACK){
	    QP->color[q] = BLACK;
	    label[q] = label[p];
	    Qtop++;
	    Q[Qtop] = q;
	  }
	}
      }
	
      gft::FreeFloatArray(&value);
      HeapPair::Destroy(&QP);
      gft::FreeIntArray(&Q);
    }
    

    
    void EOIFT_Heap_2(sGraph *graph,
		      sGraph *transpose,
		      int *S,
		      int *label,
		      float e_max){
      sHeapPair *QP=NULL;
      sGraph *g;
      int i,j,p,p_obj,p_bkg,q,n,lp;
      int l_ant;
      float e_obj,e_bkg;
      int wi;
      float w;
      float *value;
      int *Q=NULL;
      int Qtop = -1;
      int *T=NULL;
      float *Tv = NULL;
      int Ttop = -1;
      n = graph->nnodes;
      value = gft::AllocFloatArray(n);
      QP = HeapPair::Create(n, value);
      Q	= gft::AllocIntArray(n);
      T	= gft::AllocIntArray(n);
      Tv = gft::AllocFloatArray(n);
      
      if(S != NULL){
	for(i=1; i<=S[0]; i++){
	  Tv[S[i]] = value[S[i]] = 0.0;
	  if(label[S[i]] == 0)
	    HeapPair::Insert_MinPolicy_0(QP, S[i]);
	  else if(label[S[i]] != NIL)
	    HeapPair::Insert_MinPolicy_1(QP, S[i]);
	}
	for(p=0; p<n; p++){
	  if(label[p]==NIL){
	    Tv[p] = value[p] = FLT_MAX;
	    label[p] = 0;
	  }
	  else
	    Tv[p] = value[p] = 0.0;
	}
      }
      else{
	for(p=0; p<n; p++){
	  if(label[p] == 0){
	    Tv[p] = value[p] = 0.0;
	    HeapPair::Insert_MinPolicy_0(QP, p);
	  }
	  else if(label[p] != NIL){
	    Tv[p] = value[p] = 0.0;
	    HeapPair::Insert_MinPolicy_1(QP, p);
	  }
	  else{
	    Tv[p] = value[p] = FLT_MAX;
	    label[p] = 0;
	  }
	}
      }

      l_ant = 0;
      while(QP->last_0 > 0 && QP->last_1 <= QP->n){
	p_obj = QP->pixel[QP->n];
	p_bkg = QP->pixel[1];
	
	e_obj = value[p_obj];
	e_bkg = value[p_bkg];
	if(e_obj < e_bkg){
	  e_max = MAX(e_max, e_bkg);
	  HeapPair::Remove_MinPolicy_1(QP, &p);
	}
	else if(e_obj > e_bkg){
	  e_max = MAX(e_max, e_obj);
	  HeapPair::Remove_MinPolicy_0(QP, &p);
	}
	else{
	  e_max = MAX(e_max, e_obj);
          if(l_ant == 0) HeapPair::Remove_MinPolicy_1(QP, &p);
          else   	 HeapPair::Remove_MinPolicy_0(QP, &p);
          l_ant = 1 - l_ant;
	}

	//QP->color[p] = BLACK;

	lp = label[p];
	
	if(lp==0) g = transpose;
	else   	  g = graph;
	
	goto label02;
	  
	Qtop++;
	Q[Qtop] = p;
	while(Qtop > -1){
	  p = Q[Qtop];
	  Qtop--;
	  
	label02:
	  
	  for(i = 0; i < g->nodes[p].outdegree; i++){
	    q = g->nodes[p].adjList[i];
	    if(QP->color[q] != BLACK){
	      
	      wi = g->nodes[p].Warcs[i];

	      if(wi == INT_MAX) continue;

	      w = (float)wi;

	      if(w < e_max){
		if(QP->color[q] == GRAY)
		  HeapPair::Delete_MinPolicy(QP, q);
		QP->color[q] = BLACK;
		label[q] = lp;
		Qtop++;
		Q[Qtop] = q;
	      }
	      else if(w < Tv[q]){
		if(Tv[q] >= value[q]){
		  label[q] = lp;
		  Ttop++;
		  T[Ttop] = q;
		}
		Tv[q] = w;
	      }
	    }
	  }
	}
	while(Ttop > -1){
	  p = T[Ttop];
	  w = Tv[p];
	  Ttop--;
	  if(QP->color[p] != BLACK){
	    if(label[p] > 0)
	      HeapPair::Update_MinPolicy_1(QP, p, w);
	    else
	      HeapPair::Update_MinPolicy_0(QP, p, w);
	  }
	}	
      }

      for(i = QP->last_1; i <= QP->n; i++){
	p = QP->pixel[i];
	Qtop++;
	Q[Qtop] = p;
	QP->color[p] = BLACK;
      }

      while(Qtop > -1){
	p = Q[Qtop];
	Qtop--;      
	for(i = 0; i < graph->nodes[p].outdegree; i++){
	  q = graph->nodes[p].adjList[i];
	  if(QP->color[q] != BLACK){
	    QP->color[q] = BLACK;
	    label[q] = label[p];
	    Qtop++;
	    Q[Qtop] = q;
	  }
	}
      }
	
      gft::FreeFloatArray(&value);
      HeapPair::Destroy(&QP);
      gft::FreeIntArray(&Q);
      gft::FreeIntArray(&T);
      gft::FreeFloatArray(&Tv);
    }
    
    
    
    
    void IFT_fmax_Heap(sGraph *graph,
		       int *S,
		       int *label,
		       float *cost){
      sHeap *Q;
      float tmp, w;
      int n,p,q,i;
      n = graph->nnodes;
      Q = gft::Heap::Create(n, cost);

      for(p = 0; p < n; p++){
	if(label[p]==NIL) cost[p] = FLT_MAX;
	else              cost[p] = 0.0;
      }
      
      if(S != NULL){
	for(i = 1; i <= S[0]; i++)
	  gft::Heap::Insert_MinPolicy(Q, S[i]);
      }
      else{
	for(p = 0; p < n; p++)
	  if(label[p]!=NIL)
	    gft::Heap::Insert_MinPolicy(Q, p);
      }
      
      while(!gft::Heap::IsEmpty(Q)){
	gft::Heap::Remove_MinPolicy(Q, &p);
	
	for(i = 0; i < graph->nodes[p].outdegree; i++){
	  q = graph->nodes[p].adjList[i];
	  if(Q->color[q] != BLACK){
	    w = graph->nodes[p].Warcs[i];
	    tmp = MAX(cost[p], w);

	    if(tmp < cost[q]){
	      gft::Heap::Update_MinPolicy(Q, q, tmp);
	      label[q] = label[p];
	    }
	  }
	}
      }
      gft::Heap::Destroy(&Q);
    }




    void IFT_fmax(sGraph *graph,
		  int *S,
		  int *label,
		  int *cost){
      sPQueue32 *Q;
      int tmp, w, Wmax;
      int n,p,q,i;
      n = graph->nnodes;
      Wmax = Graph::GetMaximumArc(graph);
      Q = PQueue32::Create(Wmax+2, n, cost);

      for(p = 0; p < n; p++){
	if(label[p]==NIL) cost[p] = INT_MAX;
	else              cost[p] = 0;
      }
      
      if(S != NULL){
	for(i = 1; i <= S[0]; i++)
	  PQueue32::FastInsertElem(Q, S[i]);
      }
      else{
	for(p = 0; p < n; p++)
	  if(label[p]!=NIL)
	    PQueue32::FastInsertElem(Q, p);
      }
      
      while(!PQueue32::IsEmpty(Q)){
	p = PQueue32::FastRemoveMinFIFO(Q);
	
	for(i = 0; i < graph->nodes[p].outdegree; i++){
	  q = graph->nodes[p].adjList[i];
	  if(Q->L.elem[q].color != BLACK){
	    w = graph->nodes[p].Warcs[i];
	    tmp = MAX(cost[p], w);

	    if(tmp < cost[q]){
	      if(Q->L.elem[q].color == GRAY)
		PQueue32::FastRemoveElem(Q, q);
	      cost[q] = tmp;
	      label[q] = label[p];
	      PQueue32::FastInsertElem(Q, q);
	    }
	  }
	}
      }
      PQueue32::Destroy(&Q);
    }



    void IFT_fw(sGraph *graph,
		int *S,
		int *label,
		int *cost){
      sPQueue32 *Q;
      int tmp, w, Wmax;
      int n,p,q,i;
      n = graph->nnodes;
      Wmax = Graph::GetMaximumArc(graph);
      Q = PQueue32::Create(Wmax+2, n, cost);

      for(p = 0; p < n; p++){
	if(label[p]==NIL) cost[p] = INT_MAX;
	else              cost[p] = 0;
      }
      
      if(S != NULL){
	for(i = 1; i <= S[0]; i++)
	  PQueue32::FastInsertElem(Q, S[i]);
      }
      else{
	for(p = 0; p < n; p++)
	  if(label[p]!=NIL)
	    PQueue32::FastInsertElem(Q, p);
      }
      
      while(!PQueue32::IsEmpty(Q)){
	p = PQueue32::FastRemoveMinFIFO(Q);
	
	for(i = 0; i < graph->nodes[p].outdegree; i++){
	  q = graph->nodes[p].adjList[i];
	  if(Q->L.elem[q].color != BLACK){
	    w = graph->nodes[p].Warcs[i];
	    tmp = w;

	    if(tmp < cost[q]){
	      if(Q->L.elem[q].color == GRAY)
		PQueue32::FastRemoveElem(Q, q);
	      cost[q] = tmp;
	      label[q] = label[p];
	      PQueue32::FastInsertElem(Q, q);
	    }
	  }
	}
      }
      PQueue32::Destroy(&Q);
    }
    
    

    
    void IFT_fw_Heap(sGraph *graph,
		     int *S,
		     int *label,
		     float *cost,
		     int *pred){
      sHeap *Q;
      float tmp, w;
      int n,p,q,i;
      n = graph->nnodes;
      Q = gft::Heap::Create(n, cost);

      for(p = 0; p < n; p++){
	pred[p] = NIL;
	if(label[p]==NIL) cost[p] = FLT_MAX;
	else              cost[p] = 0.0;
      }
      
      if(S != NULL){
	for(i = 1; i <= S[0]; i++)
	  gft::Heap::Insert_MinPolicy(Q, S[i]);
      }
      else{
	for(p = 0; p < n; p++)
	  if(label[p]!=NIL)
	    gft::Heap::Insert_MinPolicy(Q, p);
      }
      
      while(!gft::Heap::IsEmpty(Q)){
	gft::Heap::Remove_MinPolicy(Q, &p);
	
	for(i = 0; i < graph->nodes[p].outdegree; i++){
	  q = graph->nodes[p].adjList[i];
	  if(Q->color[q] != BLACK){
	    w = graph->nodes[p].Warcs[i];
	    tmp = w;

	    if(tmp < cost[q]){
	      gft::Heap::Update_MinPolicy(Q, q, tmp);
	      label[q] = label[p];
	      pred[q] = p;
	    }
	  }
	}
      }
      gft::Heap::Destroy(&Q);
    }
   


    void IFT_fw(sGraph *graph,
		int *S,
		int *label,
		int *cost,
		int *pred){
      sPQueue32 *Q;
      int tmp, w, Wmax;
      int n,p,q,i;
      n = graph->nnodes;
      Wmax = Graph::GetMaximumArc(graph);
      Q = PQueue32::Create(Wmax+2, n, cost);
      
      for(p = 0; p < n; p++){
	pred[p] = NIL;
	if(label[p]==NIL) cost[p] = INT_MAX;
	else              cost[p] = 0;
      }
      
      if(S != NULL){
	for(i = 1; i <= S[0]; i++)
	  PQueue32::FastInsertElem(Q, S[i]);
      }
      else{
	for(p = 0; p < n; p++)
	  if(label[p]!=NIL)
	    PQueue32::FastInsertElem(Q, p);
      }
      
      while(!PQueue32::IsEmpty(Q)){
	p = PQueue32::FastRemoveMinFIFO(Q);
	
	for(i = 0; i < graph->nodes[p].outdegree; i++){
	  q = graph->nodes[p].adjList[i];
	  if(Q->L.elem[q].color != BLACK){
	    w = graph->nodes[p].Warcs[i];
	    tmp = w;

	    if(tmp < cost[q]){
	      if(Q->L.elem[q].color == GRAY)
		PQueue32::FastRemoveElem(Q, q);
	      cost[q] = tmp;
	      label[q] = label[p];
	      pred[q] = p;
	      PQueue32::FastInsertElem(Q, q);
	    }
	  }
	}
      }
      PQueue32::Destroy(&Q);
    }
   
    
    //---------------------------------------
    // Convex IFT:

    /* P_sum = Predecessor map obtained by the IFT fsum.*/
    void SC_conquer_path(int p,
			 sImageGraph *sg,
			 sImage32 *P_sum, 
			 sImage32 *V,
			 sPQueue32 *Q,
			 sImage32 *label){
      int i,q,edge;
      Pixel u,v;
      sAdjRel *A;
      
      A = sg->A;
      do{
	if(Q->L.elem[p].color == GRAY)
	  PQueue32::FastRemoveElem(Q, p);
	Q->L.elem[p].color = BLACK;
	
	label->data[p] = 1;
	
	u.x = p%label->ncols;
	u.y = p/label->ncols;
	
	for(i=1; i<A->n; i++){
	  v.x = u.x + A->dx[i];
	  v.y = u.y + A->dy[i];
	  //if(Image32::IsValidPixel(label, v.x, v.y)){
	  if(v.x >= 0 && v.x < label->ncols &&
	     v.y >= 0 && v.y < label->nrows){
	    q = v.x + label->ncols*v.y;
	    if(Q->L.elem[q].color != BLACK){
	      
	      edge = (sg->n_link[p])[i];
	      
	      if(edge < V->data[q] && q != P_sum->data[p]){
		if(Q->L.elem[q].color == GRAY)
		  PQueue32::FastRemoveElem(Q, q);
		V->data[q] = edge;
		label->data[q] = 1; //label->data[p];
		PQueue32::FastInsertElem(Q, q);
	      }
	    }
	  }
	}
	p = P_sum->data[p];
	if(p == NIL) break;
      }while(Q->L.elem[p].color != BLACK);
    }

    
    /* P_sum = Predecessor map obtained by the IFT fsum.*/
    void SC_prune_tree(int p,
		       sImageGraph *sg,
		       sImage32 *P_sum, 
		       sImage32 *V,
		       sPQueue32 *Q,
		       sQueue *Qfifo,
		       sImage32 *label){
      Pixel u,v;
      int i,q,edge;
      sAdjRel *A = sg->A;
      
      if(Q->L.elem[p].color == GRAY)
	PQueue32::FastRemoveElem(Q, p);
      Q->L.elem[p].color = BLACK;
      
      label->data[p] = 0;
      
      Queue::Push(Qfifo, p);
      
      //printf("Prune tree\n");
      
      while(!Queue::IsEmpty(Qfifo)){
	p = Queue::Pop(Qfifo);
	u.x = p%label->ncols; 
	u.y = p/label->ncols; 
	
	for(i=1; i<A->n; i++){
	  v.x = u.x + A->dx[i];
	  v.y = u.y + A->dy[i];
	  //if(Image32::IsValidPixel(label,v.x,v.y)){
	  if(v.x >= 0 && v.x < label->ncols &&
	     v.y >= 0 && v.y < label->nrows){
	    q = v.x + label->ncols*v.y;
	    if(Q->L.elem[q].color != BLACK){
	      
	      if(P_sum->data[q] == p){
		Queue::Push(Qfifo, q);
		
		if(Q->L.elem[q].color == GRAY)
		  PQueue32::FastRemoveElem(Q, q);
		Q->L.elem[q].color = BLACK;
		
		label->data[q] = 0;
	      }
	      else{
		edge = (sg->n_link[p])[i];
		if(edge < V->data[q]){
		  if(Q->L.elem[q].color == GRAY)
		    PQueue32::FastRemoveElem(Q, q);
		  V->data[q] = edge;
		  label->data[q] = 0;
		  PQueue32::FastInsertElem(Q, q);
		}
	      }
	      
	    }
	  }
	}
	
      }
    }

    
    void SC_IFT(sImageGraph *sg,
		int *S,
		sImage32 *label,
		sImage32 *P_sum){
      sPQueue32 *Q=NULL;
      sQueue *Qfifo=NULL;
      sImage32 *V;
      int p,n,i;
      
      n = sg->ncols*sg->nrows;
      V = Image32::Create(sg->ncols, sg->nrows);
      Q = PQueue32::Create(sg->Wmax+2, n, V->data);
      Qfifo = Queue::Create(n);
      
      for(p = 0; p < n; p++){
	if(label->data[p]==NIL) V->data[p] = INT_MAX;
	else                    V->data[p] = 0;
      }

      for(i = 1; i <= S[0]; i++)
	PQueue32::FastInsertElem(Q, S[i]);
      
      while(!PQueue32::IsEmpty(Q)){
	p = PQueue32::FastRemoveMinFIFO(Q);
	
	if(label->data[p] > 0)
	  SC_conquer_path(p, sg, P_sum, V, Q, label);
	else if(label->data[p] == 0)
	  SC_prune_tree(p, sg, P_sum, V, Q, Qfifo, label);
      }
      
      Image32::Destroy(&V);
      Queue::Destroy(&Qfifo);
      PQueue32::Destroy(&Q);
    }


    //----------------------------------------

    sImage32 *Cost_fmin(sImageGraph *sg,
			int *S, int lb,
			sImage32 *label){
      sPQueue32 *Q = NULL; 
      sImage32 *V;
      int i,p,q,n, edge,tmp;
      Pixel u,v;
      sAdjRel *A;
      
      n = sg->ncols*sg->nrows;
      V = Image32::Create(sg->ncols, sg->nrows);
      Q = PQueue32::Create(sg->Wmax+2, n, V->data);
      A = sg->A;
      
      for(p=0; p<n; p++){
	if(label->data[p] == lb) V->data[p] = sg->Wmax+1;
	else                     V->data[p] = INT_MIN;
      }

      if(S != NULL){
	for(i = 1; i <= S[0]; i++){
	  if(label->data[S[i]] == lb)
	    PQueue32::FastInsertElem(Q, S[i]);
	}
      }
      else{
	for(p=0; p<n; p++)
	  if(label->data[p] == lb)
	    PQueue32::FastInsertElem(Q, p);
      }
      
      while(!PQueue32::IsEmpty(Q)){
	p = PQueue32::FastRemoveMaxFIFO(Q);
	u.x = p%label->ncols; 
	u.y = p/label->ncols; 
	
	for(i=1; i<A->n; i++){
	  v.x = u.x + A->dx[i];
	  v.y = u.y + A->dy[i];
	  //if(Image32::IsValidPixel(label, v.x, v.y)){
	  if(v.x >= 0 && v.x < label->ncols &&
	     v.y >= 0 && v.y < label->nrows){
	    q = v.x + label->ncols*v.y;
	    if(Q->L.elem[q].color != BLACK){
	      edge = (sg->n_link[p])[i];          
	      tmp  = MIN(V->data[p], edge);
	      if(tmp > V->data[q]){
		if(Q->L.elem[q].color == GRAY)
		  PQueue32::FastRemoveElem(Q, q);
		V->data[q] = tmp; //mapa de custos
		PQueue32::FastInsertElem(Q, q);
	      }
	    }
	  }
	}
      }
      PQueue32::Destroy(&Q);
      return V;
    }

    //----------------------------------------
    
    sImage32 *Cost_fmax(sImageGraph *sg,
			int *S, int lb,
			sImage32 *label){
      sPQueue32 *Q = NULL; 
      sImage32 *V;
      int i,p,q,n, edge,tmp;
      Pixel u,v;
      sAdjRel *A;
      
      n = sg->ncols*sg->nrows;
      V = Image32::Create(sg->ncols, sg->nrows);
      Q = PQueue32::Create(sg->Wmax+2, n, V->data);
      A = sg->A;
      
      for(p=0; p<n; p++){
	if(label->data[p] == lb) V->data[p] = 0;
	else                     V->data[p] = INT_MAX;
      }

      if(S != NULL){
	for(i = 1; i <= S[0]; i++){
	  if(label->data[S[i]] == lb)
	    PQueue32::FastInsertElem(Q, S[i]);
	}
      }
      else{
	for(p=0; p<n; p++)
	  if(label->data[p] == lb)
	    PQueue32::FastInsertElem(Q, p);
      }
      
      while(!PQueue32::IsEmpty(Q)){
	p = PQueue32::FastRemoveMinFIFO(Q);
	u.x = p%label->ncols; 
	u.y = p/label->ncols; 
	
	for(i=1; i<A->n; i++){
	  v.x = u.x + A->dx[i];
	  v.y = u.y + A->dy[i];
	  //if(Image32::IsValidPixel(label, v.x, v.y)){
	  if(v.x >= 0 && v.x < label->ncols &&
	     v.y >= 0 && v.y < label->nrows){
	    q = v.x + label->ncols*v.y;
	    if(Q->L.elem[q].color != BLACK){
	      edge = (sg->n_link[p])[i];          
	      tmp  = MAX(V->data[p], edge);
	      if(tmp < V->data[q]){
		if(Q->L.elem[q].color == GRAY)
		  PQueue32::FastRemoveElem(Q, q);
		V->data[q] = tmp; //mapa de custos
		PQueue32::FastInsertElem(Q, q);
	      }
	    }
	  }
	}
      }
      PQueue32::Destroy(&Q);
      return V;
    }
    
    
    sSet *COIFT_new_seeds(sImageGraph *sg,
			  int *S,
			  sImage32 *E,
			  sImage32 *label_oift){
      sPQueue32 *Q = NULL;
      sImage32 *V, *pred, *label;
      Pixel u,v;
      sAdjRel *A;
      sSet *newSi = NULL;
      int i,j,p,q,n, tmp, Emax, Emin, e_q, e_p, w_pq; //Wmin;
      
      //Wmin = MinimumWeight(sg);
      Emax = Image32::GetMaxVal(E);
      Emin = E->data[S[1]];
      
      V = Image32::Create(sg->ncols, sg->nrows);
      label = Image32::Create(sg->ncols, sg->nrows);
      pred = Image32::Create(sg->ncols, sg->nrows);
      n = label->ncols*label->nrows;
      Q = PQueue32::Create((Emax+1)+2, n, V->data);
      A = sg->A;
      
      Image32::Set(pred,  NIL);
      Image32::Set(V, INT_MAX);
      Image32::Set(label, NIL);

      for(i = 1; i <= S[0]; i++){
	if(label_oift->data[S[i]] == 1)
	  label->data[S[i]] = 1;
      }

      V->data[S[1]] = 0;
      PQueue32::FastInsertElem(Q, S[1]);
      
      while(!PQueue32::IsEmpty(Q)){
	p = PQueue32::FastRemoveMinFIFO(Q);
	u.x = p%label->ncols;
	u.y = p/label->ncols;
	
	if(E->data[p] < Emin)
	  Emin = E->data[p];
	
	for(i = 1; i < A->n; i++){
          v.x = u.x + A->dx[i];
          v.y = u.y + A->dy[i];
	  
          //if(Image32::IsValidPixel(label, v.x, v.y)){
	  if(v.x >= 0 && v.x < label->ncols &&
	     v.y >= 0 && v.y < label->nrows){
	    q = v.x + label->ncols * v.y;
            
	    if(Q->L.elem[q].color != BLACK){
	      //------------------------------------
	      e_p = E->data[p];
	      e_q = E->data[q];
	      w_pq = (sg->n_link[p])[i];
	      if(label->data[q] == 1){ //Seed pixel in Si
		tmp = 0;
	      }
	      else if(e_q >= Emin && label_oift->data[q] > 0){
		tmp = 1;
	      }
	      /*
		else if(e_q == e_p){
		tmp = 1;
		}
	      */
	      /*
		else if(e_q > w_pq){
		tmp = 0;
		}
		else if(e_q == w_pq && w_pq > Wmin){
		tmp = 1;
		}
	      */
	      else{
		tmp = Emax - e_q + 2;
	      }
	      //------------------------------------
	      if(tmp < V->data[q]){
		if(Q->L.elem[q].color == GRAY) 
		  PQueue32::FastRemoveElem(Q, q);                     
		V->data[q] = tmp;
		pred->data[q] = p;
		PQueue32::FastInsertElem(Q, q);
	      }
	    }
          }
	}
      }
      
      // calculando os predecessores das sementes Si2, para obter as novas sementes internas
      for(i = 1; i <= S[0]; i++){
	if(label->data[S[i]] == 1){
	  p = S[i];
	  q = pred->data[p];
	  while(q != NIL && (label->data[q] != 1)){ //enquanto o predecessor nao for raiz(i.e.,Si1) e seu rótulo nao for marcado antes(i.e., NIL)
	    Set::Insert(&newSi, q);
	    //printf("newSi->elem: %i\n",newSi->elem);
	    //printf("p: %i\n",q);
	    label->data[q] = 1;
	    q = pred->data[q];
	  }
	}
      }
      
      Image32::Destroy(&V);
      Image32::Destroy(&label);
      Image32::Destroy(&pred);
      PQueue32::Destroy(&Q);
      return newSi;
    }
   

    void COIFT(sImageGraph *sg,
	       int *S,
	       sImage32 *label){
      sImage32 *V_bkg, *Vc_bkg, *ero_V_bkg, *label_oift;
      sSet *Snew=NULL;
      int i_worst,i,n,p,energy,e;
      sAdjRel *A = AdjRel::Circular(2.5);
      
      if(S == NULL) return;
      
      label_oift = Image32::Clone(label);
      OIFT(sg, S, label_oift);
      
      /*1. IFT_max com semente do fundo*/
      ImageGraph::Transpose(sg);

      // mapa de custos via IFT_MAX usando só sementes do fundo
      V_bkg = Cost_fmax(sg, S, 0, label);
      ImageGraph::Transpose(sg);
      
      ero_V_bkg = Image32::Erode(V_bkg, A);
      //Vc_bkg = Image32::Complement(ero_V_bkg);
      
      //Image32::Write(V_bkg, "./out/V_bkg.pgm");
      //Image32::Write(ero_V_bkg, "./out/ero_V_bkg.pgm");
      //Image32::Write(Vc_bkg, "./out/Vc_bkg.pgm");
      
      /*2. IFT com f(pi_s.<s,t>)= Vc_bkg(t) */
      /*
	Si1 = Si;
	Si2 = Si->next;
	Si->next = NULL;
      */
      
      energy = INT_MAX;
      i_worst = 1;
      for(i = 1; i <= S[0]; i++){
	e = ero_V_bkg->data[S[i]];
	if(label->data[S[i]] == 1 && e < energy){
	  energy = e;
	  i_worst = i;
	}
      }

      p = S[1];
      S[1] = S[i_worst];
      S[i_worst] = p;
      
      Snew = COIFT_new_seeds(sg, S, ero_V_bkg, label_oift);
      
      /*3. IFT com novas sementes*/
      Image32::DrawSet(label, Snew,   1);

      OIFT(sg, NULL, label);
      
      AdjRel::Destroy(&A);
      Image32::Destroy(&ero_V_bkg);
      Image32::Destroy(&V_bkg);
      Image32::Destroy(&label_oift);
      Set::Destroy(&Snew);
    }
    

    //----------------------------------------

    sImage32 *BB_Geodesic_Cost(sImage32 *pred,
			       sAdjRel *A){
      sQueue *Qfifo=NULL;
      sImage32 *cost;
      Pixel u,v;
      int n,p,q,i;
      int *Dpq;
      
      Dpq = gft::AllocIntArray(A->n);
      for(i = 1; i < A->n; i++){
	Dpq[i] = ROUND(10*sqrtf(A->dx[i]*A->dx[i] + A->dy[i]*A->dy[i]));
      }
      
      n = pred->ncols*pred->nrows;
      cost  = Image32::Create(pred->ncols, pred->nrows);
      Image32::Set(cost, INT_MAX);
      Qfifo = Queue::Create(n);
      
      for(p = 0; p < n; p++){
	if(pred->data[p] == NIL){
	  cost->data[p] = 0;
	  Queue::Push(Qfifo, p);
	}
      }
      
      while(!Queue::IsEmpty(Qfifo)){
	p = Queue::Pop(Qfifo);
	u.x = p%pred->ncols; 
	u.y = p/pred->ncols; 
	
	for(i = 1; i < A->n; i++){
	  v.x = u.x + A->dx[i];
	  v.y = u.y + A->dy[i];
	  //if(Image32::IsValidPixel(pred, v.x, v.y)){
	  if(v.x >= 0 && v.x < pred->ncols &&
	     v.y >= 0 && v.y < pred->nrows){
	    q = v.x + pred->ncols*v.y;
	    if(p == pred->data[q]){
	      cost->data[q] = cost->data[p] + Dpq[i];
	      Queue::Push(Qfifo, q);
	    }
	  }
	}
      }
      gft::FreeIntArray(&Dpq);
      Queue::Destroy(&Qfifo);
      return cost;
    }



    sScene32 *BB_Geodesic_Cost(sScene32 *pred,
			       sAdjRel3 *A){
      sQueue *Qfifo=NULL;
      sScene32 *cost;
      Voxel u,v;
      int n,p,q,i;
      int *Dpq;
      
      Dpq = gft::AllocIntArray(A->n);
      for(i = 1; i < A->n; i++){
	Dpq[i] = ROUND(10*sqrtf(A->d[i].axis.x*A->d[i].axis.x +
				A->d[i].axis.y*A->d[i].axis.y +
				A->d[i].axis.z*A->d[i].axis.z));
      }
      
      n = pred->n;
      cost  = Scene32::Create(pred);
      Scene32::Fill(cost, INT_MAX);
      Qfifo = Queue::Create(n);
      
      for(p = 0; p < n; p++){
	if(pred->data[p] == NIL){
	  cost->data[p] = 0;
	  Queue::Push(Qfifo, p);
	}
      }
      
      while(!Queue::IsEmpty(Qfifo)){
	p = Queue::Pop(Qfifo);
	u.c.x = gft::Scene32::GetAddressX(pred, p);
	u.c.y = gft::Scene32::GetAddressY(pred, p);
	u.c.z = gft::Scene32::GetAddressZ(pred, p);	
	
	for(i = 1; i < A->n; i++){
	  v.c.x = u.c.x + A->d[i].axis.x;
	  v.c.y = u.c.y + A->d[i].axis.y;
	  v.c.z = u.c.z + A->d[i].axis.z;
	  if(Scene32::IsValidVoxel(pred, v)){
	    q = Scene32::GetVoxelAddress(pred, v);
	    if(p == pred->data[q]){
	      cost->data[q] = cost->data[p] + Dpq[i];
	      Queue::Push(Qfifo, q);
	    }
	  }
	}
      }
      gft::FreeIntArray(&Dpq);
      Queue::Destroy(&Qfifo);
      return cost;
    }


    
    //Always assuming that raw_map is at least double
    //the size on each dimension
    sImage32 *BB_CropTemplate(sImage32 *cost_template,
			      int *S,
			      sImage32 *label,
			      int lb){
      sImage32 *cropped;
      int x_center = 0, y_center = 0, tx_center, ty_center;
      int start_x, start_y, x, y, Imin;
      int i, p,q, ns = 0;
      //Calcua centro de massa das sementes internas:
      for(i = 1; i <= S[0]; i++){
	p = S[i];
	if(label->data[p] == lb){
	  x_center += p%label->ncols;
	  y_center += p/label->ncols;
	  ns++;
	}
      }
      if(ns != 0){
	x_center /= ns;
	y_center /= ns;
      }
      //------------------------
      cropped = Image32::Create(label->ncols, label->nrows);
      tx_center = cost_template->ncols/2;
      ty_center = cost_template->nrows/2;
      if(cost_template->array[ty_center][tx_center] != 0){
	Imin = gft::Image32::GetMinVal(cost_template, &p);
	tx_center = p%cost_template->ncols;
	ty_center = p/cost_template->ncols;
      }
      start_x = tx_center - x_center;
      start_y = ty_center - y_center;
      p = 0;
      for (y = start_y; y < start_y + label->nrows; y++) {  
	for (x = start_x; x < start_x + label->ncols; x++) {      
	  if(x >= 0 && x < cost_template->ncols &&
	     y >= 0 && y < cost_template->nrows){
	    q = y * cost_template->ncols + x;
	    cropped->data[p] = cost_template->data[q];
	    p++;
	  }
	  else
	    cropped->data[p] = INT_MAX;
	}
      }
      return cropped;
    }

    /*
    sImage32 *BB_OIFT_View(sImageGraph *G,
                           sImage32 *L,
			   sImage32 *C,
			   float delta,
			   int band_id) {
      sAdjRel *A = G->A;
      int p, q, i, n, Cmin;
      Pixel u, v;
      sImage32 *B;

      n = G->ncols * G->nrows;
      B = Image32::Create(G->ncols, G->nrows);
      
      Cmin = INT_MAX;
      for (p = 0; p < n; p++) {
	if (L->data[p] == 0) continue;
	u.x = p % G->ncols;  
	u.y = p / G->ncols;
	for (i = 1; i < A->n; i++) {
	  v.x = u.x + A->dx[i];
	  v.y = u.y + A->dy[i];
	  if (Image32::IsValidPixel(L, v.x, v.y)) {
	    q = v.x + G->ncols * v.y;
	    if (L->data[q] == band_id) 
	      if (C->data[p] < Cmin) Cmin = C->data[p];
	  }
	}
      }
      
      for (p = 0; p < n; p++) {
	if (C->data[p] <= Cmin + delta && C->data[p] >= Cmin)
	  B->data[p] = 1;
      }
      return B;
    }
    */
    
    void BB_OIFT_Propagate(int s,
			   sImageGraph *G,
			   sPQueue32 *Q,
			   sImage32 *V,
			   sImage32 *L,
			   int *i_inv) {
      Pixel u, v;
      sAdjRel *A = G->A;
      int j, i, t, tmp;
      
      u.x = s % L->ncols;  
      u.y = s / L->ncols;  
      
      for (i = 1; i < A->n; i++) {
	v.x = u.x + A->dx[i];
	v.y = u.y + A->dy[i];
	//if (Image32::IsValidPixel(L, v.x, v.y)) {
	if(v.x >= 0 && v.x < L->ncols &&
	   v.y >= 0 && v.y < L->nrows){
	  t = v.x + L->ncols * v.y;
	  if (Q->L.elem[t].color != BLACK) {

	    if(L->data[s]==0){
	      j = i_inv[i];
	      tmp = (G->n_link[t])[j];
	    }
	    else
	      tmp = (G->n_link[s])[i];
	    
	    if (tmp < V->data[t]) {
	      if (Q->L.elem[t].color == GRAY)
		PQueue32::FastRemoveElem(Q, t);
	      V->data[t] = tmp;
	      L->data[t] = L->data[s];
	      PQueue32::FastInsertElem(Q, t);
	    }
	  }
	}
      }
    }

    
    void BB_OIFT_Propagate_bkg(int s,
			       sImageGraph *G,
			       sPQueue32 *Q,
			       sPQueue32 *Qe,
			       sImage32 *V,
			       sImage32 *L,
			       int *i_inv) {
      sAdjRel *A = G->A;
      Pixel u, v;
      int i, t;
      
      BB_OIFT_Propagate(s, G, Q, V, L, i_inv);

      u.x = s % L->ncols;  
      u.y = s / L->ncols;  
      for (i = 1; i < A->n; i++) {
	v.x = u.x + A->dx[i];
	v.y = u.y + A->dy[i];
	//if (Image32::IsValidPixel(L, v.x, v.y)) {
	if(v.x >= 0 && v.x < L->ncols &&
	   v.y >= 0 && v.y < L->nrows){
	  t = v.x + L->ncols * v.y;
	  if (Qe->L.elem[t].color == WHITE) {
	    if ((L->data[t] != 0) || (Q->L.elem[t].color != BLACK))
	      PQueue32::FastInsertElem(Qe, t);
	  }
	}
      }
    }

    
    void BB_OIFT_Propagate_obj(int s,
			       sImageGraph *G,
			       sPQueue32 *Q,
			       sPQueue32 *Qi,
			       sImage32 *V,
			       sImage32 *L,
			       int *i_inv) {
      sAdjRel *A = G->A;
      Pixel u, v;
      int i, t;
      int white, black;
      white = black = 0;
      
      u.x = s % L->ncols;
      u.y = s / L->ncols;
      for (i = 1; i < A->n; i++) {
	v.x = u.x + A->dx[i];
	v.y = u.y + A->dy[i];
	//if (Image32::IsValidPixel(L, v.x, v.y)) {
	if(v.x >= 0 && v.x < L->ncols &&
	   v.y >= 0 && v.y < L->nrows){
	  t = v.x + L->ncols * v.y;
	  if (Qi->L.elem[t].color == WHITE) {
	    if (Q->L.elem[t].color == GRAY) { // mais um daqueles erros que não acontece nunca...
	      L->data[t] = 1;
	      PQueue32::FastRemoveElem(Q, t);
	      PQueue32::FastInsertElem(Qi, t);
	      BB_OIFT_Propagate(t, G, Q, V, L, i_inv);
	    }
	    else {
	      //printf("Erro 2 - %d not GRAY on queue Q - Color = %d\n", t, Q->L.elem[t].color);
	      if (Q->L.elem[t].color == BLACK) black++;
	      else if (Q->L.elem[t].color == WHITE) white++;
	      else printf("Error!\n");
	    }
	  }
	}
      }
    }

    
    sImage32 *BB_OIFT_GetLeafNodes(sImage32 *pred,
				   sAdjRel *A) {
      int n, p, q, i, isleaf;
      Pixel u, v;
      sImage32 *bin;
      
      bin = Image32::Create(pred->ncols, pred->nrows);
      n = pred->ncols * pred->nrows;
      for (p = 0; p < n; p++) {
	isleaf = 1;
	u.x = p % pred->ncols;  
	u.y = p / pred->ncols;  
	for (i = 1; i < A->n; i++) {
	  v.x = u.x + A->dx[i];
	  v.y = u.y + A->dy[i];
	  //if (Image32::IsValidPixel(pred, v.x, v.y)) {
	  if(v.x >= 0 && v.x < pred->ncols &&
	     v.y >= 0 && v.y < pred->nrows){
	    q = v.x + pred->ncols*v.y;
	    if(pred->data[q] == p) {
	      isleaf = 0;
	      break;
	    }
	  }
	}
	if (isleaf) bin->data[p] = 1;
	else        bin->data[p] = 0;
      }
      return bin;
    }


    
    void BB_OIFT(sImageGraph *G,
		 int *S,
		 sImage32 *L,
		 sImage32 *C,
		 sImage32 *P,
		 float delta) {
      sPQueue32 *Q  = NULL; // Fila "principal"
      sPQueue32 *Qi = NULL; // Candidatos a borda do objeto (interna - OBJ)
      sPQueue32 *Qe = NULL; // Candidatos a borda do fundo  (externa - BKG)
      sPQueue32 *Qt = NULL; // Leaves
      sImage32 *V, *Leaf, *B;
      int i, j, s, t, n, Cmin, Cmax, p, q, change_Qt;
      sSet *seed = NULL;
      sAdjRel *A = G->A;
      float tmp;
      Pixel u, v;
      int *i_inv = NULL;

      i_inv = gft::AdjRel::InverseIndexes(A);
      
      n = G->ncols * G->nrows;
      V = Image32::Create(G->ncols, G->nrows);
      
      Cmax = Image32::GetMaxVal(C);
     
      Qt = PQueue32::Create(Cmax + 2, n, C->data);
      if(P != NULL){
	Leaf = BB_OIFT_GetLeafNodes(P, G->A);
	for (t = 0; t < n; t++)
	  if (Leaf->data[t])
	    PQueue32::FastInsertElem(Qt, t);
	Image32::Destroy(&Leaf);
      }
      
      Q  = PQueue32::Create(G->Wmax + 2, n, V->data);
      Qi = PQueue32::Create(Cmax + 2, n, C->data);
      Qe = PQueue32::Create(Cmax + 2, n, C->data);
      
      for (t = 0; t < n; t++)
	V->data[t] = INT_MAX;

      Cmin = INT_MAX;
      for(i = 1; i <= S[0]; i++){
	t = S[i];
	if(L->data[t] == 1){
	  V->data[t] = 0;
	  PQueue32::FastInsertElem(Q, t);
	}
	else if(L->data[t] == 0){
	  V->data[t] = 0;
	  PQueue32::FastInsertElem(Q, t);
	  PQueue32::FastInsertElem(Qe, t);
	  if (C->data[t] < Cmin)
	    Cmin = C->data[t];	  
	}
      }
      
      for (t = 0; t < n; t++) {
	if ((C->data[t] - Cmin >= delta) && (V->data[t] != 0)) {
	  V->data[t] = 0;
	  PQueue32::FastInsertElem(Q, t);
	  PQueue32::FastInsertElem(Qe, t);
	}
      }
      
      while (!PQueue32::IsEmpty(Q)) {
	s = PQueue32::FastRemoveMinFIFO(Q);
	if (L->data[s] == 0) { 
	  if (Qe->L.elem[s].color == GRAY)
	    PQueue32::FastRemoveElem(Qe, s);
	  BB_OIFT_Propagate_bkg(s, G, Q, Qe, V, L, i_inv);
	  
	  change_Qt = true;
	  while(change_Qt == true) {
	    change_Qt = false;
	    while (!PQueue32::IsEmpty(Qe) &&
		   (PQueue32::FastGetMaxVal(Qe) - PQueue32::FastGetMinVal(Qe) >=
		    delta)) {
	      t = PQueue32::FastRemoveMaxFIFO(Qe);
	      if (Q->L.elem[t].color == GRAY) { 
		L->data[t] = 0;
		PQueue32::FastRemoveElem(Q, t);
		BB_OIFT_Propagate_bkg(t, G, Q, Qe, V, L, i_inv);
	      }
	      else {
		printf("Erro 1: C(t): %d\n", C->data[t]);
		/*
		printf("Bmin: %d, Bmax: %d\n", Qe->C.minvalue, Qe->C.maxvalue);
		printf("Omin: %d, Omax: %d\n", Qi->C.minvalue, Qi->C.maxvalue);
		printf("Ocolor: %d\n", Qi->L.elem[t].color);
		printf("q: %d, (x,y)=(%d,%d), Qcolor: %d, ", t, t % L->ncols, t/L->ncols, Q->L.elem[t].color);
		printf("L: %d\n", L->data[t]);
		*/
	      }
	    }

	    while (!PQueue32::IsEmpty(Qt) &&
		   (PQueue32::FastGetMaxVal(Qt) - PQueue32::FastGetMinVal(Qe) >=
		    delta)) {
	      t = PQueue32::FastRemoveMaxFIFO(Qt);
	      change_Qt = true;
	      if (Q->L.elem[t].color != BLACK) {
		L->data[t] = 0;
		if(Q->L.elem[t].color == GRAY)  PQueue32::FastRemoveElem(Q, t);
		if(Qe->L.elem[t].color == GRAY) PQueue32::FastRemoveElem(Qe, t);
		BB_OIFT_Propagate_bkg(t, G, Q, Qe, V, L, i_inv);
	      }
	      else if (L->data[t] == 1) printf("Erro 3:\n");
	    }
	  }
	}
	else if (L->data[s] == 1) {
	  PQueue32::FastInsertElem(Qi, s);
	  BB_OIFT_Propagate(s, G, Q, V, L, i_inv);
	  while (PQueue32::FastGetMaxVal(Qi) - PQueue32::FastGetMinVal(Qi) >=
		 delta) {
	    t = PQueue32::FastRemoveMinFIFO(Qi);
	    BB_OIFT_Propagate_obj(t, G, Q, Qi, V, L, i_inv);
	  }
	}
      }
      Image32::Destroy(&V);
      PQueue32::Destroy(&Qt);
      PQueue32::Destroy(&Qi);
      PQueue32::Destroy(&Qe);
      PQueue32::Destroy(&Q);
      free(i_inv);
    }

    
    
    void RBB_OIFT(sImageGraph *G,
		  int *S,
		  sImage32 *L,
		  sImage32 *C,
		  sImage32 *P,
		  float delta) {
      sPQueue32 *Q  = NULL; // Fila "principal"
      sPQueue32 *Qi = NULL; // Candidatos a borda do objeto (interna - OBJ)
      sPQueue32 *Qe = NULL; // Candidatos a borda do fundo  (externa - BKG)
      sPQueue32 *Qt = NULL; // Leaves
      sImage32 *V, *Leaf, *B;
      int i, j, s, t, n, Cmin, Cmax, p, q, change_Qt;
      sSet *seed = NULL;
      sAdjRel *A = G->A;
      float tmp;
      Pixel u, v;
      int *i_inv = NULL;

      i_inv = gft::AdjRel::InverseIndexes(A);
      
      n = G->ncols * G->nrows;
      V = Image32::Create(G->ncols, G->nrows);
      
      Cmax = Image32::GetMaxVal(C);
     
      Qt = PQueue32::Create(Cmax + 2, n, C->data);
      if(P != NULL){
	Leaf = BB_OIFT_GetLeafNodes(P, G->A);
	for (t = 0; t < n; t++)
	  if (Leaf->data[t])
	    PQueue32::FastInsertElem(Qt, t);
	Image32::Destroy(&Leaf);
      }
      
      Q  = PQueue32::Create(G->Wmax + 2, n, V->data);
      Qi = PQueue32::Create(Cmax + 2, n, C->data);
      Qe = PQueue32::Create(Cmax + 2, n, C->data);
      
      for (t = 0; t < n; t++)
	V->data[t] = INT_MAX;

      Cmin = INT_MAX;
      for(i = 1; i <= S[0]; i++){
	t = S[i];
	if(L->data[t] == 1){
	  V->data[t] = 0;
	  PQueue32::FastInsertElem(Q, t);
	}
	else if(L->data[t] == 0){
	  V->data[t] = 0;
	  PQueue32::FastInsertElem(Q, t);
	  PQueue32::FastInsertElem(Qe, t);
	  if (C->data[t] < Cmin)
	    Cmin = C->data[t];	  
	}
      }
      
      for (t = 0; t < n; t++) {
	if ((C->data[t] - Cmin >= delta * Cmin) && (V->data[t] != 0)) {
	  V->data[t] = 0;
	  PQueue32::FastInsertElem(Q, t);
	  PQueue32::FastInsertElem(Qe, t);
	}
      }
      
      while (!PQueue32::IsEmpty(Q)) {
	s = PQueue32::FastRemoveMinFIFO(Q);
	if (L->data[s] == 0) { 
	  if (Qe->L.elem[s].color == GRAY)
	    PQueue32::FastRemoveElem(Qe, s);
	  BB_OIFT_Propagate_bkg(s, G, Q, Qe, V, L, i_inv);
	  
	  change_Qt = true;
	  while(change_Qt == true) {
	    change_Qt = false;
	    while (!PQueue32::IsEmpty(Qe) &&
		   (PQueue32::FastGetMaxVal(Qe) - PQueue32::FastGetMinVal(Qe) >=
		    delta * PQueue32::FastGetMinVal(Qe))) {
	      t = PQueue32::FastRemoveMaxFIFO(Qe);
	      if (Q->L.elem[t].color == GRAY) { 
		L->data[t] = 0;
		PQueue32::FastRemoveElem(Q, t);
		BB_OIFT_Propagate_bkg(t, G, Q, Qe, V, L, i_inv);
	      }
	      else {
		printf("Erro 1: C(t): %d\n", C->data[t]);
		/*
		printf("Bmin: %d, Bmax: %d\n", Qe->C.minvalue, Qe->C.maxvalue);
		printf("Omin: %d, Omax: %d\n", Qi->C.minvalue, Qi->C.maxvalue);
		printf("Ocolor: %d\n", Qi->L.elem[t].color);
		printf("q: %d, (x,y)=(%d,%d), Qcolor: %d, ", t, t % L->ncols, t/L->ncols, Q->L.elem[t].color);
		printf("L: %d\n", L->data[t]);
		*/
	      }
	    }
	    
	    while (!PQueue32::IsEmpty(Qt) &&
		   (PQueue32::FastGetMaxVal(Qt) - PQueue32::FastGetMinVal(Qe) >=
		    delta * PQueue32::FastGetMinVal(Qe))) {
	      t = PQueue32::FastRemoveMaxFIFO(Qt);
	      change_Qt = true;
	      if (Q->L.elem[t].color != BLACK) {
		L->data[t] = 0;
		if(Q->L.elem[t].color == GRAY)  PQueue32::FastRemoveElem(Q, t);
		if(Qe->L.elem[t].color == GRAY) PQueue32::FastRemoveElem(Qe, t);
		BB_OIFT_Propagate_bkg(t, G, Q, Qe, V, L, i_inv);
	      }
	      else if (L->data[t] == 1) printf("Erro 3:\n");
	    }
	  }
	}
	else if (L->data[s] == 1) {
	  PQueue32::FastInsertElem(Qi, s);
	  BB_OIFT_Propagate(s, G, Q, V, L, i_inv);
	  while (PQueue32::FastGetMaxVal(Qi) - PQueue32::FastGetMinVal(Qi) >=
		 delta * PQueue32::FastGetMinVal(Qi)) {
	    t = PQueue32::FastRemoveMinFIFO(Qi);
	    BB_OIFT_Propagate_obj(t, G, Q, Qi, V, L, i_inv);
	  }
	}
      }
      Image32::Destroy(&V);
      PQueue32::Destroy(&Qt);
      PQueue32::Destroy(&Qi);
      PQueue32::Destroy(&Qe);
      PQueue32::Destroy(&Q);
      free(i_inv);
    }

    //----------------------------------------

    void B_OIFT(sImageGraph *G,
		int *S,
		sImage32 *L,
		sImage32 *C,
		float delta) {
      sPQueue32 *Q  = NULL; // Main queue.
      sPQueue32 *Qo = NULL; // Object pixels.
      sPQueue32 *Qb = NULL; // Background pixels.
      sPQueue32 *Qx = NULL; // Undefined pixels.
      sImage32 *V;
      int i, j, s, t, n, Cmin, Cmax, p, q;
      sSet *seed = NULL;
      sAdjRel *A = G->A;
      float tmp;
      Pixel u, v;
      int *i_inv = NULL;

      i_inv = gft::AdjRel::InverseIndexes(A);
      
      n = G->ncols * G->nrows;
      V = Image32::Create(G->ncols, G->nrows);
      
      Cmax = Image32::GetMaxVal(C);
     
      Q  = PQueue32::Create(G->Wmax + 2, n, V->data);
      Qo = PQueue32::Create(Cmax + 2, n, C->data);
      Qb = PQueue32::Create(Cmax + 2, n, C->data);
      Qx = PQueue32::Create(Cmax + 2, n, C->data);
      
      for (t = 0; t < n; t++){
	V->data[t] = INT_MAX;
	PQueue32::FastInsertElem(Qx, t);
      }

      for(i = 1; i <= S[0]; i++){
	t = S[i];
	if(L->data[t] == 1){
	  V->data[t] = 0;
	  PQueue32::FastInsertElem(Q, t);
	}
	else if(L->data[t] == 0){
	  V->data[t] = 0;
	  PQueue32::FastInsertElem(Q, t);
	}
      }
      
      while (!PQueue32::IsEmpty(Q)) {
	s = PQueue32::FastRemoveMinFIFO(Q);
	BB_OIFT_Propagate(s, G, Q, V, L, i_inv);
	PQueue32::FastRemoveElem(Qx, s);
	
	if(L->data[s] == 0){ 
	  PQueue32::FastInsertElem(Qb, s);

	  while (!PQueue32::IsEmpty(Qx) &&
		 (PQueue32::FastGetMaxVal(Qx) - PQueue32::FastGetMinVal(Qb) >= delta)){
	    t = PQueue32::FastRemoveMaxFIFO(Qx);
	    PQueue32::FastInsertElem(Qb, t);
	    L->data[t] = 0;
	    if (Q->L.elem[t].color == GRAY)
	      PQueue32::FastRemoveElem(Q, t);
	    else
	      Q->L.elem[t].color = BLACK;
	    BB_OIFT_Propagate(t, G, Q, V, L, i_inv);
	  }
	}
	else if(L->data[s] == 1){
	  PQueue32::FastInsertElem(Qo, s);
	  
	  while (!PQueue32::IsEmpty(Qx) &&
		 (PQueue32::FastGetMaxVal(Qo) - PQueue32::FastGetMinVal(Qx) >= delta)){
	    t = PQueue32::FastRemoveMinFIFO(Qx);
	    PQueue32::FastInsertElem(Qo, t);
	    L->data[t] = 1;
	    if (Q->L.elem[t].color == GRAY)
	      PQueue32::FastRemoveElem(Q, t);
	    else
	      Q->L.elem[t].color = BLACK;
	    BB_OIFT_Propagate(t, G, Q, V, L, i_inv);
	  }
	}
      }
      Image32::Destroy(&V);
      PQueue32::Destroy(&Qo);
      PQueue32::Destroy(&Qb);
      PQueue32::Destroy(&Qx);
      PQueue32::Destroy(&Q);
      free(i_inv);
    }



    void B_OIFT_Propagate(int s,
			  sGraph *graph,
			  sGraph *transpose,
			  sPQueue32 *Q,
			  int *V,
			  int *L){
      int i, t, tmp;
      sGraph *g;

      if(L[s]==0) g = transpose;
      else   	  g = graph;
      
      for(i = 0; i < g->nodes[s].outdegree; i++){
	t = g->nodes[s].adjList[i];

	if (Q->L.elem[t].color != BLACK) {
	  tmp = g->nodes[s].Warcs[i];
  
	  if (tmp < V[t]) {
	    if (Q->L.elem[t].color == GRAY)
	      PQueue32::FastRemoveElem(Q, t);
	    V[t] = tmp;
	    L[t] = L[s];
	    PQueue32::FastInsertElem(Q, t);
	  }
	}
      }
    }
    


    void B_OIFT(sGraph *graph,
		sGraph *transpose,
		int *S,
		int *L,
		int *C,
		float delta) {
      sPQueue32 *Q  = NULL; // Main queue.
      sPQueue32 *Qo = NULL; // Object pixels.
      sPQueue32 *Qb = NULL; // Background pixels.
      sPQueue32 *Qx = NULL; // Undefined pixels.
      int *V;
      int i, j, s, t, n, Cmin, Cmax, p, q, Wmax;
      sSet *seed = NULL;
      float tmp;
      
      n = graph->nnodes;
      V = gft::AllocIntArray(n);

      Cmax = INT_MIN;
      for(p=0; p<n; p++)
	if(C[p] > Cmax)
	  Cmax = C[p];

      Wmax = gft::Graph::GetMaximumArc(graph);
      
      Q  = PQueue32::Create(Wmax + 2, n, V);
      Qo = PQueue32::Create(Cmax + 2, n, C);
      Qb = PQueue32::Create(Cmax + 2, n, C);
      Qx = PQueue32::Create(Cmax + 2, n, C);
      
      for (t = 0; t < n; t++){
	V[t] = INT_MAX;
	PQueue32::FastInsertElem(Qx, t);
      }

      for(i = 1; i <= S[0]; i++){
	t = S[i];
	if(L[t] == 1){
	  V[t] = 0;
	  PQueue32::FastInsertElem(Q, t);
	}
	else if(L[t] == 0){
	  V[t] = 0;
	  PQueue32::FastInsertElem(Q, t);
	}
      }
      
      while (!PQueue32::IsEmpty(Q)) {
	s = PQueue32::FastRemoveMinFIFO(Q);
	B_OIFT_Propagate(s, graph, transpose, Q, V, L);
	PQueue32::FastRemoveElem(Qx, s);
	
	if(L[s] == 0){ 
	  PQueue32::FastInsertElem(Qb, s);

	  while (!PQueue32::IsEmpty(Qx) &&
		 (PQueue32::FastGetMaxVal(Qx) - PQueue32::FastGetMinVal(Qb) >= delta)){
	    t = PQueue32::FastRemoveMaxFIFO(Qx);
	    PQueue32::FastInsertElem(Qb, t);
	    L[t] = 0;
	    if (Q->L.elem[t].color == GRAY)
	      PQueue32::FastRemoveElem(Q, t);
	    else
	      Q->L.elem[t].color = BLACK;
	    B_OIFT_Propagate(t, graph, transpose, Q, V, L);
	  }
	}
	else if(L[s] == 1){
	  PQueue32::FastInsertElem(Qo, s);
	  
	  while (!PQueue32::IsEmpty(Qx) &&
		 (PQueue32::FastGetMaxVal(Qo) - PQueue32::FastGetMinVal(Qx) >= delta)){
	    t = PQueue32::FastRemoveMinFIFO(Qx);
	    PQueue32::FastInsertElem(Qo, t);
	    L[t] = 1;
	    if (Q->L.elem[t].color == GRAY)
	      PQueue32::FastRemoveElem(Q, t);
	    else
	      Q->L.elem[t].color = BLACK;
	    B_OIFT_Propagate(t, graph, transpose, Q, V, L);
	  }
	}
      }
      gft::FreeIntArray(&V);
      PQueue32::Destroy(&Qo);
      PQueue32::Destroy(&Qb);
      PQueue32::Destroy(&Qx);
      PQueue32::Destroy(&Q);
    }

    
    
    //----------------------------------------

    int *GetSeedsByLabel(int *S,
			 sImage32 *label,
			 int lb){
      int *Slb = NULL;
      int l, s, i, j, size = 0;
      for(i = 1; i <= S[0]; i++){
	s = S[i];
	l = label->data[s];
	if(l == lb) size++;
      }
      Slb = gft::AllocIntArray(size+1);
      Slb[0] = size;
      i = 1;
      for(j = 1; j <= S[0]; j++){
	s = S[j];
	l = label->data[s];
	if(l == lb){
	  Slb[i] = s;
	  i++;
	}
      }
      return Slb;
    }
    
    
    int *GetAllInternalSeedsByLabel(int *S,
				    sImage32 *label,
				    int lb,
				    int *hr,
				    int nlayers){
      int *Sall, *Sl, *Stmp;
      int l,size,i,j;
      Sall = GetSeedsByLabel(S, label, lb);
      for(l=0; l < nlayers; l++){
	if(hr[l] == lb-1){
	  //Sl = GetSeedsByLabel(S, label, l+1);
	  Sl = GetAllInternalSeedsByLabel(S, label, l+1,
					  hr, nlayers);
	  size = Sl[0] + Sall[0];
	  Stmp = gft::AllocIntArray(size+1);
	  Stmp[0] = size;
	  j =  1;
	  for(i = 1; i <= Sall[0]; i++){
	    Stmp[j] = Sall[i];
	    j++;
	  }
	  for(i = 1; i <= Sl[0]; i++){
	    Stmp[j] = Sl[i];
	    j++;
	  }
	  gft::FreeIntArray(&Sall);
	  gft::FreeIntArray(&Sl);
	  Sall = Stmp;
	}
      }
      return Sall;
    }
    

    sImageGraph *GetPolarityGraph(sImageGraph *graph,
				  sCImage *cimg,
				  sImage32 *img,
				  int pol){
      sImageGraph *sg_oriented;
      sg_oriented = gft::ImageGraph::Clone(graph);
      if(cimg != NULL){
	gft::sImage32 *lumi;
	lumi = gft::CImage::Luminosity(cimg);
	gft::ImageGraph::Orient2Digraph(sg_oriented, lumi, pol);
	gft::Image32::Destroy(&lumi);
      }
      else
	gft::ImageGraph::Orient2Digraph(sg_oriented, img, pol);
      return sg_oriented;
    }
    

    
    void HL_OIFT(sImageGraph *graph,
		 sCImage *cimg,
		 sImage32 *img,
		 float radius,
		 char *hierarchy,
		 int *S,
		 sImage32 *label){
      int p;
      //-----------------------
      gft::sLayeredGraph *lg;
      gft::sImageGraph *sg, *sg_oriented;
      //------------------------
      FILE *fpHierarchy;
      int nlayers, typeL, nRelations, relationT, layer1, layer2;
      int i, j, n, nSeeds, k, s, Wmax, pol;
      int x, y, lb, ncols, nrows;
      int *lg_label;
      
      ncols = label->ncols;
      nrows = label->nrows;
      gft::ImageGraph::ChangeType(graph, DISSIMILARITY);
      
      fpHierarchy = fopen((char *)hierarchy, "r");
      if (!fpHierarchy){
	gft::Error((char *)MSG1,(char *)"Error opening file of Hierarchy Relations!");
      }
      
      /*------------- Set number of layers -------------*/
      (void)fscanf(fpHierarchy,"%d\n",&nlayers);
      //printf("NLayers =  %d\n",nlayers);
      
      /*-------------GIVEN BY USER: TYPE OF LAYER -------------*/
      /*typeLayer:
	0 = normal,
	1 = GSC,
      */
      int *typeLayer = gft::AllocIntArray(nlayers);
      int *polLayer = gft::AllocIntArray(nlayers);
      
      for(i = 0; i < nlayers; i++) {
	(void)fscanf(fpHierarchy,"%d %d\n",&typeL, &pol);
	typeLayer[i] = typeL;
	polLayer[i] = pol;
      }
      
      /*------------- Creating HIERARCHY -------------*/
      /*Relations is given by inclusion(=1) or exclusion(=2) binary relation between: (son,father) and (brothers) that defines a hierarchy, in a .txt*/
      /*Example: 1|0|1 = inclusion of layer 0 in layer 1 
	or 2|0|1 = exclusion between layer 0 and layer 1*/
      
      /*Initialization of Hierarchy*/
      int *hr = gft::AllocIntArray(nlayers);
      for(i = 0; i < nlayers; i++){
	hr[i]= -1; 
      }
      
      /*Reading a .txt fill, and set in "hr" with only inclusion cases*/
      // hr[0] = 1;  // i(0,1) = 1 0 1  // means: "1 is father of 0" or "0 is included in 1"
      (void)fscanf(fpHierarchy,"%d\n",&nRelations);
      for(i = 0; i < nRelations; i++) {
	(void)fscanf(fpHierarchy,"%d %d %d\n",&relationT,&layer1,&layer2);
	if(relationT == 1){ // it's inclusion relation
	  hr[layer1] = layer2;
	}
      }
      
      fclose(fpHierarchy);
      
      /* ------------- Create label image -------------*/
      n = ncols*nrows*nlayers;
      lg_label = gft::AllocIntArray(n);
      /*Set NIL (-1) in all matrix positions*/
      for(i = 0; i < n; i++){
	lg_label[i] = NIL;
      }
      
      nSeeds = S[0];
      int* Seeds = gft::AllocIntArray(nSeeds*nlayers + 1);
      
      k = 1;
      for(i = 1; i <= S[0]; i++){
	s = S[i];
	x = s%ncols;
	y = s/ncols;
	lb = label->data[s];
	
	if(lb != 0){
	  p = x + y*ncols + nrows*ncols*(lb-1); 
	  lg_label[p] = 1;
	  Seeds[k] = p;
	  k++;
	  /*Essa parte faltou no algoritmo do artigo:*/
	  for(j = 0; j < nlayers ; j++){
	    /*
	    if(j != lb-1 && hr[lb-1] == j){
	      p = x + y*ncols + nrows*ncols*(j);
	      lg_label[p] = 1;
	      Seeds[k] = p;
	      k++;
	    }
	    else if(j != lb-1){
	      p = x + y*ncols + nrows*ncols*(j);
	      lg_label[p] = 0;
	      Seeds[k] = p;
	      k++;
	    }
	    */
	    if(j != lb-1 && hr[j] == lb-1){
	      p = x + y*ncols + nrows*ncols*(j);
	      lg_label[p] = 0;
	      Seeds[k] = p;
	      k++;
	    }
	  }
	}
	else{ // it is background seed, then transfer for all layers
	  for(j = 0; j < nlayers ; j++){
	    p = x + y*ncols + nrows*ncols*(j);
	    lg_label[p] = 0;
	    Seeds[k] = p;
	    k++;
	  }
	}
      }
      Seeds[0] = k-1;

      /* ------------- Create LAYERED GRAPH -------------*/
      lg = gft::LayeredGraph::Create(nlayers, ncols*nrows);
      
      /*Create a graph image with 8-neighborhood*/
      sg = graph;
      //gft::SparseGraph::SuppressZeroWeightedArcs(sg);
      Wmax = sg->Wmax;

      /*Set cost of arcs for each layer*/
      for(i = 0; i < nlayers; i++){
	if(typeLayer[i] == 0){ /* 0 = normal*/
	  if(polLayer[i] == 0)
	    gft::LayeredGraph::SetArcs(lg, sg, i);
	  else{
	    sg_oriented = GetPolarityGraph(sg, cimg, img, polLayer[i]);
	    if(sg_oriented->Wmax > Wmax) Wmax = sg_oriented->Wmax;
	    gft::LayeredGraph::SetArcs(lg, sg_oriented, i);
	    gft::ImageGraph::Destroy(&sg_oriented);
	  }
	}
	else if(typeLayer[i] == 1){ /* 1 = GSC */
	  gft::sImageGraph *sg_GSC;
	  gft::sImage32 *P_sum;
	  int *Slb = NULL;
	  if(polLayer[i] == 0)
	    sg_GSC = gft::ImageGraph::Clone(sg);
	  else{
	    sg_GSC = GetPolarityGraph(sg, cimg, img, polLayer[i]);
	    if(sg_GSC->Wmax > Wmax) Wmax = sg_GSC->Wmax;
	  }
	  Slb = GetAllInternalSeedsByLabel(S, label, i+1, hr, nlayers);
	  //gft::SparseGraph::SuppressZeroWeightedArcs(sg_GSC);
	  P_sum = gft::ift::SC_Pred_fsum(sg_GSC, Slb, 0.1);
	  gft::ImageGraph::Orient2DigraphOuter(sg_GSC, P_sum);
	  gft::LayeredGraph::SetArcs(lg, sg_GSC, i);  
	  gft::Image32::Destroy(&P_sum);
	  gft::ImageGraph::Destroy(&sg_GSC);
	  gft::FreeIntArray(&Slb);
	}
      }

      
      for(i = 0; i < nlayers ; i++){
	if(hr[i] != -1){ //son and father
	  gft::LayeredGraph::SetArcs(lg, i, hr[i], ncols, 0.0, radius);
	  gft::LayeredGraph::SetArcs(lg, hr[i], i, ncols, (float)(Wmax+1), radius);
	}
      }
      

      //Second exclusion
      for(i = 0; i < nlayers-1; i++){
	for(j = i+1; j < nlayers; j++){
	  if( hr[i] == hr[j] ){ //same father -> brothers -> exclusion
	    gft::LayeredGraph::SetArcs(lg, i, j, ncols, (float)(Wmax+1), radius);
	    gft::LayeredGraph::SetArcs(lg, j, i, ncols, 0.0, radius);
	  }
	}
      }
      

      /* ------------- Executa a Hierarchical Layered OIFT -------------*/ 
      gft::ift::HL_OIFT(lg, Wmax, Seeds, lg_label, hr);

      
      /* ------------- OUTPUT RESULT-------------*/
      
      gft::sQueue *FIFO;
      FIFO = gft::Queue::Create(nlayers);
      int *depth = gft::AllocIntArray(nlayers);
      int dmax = 0;
      for(i = 0; i < nlayers ; i++){
	if(hr[i] == -1){
	  depth[i] = 0;
	  gft::Queue::Push(FIFO, i);
	}
	else
	  depth[i] = -1;
      }
      while(!gft::Queue::IsEmpty(FIFO)){
	j = gft::Queue::Pop(FIFO);
	for(i = 0; i < nlayers ; i++){
	  if(hr[i] == j){
	    depth[i] = depth[j] + 1;
	    if(depth[i] > dmax) dmax = depth[i];
	    gft::Queue::Push(FIFO, i);
	  }
	}
      }
      
      int sizeImg = ncols*nrows;
      
      gft::Image32::Set(label, 0);
      int d;
      for(d = 0; d <= dmax; d++){
	for(i = 0; i < nlayers; i++){
	  if(depth[i] == d){
	    for(j = 0; j < sizeImg; j++){
	      if(lg_label[j+ sizeImg*i] != 0)
		label->data[j] = i+1;
	    }
	  }
	}
      }

      gft::LayeredGraph::Destroy(&lg);
      
      gft::FreeIntArray(&typeLayer);
      gft::FreeIntArray(&polLayer);
      gft::FreeIntArray(&hr);
      gft::FreeIntArray(&lg_label);
      gft::FreeIntArray(&Seeds);
      gft::FreeIntArray(&depth);
      gft::Queue::Destroy(&FIFO);
    }


    //----------------------------------------

    
    void HL_OIFT(sLayeredGraph *lg,
		 int Wmax, int *S, int *L, int *hr){
      sPQueue32 *Q=NULL;
      int i,j,p,q,n,layer_p,layer_q,size_layer,tmp2=0;
      int w_pq=0, w_qp=0, tmp=0;
      int* value;
      sGraphNode *A;
      int exclusionCase = 0;
      int *i_inv = NULL;
      int i_inv_size;
      bool flag;

      A = lg->graph->nodes;
      p = 0;
      i_inv_size = A[p].outdegree;
      i_inv = gft::AllocIntArray(i_inv_size);
      flag = true;
      while(flag){
	flag = false;
	for(i = 0; i < A[p].outdegree && i < i_inv_size; i++){
	  q = A[p].adjList[i];
	  if(p == q){ flag = true; break; }
	  for(j = 0; j < A[q].outdegree; j++)
	    if(A[q].adjList[j] == p)
	      i_inv[i] = j;
	}
	p++;
      }
      
      size_layer = lg->nnodesperlayer;
      n = size_layer*lg->nlayers;
      
      //printf("Size Layer: %d\n",size_layer);
      //printf("N pixels in the Graph: %d\n",n);
      
      /*Initialization*/
      value = gft::AllocIntArray(n);
      //Q = gft::PQueue32::Create(Wmax*lg->nlayers+2,n,value);
      Q = gft::PQueue32::Create(Wmax+3,n,value);
      
      /*Insert in value*/
      for(p=0; p<n; p++){
	if(L[p]==NIL) value[p] = INT_MAX; //(int)floor(FLT_MAX+0.5);
	else          value[p] = 0;
      }
      
      /*Insert in PQueue*/
      if(S != NULL){
	for(i=1; i<=S[0]; i++){
	  gft::PQueue32::FastInsertElem(Q, S[i]);
	}
      }
      else{
	for(p=0; p<n; p++)
	  if(L[p]!=NIL)
	    gft::PQueue32::FastInsertElem(Q, p);	    
      }
      
      /*
	clock_t t;
	t = clock();
      */
      
      /*Starting OIFT*/
      while(!gft::PQueue32::IsEmpty(Q)){
	p = gft::PQueue32::FastRemoveMinFIFO(Q);
	
	for(i=0; i<A[p].outdegree; i++){
	  q = A[p].adjList[i];
	  //if(q != NIL){
	  if(Q->L.elem[q].color != BLACK){
	    //w_pq = ROUND(gft::Graph::GetArcWeight(lg->graph,p,q));
	    w_pq = ROUND(A[p].Warcs[i]);
	    //if(w_pq < 0) continue;
	    
	    layer_p = p/size_layer;
	    layer_q = q/size_layer;
	    exclusionCase = 0;
	    /* Analize each relation type*/
	    if (layer_p == layer_q){  /*SAME LAYER*/
	      // Get w_qp
	      //w_qp = ROUND(gft::Graph::GetArcWeight(lg->graph,q,p));
	      w_qp = ROUND(A[q].Warcs[i_inv[i]]);
	      exclusionCase = 0;
	      if(L[p] != 0)
		tmp = w_pq;
	      else 
		tmp = w_qp; 
	    }
	    /*INTER LAYERS: relations where defined*/
	    else if(layer_p != layer_q){  
	      // Get w_qp
	      if(w_pq == 0)
		w_qp = Wmax+1; //INT_MAX;
	      else
		w_qp = 0;
	      // inclusion
	      /*layer_p is "INCLUDED" in layer_q /OR/ layer_q is "INCLUDED" in layer_p */
	      if((hr[layer_p] == layer_q) || (hr[layer_q] == layer_p)){  
		exclusionCase = 0;
		if(L[p] != 0)
		  tmp = w_pq;
		else 
		  tmp = w_qp; 
	      }
	      // exclusion
	      /*layer_q is EXCLUDED from layer_p */
	      else if(hr[layer_p] == hr[layer_q]){  
		exclusionCase = 1;
		/* It was defined to change "inputs"/"outputs", in/from the higher layer for exclusion case*/
		if (layer_p < layer_q){ 
		  if(L[p] != 0)
		    tmp = w_qp;
		  else 
		    tmp = w_pq;
		}
		else{
		  if(L[p] != 0)
		    tmp = w_pq;
		  else 
		    tmp = w_qp;
		}
	      }
	    }
	    
	    if(tmp < value[q]){
	      if(Q->L.elem[q].color == GRAY)
		gft::PQueue32::FastRemoveElem(Q, q);
	      
	      value[q] = tmp;
	      
	      if(exclusionCase == 1){ /*Treat exclusion case*/
		if(L[p] != 0)
		  L[q] = 0;
	      }
	      else{ /*exclusion == 0*/
		L[q] = L[p];
	      }
	      
	      gft::PQueue32::FastInsertElem(Q, q);
	    } 
	  }
	  //}
	}
      }

      /*
	t = clock() - t;
	double time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds
	printf("HLOIFT took %d seconds(t) to execute \n", ((int)t));
	printf("HLOIFT took %f seconds(t/C_P_S) to execute \n", time_taken);      
      */
      
      gft::PQueue32::Destroy(&Q);
      gft::FreeIntArray(&value);
      gft::FreeIntArray(&i_inv);
    }


    //--------------------------------------------------------
 
    void HL_OIFT_2(sLayeredGraph *lg,
		   int Wmax, int *S, int *L, int *hr){
      sPQueue32 *Q=NULL;
      int i,j,p,q,n,layer_p,layer_q,size_layer,tmp2=0;
      int w_pq=0, w_qp=0, tmp=0;
      int* value;
      sGraphNode *A;
      int exclusionCase = 0;

      A = lg->graph->nodes;
      size_layer = lg->nnodesperlayer;
      n = size_layer*lg->nlayers;

      /*Initialization*/
      value = gft::AllocIntArray(n);
      Q = gft::PQueue32::Create(Wmax+3,n,value);
      
      /*Insert in value*/
      for(p=0; p<n; p++){
	if(L[p]==NIL) value[p] = INT_MAX; //(int)floor(FLT_MAX+0.5);
	else          value[p] = 0;
      }
      
      /*Insert in PQueue*/
      if(S != NULL){
	for(i=1; i<=S[0]; i++){
	  gft::PQueue32::FastInsertElem(Q, S[i]);
	}
      }
      else{
	for(p=0; p<n; p++)
	  if(L[p]!=NIL)
	    gft::PQueue32::FastInsertElem(Q, p);	    
      }
      
      /*Starting OIFT*/
      while(!gft::PQueue32::IsEmpty(Q)){
	p = gft::PQueue32::FastRemoveMinFIFO(Q);
	
	for(i=0; i<A[p].outdegree; i++){
	  q = A[p].adjList[i];
	  //if(q != NIL){
	  if(Q->L.elem[q].color != BLACK){
	    //w_pq = ROUND(gft::Graph::GetArcWeight(lg->graph,p,q));
	    w_pq = ROUND(A[p].Warcs[i]);
	    //if(w_pq < 0) continue;
	    
	    layer_p = p/size_layer;
	    layer_q = q/size_layer;
	    exclusionCase = 0;
	    /* Analize each relation type*/
	    if (layer_p == layer_q){  /*SAME LAYER*/
	      // Get w_qp
	      w_qp = ROUND(gft::Graph::GetArcWeight(lg->graph,q,p));
	      exclusionCase = 0;
	      if(L[p] != 0)
		tmp = w_pq;
	      else 
		tmp = w_qp; 
	    }
	    /*INTER LAYERS: relations where defined*/
	    else if(layer_p != layer_q){  
	      // Get w_qp
	      if(w_pq == 0)
		w_qp = Wmax+1; //INT_MAX;
	      else
		w_qp = 0;
	      // inclusion
	      /*layer_p is "INCLUDED" in layer_q /OR/ layer_q is "INCLUDED" in layer_p */
	      if((hr[layer_p] == layer_q) || (hr[layer_q] == layer_p)){  
		exclusionCase = 0;
		if(L[p] != 0)
		  tmp = w_pq;
		else 
		  tmp = w_qp; 
	      }
	      // exclusion
	      /*layer_q is EXCLUDED from layer_p */
	      else if(hr[layer_p] == hr[layer_q]){  
		exclusionCase = 1;
		/* It was defined to change "inputs"/"outputs", in/from the higher layer for exclusion case*/
		if (layer_p < layer_q){ 
		  if(L[p] != 0)
		    tmp = w_qp;
		  else 
		    tmp = w_pq;
		}
		else{
		  if(L[p] != 0)
		    tmp = w_pq;
		  else 
		    tmp = w_qp;
		}
	      }
	    }
	    
	    if(tmp < value[q]){
	      if(Q->L.elem[q].color == GRAY)
		gft::PQueue32::FastRemoveElem(Q, q);
	      
	      value[q] = tmp;
	      
	      if(exclusionCase == 1){ /*Treat exclusion case*/
		if(L[p] != 0)
		  L[q] = 0;
	      }
	      else{ /*exclusion == 0*/
		L[q] = L[p];
	      }
	      
	      gft::PQueue32::FastInsertElem(Q, q);
	    } 
	  }
	}
      }

      gft::PQueue32::Destroy(&Q);
      gft::FreeIntArray(&value);
    }
   


    //Compute watershed by fpeak from markers
    void IFT_fpeak(sImage32 *grad,
		   sAdjRel *A,
		   sImage32 *label){
      sPQueue32 *Q;
      sImage32 *cost; //*Rmin;
      int tmp, w, Wmax;
      int n,p,q,i,px,py,qx,qy;
      n = grad->n;
      Wmax = gft::Image32::GetMaxVal(grad);
      //Rmin = gft::Image32::RegMin(grad, A);
      cost = gft::Image32::Create(grad);
      Q = gft::PQueue32::Create(Wmax+2, n, cost->data);

      for(p = 0; p < n; p++){
	//if(Rmin->data[p] == 0){
	if(label->data[p] == NIL){
	  cost->data[p] = INT_MAX;
	  //label->data[p] = NIL;
	}
	else{
	  cost->data[p] = grad->data[p];
	  //label->data[p] = Rmin->data[p]-1;
	  PQueue32::FastInsertElem(Q, p);
	}
      }
      
      while(!PQueue32::IsEmpty(Q)){
	p = PQueue32::FastRemoveMinFIFO(Q);
	px = p%grad->ncols;
	py = p/grad->ncols;
	for(i=1; i < A->n; i++){
	  qx = px + A->dx[i];
	  qy = py + A->dy[i];
	  //if(gft::Image32::IsValidPixel(grad, qx,qy)){
	  if(qx >= 0 && qx < grad->ncols &&
	     qy >= 0 && qy < grad->nrows){
	    q = qx + qy*grad->ncols;
	      
	    if(Q->L.elem[q].color != BLACK){
	      w = grad->data[q];
	      tmp = MAX(cost->data[p], w);
	      
	      if(tmp < cost->data[q]){
		if(Q->L.elem[q].color == GRAY)
		  PQueue32::FastRemoveElem(Q, q);
		cost->data[q] = tmp;
		label->data[q] = label->data[p];
		PQueue32::FastInsertElem(Q, q);
	      }
	    }
	  }
	}
      }
      gft::PQueue32::Destroy(&Q);
      gft::Image32::Destroy(&cost);
      //gft::Image32::Destroy(&Rmin);
    }

    
    //Compute watershed by fwv from markers
    void IFT_fwv(sImage32 *grad,
		 sAdjRel *A,
		 sImage32 *label){
      sPQueue32 *Q;
      sImage32 *cost; //*Rmin;
      int tmp, Wmax;
      int n,p,q,i,px,py,qx,qy;
      n = grad->n;
      Wmax = gft::Image32::GetMaxVal(grad);
      //Rmin = gft::Image32::RegMin(grad, A);
      cost = gft::Image32::Create(grad);
      Q = gft::PQueue32::Create(Wmax+2, n, cost->data);

      for(p = 0; p < n; p++){
	//if(Rmin->data[p] == 0){
	if(label->data[p] == NIL){
	  cost->data[p] = INT_MAX;
	  //label->data[p] = NIL;
	}
	else{
	  cost->data[p] = grad->data[p];
	  //label->data[p] = Rmin->data[p]-1;
	  PQueue32::FastInsertElem(Q, p);
	}
      }
      
      while(!PQueue32::IsEmpty(Q)){
	p = PQueue32::FastRemoveMinFIFO(Q);
	px = p%grad->ncols;
	py = p/grad->ncols;
	for(i=1; i < A->n; i++){
	  qx = px + A->dx[i];
	  qy = py + A->dy[i];
	  //if(gft::Image32::IsValidPixel(grad, qx,qy)){
	  if(qx >= 0 && qx < grad->ncols &&
	     qy >= 0 && qy < grad->nrows){
	      q = qx + qy*grad->ncols;
	      
	    if(Q->L.elem[q].color != BLACK){
	      tmp = grad->data[q];
	      
	      if(tmp < cost->data[q]){
		if(Q->L.elem[q].color == GRAY)
		  PQueue32::FastRemoveElem(Q, q);
		cost->data[q] = tmp;
		label->data[q] = label->data[p];
		PQueue32::FastInsertElem(Q, q);
	      }
	    }
	  }
	}
      }
      gft::PQueue32::Destroy(&Q);
      gft::Image32::Destroy(&cost);
      //gft::Image32::Destroy(&Rmin);
    }




    void RelaxMobj(sImageGraph *sg,
		   int *S,
		   sImage32 *label, int ntimes){
      sImage32 *mask;
      float *flabel_1[100],*flabel_2[100],*tmp;
      sAdjRel *A;
      int *mask_nodes;
      int n,p,q,i,j,k,px,py,qx,qy,nlast,ninic,Lmax,l,lb,lm;
      float sw,w;

      Lmax = gft::Image32::GetMaxVal(label);
      gft::ImageGraph::ChangeType(sg, DISSIMILARITY);
      
      A = sg->A;
      ninic = 1;
      n = label->n;
      for(l = 0; l <= Lmax; l++){
	flabel_1[l] = (float *)calloc(n, sizeof(float));
	flabel_2[l] = (float *)calloc(n, sizeof(float));
      }
      mask_nodes = (int *)malloc(sizeof(int)*(n+1));
      mask_nodes[0] = 0;
      mask = gft::Image32::Create(label);
      for(p = 0; p < n; p++){
	lb = label->data[p];
	flabel_1[lb][p] = 1.0;
	px = p%label->ncols;
	py = p/label->ncols;
	for(i=1; i < A->n; i++){
	  qx = px + A->dx[i];
	  qy = py + A->dy[i];
	  //if(gft::Image32::IsValidPixel(label, qx,qy)){
	  if(qx >= 0 && qx < label->ncols &&
	     qy >= 0 && qy < label->nrows){
	    q = qx + qy*label->ncols;
	    if(label->data[p] != label->data[q]){
	      mask->data[p] = 1;
	      mask_nodes[0]++;
	      mask_nodes[mask_nodes[0]] = p;
	      break;
	    }
	  }
	}
      }

      while(ntimes > 0){
	//Update flabel:
	for(l = 0; l <= Lmax; l++)
	  memcpy(flabel_2[l], flabel_1[l], n*sizeof(float));

	for(j = 1; j <= mask_nodes[0]; j++){
	  p = mask_nodes[j];
	  px = p%label->ncols;
	  py = p/label->ncols;
	  for(l = 0; l <= Lmax; l++)
	    flabel_2[l][p] = 0.0;
	  sw = 0.0;
	  for(i=1; i < A->n; i++){
	    qx = px + A->dx[i];
	    qy = py + A->dy[i];
	    //if(gft::Image32::IsValidPixel(label, qx,qy)){
	    if(qx >= 0 && qx < label->ncols &&
	       qy >= 0 && qy < label->nrows){
	      q = qx + qy*label->ncols;
	      w = (sg->n_link[p])[i];
	      w = sg->Wmax - w;
	      w = w*w;
	      w = w*w;
	      w = w*w;

	      sw += w;
	      for(l = 0; l <= Lmax; l++)
		flabel_2[l][p] += w*flabel_1[l][q];
	    }
	  }
	  for(l = 0; l <= Lmax; l++)
	    flabel_2[l][p] /= sw;
	}
	for(l = 0; l <= Lmax; l++){
	  tmp = flabel_1[l];
	  flabel_1[l] = flabel_2[l];
	  flabel_2[l] = tmp;
	}
	
	for(i = 1; i <= S[0]; i++){
	  p = S[i];
	  lb = label->data[p];
	  for(l = 0; l <= Lmax; l++){
	    if(lb == l)
	      flabel_1[l][p] = 1.0;
	    else
	      flabel_1[l][p] = 0.0;
	  }
	}
	
	ntimes--;
	
	if(ntimes > 0){
	  //Dilate mask:
	  nlast = mask_nodes[0];
	  for(j = ninic; j <= mask_nodes[0]; j++){
	    p = mask_nodes[j];
	    px = p%label->ncols;
	    py = p/label->ncols;
	    for(i=1; i < A->n; i++){
	      qx = px + A->dx[i];
	      qy = py + A->dy[i];
	      //if(gft::Image32::IsValidPixel(label, qx,qy)){
	      if(qx >= 0 && qx < label->ncols &&
		 qy >= 0 && qy < label->nrows){
		q = qx + qy*label->ncols;
		
		if(mask->data[q] == 0){
		  mask->data[q] = 1;
		  nlast++;
		  mask_nodes[nlast] = q;
		}
	      }
	    }
	  }
	  ninic = mask_nodes[0] + 1;
	  mask_nodes[0] = nlast;
	}
      }


      for(p = 0; p < n; p++){
	lm = 0;
	for(l = 1; l <= Lmax; l++){
	  if(flabel_1[l][p] > flabel_1[lm][p])
	    lm = l;
	}
	label->data[p] = lm;
      }

      for(l = 0; l <= Lmax; l++){
	free(flabel_1[l]);
	free(flabel_2[l]);
      }
      free(mask_nodes); 
      gft::Image32::Destroy(&mask);
    }
    


    float *Relax_dual(sImageGraph *sg,
		      int *S,
		      sImage32 *label, int ntimes){
      sImage32 *mask;
      float *flabel_1,*flabel_2,*tmp;
      sAdjRel *A;
      int *mask_nodes;
      int n,p,q,i,j,k,px,py,qx,qy,nlast,ninic;
      float sw,w;

      A = sg->A;
      ninic = 1;
      n = label->n;
      flabel_1 = (float *)malloc(sizeof(float)*n);
      flabel_2 = (float *)malloc(sizeof(float)*n);
      mask_nodes = (int *)malloc(sizeof(int)*(n+1));
      mask_nodes[0] = 0;
      mask = gft::Image32::Create(label);
      for(p = 0; p < n; p++){
	if(label->data[p] > 0)
	  flabel_1[p] = 1.0;
	else
	  flabel_1[p] = 0.0;
	px = p%label->ncols;
	py = p/label->ncols;
	for(i=1; i < A->n; i++){
	  qx = px + A->dx[i];
	  qy = py + A->dy[i];
	  //if(gft::Image32::IsValidPixel(label, qx,qy)){
	  if(qx >= 0 && qx < label->ncols &&
	     qy >= 0 && qy < label->nrows){
	    q = qx + qy*label->ncols;
	    if(label->data[p] != label->data[q]){
	      mask->data[p] = 1;
	      mask_nodes[0]++;
	      mask_nodes[mask_nodes[0]] = p;
	      break;
	    }
	  }
	}
      }

      while(ntimes > 0){
	//Update flabel:
	memcpy(flabel_2, flabel_1, n*sizeof(float));
	for(j = 1; j <= mask_nodes[0]; j++){
	  p = mask_nodes[j];
	  px = p%label->ncols;
	  py = p/label->ncols;
	  flabel_2[p] = 0.0;
	  sw = 0.0;
	  for(i=1; i < A->n; i++){
	    qx = px + A->dx[i];
	    qy = py + A->dy[i];
	    //if(gft::Image32::IsValidPixel(label, qx,qy)){
	    if(qx >= 0 && qx < label->ncols &&
	       qy >= 0 && qy < label->nrows){
	      q = qx + qy*label->ncols;
	      w = (sg->n_link[p])[i];

	      //w = w*w;
	      //w = w*w;
	      //w = w*w;

	      sw += w;
	      flabel_2[p] += w*flabel_1[q];
	    }
	  }
	  flabel_2[p] /= sw;
	}
	tmp = flabel_1;
	flabel_1 = flabel_2;
	flabel_2 = tmp;

	for(i = 1; i <= S[0]; i++){
	  p = S[i];
	  if(label->data[p] > 0)
	    flabel_1[p] = 1.0;
	  else
	    flabel_1[p] = 0.0;
	}
	
	ntimes--;
	
	if(ntimes > 0){
	  //Dilate mask:
	  nlast = mask_nodes[0];
	  for(j = ninic; j <= mask_nodes[0]; j++){
	    p = mask_nodes[j];
	    px = p%label->ncols;
	    py = p/label->ncols;
	    for(i=1; i < A->n; i++){
	      qx = px + A->dx[i];
	      qy = py + A->dy[i];
	      //if(gft::Image32::IsValidPixel(label, qx,qy)){
	      if(qx >= 0 && qx < label->ncols &&
		 qy >= 0 && qy < label->nrows){
		q = qx + qy*label->ncols;
		
		if(mask->data[q] == 0){
		  mask->data[q] = 1;
		  nlast++;
		  mask_nodes[nlast] = q;
		}
	      }
	    }
	  }
	  ninic = mask_nodes[0] + 1;
	  mask_nodes[0] = nlast;
	}
      }


      for(p = 0; p < n; p++){
	if(flabel_1[p] < 0.5)
	  label->data[p] = 0;
	else
	  label->data[p] = 1;
      }
      
      //free(flabel_1);
      free(flabel_2);
      free(mask_nodes); 
      gft::Image32::Destroy(&mask);

      return flabel_1;
    }
    

    void Relax(sImageGraph *sg,
	       int *S,
	       sImage32 *label, int ntimes){
      sImage32 *mask;
      float *flabel_1,*flabel_2,*tmp;
      sAdjRel *A;
      int *mask_nodes;
      int n,p,q,i,j,k,px,py,qx,qy,nlast,ninic;
      float sw,w;

      gft::ImageGraph::ChangeType(sg, DISSIMILARITY);
      
      A = sg->A;
      ninic = 1;
      n = label->n;
      flabel_1 = (float *)malloc(sizeof(float)*n);
      flabel_2 = (float *)malloc(sizeof(float)*n);
      mask_nodes = (int *)malloc(sizeof(int)*(n+1));
      mask_nodes[0] = 0;
      mask = gft::Image32::Create(label);
      for(p = 0; p < n; p++){
	if(label->data[p] > 0)
	  flabel_1[p] = 1.0;
	else
	  flabel_1[p] = 0.0;
	px = p%label->ncols;
	py = p/label->ncols;
	for(i=1; i < A->n; i++){
	  qx = px + A->dx[i];
	  qy = py + A->dy[i];
	  //if(gft::Image32::IsValidPixel(label, qx,qy)){
	  if(qx >= 0 && qx < label->ncols &&
	     qy >= 0 && qy < label->nrows){
	    q = qx + qy*label->ncols;
	    if(label->data[p] != label->data[q]){
	      mask->data[p] = 1;
	      mask_nodes[0]++;
	      mask_nodes[mask_nodes[0]] = p;
	      break;
	    }
	  }
	}
      }

      while(ntimes > 0){
	//Update flabel:
	memcpy(flabel_2, flabel_1, n*sizeof(float));
	for(j = 1; j <= mask_nodes[0]; j++){
	  p = mask_nodes[j];
	  px = p%label->ncols;
	  py = p/label->ncols;
	  flabel_2[p] = 0.0;
	  sw = 0.0;
	  for(i=1; i < A->n; i++){
	    qx = px + A->dx[i];
	    qy = py + A->dy[i];
	    //if(gft::Image32::IsValidPixel(label, qx,qy)){
	    if(qx >= 0 && qx < label->ncols &&
	       qy >= 0 && qy < label->nrows){
	      q = qx + qy*label->ncols;
	      w = (sg->n_link[p])[i];
	      w = sg->Wmax - w;
	      w = w*w;
	      w = w*w;
	      w = w*w;

	      sw += w;
	      flabel_2[p] += w*flabel_1[q];
	    }
	  }
	  flabel_2[p] /= sw;
	}
	tmp = flabel_1;
	flabel_1 = flabel_2;
	flabel_2 = tmp;

	for(i = 1; i <= S[0]; i++){
	  p = S[i];
	  if(label->data[p] > 0)
	    flabel_1[p] = 1.0;
	  else
	    flabel_1[p] = 0.0;
	}
	
	ntimes--;
	
	if(ntimes > 0){
	  //Dilate mask:
	  nlast = mask_nodes[0];
	  for(j = ninic; j <= mask_nodes[0]; j++){
	    p = mask_nodes[j];
	    px = p%label->ncols;
	    py = p/label->ncols;
	    for(i=1; i < A->n; i++){
	      qx = px + A->dx[i];
	      qy = py + A->dy[i];
	      //if(gft::Image32::IsValidPixel(label, qx,qy)){
	      if(qx >= 0 && qx < label->ncols &&
		 qy >= 0 && qy < label->nrows){
		q = qx + qy*label->ncols;
		
		if(mask->data[q] == 0){
		  mask->data[q] = 1;
		  nlast++;
		  mask_nodes[nlast] = q;
		}
	      }
	    }
	  }
	  ninic = mask_nodes[0] + 1;
	  mask_nodes[0] = nlast;
	}
      }


      for(p = 0; p < n; p++){
	if(flabel_1[p] < 0.5)
	  label->data[p] = 0;
	else
	  label->data[p] = 1;
      }
      
      free(flabel_1);
      free(flabel_2);
      free(mask_nodes); 
      gft::Image32::Destroy(&mask);
    }


    void ORelax_1(sImageGraph *sg,
		  int *S,
		  sImage32 *label, int ntimes){
      sImage32 *mask;
      float *flabel_1,*flabel_2,*tmp;
      sAdjRel *A;
      int *mask_nodes;
      int n,p,q,i,j,k,px,py,qx,qy,nlast,ninic;
      float sw,w;
      int *i_inv = NULL;
      int Imax;
      Imax = gft::Image32::GetMaxVal(label);
      gft::ImageGraph::ChangeType(sg, DISSIMILARITY);
      //printf("Wmax: %d\n", sg->Wmax);
      
      A = sg->A;

      i_inv = gft::AdjRel::InverseIndexes(A);
      
      ninic = 1;
      n = label->n;
      flabel_1 = (float *)malloc(sizeof(float)*n);
      flabel_2 = (float *)malloc(sizeof(float)*n);
      mask_nodes = (int *)malloc(sizeof(int)*(n+1));
      mask_nodes[0] = 0;
      mask = gft::Image32::Create(label);
      for(p = 0; p < n; p++){
	/*
	if(label->data[p] > 0)
	  flabel_1[p] = 1.0;
	else
	  flabel_1[p] = 0.0;
	*/
	flabel_1[p] = (float)label->data[p]/(float)Imax;
	px = p%label->ncols;
	py = p/label->ncols;
	for(i=1; i < A->n; i++){
	  qx = px + A->dx[i];
	  qy = py + A->dy[i];
	  //if(gft::Image32::IsValidPixel(label, qx,qy)){
	  if(qx >= 0 && qx < label->ncols &&
	     qy >= 0 && qy < label->nrows){
	    q = qx + qy*label->ncols;
	    if(label->data[p] != label->data[q]){
	      mask->data[p] = 1;
	      mask_nodes[0]++;
	      mask_nodes[mask_nodes[0]] = p;
	      break;
	    }
	  }
	}
      }

      while(ntimes > 0){
	//Update flabel:
	memcpy(flabel_2, flabel_1, n*sizeof(float));
	for(j = 1; j <= mask_nodes[0]; j++){
	  p = mask_nodes[j];
	  px = p%label->ncols;
	  py = p/label->ncols;
	  flabel_2[p] = 0.0;
	  sw = 0.0;
	  for(i=1; i < A->n; i++){
	    qx = px + A->dx[i];
	    qy = py + A->dy[i];
	    //if(gft::Image32::IsValidPixel(label, qx,qy)){
	    if(qx >= 0 && qx < label->ncols &&
	       qy >= 0 && qy < label->nrows){
	      q = qx + qy*label->ncols;
	      
	      if(flabel_1[q] < 0.5){
		w = (sg->n_link[p])[i];
	      }
	      else{
		k = i_inv[i];
		w = (sg->n_link[q])[k];
	      }
	      
	      w = sg->Wmax - w;
	      w = w*w;
	      w = w*w;
	      w = w*w;

	      sw += w;
	      flabel_2[p] += w*flabel_1[q];
	    }
	  }
	  flabel_2[p] /= sw;
	}
	tmp = flabel_1;
	flabel_1 = flabel_2;
	flabel_2 = tmp;

	for(i = 1; i <= S[0]; i++){
	  p = S[i];
	  flabel_1[p] = (float)label->data[p]/(float)Imax;
	  /*
	  if(label->data[p] > 0)
	    flabel_1[p] = 1.0;
	  else
	    flabel_1[p] = 0.0;
	  */
	}
	
	ntimes--;
	
	if(ntimes > 0){
	  //Dilate mask:
	  nlast = mask_nodes[0];
	  for(j = ninic; j <= mask_nodes[0]; j++){
	    p = mask_nodes[j];
	    px = p%label->ncols;
	    py = p/label->ncols;
	    for(i=1; i < A->n; i++){
	      qx = px + A->dx[i];
	      qy = py + A->dy[i];
	      //if(gft::Image32::IsValidPixel(label, qx,qy)){
	      if(qx >= 0 && qx < label->ncols &&
		 qy >= 0 && qy < label->nrows){
		q = qx + qy*label->ncols;
		
		if(mask->data[q] == 0){
		  mask->data[q] = 1;
		  nlast++;
		  mask_nodes[nlast] = q;
		}
	      }
	    }
	  }
	  ninic = mask_nodes[0] + 1;
	  mask_nodes[0] = nlast;
	}
      }


      for(p = 0; p < n; p++){
	if(flabel_1[p] < 0.5)
	  label->data[p] = 0;
	else
	  label->data[p] = 1;
      }
      
      free(flabel_1);
      free(flabel_2);
      free(mask_nodes); 
      gft::Image32::Destroy(&mask);
      free(i_inv);
    }



    float *ORelax_dual(sImageGraph *sg,
		       int *S,
		       sImage32 *label, int ntimes){
      sImage32 *mask;
      float *flabel_1,*flabel_2,*tmp;
      sAdjRel *A;
      int *mask_nodes;
      int n,p,q,i,j,k,px,py,qx,qy,nlast,ninic;
      float sw,w;
      int *i_inv = NULL;

      A = sg->A;
      i_inv = gft::AdjRel::InverseIndexes(A);
      
      ninic = 1;
      n = label->n;
      flabel_1 = (float *)malloc(sizeof(float)*n);
      flabel_2 = (float *)malloc(sizeof(float)*n);
      mask_nodes = (int *)malloc(sizeof(int)*(n+1));
      mask_nodes[0] = 0;
      mask = gft::Image32::Create(label);
      for(p = 0; p < n; p++){
	if(label->data[p] > 0)
	  flabel_1[p] = 1.0;
	else
	  flabel_1[p] = 0.0;
	px = p%label->ncols;
	py = p/label->ncols;
	for(i=1; i < A->n; i++){
	  qx = px + A->dx[i];
	  qy = py + A->dy[i];
	  //if(gft::Image32::IsValidPixel(label, qx,qy)){
	  if(qx >= 0 && qx < label->ncols &&
	     qy >= 0 && qy < label->nrows){
	    q = qx + qy*label->ncols;
	    if(label->data[p] != label->data[q]){
	      mask->data[p] = 1;
	      mask_nodes[0]++;
	      mask_nodes[mask_nodes[0]] = p;
	      break;
	    }
	  }
	}
      }

      while(ntimes > 0){
	//Update flabel:
	memcpy(flabel_2, flabel_1, n*sizeof(float));
	for(j = 1; j <= mask_nodes[0]; j++){
	  p = mask_nodes[j];
	  px = p%label->ncols;
	  py = p/label->ncols;
	  flabel_2[p] = 0.0;
	  sw = 0.0;
	  for(i=1; i < A->n; i++){
	    qx = px + A->dx[i];
	    qy = py + A->dy[i];
	    //if(gft::Image32::IsValidPixel(label, qx,qy)){
	    if(qx >= 0 && qx < label->ncols &&
	       qy >= 0 && qy < label->nrows){
	      q = qx + qy*label->ncols;

	      /*
	      if(flabel_1[p] >= flabel_1[q]){
		w = (sg->n_link[p])[i];
	      }
	      else{
		k = i_inv[i];
		w = (sg->n_link[q])[k];
	      }
	      */
	      if(flabel_1[q] < 0.5){
		w = (sg->n_link[p])[i];
	      }
	      else{
		k = i_inv[i];
		w = (sg->n_link[q])[k];
	      }
	      
	      //w = w*w;
	      //w = w*w;
	      //w = w*w;

	      sw += w;
	      flabel_2[p] += w*flabel_1[q];
	    }
	  }
	  flabel_2[p] /= sw;
	}
	tmp = flabel_1;
	flabel_1 = flabel_2;
	flabel_2 = tmp;

	for(i = 1; i <= S[0]; i++){
	  p = S[i];
	  if(label->data[p] > 0)
	    flabel_1[p] = 1.0;
	  else
	    flabel_1[p] = 0.0;
	}
	
	ntimes--;
	
	if(ntimes > 0){
	  //Dilate mask:
	  nlast = mask_nodes[0];
	  for(j = ninic; j <= mask_nodes[0]; j++){
	    p = mask_nodes[j];
	    px = p%label->ncols;
	    py = p/label->ncols;
	    for(i=1; i < A->n; i++){
	      qx = px + A->dx[i];
	      qy = py + A->dy[i];
	      //if(gft::Image32::IsValidPixel(label, qx,qy)){
	      if(qx >= 0 && qx < label->ncols &&
		 qy >= 0 && qy < label->nrows){
		q = qx + qy*label->ncols;
		
		if(mask->data[q] == 0){
		  mask->data[q] = 1;
		  nlast++;
		  mask_nodes[nlast] = q;
		}
	      }
	    }
	  }
	  ninic = mask_nodes[0] + 1;
	  mask_nodes[0] = nlast;
	}
      }


      for(p = 0; p < n; p++){
	if(flabel_1[p] < 0.5)
	  label->data[p] = 0;
	else
	  label->data[p] = 1;
      }
      
      //free(flabel_1);
      free(flabel_2);
      free(mask_nodes); 
      gft::Image32::Destroy(&mask);
      free(i_inv);

      return flabel_1;
    }


    

    void ORelax(sImageGraph *sg,
		int *S,
		sImage32 *label, int ntimes){
      sImage32 *mask;
      float *flabel_1,*flabel_2,*tmp;
      sAdjRel *A;
      int *mask_nodes;
      int n,p,q,i,j,k,px,py,qx,qy,nlast,ninic;
      float sw,w;
      int *i_inv = NULL;
      int Imax;
      Imax = gft::Image32::GetMaxVal(label);
      gft::ImageGraph::ChangeType(sg, DISSIMILARITY);
      
      A = sg->A;
      i_inv = gft::AdjRel::InverseIndexes(A);
      
      ninic = 1;
      n = label->n;
      flabel_1 = (float *)malloc(sizeof(float)*n);
      flabel_2 = (float *)malloc(sizeof(float)*n);
      mask_nodes = (int *)malloc(sizeof(int)*(n+1));
      mask_nodes[0] = 0;
      mask = gft::Image32::Create(label);
      for(p = 0; p < n; p++){
	/*
	if(label->data[p] > 0)
	  flabel_1[p] = 1.0;
	else
	  flabel_1[p] = 0.0;
	*/
	flabel_1[p] = (float)label->data[p]/(float)Imax;
	px = p%label->ncols;
	py = p/label->ncols;
	for(i=1; i < A->n; i++){
	  qx = px + A->dx[i];
	  qy = py + A->dy[i];
	  //if(gft::Image32::IsValidPixel(label, qx,qy)){
	  if(qx >= 0 && qx < label->ncols &&
	     qy >= 0 && qy < label->nrows){
	    q = qx + qy*label->ncols;
	    if(label->data[p] != label->data[q]){
	      mask->data[p] = 1;
	      mask_nodes[0]++;
	      mask_nodes[mask_nodes[0]] = p;
	      break;
	    }
	  }
	}
      }

      while(ntimes > 0){
	//Update flabel:
	memcpy(flabel_2, flabel_1, n*sizeof(float));
	for(j = 1; j <= mask_nodes[0]; j++){
	  p = mask_nodes[j];
	  px = p%label->ncols;
	  py = p/label->ncols;
	  flabel_2[p] = 0.0;
	  sw = 0.0;
	  for(i=1; i < A->n; i++){
	    qx = px + A->dx[i];
	    qy = py + A->dy[i];
	    //if(gft::Image32::IsValidPixel(label, qx,qy)){
	    if(qx >= 0 && qx < label->ncols &&
	       qy >= 0 && qy < label->nrows){
	      q = qx + qy*label->ncols;
	      
	      if(flabel_1[p] >= flabel_1[q]){
		w = (sg->n_link[p])[i];
	      }
	      else{
		k = i_inv[i];
		w = (sg->n_link[q])[k];
	      }
	      
	      w = sg->Wmax - w;
	      w = w*w;
	      w = w*w;
	      w = w*w;

	      sw += w;
	      flabel_2[p] += w*flabel_1[q];
	    }
	  }
	  flabel_2[p] /= sw;
	}
	tmp = flabel_1;
	flabel_1 = flabel_2;
	flabel_2 = tmp;

	for(i = 1; i <= S[0]; i++){
	  p = S[i];
	  flabel_1[p] = (float)label->data[p]/(float)Imax;
	  /*
	  if(label->data[p] > 0)
	    flabel_1[p] = 1.0;
	  else
	    flabel_1[p] = 0.0;
	  */
	}
	
	ntimes--;
	
	if(ntimes > 0){
	  //Dilate mask:
	  nlast = mask_nodes[0];
	  for(j = ninic; j <= mask_nodes[0]; j++){
	    p = mask_nodes[j];
	    px = p%label->ncols;
	    py = p/label->ncols;
	    for(i=1; i < A->n; i++){
	      qx = px + A->dx[i];
	      qy = py + A->dy[i];
	      //if(gft::Image32::IsValidPixel(label, qx,qy)){
	      if(qx >= 0 && qx < label->ncols &&
		 qy >= 0 && qy < label->nrows){
		q = qx + qy*label->ncols;
		
		if(mask->data[q] == 0){
		  mask->data[q] = 1;
		  nlast++;
		  mask_nodes[nlast] = q;
		}
	      }
	    }
	  }
	  ninic = mask_nodes[0] + 1;
	  mask_nodes[0] = nlast;
	}
      }


      for(p = 0; p < n; p++){
	if(flabel_1[p] < 0.5)
	  label->data[p] = 0;
	else
	  label->data[p] = 1;
      }
      
      free(flabel_1);
      free(flabel_2);
      free(mask_nodes); 
      gft::Image32::Destroy(&mask);
      free(i_inv);
    }




    void ORelax_i(sImageGraph *sg,
		  int *S,
		  sImage32 *label, int ntimes){
      sImage32 *mask;
      float *flabel_1,*flabel_2,*tmp;
      sAdjRel *A;
      int *mask_nodes;
      int n,p,q,i,j,k,px,py,qx,qy,nlast,ninic;
      float sw,w;
      int *i_inv = NULL;

      gft::ImageGraph::ChangeType(sg, DISSIMILARITY);
      
      A = sg->A;
      i_inv = gft::AdjRel::InverseIndexes(A);
      
      ninic = 1;
      n = label->n;
      flabel_1 = (float *)malloc(sizeof(float)*n);
      flabel_2 = (float *)malloc(sizeof(float)*n);
      mask_nodes = (int *)malloc(sizeof(int)*(n+1));
      mask_nodes[0] = 0;
      mask = gft::Image32::Create(label);
      for(p = 0; p < n; p++){
	if(label->data[p] > 0)
	  flabel_1[p] = 1.0;
	else
	  flabel_1[p] = 0.0;
	px = p%label->ncols;
	py = p/label->ncols;
	for(i=1; i < A->n; i++){
	  qx = px + A->dx[i];
	  qy = py + A->dy[i];
	  //if(gft::Image32::IsValidPixel(label, qx,qy)){
	  if(qx >= 0 && qx < label->ncols &&
	     qy >= 0 && qy < label->nrows){
	    q = qx + qy*label->ncols;
	    if(label->data[p] != label->data[q]){
	      mask->data[p] = 1;
	      mask_nodes[0]++;
	      mask_nodes[mask_nodes[0]] = p;
	      break;
	    }
	  }
	}
      }

      while(ntimes > 0){
	//Update flabel:
	memcpy(flabel_2, flabel_1, n*sizeof(float));
	for(j = 1; j <= mask_nodes[0]; j++){
	  p = mask_nodes[j];
	  px = p%label->ncols;
	  py = p/label->ncols;
	  flabel_2[p] = 0.0;
	  sw = 0.0;
	  for(i=1; i < A->n; i++){
	    qx = px + A->dx[i];
	    qy = py + A->dy[i];
	    //if(gft::Image32::IsValidPixel(label, qx,qy)){
	    if(qx >= 0 && qx < label->ncols &&
	       qy >= 0 && qy < label->nrows){
	      q = qx + qy*label->ncols;
	      
	      if(flabel_1[p] > flabel_1[q]){
		w = (sg->n_link[p])[i];
	      }
	      else{
		k = i_inv[i];
		w = (sg->n_link[q])[k];
	      }
	      
	      w = sg->Wmax - w;
	      w = w*w;
	      w = w*w;
	      w = w*w;

	      sw += w;
	      flabel_2[p] += w*flabel_1[q];
	    }
	  }
	  flabel_2[p] /= sw;
	}
	tmp = flabel_1;
	flabel_1 = flabel_2;
	flabel_2 = tmp;

	for(i = 1; i <= S[0]; i++){
	  p = S[i];
	  if(label->data[p] > 0)
	    flabel_1[p] = 1.0;
	  else
	    flabel_1[p] = 0.0;
	}
	
	ntimes--;
	
	if(ntimes > 0){
	  //Dilate mask:
	  nlast = mask_nodes[0];
	  for(j = ninic; j <= mask_nodes[0]; j++){
	    p = mask_nodes[j];
	    px = p%label->ncols;
	    py = p/label->ncols;
	    for(i=1; i < A->n; i++){
	      qx = px + A->dx[i];
	      qy = py + A->dy[i];
	      //if(gft::Image32::IsValidPixel(label, qx,qy)){
	      if(qx >= 0 && qx < label->ncols &&
		 qy >= 0 && qy < label->nrows){
		q = qx + qy*label->ncols;
		
		if(mask->data[q] == 0){
		  mask->data[q] = 1;
		  nlast++;
		  mask_nodes[nlast] = q;
		}
	      }
	    }
	  }
	  ninic = mask_nodes[0] + 1;
	  mask_nodes[0] = nlast;
	}
      }


      for(p = 0; p < n; p++){
	if(flabel_1[p] < 0.5)
	  label->data[p] = 0;
	else
	  label->data[p] = 1;
      }
      
      free(flabel_1);
      free(flabel_2);
      free(mask_nodes); 
      gft::Image32::Destroy(&mask);
      free(i_inv);
    }



    void ORelax_s(sImageGraph *sg,
		  int *S,
		  sImage32 *label, int ntimes){
      sImage32 *mask;
      float *flabel_1,*flabel_2,*tmp;
      sAdjRel *A;
      int *mask_nodes;
      int n,p,q,i,j,k,px,py,qx,qy,nlast,ninic;
      float sw,w;
      int *i_inv = NULL;
      //printf("Entrou\n");
      
      gft::ImageGraph::ChangeType(sg, DISSIMILARITY);
      
      A = sg->A;
      i_inv = gft::AdjRel::InverseIndexes(A);
      
      ninic = 1;
      n = label->n;
      flabel_1 = (float *)malloc(sizeof(float)*n);
      flabel_2 = (float *)malloc(sizeof(float)*n);
      mask_nodes = (int *)malloc(sizeof(int)*(n+1));
      mask_nodes[0] = 0;
      mask = gft::Image32::Create(label);
      for(p = 0; p < n; p++){
	if(label->data[p] > 0)
	  flabel_1[p] = 1.0;
	else
	  flabel_1[p] = 0.0;
	px = p%label->ncols;
	py = p/label->ncols;
	for(i=1; i < A->n; i++){
	  qx = px + A->dx[i];
	  qy = py + A->dy[i];
	  //if(gft::Image32::IsValidPixel(label, qx,qy)){
	  if(qx >= 0 && qx < label->ncols &&
	     qy >= 0 && qy < label->nrows){
	    q = qx + qy*label->ncols;
	    if(label->data[p] != label->data[q]){
	      mask->data[p] = 1;
	      mask_nodes[0]++;
	      mask_nodes[mask_nodes[0]] = p;
	      break;
	    }
	  }
	}
      }

      while(ntimes > 0){
	//Update flabel:
	memcpy(flabel_2, flabel_1, n*sizeof(float));
	for(j = 1; j <= mask_nodes[0]; j++){
	  p = mask_nodes[j];
	  px = p%label->ncols;
	  py = p/label->ncols;
	  flabel_2[p] = 0.0;
	  sw = 0.0;
	  for(i=1; i < A->n; i++){
	    qx = px + A->dx[i];
	    qy = py + A->dy[i];
	    //if(gft::Image32::IsValidPixel(label, qx,qy)){
	    if(qx >= 0 && qx < label->ncols &&
	       qy >= 0 && qy < label->nrows){
	      q = qx + qy*label->ncols;

	      if(flabel_1[p] == flabel_1[q]){
		k = i_inv[i];
		w = ((sg->n_link[p])[i] + (sg->n_link[q])[k])/2.;
	      }
	      else if(flabel_1[p] > flabel_1[q]){
		w = (sg->n_link[p])[i];
	      }
	      else{
		k = i_inv[i];
		w = (sg->n_link[q])[k];
	      }
	      
	      w = sg->Wmax - w;
	      w = w*w;
	      w = w*w;
	      w = w*w;

	      sw += w;
	      flabel_2[p] += w*flabel_1[q];
	    }
	  }
	  flabel_2[p] /= sw;
	}
	tmp = flabel_1;
	flabel_1 = flabel_2;
	flabel_2 = tmp;

	for(i = 1; i <= S[0]; i++){
	  p = S[i];
	  if(label->data[p] > 0)
	    flabel_1[p] = 1.0;
	  else
	    flabel_1[p] = 0.0;
	}
	
	ntimes--;
	
	if(ntimes > 0){
	  //Dilate mask:
	  nlast = mask_nodes[0];
	  for(j = ninic; j <= mask_nodes[0]; j++){
	    p = mask_nodes[j];
	    px = p%label->ncols;
	    py = p/label->ncols;
	    for(i=1; i < A->n; i++){
	      qx = px + A->dx[i];
	      qy = py + A->dy[i];
	      //if(gft::Image32::IsValidPixel(label, qx,qy)){
	      if(qx >= 0 && qx < label->ncols &&
		 qy >= 0 && qy < label->nrows){
		q = qx + qy*label->ncols;
		
		if(mask->data[q] == 0){
		  mask->data[q] = 1;
		  nlast++;
		  mask_nodes[nlast] = q;
		}
	      }
	    }
	  }
	  ninic = mask_nodes[0] + 1;
	  mask_nodes[0] = nlast;
	}
      }


      for(p = 0; p < n; p++){
	if(flabel_1[p] < 0.5)
	  label->data[p] = 0;
	else
	  label->data[p] = 1;
      }
      
      free(flabel_1);
      free(flabel_2);
      free(mask_nodes); 
      gft::Image32::Destroy(&mask);
      free(i_inv);
    }




    void ORelax_1(sAdjRel3 *A,
		  sScene32 *scn,
		  float per,
		  int *S,
		  sScene32 *label,
		  int ntimes){
      sScene32 *mask;
      float *flabel_1,*flabel_2,*tmp;
      int *mask_nodes;
      int n,p,q,i,j,k,nlast,ninic;
      Voxel u,v;
      float sw,w,per_pq,dmin;
      int Imax,Wmax;
      float *Dpq;

      //--------------------
      dmin = MIN(scn->dx, MIN(scn->dy, scn->dz));
      Dpq = (float *)malloc(A->n*sizeof(float));
      for(i=1; i<A->n; i++){
	Dpq[i] = sqrtf(A->d[i].axis.x*A->d[i].axis.x*scn->dx*scn->dx + 
		       A->d[i].axis.y*A->d[i].axis.y*scn->dy*scn->dy + 
		       A->d[i].axis.z*A->d[i].axis.z*scn->dz*scn->dz)/dmin;
      }
      //--------------------
      Imax = gft::Scene32::GetMaximumValue(label);
      Wmax = gft::Scene32::GetMaximumValue(scn);
      Wmax *= (1.0 + fabsf(per)/100.0);
      ninic = 1;
      n = label->n;
      flabel_1 = (float *)malloc(sizeof(float)*n);
      flabel_2 = (float *)malloc(sizeof(float)*n);
      mask_nodes = (int *)malloc(sizeof(int)*(n+1));
      mask_nodes[0] = 0;
      mask = gft::Scene32::Create(label);
      for(p = 0; p < n; p++){
	flabel_1[p] = (float)label->data[p]/(float)Imax;
	u.c.x = gft::Scene32::GetAddressX(label,p);
	u.c.y = gft::Scene32::GetAddressY(label,p);
	u.c.z = gft::Scene32::GetAddressZ(label,p);
	for(i=1; i < A->n; i++){
	  v.v = u.v + A->d[i].v;
	  if(gft::Scene32::IsValidVoxel(label, v)){
	    q = gft::Scene32::GetVoxelAddress(label,v);
	    if(label->data[p] != label->data[q]){
	      mask->data[p] = 1;
	      mask_nodes[0]++;
	      mask_nodes[mask_nodes[0]] = p;
	      break;
	    }
	  }
	}
      }

      while(ntimes > 0){
	//Update flabel:
	memcpy(flabel_2, flabel_1, n*sizeof(float));
	for(j = 1; j <= mask_nodes[0]; j++){
	  p = mask_nodes[j];
	  u.c.x = gft::Scene32::GetAddressX(label,p);
	  u.c.y = gft::Scene32::GetAddressY(label,p);
	  u.c.z = gft::Scene32::GetAddressZ(label,p);
	  flabel_2[p] = 0.0;
	  sw = 0.0;
	  for(i=1; i < A->n; i++){
	    v.v = u.v + A->d[i].v;
	    if(gft::Scene32::IsValidVoxel(label, v)){
	      q = gft::Scene32::GetVoxelAddress(label,v);
	      w = abs(scn->data[p] - scn->data[q]);
	      if(flabel_1[q] < 0.5)
		per_pq = per;
	      else
		per_pq = -per;

	      if(scn->data[p] > scn->data[q])
		w *= (1.0 + per_pq/100.0);
	      else if(scn->data[p] < scn->data[q])
		w *= (1.0 - per_pq/100.0);
		
	      w = Wmax - w;
	      w = w*w;
	      w = w*w;
	      w = w*w;
	      w /= Dpq[i];
	      
	      sw += w;
	      flabel_2[p] += w*flabel_1[q];
	    }
	  }
	  flabel_2[p] /= sw;
	}
	tmp = flabel_1;
	flabel_1 = flabel_2;
	flabel_2 = tmp;

	for(i = 1; i <= S[0]; i++){
	  p = S[i];
	  flabel_1[p] = (float)label->data[p]/(float)Imax;
	}
	
	ntimes--;
	
	if(ntimes > 0){
	  //Dilate mask:
	  nlast = mask_nodes[0];
	  for(j = ninic; j <= mask_nodes[0]; j++){
	    p = mask_nodes[j];
	    u.c.x = gft::Scene32::GetAddressX(label,p);
	    u.c.y = gft::Scene32::GetAddressY(label,p);
	    u.c.z = gft::Scene32::GetAddressZ(label,p);
	    for(i=1; i < A->n; i++){
	      v.v = u.v + A->d[i].v;
	      if(gft::Scene32::IsValidVoxel(label, v)){
		q = gft::Scene32::GetVoxelAddress(label,v);
		if(mask->data[q] == 0){
		  mask->data[q] = 1;
		  nlast++;
		  mask_nodes[nlast] = q;
		}
	      }
	    }
	  }
	  ninic = mask_nodes[0] + 1;
	  mask_nodes[0] = nlast;
	}
      }

      for(p = 0; p < n; p++){
	if(flabel_1[p] < 0.5)
	  label->data[p] = 0;
	else
	  label->data[p] = 1;
      }
      
      free(flabel_1);
      free(flabel_2);
      free(mask_nodes); 
      gft::Scene32::Destroy(&mask);
    }

    


    void ORelax_1(sScene32 *W,
		  sAdjRel3 *A,
		  sScene32 *scn,
		  float per,
		  int *S,
		  sScene32 *label,
		  int ntimes){
      sScene32 *mask;
      float *flabel_1,*flabel_2,*tmp;
      int *mask_nodes;
      int n,p,q,i,j,k,nlast,ninic;
      Voxel u,v;
      float sw,w,per_pq;
      int Imax,Wmax;
      Imax = gft::Scene32::GetMaximumValue(label);
      Wmax = gft::Scene32::GetMaximumValue(W)*2;
      Wmax *= (1.0 + fabsf(per)/100.0);
      ninic = 1;
      n = label->n;
      flabel_1 = (float *)malloc(sizeof(float)*n);
      flabel_2 = (float *)malloc(sizeof(float)*n);
      mask_nodes = (int *)malloc(sizeof(int)*(n+1));
      mask_nodes[0] = 0;
      mask = gft::Scene32::Create(label);
      for(p = 0; p < n; p++){
	flabel_1[p] = (float)label->data[p]/(float)Imax;
	u.c.x = gft::Scene32::GetAddressX(label,p);
	u.c.y = gft::Scene32::GetAddressY(label,p);
	u.c.z = gft::Scene32::GetAddressZ(label,p);
	for(i=1; i < A->n; i++){
	  v.v = u.v + A->d[i].v;
	  if(gft::Scene32::IsValidVoxel(label, v)){
	    q = gft::Scene32::GetVoxelAddress(label,v);
	    if(label->data[p] != label->data[q]){
	      mask->data[p] = 1;
	      mask_nodes[0]++;
	      mask_nodes[mask_nodes[0]] = p;
	      break;
	    }
	  }
	}
      }

      while(ntimes > 0){
	//Update flabel:
	memcpy(flabel_2, flabel_1, n*sizeof(float));
	for(j = 1; j <= mask_nodes[0]; j++){
	  p = mask_nodes[j];
	  u.c.x = gft::Scene32::GetAddressX(label,p);
	  u.c.y = gft::Scene32::GetAddressY(label,p);
	  u.c.z = gft::Scene32::GetAddressZ(label,p);
	  flabel_2[p] = 0.0;
	  sw = 0.0;
	  for(i=1; i < A->n; i++){
	    v.v = u.v + A->d[i].v;
	    if(gft::Scene32::IsValidVoxel(label, v)){
	      q = gft::Scene32::GetVoxelAddress(label,v);
	      w = W->data[p] + W->data[q];
	      if(flabel_1[q] < 0.5)
		per_pq = per;
	      else
		per_pq = -per;

	      if(scn->data[p] > scn->data[q])
		w *= (1.0 + per_pq/100.0);
	      else if(scn->data[p] < scn->data[q])
		w *= (1.0 - per_pq/100.0);
		
	      w = Wmax - w;
	      w = w*w;
	      w = w*w;
	      w = w*w;

	      sw += w;
	      flabel_2[p] += w*flabel_1[q];
	    }
	  }
	  flabel_2[p] /= sw;
	}
	tmp = flabel_1;
	flabel_1 = flabel_2;
	flabel_2 = tmp;

	for(i = 1; i <= S[0]; i++){
	  p = S[i];
	  flabel_1[p] = (float)label->data[p]/(float)Imax;
	}
	
	ntimes--;
	
	if(ntimes > 0){
	  //Dilate mask:
	  nlast = mask_nodes[0];
	  for(j = ninic; j <= mask_nodes[0]; j++){
	    p = mask_nodes[j];
	    u.c.x = gft::Scene32::GetAddressX(label,p);
	    u.c.y = gft::Scene32::GetAddressY(label,p);
	    u.c.z = gft::Scene32::GetAddressZ(label,p);
	    for(i=1; i < A->n; i++){
	      v.v = u.v + A->d[i].v;
	      if(gft::Scene32::IsValidVoxel(label, v)){
		q = gft::Scene32::GetVoxelAddress(label,v);
		if(mask->data[q] == 0){
		  mask->data[q] = 1;
		  nlast++;
		  mask_nodes[nlast] = q;
		}
	      }
	    }
	  }
	  ninic = mask_nodes[0] + 1;
	  mask_nodes[0] = nlast;
	}
      }

      for(p = 0; p < n; p++){
	if(flabel_1[p] < 0.5)
	  label->data[p] = 0;
	else
	  label->data[p] = 1;
      }
      
      free(flabel_1);
      free(flabel_2);
      free(mask_nodes); 
      gft::Scene32::Destroy(&mask);
    }



    

    void ORelax_1(sScene32 *Wx,
		  sScene32 *Wy,
		  sScene32 *Wz,
		  sScene32 *scn,
		  float per,
		  int *S,
		  sScene32 *label,
		  int ntimes){
      sScene32 *mask;
      float *flabel_1,*flabel_2,*tmp;
      int *mask_nodes;
      int n,p,q,i,j,k,nlast,ninic,t;
      Voxel u,v;
      float sw,w,per_pq;
      int Imax,Wmax;
      sAdjRel3 *A;
      sScene32 *W[3];
      int T[7];
      W[0] = Wx;
      W[1] = Wy;
      W[2] = Wz;
      A = gft::AdjRel3::Spheric(1.0);
      for(i=1; i<A->n; i++){
	if(A->d[i].axis.x != 0) T[i] = 0;
	if(A->d[i].axis.y != 0) T[i] = 1;
	if(A->d[i].axis.z != 0) T[i] = 2;
      }
      Wmax = MAX(gft::Scene32::GetMaximumValue(Wx),
		 MAX(gft::Scene32::GetMaximumValue(Wy),
		     gft::Scene32::GetMaximumValue(Wz)))*2;
      Imax = gft::Scene32::GetMaximumValue(label);

      Wmax *= (1.0 + fabsf(per)/100.0);
      ninic = 1;
      n = label->n;
      flabel_1 = (float *)malloc(sizeof(float)*n);
      flabel_2 = (float *)malloc(sizeof(float)*n);
      mask_nodes = (int *)malloc(sizeof(int)*(n+1));
      mask_nodes[0] = 0;
      mask = gft::Scene32::Create(label);
      for(p = 0; p < n; p++){
	flabel_1[p] = (float)label->data[p]/(float)Imax;
	u.c.x = gft::Scene32::GetAddressX(label,p);
	u.c.y = gft::Scene32::GetAddressY(label,p);
	u.c.z = gft::Scene32::GetAddressZ(label,p);
	for(i=1; i < A->n; i++){
	  v.v = u.v + A->d[i].v;
	  if(gft::Scene32::IsValidVoxel(label, v)){
	    q = gft::Scene32::GetVoxelAddress(label,v);
	    if(label->data[p] != label->data[q]){
	      mask->data[p] = 1;
	      mask_nodes[0]++;
	      mask_nodes[mask_nodes[0]] = p;
	      break;
	    }
	  }
	}
      }

      while(ntimes > 0){
	//Update flabel:
	memcpy(flabel_2, flabel_1, n*sizeof(float));
	for(j = 1; j <= mask_nodes[0]; j++){
	  p = mask_nodes[j];
	  u.c.x = gft::Scene32::GetAddressX(label,p);
	  u.c.y = gft::Scene32::GetAddressY(label,p);
	  u.c.z = gft::Scene32::GetAddressZ(label,p);
	  flabel_2[p] = 0.0;
	  sw = 0.0;
	  for(i=1; i < A->n; i++){
	    v.v = u.v + A->d[i].v;
	    if(gft::Scene32::IsValidVoxel(label, v)){
	      q = gft::Scene32::GetVoxelAddress(label,v);
	      t = T[i];
	      w = W[t]->data[p] + W[t]->data[q];
	      if(flabel_1[q] < 0.5)
		per_pq = per;
	      else
		per_pq = -per;

	      if(scn->data[p] > scn->data[q])
		w *= (1.0 + per_pq/100.0);
	      else if(scn->data[p] < scn->data[q])
		w *= (1.0 - per_pq/100.0);
		
	      w = Wmax - w;
	      w = w*w;
	      w = w*w;
	      w = w*w;

	      sw += w;
	      flabel_2[p] += w*flabel_1[q];
	    }
	  }
	  flabel_2[p] /= sw;
	}
	tmp = flabel_1;
	flabel_1 = flabel_2;
	flabel_2 = tmp;

	for(i = 1; i <= S[0]; i++){
	  p = S[i];
	  flabel_1[p] = (float)label->data[p]/(float)Imax;
	}
	
	ntimes--;
	
	if(ntimes > 0){
	  //Dilate mask:
	  nlast = mask_nodes[0];
	  for(j = ninic; j <= mask_nodes[0]; j++){
	    p = mask_nodes[j];
	    u.c.x = gft::Scene32::GetAddressX(label,p);
	    u.c.y = gft::Scene32::GetAddressY(label,p);
	    u.c.z = gft::Scene32::GetAddressZ(label,p);
	    for(i=1; i < A->n; i++){
	      v.v = u.v + A->d[i].v;
	      if(gft::Scene32::IsValidVoxel(label, v)){
		q = gft::Scene32::GetVoxelAddress(label,v);
		if(mask->data[q] == 0){
		  mask->data[q] = 1;
		  nlast++;
		  mask_nodes[nlast] = q;
		}
	      }
	    }
	  }
	  ninic = mask_nodes[0] + 1;
	  mask_nodes[0] = nlast;
	}
      }

      for(p = 0; p < n; p++){
	if(flabel_1[p] < 0.5)
	  label->data[p] = 0;
	else
	  label->data[p] = 1;
      }
      
      free(flabel_1);
      free(flabel_2);
      free(mask_nodes); 
      gft::Scene32::Destroy(&mask);
      gft::AdjRel3::Destroy(&A);
    }


    
    
    void GGC_maxmin(sImageGraph *graph,
		    sCImage *cimg,
		    sImage32 *img,
		    int method,
		    float power_geodesic,
		    float delta,
		    float theta_hedgehog,
		    float R,
		    int postproc,
		    int niterations,
		    int conn,
		    int pol,
		    int shapepriors,
		    int costtemplate,
		    int *S,
		    sImage32 *label){
      int *S1=NULL;
      sGQueue *Q=NULL;
      sImageGraph *sg = NULL;
      sGraph *g = NULL, *transpose = NULL;
      sImage32 *P_sum = NULL, *C_sum = NULL, *raw_map = NULL;
      sImage32 *pred, *cost, *Slabel;
      int i,p,n,ns1;
      //-------------------------
      gft::timer tic,toc;
      float totaltime;
      FILE *fp;
      int reduction = 0;
      //-------------------------      
      float A=750.0, B=75.0, n1=2.0, n2=1.0, n3=1.0;
      int m=2;
      
      Slabel = gft::Image32::Clone(label);
      pred = Image32::Create(label);
      cost = Image32::Create(label);
      Image32::Set(pred, NIL);
      Image32::Set(cost, INT_MAX);
      
      S1 = (int *)malloc((S[0]+1)*sizeof(int));
      ns1 = 0;
      for(i = 1; i <= S[0]; i++){
	p = S[i];
	if(label->data[p] == 1){
	  ns1++;
	  S1[ns1] = p;
	}
      }
      S1[0] = ns1;
      S1 = (int *)realloc(S1, (S1[0]+1)*sizeof(int));

      ImageGraph::ChangeType(graph, DISSIMILARITY);
      
      //------------------------
      if(pol == 0 && shapepriors <= 1 && !conn){
	switch(shapepriors){
	case 0:
	  IFT_fw(graph, S, label, cost, pred);
	  break;
	case 1:
	  P_sum = SC_Pred_fsum(graph, S1, power_geodesic);
	  SC_IFT(graph, S, label, P_sum);
	  Image32::Destroy(&P_sum);
	  break;
	}
      }
      else{
	sg = ImageGraph::Clone(graph);
	if(cimg != NULL){ //COLOR_IMAGE
	  sImage32 *lumi;
	  lumi = CImage::Luminosity(cimg);
	  ImageGraph::Orient2Digraph(sg, lumi, pol);
	  Image32::Destroy(&lumi);
	}
	else
	  ImageGraph::Orient2Digraph(sg, img, pol);
	
	if(shapepriors != 0){
	  if(costtemplate == 0){
	    P_sum = SC_Pred_fsum(sg, S1, power_geodesic);
	    C_sum = BB_Geodesic_Cost(P_sum, graph->A);
	  }
	  else{
	    switch(costtemplate){
	    case 1: //Circle
	      raw_map = Image32::Read("./templates/circle.pgm"); break;
	    case 2: //Square
	      raw_map = Image32::Read("./templates/square.pgm"); break;
	    case 3:
	      raw_map = Image32::Read("./templates/ellipse2:1.pgm"); break;
	    case 4:
	      raw_map = Image32::Read("./templates/ellipse3:1.pgm"); break;
	    case 5:
	      raw_map = Image32::Read("./templates/rectangle2:1.pgm"); break;
	    case 6:
	      raw_map = Image32::Read("./templates/rectangle3:1.pgm"); break;
	    case 7:
	      fp = fopen("gielis.txt", "r");
	      if(fp != NULL){
		fscanf(fp," %f %f %f %f %f %d", &A, &B, &n1, &n2, &n3, &m);
		fclose(fp);
	      }
	      raw_map = gft::Band::GielisEquationTemplate(2001, A, B,
							  n1, n2, n3, m);
	      break;	      
	    }
	    
	    C_sum = BB_CropTemplate(raw_map, S, label, 1);
	    P_sum = NULL;
	    Image32::Destroy(&raw_map);
	  }
	}
	
	if(shapepriors == 1 && costtemplate == 0){ //Star Convexity
	  ImageGraph::Orient2DigraphOuter(sg, P_sum);
	}
	else if(shapepriors == 2){ //Hedgehog
	  if(R <= 1.6 && sg->A->n == 9)
	    ImageGraph::Convert2HedgehogDigraph(sg, C_sum, theta_hedgehog);
	  else{
	    g = Graph::Clone(sg);
	    Graph::HedgehogDigraph(g, C_sum, theta_hedgehog, R);
	    transpose = Graph::Transpose(g);
	  }
	}
	else if(shapepriors == 3){ //Local Band
	  g = Graph::Clone(sg);
	  printf("NumberOfArcs: %d (image graph)\n", Graph::GetNumberOfArcs(g));
	  Graph::LocalBandConstraint(g, C_sum, ROUND(delta), R);
	  printf("NumberOfArcs: %d (local band)\n", Graph::GetNumberOfArcs(g));
	  fp = fopen("tmp.txt", "r");
	  if(fp !=NULL){
	    fscanf(fp, "%d", &reduction);
	    fclose(fp);
	  }
	  if(reduction){
	    gettimeofday(&tic,NULL);
	    Graph::LocalBandReduction(&g, C_sum, ROUND(delta));
	    gettimeofday(&toc,NULL);
	    totaltime = (toc.tv_sec-tic.tv_sec)*1000.0 + (toc.tv_usec-tic.tv_usec)*0.001;
	    printf("\nReduction Time: %f ms\n",totaltime);	    
	    printf("NumberOfArcs: %d (reduction)\n", Graph::GetNumberOfArcs(g));
	  }
	  transpose = Graph::Transpose(g);
	}
	
	if(conn){
	  COIFT(sg, S, label);
	}
	else if(shapepriors < 4 && g == NULL){
	  gettimeofday(&tic,NULL);
	  switch(method){
	  case 0:
	    OIFT(sg, NULL, label); break;
	  case 1:
	    EOIFT(sg, NULL, label, 0); break;
	  case 2:
	    OIFT_Heap(sg, NULL, label);	break;
	  case 3:
	    EOIFT_Heap(sg, NULL, label, 0.0); break;
	  case 4:
	    OIFT_TZ2Bkg(sg, NULL, label); break;
	  case 5:
	    OIFT_TZ2Obj(sg, NULL, label); break;	    
	  case 6:
	    EOIFT_Heap_2(sg, NULL, label, 0.0); break;
	  }
	  gettimeofday(&toc,NULL);
	  totaltime = (toc.tv_sec-tic.tv_sec)*1000.0 + (toc.tv_usec-tic.tv_usec)*0.001;
	  printf("\nInner Time: %f ms\n",totaltime);

	  /****** REMOVER *****/

	  printf("energy: %d\n", GetEnergy_Min(sg, label, 1));
	  
	  if( isOIFT_Segmentation(sg, S, Slabel, label) )
	    printf("isOIFT: OK\n");
	  else
	    printf("isOIFT: Error\n");
	  
	  /********************/
	}
	else if(shapepriors < 4 && g != NULL){
	  gettimeofday(&tic,NULL);
	  switch(method){
	  case 0:
	    OIFT(g, transpose, NULL, label->data); break;
	  case 1:
	    EOIFT(g, transpose, NULL, label->data, 0); break;
	  case 2:
	    OIFT_Heap(g, transpose, NULL, label->data); break;
	  case 3:
	    EOIFT_Heap(g, transpose, NULL, label->data, 0.0); break;
	  case 6:
	    EOIFT_Heap_2(g, transpose, NULL, label->data, 0.0); break;
	  }
	  gettimeofday(&toc,NULL);
	  totaltime = (toc.tv_sec-tic.tv_sec)*1000.0 + (toc.tv_usec-tic.tv_usec)*0.001;
	  printf("\nInner Time: %f ms\n",totaltime);
	  printf("Energy: %d\n", GetEnergy_Min(g, label->data, 1));
	}
	else{
	  if(shapepriors == 4)
	    BB_OIFT(sg, S, label,
		    C_sum, P_sum, delta);
	  else if(shapepriors == 5){
	    //---code included to fix errors with zero costs---
	    //---because delta*C(t) cannot be zero-------------
	    for(p = 0; p < C_sum->n; p++) 
	      C_sum->data[p] += 1;
	    //-------------------------------------------------
	    RBB_OIFT(sg, S, label,
		     C_sum, P_sum, delta);
	  }
	  else if(shapepriors == 6){
	    B_OIFT(sg, S, label, C_sum, delta);
	  }
	}
	//printf("energy: %d\n", gft::ift::GetEnergy_Min(sg, APP->Data.label, 1));
      }
      
      //--------REMOVER----------------
      /*      
      if(sg != NULL)
	ImageGraph::Destroy(&sg);
      sg = ImageGraph::Clone(graph);
      ImageGraph::Orient2Digraph(sg, img, pol);
      */
      //-------------------------------
      
      switch(postproc){
      case 1:
	if(sg != NULL)
	  Relax(sg, S, label, niterations);
	else
	  Relax(graph, S, label, niterations);
	break;
      case 2:
	if(sg != NULL)
	  ORelax(sg, S, label, niterations);
	else
	  ORelax(graph, S, label, niterations);
	break;
      case 3:
	if(sg != NULL)
	  ORelax_1(sg, S, label, niterations);
	else
	  ORelax_1(graph, S, label, niterations);
	break;
      case 4:
	Image32::ModeFilterLabel(label, (float)niterations);
	break;
      case 5:
	if(sg != NULL)
	  ORelax_i(sg, S, label, niterations);
	else
	  ORelax_i(graph, S, label, niterations);
	break;
      case 6:
	if(sg != NULL)
	  ORelax_s(sg, S, label, niterations);
	else
	  ORelax_s(graph, S, label, niterations);
	break;
      case 7:
	Relax(graph, S, label, niterations);
	break;
      }

      /****** REMOVER *****/
      /*
      if(sg != NULL){
	if( isOIFT(sg, S, Slabel, label) )
	  printf("pos isOIFT: OK\n");
	else
	  printf("pos isOIFT: Error\n");
      }
      */
      /********************/
      
      //------------------------

      if(Slabel != NULL) Image32::Destroy(&Slabel);
      
      if(sg != NULL) ImageGraph::Destroy(&sg);
      if(g  != NULL) Graph::Destroy(&g);
      if(transpose != NULL) Graph::Destroy(&transpose);
      if(C_sum != NULL) Image32::Destroy(&C_sum);
      if(P_sum != NULL) Image32::Destroy(&P_sum);
      if(pred != NULL) Image32::Destroy(&pred);
      if(cost != NULL) Image32::Destroy(&cost);
      free(S1);
    }



    //--------------------------------------

    
    void DOIFT_removeSubTree(int q_in,
			     sHeap *Q,
			     sImage32 *pred,
			     sImage32 *root,
			     sAdjRel *A){
      std::queue<int> path, frontier_path;
      int i, p, q;
      int px,py,qx,qy;
      
      path.push(q_in);
      
      while (!path.empty()) {
	p = path.front();
	path.pop();
	px = p%root->ncols;
	py = p/root->ncols;
	
	pred->data[p] = NIL;
	root->data[p] = p;
	
	if (Q->color[p] == GRAY)
	  gft::Heap::Delete_MinPolicy(Q, p);
	Q->color[p] = WHITE;
	Q->cost[p] = FLT_MAX;
	
	for (i = 1; i < A->n; i++) {
	  qx = px + A->dx[i];
	  qy = py + A->dy[i];
	  if(gft::Image32::IsValidPixel(root, qx, qy)){
	    q = qx + qy*root->ncols;
	    if (p == pred->data[q])
	      path.push(q);
	    else if(Q->cost[q] != FLT_MAX)
	      frontier_path.push(q);
	  }
	}
      }
      
      while (!frontier_path.empty()) {
	p = frontier_path.front();
	frontier_path.pop();
	if (Q->cost[p] != FLT_MAX) {
	  if (Q->color[p] != BLACK) {
	    gft::Heap::Update_MinPolicy(Q, p, Q->cost[p]);
	  } else {
	    Q->color[p] = WHITE;
	    gft::Heap::Insert_MinPolicy(Q, p);
	  }
	}
      }
    }




    void DOIFT_removeSubTree(int q_in,
			     sHeap *Q,
			     int *pred,
			     int *root,
			     sGraph *graph){
      std::queue<int> path, frontier_path;
      int i, p, q;
      
      path.push(q_in);
      
      while (!path.empty()) {
	p = path.front();
	path.pop();
	
	pred[p] = NIL;
	root[p] = p;
	
	if (Q->color[p] == GRAY)
	  gft::Heap::Delete_MinPolicy(Q, p);
	Q->color[p] = WHITE;
	Q->cost[p] = FLT_MAX;

	for(i = 0; i < graph->nodes[p].outdegree; i++){
	  q = graph->nodes[p].adjList[i];
	  if (p == pred[q])
	    path.push(q);
	  else if(Q->cost[q] != FLT_MAX)
	    frontier_path.push(q);
	}
      }
      
      while (!frontier_path.empty()) {
	p = frontier_path.front();
	frontier_path.pop();
	if (Q->cost[p] != FLT_MAX) {
	  if (Q->color[p] != BLACK) {
	    gft::Heap::Update_MinPolicy(Q, p, Q->cost[p]);
	  } else {
	    Q->color[p] = WHITE;
	    gft::Heap::Insert_MinPolicy(Q, p);
	  }
	}
      }
    }

    


    void DOIFT(sImageGraph *sg,
	       int *S, //new seeds only.
	       int *R, //new removal only.
	       sImage32 *label,
	       sHeap *Q,
	       sImage32 *pred,
	       sImage32 *root){
      //gft::sQueue *frontier = NULL;
      int *i_inv;
      int u_x,u_y,v_x,v_y;
      int p,q,i,j,rq,rp;
      float w,tmp;
      
      i_inv = gft::AdjRel::InverseIndexes(sg->A);

      if(R != NULL && R[0] > 0){
	DOIFT_TreeRemoval(R, sg, Q,
			  pred, root);

	/*
	while(!gft::Queue::IsEmpty(frontier)){
	  p = gft::Queue::Pop(frontier);
	  if (Q->color[p] != BLACK) {
	    gft::Heap::Update_MinPolicy(Q, p, Q->cost[p]);
	  } else {
	    Q->color[p] = WHITE;
	    gft::Heap::Insert_MinPolicy(Q, p);
	  }
	  Q->color[p] = GRAY;
	}
	*/
      }

      for(i = 1; i <= S[0]; i++){
	p = S[i];
	if(Q->cost[p] > 0.0){
	  pred->data[p] = NIL;
	  root->data[p] = p;
	  
	  if (Q->color[p] != BLACK) {
	    gft::Heap::Update_MinPolicy(Q, p, 0.0);
	  } else {
	    Q->color[p] = WHITE;
	    Q->cost[p] = 0.0;
	    gft::Heap::Insert_MinPolicy(Q, p);
	  }
	  Q->color[p] = GRAY;
	}
      }

      //----------
      while(!Heap::IsEmpty(Q)) {
	Heap::Remove_MinPolicy(Q, &p);
	u_x = p%label->ncols; //PixelX(label, p);
	u_y = p/label->ncols; //PixelY(label, p);

	for(i=1; i<sg->A->n; i++){
	  v_x = u_x + sg->A->dx[i];
	  v_y = u_y + sg->A->dy[i];
	  if(Image32::IsValidPixel(label,v_x,v_y)){
	    q = v_x + label->ncols*v_y;
	    
	    //if (Q->color[q] != BLACK){
	      
	      if(label->data[p]==0){
		j = i_inv[i];
		w = 2.0*(sg->n_link[q])[j];
	      }
	      else
		w = 2.0*(sg->n_link[p])[i] + 1.0;
	      
	      tmp = MAX(w, Q->cost[p]);
	      
	      rq = root->data[q];
	      rp = root->data[p];
	      
	      if(tmp < Q->cost[q]){
		label->data[q] = label->data[p];
		root->data[q] = root->data[p];
		pred->data[q] = p;
		
		if (Q->color[q] != BLACK) {
		  gft::Heap::Update_MinPolicy(Q, q, tmp);
		} else {
		  Q->color[q] = WHITE;
		  Q->cost[q] = tmp;
		  gft::Heap::Insert_MinPolicy(Q, q);
		}
	      }
	      else if(p == pred->data[q]){
		if(rp != rq || tmp > Q->cost[q]){
		  DOIFT_removeSubTree(q,
				      Q,
				      pred,
				      root,
				      sg->A);
		  break;
		}
	      }
	   //}
	  }
	}
      }
      //----------
      //gft::Queue::Destroy(&frontier);
      free(i_inv);
    }





    void DOIFT(sGraph *graph,
	       sGraph *transpose,
	       int *S, //new seeds only.
	       int *R, //new removal only.
	       int *label,
	       sHeap *Q,
	       int *pred,
	       int *root){
      int p,q,i,j,rq,rp;
      float w,tmp;
      sGraph *g;
      
      if(R != NULL && R[0] > 0){
	DOIFT_TreeRemoval(R, graph, transpose, Q,
			  pred, root);
      }

      for(i = 1; i <= S[0]; i++){
	p = S[i];
	if(Q->cost[p] > 0.0){
	  pred[p] = NIL;
	  root[p] = p;
	  
	  if (Q->color[p] != BLACK) {
	    gft::Heap::Update_MinPolicy(Q, p, 0.0);
	  } else {
	    Q->color[p] = WHITE;
	    Q->cost[p] = 0.0;
	    gft::Heap::Insert_MinPolicy(Q, p);
	  }
	  Q->color[p] = GRAY;
	}
      }

      //----------
      while(!Heap::IsEmpty(Q)) {
	Heap::Remove_MinPolicy(Q, &p);

	if(label[p]==0) g = transpose;
	else   	        g = graph;
	
	for(i = 0; i < g->nodes[p].outdegree; i++){
	  q = g->nodes[p].adjList[i];

	  w = g->nodes[p].Warcs[i];
	  if(label[p] == 0) w = 2.0*w;
	  else              w = 2.0*w + 1.0;
	    
	  tmp = MAX(w, Q->cost[p]);
	  
	  rq = root[q];
	  rp = root[p];
	    
	  if(tmp < Q->cost[q]){
	    label[q] = label[p];
	    root[q] = root[p];
	    pred[q] = p;
	    
	    if (Q->color[q] != BLACK) {
	      gft::Heap::Update_MinPolicy(Q, q, tmp);
	    } else {
	      Q->color[q] = WHITE;
	      Q->cost[q] = tmp;
	      gft::Heap::Insert_MinPolicy(Q, q);
	    }
	  }
	  else if(p == pred[q]){
	    if(rp != rq || tmp > Q->cost[q]){
	      DOIFT_removeSubTree(q,
				  Q,
				  pred,
				  root,
				  graph);
	      break;
	    }
	  }
	}
      }
    }
    

    
    

    void DOIFT_PaintSubTree(int r,
			    sImageGraph *sg,
			    sImage32 *pred,
			    sImage32 *img,
			    int val){
      gft::sQueue *Q;
      int i,p,q;
      int px,py,qx,qy;
      Q = gft::Queue::Create(pred->n);
      gft::Queue::Push(Q, r);
      while(!gft::Queue::IsEmpty(Q)){
	p = gft::Queue::Pop(Q);
	img->data[p] = val;

	px = p%pred->ncols;
	py = p/pred->ncols;
	for(i = 1; i < (sg->A)->n; i++){
	  qx = px + (sg->A)->dx[i];
	  qy = py + (sg->A)->dy[i];
	  //if(gft::Image32::IsValidPixel(pred, qx, qy)){
	  if(qx >= 0 && qx < pred->ncols &&
	     qy >= 0 && qy < pred->nrows){
	    q = qx + qy*pred->ncols;
	    if(p == pred->data[q])
	      gft::Queue::Push(Q, q);
	  }
	}
      }
      gft::Queue::Destroy(&Q);
    }




    void DOIFT_SubTrees2Trees(int *R,
			      sImageGraph *sg,
			      sImage32 *pred,
			      sImage32 *root){
      gft::sQueue *Q;
      int i,p,q,r;
      int px,py,qx,qy;
      Q = gft::Queue::Create(pred->n);

      for(i = 1; i <= R[0]; i++){
	r = R[i];
	pred->data[r] = NIL;
	root->data[r] = r;
	gft::Queue::Push(Q, r);
      }
      while(!gft::Queue::IsEmpty(Q)){
	p = gft::Queue::Pop(Q);
	r = root->data[p];

	px = p%pred->ncols;
	py = p/pred->ncols;
	for(i = 1; i < (sg->A)->n; i++){
	  qx = px + (sg->A)->dx[i];
	  qy = py + (sg->A)->dy[i];
	  //if(gft::Image32::IsValidPixel(pred, qx, qy)){
	  if(qx >= 0 && qx < pred->ncols &&
	     qy >= 0 && qy < pred->nrows){
	    q = qx + qy*pred->ncols;
	    if(p == pred->data[q]){
	      gft::Queue::Push(Q, q);
	      root->data[q] = r;
	    }
	  }
	}
      }
      gft::Queue::Destroy(&Q);
    }

    


    void DOIFT_TreeRemoval(int *R,
			   sImageGraph *sg,
			   sHeap *Q,
			   sImage32 *pred,
			   sImage32 *root){
      gft::sQueue *T;
      int i,k,r,p,q;
      int px,py,qx,qy;
      
      T = gft::Queue::Create(root->n);

      for(i = 1; i <= R[0]; i++){
	r = root->data[R[i]];
	if(Q->cost[r] != FLT_MAX){
	  if (Q->color[r] == GRAY)
	    gft::Heap::Delete_MinPolicy(Q, r);
	  Q->cost[r] = FLT_MAX;
	  Q->color[r] = WHITE;
	  pred->data[r] = NIL;
	  root->data[r] = r;
	  gft::Queue::Push(T, r);
	}
      }
      while (!gft::Queue::IsEmpty(T)) {
	p = gft::Queue::Pop(T);
	px = p%root->ncols;
	py = p/root->ncols;

	for (i = 1; i < (sg->A)->n; i++) {
	  qx = px + (sg->A)->dx[i];
	  qy = py + (sg->A)->dy[i];
	  if(qx >= 0 && qx < root->ncols &&
	     qy >= 0 && qy < root->nrows){
	    q = qx + qy*root->ncols;
	    if(p == pred->data[q]){
	      if (Q->color[q] == GRAY)
		gft::Heap::Delete_MinPolicy(Q, q);
	      Q->cost[q] = FLT_MAX;
	      Q->color[q] = WHITE;
	      pred->data[q] = NIL;
	      root->data[q] = q;
	      gft::Queue::Push(T, q);
	    }
	    else if(Q->cost[root->data[q]] != FLT_MAX){
	      if(Q->color[q] != GRAY){
		Q->color[q] = WHITE;
		gft::Heap::Insert_MinPolicy(Q, q);
	      }
	    }
	  }
	}
      }
      gft::Queue::Destroy(&T);
    }




    void DOIFT_TreeRemoval(int *R,
			   sGraph *graph,
			   sGraph *transpose,
			   sHeap *Q,
			   int *pred,
			   int *root){
      gft::sQueue *T;
      int i,k,r,p,q;
      
      T = gft::Queue::Create(graph->nnodes);

      for(i = 1; i <= R[0]; i++){
	r = root[R[i]];
	if(Q->cost[r] != FLT_MAX){
	  if (Q->color[r] == GRAY)
	    gft::Heap::Delete_MinPolicy(Q, r);
	  Q->cost[r] = FLT_MAX;
	  Q->color[r] = WHITE;
	  pred[r] = NIL;
	  root[r] = r;
	  gft::Queue::Push(T, r);
	}
      }
      while (!gft::Queue::IsEmpty(T)) {
	p = gft::Queue::Pop(T);

	for(i = 0; i < graph->nodes[p].outdegree; i++){
	  q = graph->nodes[p].adjList[i];

	  if(p == pred[q]){
	    if (Q->color[q] == GRAY)
	      gft::Heap::Delete_MinPolicy(Q, q);
	    Q->cost[q] = FLT_MAX;
	    Q->color[q] = WHITE;
	    pred[q] = NIL;
	    root[q] = q;
	    gft::Queue::Push(T, q);
	  }
	  else if(Q->cost[root[q]] != FLT_MAX){
	    if(Q->color[q] != GRAY){
	      Q->color[q] = WHITE;
	      gft::Heap::Insert_MinPolicy(Q, q);
	    }
	  }
	}
      }
      gft::Queue::Destroy(&T);
    }

    
    

    void DOIFT_TreeRemoval(int *R,
			   sImageGraph *sg,
			   sHeap32fiif_lex *Q,
			   sImage32 *pred,
			   sImage32 *root,
			   sImage32 *maxorder,
			   sImage32 *order){
      gft::sQueue *T;
      gft::sBMap *Fcolor;
      int i,k,r,p,q;
      int px,py,qx,qy;
      
      //Fcolor = gft::BMap::Create(root->n);
      T = gft::Queue::Create(root->n);

      for(i = 1; i <= R[0]; i++){
	r = root->data[R[i]];
	if(Q->cost1[r] != FLT_MAX){
	  if (Q->color[r] == GRAY)
	    gft::Heap32fiif_lex::Delete_MinPolicy(Q, r);	  
	  Q->cost1[r] = FLT_MAX;
	  Q->cost4[r] = FLT_MAX;
	  Q->color[r] = WHITE;
	  pred->data[r] = NIL;
	  root->data[r] = r;
	  order->data[r] = INT_MAX;
	  maxorder->data[r] = INT_MAX;	  
	  gft::Queue::Push(T, r);
	}
      }
      while (!gft::Queue::IsEmpty(T)) {
	p = gft::Queue::Pop(T);
	px = p%root->ncols;
	py = p/root->ncols;

	for (i = 1; i < (sg->A)->n; i++) {
	  qx = px + (sg->A)->dx[i];
	  qy = py + (sg->A)->dy[i];
	  //if(gft::Image32::IsValidPixel(root, qx, qy)){
	  if(qx >= 0 && qx < root->ncols &&
	     qy >= 0 && qy < root->nrows){
	    q = qx + qy*root->ncols;
	    if(p == pred->data[q]){
	      if (Q->color[q] == GRAY)
		gft::Heap32fiif_lex::Delete_MinPolicy(Q, q);
	      Q->cost1[q] = FLT_MAX;
	      Q->cost4[q] = FLT_MAX;
	      Q->color[q] = WHITE;
	      pred->data[q] = NIL;
	      root->data[q] = q;
	      order->data[q] = INT_MAX;
	      maxorder->data[q] = INT_MAX;
	      gft::Queue::Push(T, q);
	    }
	    else if(Q->cost1[root->data[q]] != FLT_MAX &&
		    order->data[q] != INT_MAX){
	      //if(gft::BMap::Get(Fcolor, q)==0){
	      if(Q->color[q] != GRAY){
		Q->color[q] = WHITE;
		gft::Heap32fiif_lex::Insert_MinPolicy(Q, q);
	      }
	      //gft::BMap::Set1(Fcolor, q);
	      //}
	    }
	  }
	}
      }
      /*     
      int c2 = INT_MAX;
      gft::Heap32fiif_lex::Get_MinPolicy(Q, &p);
      while(order->data[p] != INT_MAX){
	gft::Heap32fiif_lex::Remove_MinPolicy(Q, &p);
	gft::Queue::Push(T, p);
	if(c2 != Q->cost2[p]){
	  c2 = Q->cost2[p];
	  
	}
	else{
	  
	}
	gft::Heap32fiif_lex::Get_MinPolicy(Q, &p);
      }
      */
      //gft::BMap::Destroy(&Fcolor);
      gft::Queue::Destroy(&T);
    }
    


    /*
    // versao anterior sem break e sem RemovalByWrongOrder.
    void DOIFT(sImageGraph *sg,
	       int *S, //new seeds only.
	       int *R, //new removal only.
	       sImage32 *label,
	       sHeap32fiif_lex *Q,
	       sImage32 *pred,
	       sImage32 *root,
	       sImage32 *maxorder,
	       sImage32 *order,
	       int *iter){
      int *i_inv;
      int u_x,u_y,v_x,v_y;
      int p,q,i,j,rq,rp;
      float w,Wm;
      bool lex, permitted, eq;
      int Removal[2];
      int maxord;
      
      i_inv = gft::AdjRel::InverseIndexes(sg->A);

      if(R != NULL && R[0] > 0){
	DOIFT_TreeRemoval(R, sg, Q, pred, root,
			  maxorder, order, i_inv);
      }

      for(i = 1; i <= S[0]; i++){
	p = S[i];
	//if(Q->cost1[p] > 0.0){
	pred->data[p] = NIL;
	root->data[p] = p;
	if (Q->color[p] != BLACK) {
	  gft::Heap32fiif_lex::Update_MinPolicy(Q, p, 0.0, INT_MAX, INT_MAX, 0.0);
	} else {
	  Q->color[p] = WHITE;
	  Q->cost1[p] = 0.0;
	  Q->cost4[p] = 0.0;
	  maxorder->data[p] = INT_MAX;
	  order->data[p] = INT_MAX;
	  gft::Heap32fiif_lex::Insert_MinPolicy(Q, p);
	}
	Q->color[p] = GRAY;
	//}
      }

      //----------
      while(!Heap32fiif_lex::IsEmpty(Q)) {
	Heap32fiif_lex::Remove_MinPolicy(Q, &p);

	if(order->data[p] == INT_MAX){
	  (*iter)++;
	  order->data[p] = *iter;
	}

	if(Q->cost1[p] == Q->cost4[p]){
	  maxorder->data[p] = order->data[p];
	}
	
	//-------------------------------------------------------
	u_x = p%label->ncols; //PixelX(label, p);
	u_y = p/label->ncols; //PixelY(label, p);
	
	for(i=1; i<sg->A->n; i++){
	  v_x = u_x + sg->A->dx[i];
	  v_y = u_y + sg->A->dy[i];
	  //if(Image32::IsValidPixel(label,v_x,v_y)){
	  if(v_x >= 0 && v_x < label->ncols &&
	     v_y >= 0 && v_y < label->nrows){
	    q = v_x + label->ncols*v_y;
	    
	    if(label->data[p]==0){
	      j = i_inv[i];
	      w = (sg->n_link[q])[j];
	    }
	    else
	      w = (sg->n_link[p])[i];
	    
	    Wm = MAX(w, Q->cost1[p]);
	    
	    rq = root->data[q];
	    rp = root->data[p];
	    
	    permitted = (order->data[q] == INT_MAX);
	    
	    eq = (Q->cost1[q] == Q->cost4[q] && Wm == w && Wm == Q->cost1[q]);

	    if(Wm != w)
	      maxord = maxorder->data[p];
	    else
	      maxord = INT_MAX;

	    lex = (Wm == Q->cost1[q] && maxord == maxorder->data[q] && w < Q->cost4[q]);
	    
	    if((Wm < Q->cost1[q] || (Wm == Q->cost1[q] && maxord < maxorder->data[q]))||
	       (lex && permitted)){
	      label->data[q] = label->data[p];
	      root->data[q] = root->data[p];
	      pred->data[q] = p;
	      
	      if (Q->color[q] != BLACK) {
		gft::Heap32fiif_lex::Update_MinPolicy(Q, q, Wm, maxord, INT_MAX, w);
	      } else {
		Q->color[q] = WHITE;
		Q->cost1[q] = Wm;
		Q->cost4[q] = w;
		maxorder->data[q] = maxord;
		order->data[q] = INT_MAX;
		gft::Heap32fiif_lex::Insert_MinPolicy(Q, q);
	      }
	    }
	    else if(eq && pred->data[q] != NIL && Q->cost1[p] < Q->cost1[pred->data[q]]){
	      pred->data[q] = NIL;
	      DOIFT_PaintSubTree(q, sg, pred,
				 root, q);
	      Removal[0] = 1; Removal[1] = q;
	      DOIFT_TreeRemoval(Removal, sg, Q, pred, root,
				maxorder, order, i_inv);
	    }
	    else if(p == pred->data[q]){
	      if(rp != rq || (order->data[p] > order->data[q])){
		pred->data[q] = NIL;
		DOIFT_PaintSubTree(q, sg, pred,
				   root, q);
		Removal[0] = 1; Removal[1] = q;
		DOIFT_TreeRemoval(Removal, sg, Q, pred, root,
				  maxorder, order, i_inv);
	      }
	    }
	  }
	}
	//-------------------------------------------------------
      }
      //----------
      free(i_inv);
    }
    */   


    
    //-----------------------------------------------------




    void DOIFT_MarkForRemovalByWrongOrder(int *R, sImageGraph *sg,
					  sImage32 *pred, sImage32 *maxorder,
					  sImage32 *order, float *cost,
					  float cst, int r){
      gft::sQueue *Q;
      int i,p,q,ord;
      int px,py,qx,qy;
      Q = gft::Queue::Create(pred->n);

      ord = order->data[r];
      while(maxorder->data[r] != order->data[r])
	r = pred->data[r];

      gft::Queue::Push(Q, r);
      while(!gft::Queue::IsEmpty(Q)){
	p = gft::Queue::Pop(Q);
	px = p%pred->ncols;
	py = p/pred->ncols;
	for(i = 1; i < (sg->A)->n; i++){
	  qx = px + (sg->A)->dx[i];
	  qy = py + (sg->A)->dy[i];
	  if(qx >= 0 && qx < pred->ncols &&
	     qy >= 0 && qy < pred->nrows){
	    q = qx + qy*pred->ncols;
	    if(p == pred->data[q]){
	      gft::Queue::Push(Q, q);
	      if(cost[q] > cst &&
		 order->data[q] > ord && order->data[q] != INT_MAX){
		R[0]++;
		R[R[0]] = q;
	      }
	    }
	  }
	}
      }
      gft::Queue::Destroy(&Q);
    }
    
    

    void DOIFT_(sImageGraph *sg,
		int *S, //new seeds only.
		int *R, //new removal only.
		sImage32 *label,
		sHeap32fiif_lex *Q,
		sImage32 *pred,
		sImage32 *root,
		sImage32 *maxorder,
		sImage32 *order,
		int *iter){
      int *i_inv;
      int u_x,u_y,v_x,v_y;
      int p,q,i,j,rq,rp,t,tm;
      float w,Wm;
      bool lex, permitted, eq;
      int Removal[2];
      int maxord;
      int Front_maxorder, Front_order, Front_p;
      float Front_cost;
      int *Front_Removal;
      bool Front,b;

      //bool debug = false;
      
      Front_Removal = (int *)malloc(sizeof(int)*(label->n+1));
      Front_Removal[0] = 0;
      Front_cost = FLT_MAX;
      Front_maxorder = INT_MAX;
      Front_order = INT_MAX;
      Front_p = NIL;
      i_inv = gft::AdjRel::InverseIndexes(sg->A);

      if(R != NULL && R[0] > 0){
	DOIFT_TreeRemoval(R, sg, Q, pred, root,
			  maxorder, order);
      }

      for(i = 1; i <= S[0]; i++){
	p = S[i];
	//if(Q->cost1[p] > 0.0){
	pred->data[p] = NIL;
	root->data[p] = p;
	if (Q->color[p] != BLACK) {
	  gft::Heap32fiif_lex::Update_MinPolicy(Q, p, 0.0, INT_MAX, INT_MAX, 0.0);
	} else {
	  Q->color[p] = WHITE;
	  Q->cost1[p] = 0.0;
	  Q->cost4[p] = 0.0;
	  maxorder->data[p] = INT_MAX;
	  order->data[p] = INT_MAX;
	  gft::Heap32fiif_lex::Insert_MinPolicy(Q, p);
	}
	Q->color[p] = GRAY;
	//}
      }

      //int k_tmp = 0;
      //----------
      while(!Heap32fiif_lex::IsEmpty(Q)) {
	Heap32fiif_lex::Remove_MinPolicy(Q, &p);

	Front = false;
	if(order->data[p] == INT_MAX){
	  (*iter)++;
	  order->data[p] = *iter;

	  Front_maxorder = INT_MAX;
	  Front_order = INT_MAX;
	  Front_cost = FLT_MAX;
	}
	else{ //frontier pixel:
	  Front = true;
	  if(Front_maxorder == maxorder->data[p]){
	    /*
	    t = pred->data[p];
	    tm = p;
	    while(Front_order < order->data[t]){
	      if(Q->cost4[t] > Q->cost4[tm])
		tm = t;
	      t = pred->data[t];
	    }

	    if(Q->cost4[tm] > Front_cost){
	      Front_Removal[0]++;
	      Front_Removal[Front_Removal[0]] = tm;
	    }
	    */
	  }
	  else{
	    Front_maxorder = maxorder->data[p];
	    Front_order = order->data[p];
	    Front_cost = FLT_MAX;
	    Front_p = p;
	  }
	}

	//----------------------------------
	//for(i = 1 ; i <= Front_Removal[0]; i++)
	//  printf(" %d", Front_Removal[i]);
	//printf("\n");
	//----------------------------------	
	
	if(Q->cost1[p] == Q->cost4[p]){
	  maxorder->data[p] = order->data[p];
	}
	
	//printf("Q: (%d,%d), iter: %d\n",p%pred->ncols, p/pred->ncols, *iter);
	
	//-------------------------------------------------------
	u_x = p%label->ncols; //PixelX(label, p);
	u_y = p/label->ncols; //PixelY(label, p);

	/*
	debug = false;
	if((u_x == 22 && u_y == 30)||
	   (u_x == 21 && u_y == 30)||
	   (u_x == 21 && u_y == 31)){
	  printf("p: %d, x: %d, y: %d, c1: %f, c2: %d, c3: %d, c4: %f, pred: %d\n", p, u_x, u_y, Q->cost1[p], maxorder->data[p], order->data[p], Q->cost4[p], pred->data[p]);
	  printf("\t1042=> c1: %f, c2: %d, c3: %d, c4: %f, pred: %d\n", Q->cost1[1042], maxorder->data[1042], order->data[1042], Q->cost4[1042], pred->data[1042]);
	  debug = true;
	}
	*/
	
	for(i=1; i<sg->A->n; i++){
	  v_x = u_x + sg->A->dx[i];
	  v_y = u_y + sg->A->dy[i];
	  //if(Image32::IsValidPixel(label,v_x,v_y)){
	  if(v_x >= 0 && v_x < label->ncols &&
	     v_y >= 0 && v_y < label->nrows){
	    q = v_x + label->ncols*v_y;
	    
	    if(label->data[p]==0){
	      j = i_inv[i];
	      w = (sg->n_link[q])[j];
	    }
	    else
	      w = (sg->n_link[p])[i];
	    
	    Wm = MAX(w, Q->cost1[p]);
	    
	    rq = root->data[q];
	    rp = root->data[p];
	    
	    permitted = (order->data[q] == INT_MAX);
	    
	    eq = (Q->cost1[q] == Q->cost4[q] && Wm == w && Wm == Q->cost1[q]);

	    if(Wm != w)
	      maxord = maxorder->data[p];
	    else
	      maxord = INT_MAX;

	    lex = (Wm == Q->cost1[q] && maxord == maxorder->data[q] && w < Q->cost4[q]);


	    /*
	    if(debug == true && q == 1042){
	      printf("\t\tpermitted: %d, lex: %d, eq: %d, maxord: %d, Wm: %f, w: %f\n",permitted,lex,eq,maxord,Wm,w);
	    }
	    */
	    
	    if((Wm < Q->cost1[q] || (Wm == Q->cost1[q] && maxord < maxorder->data[q]))||
	       (lex && permitted)){
	      /*
	      if(debug == true && q == 1042){
		printf("Entrou\n");
	      }
	      */
	      
	      label->data[q] = label->data[p];
	      root->data[q] = root->data[p];
	      pred->data[q] = p;

	      if(Front) Front_cost = MIN(Front_cost, w);
	      
	      if (Q->color[q] != BLACK) {
		gft::Heap32fiif_lex::Update_MinPolicy(Q, q, Wm, maxord, INT_MAX, w);
	      } else {
		Q->color[q] = WHITE;
		Q->cost1[q] = Wm;
		Q->cost4[q] = w;
		maxorder->data[q] = maxord;
		order->data[q] = INT_MAX;
		gft::Heap32fiif_lex::Insert_MinPolicy(Q, q);
	      }
	    }
	    else if(eq && pred->data[q] != NIL && Q->cost1[p] < Q->cost1[pred->data[q]]){
	      pred->data[q] = NIL;
	      DOIFT_PaintSubTree(q, sg, pred,
				 root, q);
	      Removal[0] = 1; Removal[1] = q;
	      DOIFT_TreeRemoval(Removal, sg, Q, pred, root,
				maxorder, order);
	      break;
	    }
	    else if(p == pred->data[q]){
	      if(rp != rq || (order->data[p] > order->data[q])){
		//printf("p: (%d,%d)\n", p%pred->ncols, p/pred->ncols);
		//printf("q: (%d,%d)\n", q%pred->ncols, q/pred->ncols);	      
		
		pred->data[q] = NIL;
		DOIFT_PaintSubTree(q, sg, pred,
				   root, q);
		Removal[0] = 1; Removal[1] = q;
		DOIFT_TreeRemoval(Removal, sg, Q, pred, root,
				  maxorder, order);
		break;
	      }
	    }
	  }
	}

	/*
	if(debug == true){
	  printf("p: %d, x: %d, y: %d, c1: %f, c2: %d, c3: %d, c4: %f, pred: %d\n", p, u_x, u_y, Q->cost1[p], maxorder->data[p], order->data[p], Q->cost4[p], pred->data[p]);
	  printf("\t1042=> c1: %f, c2: %d, c3: %d, c4: %f, pred: %d\n", Q->cost1[1042], maxorder->data[1042], order->data[1042], Q->cost4[1042], pred->data[1042]);
	}
	*/
	//-------------------------------------------------------
	/*
	if(Front_Removal[0] > 0){
	  b = false;
	  if(Heap32fiif_lex::IsEmpty(Q))
	    b = true;
	  else{
	    Heap32fiif_lex::Get_MinPolicy(Q, &p);
	    if(order->data[p] == INT_MAX || Front_maxorder != maxorder->data[p])
	      b = true;
	  }
	  
	  if(b){
	    DOIFT_SubTrees2Trees(Front_Removal, sg, pred, root);
	    DOIFT_TreeRemoval(Front_Removal, sg, Q, pred, root,
			      maxorder, order, i_inv);
	    Front_Removal[0] = 0;
	  }
	}
	*/
	b = false;
	if(Heap32fiif_lex::IsEmpty(Q))
	  b = true;
	else{
	  Heap32fiif_lex::Get_MinPolicy(Q, &p);
	  if(order->data[p] == INT_MAX || Front_maxorder != maxorder->data[p])
	    b = true;
	}

	if(b){
	  DOIFT_MarkForRemovalByWrongOrder(Front_Removal, sg,
					   pred, maxorder,
					   order, Q->cost4,
					   Front_cost, Front_p);
	  if(Front_Removal[0] > 0){
	    /*
	    k_tmp++;
	    printf("WrongOrder p: %d\n", Front_p);
	    printf("\tmaxorder: %d, order: %d, cost: %f\n",
		   Front_maxorder,
		   Front_order,
		   Front_cost);
	    printf("R: ");
	    for(i = 1; i <= Front_Removal[0]; i++)
	      printf("%d ", Front_Removal[i]);
	    printf("\n");
	    */
	    
	    DOIFT_SubTrees2Trees(Front_Removal, sg, pred, root);
	    DOIFT_TreeRemoval(Front_Removal, sg, Q, pred, root,
			      maxorder, order);
	    Front_Removal[0] = 0;
	  }
	}
	//-------------------------------------------------------
      }
      //----------
      free(i_inv);
      free(Front_Removal);
    }


    //----------------------------------------------------
    /*    
    void DOIFT_TreeRemoval_2(int *R,
			     sImageGraph *sg,
			     sHeap32fi_lex *Q,
			     sImage32 *pred,
			     sImage32 *maxorder,
			     int *i_inv,
			     int *T){
      gft::sBMap *Fcolor;
      int i,j,k,r,p,q;
      int px,py,qx,qy;

      T[0] = 0;
      for(j = 1; j <= R[0]; j++){
	r = R[j];
	T[0]++;
	T[T[0]] = r;
      }

      //---------------------------------------
      while(T[0] > 0){
	p = T[T[0]];
	T[0]--;
	
	if (Q->color[p] == GRAY)
	  gft::Heap32fi_lex::Delete_MinPolicy(Q, p);
	Q->cost1[p] = FLT_MAX;
	maxorder->data[p] = INT_MAX;
	Q->color[p] = WHITE;
	
	px = p%pred->ncols;
	py = p/pred->ncols;
	for(i = 1; i < (sg->A)->n; i++){
	  qx = px + (sg->A)->dx[i];
	  qy = py + (sg->A)->dy[i];
	  //if(gft::Image32::IsValidPixel(pred, qx, qy)){
	  if(qx >= 0 && qx < pred->ncols &&
	     qy >= 0 && qy < pred->nrows){
	    q = qx + qy*pred->ncols;
	    if(p == pred->data[q]){
	      T[0]++;
	      T[T[0]] = q;
	    }
	  }
	}
      }
      //---------------------------------------	

      for(j = 1; j <= R[0]; j++){
	r = R[j];
	T[0]++;
	T[T[0]] = r;
      }
      while(T[0] > 0){
	p = T[T[0]];
	T[0]--;

	px = p%pred->ncols;
	py = p/pred->ncols;
	pred->data[p] = NIL;
	
	for (i = 1; i < (sg->A)->n; i++) {
	  qx = px + (sg->A)->dx[i];
	  qy = py + (sg->A)->dy[i];
	  //if(gft::Image32::IsValidPixel(pred, qx, qy)){
	  if(qx >= 0 && qx < pred->ncols &&
	     qy >= 0 && qy < pred->nrows){
	    q = qx + qy*pred->ncols;
	    if(p == pred->data[q]){
	      T[0]++;
	      T[T[0]] = q;
	    }
	    else if(Q->cost1[q] != FLT_MAX){
	      if(Q->color[q] != GRAY){
		Q->color[q] = WHITE;
		gft::Heap32fi_lex::Insert_MinPolicy(Q, q);
	      }
	    }
	  }
	}
      }
    }
    */


    /*
    void DOIFT_TreeRemoval_3(int *R,
			     sImageGraph *sg,
			     sHeap32fi_lex *Q,
			     sImage32 *pred,
			     sImage32 *maxorder,
			     int *i_inv,
			     int *T){
      gft::sBMap *Fcolor;
      int i,j,k,r,p,q;
      int px,py,qx,qy;

      T[0] = 0;
      for(j = 1; j <= R[0]; j++){
	r = R[j];
	T[0]++;
	T[T[0]] = r;
      }

      //---------------------------------------
      while(T[0] > 0){
	p = T[T[0]];
	T[0]--;
	
	if (Q->color[p] == GRAY)
	  gft::Heap32fi_lex::Delete_MinPolicy(Q, p);
	Q->cost1[p] = FLT_MAX;
	maxorder->data[p] = INT_MAX;
	Q->color[p] = WHITE;
	pred->data[p] = NIL;
	
	px = p%pred->ncols;
	py = p/pred->ncols;
	for(i = 1; i < (sg->A)->n; i++){
	  qx = px + (sg->A)->dx[i];
	  qy = py + (sg->A)->dy[i];
	  //if(gft::Image32::IsValidPixel(pred, qx, qy)){
	  if(qx >= 0 && qx < pred->ncols &&
	     qy >= 0 && qy < pred->nrows){
	    q = qx + qy*pred->ncols;
	    if(p == pred->data[q]){
	      T[0]++;
	      T[T[0]] = q;
	    }
	    else if(Q->cost1[q] != FLT_MAX){
	      if(Q->color[q] != GRAY){
		Q->color[q] = WHITE;
		gft::Heap32fi_lex::Insert_MinPolicy(Q, q);
	      }
	    }
	  }
	}
      }
      //---------------------------------------	
    }
    */


    /*
    void DOIFT_TreeRemoval_4(int *R,
			     sImageGraph *sg,
			     sHeap32fi_lex *Q,
			     sImage32 *pred,
			     sImage32 *maxorder,
			     int *i_inv,
			     int *T){
      gft::sBMap *Fcolor;
      int i,j,k,r,p,q;
      int px,py,qx,qy;

      T[0] = 0;
      for(j = 1; j <= R[0]; j++){
	r = R[j];
	T[0]++;
	T[T[0]] = r;

	if (Q->color[r] == GRAY)
	  gft::Heap32fi_lex::Delete_MinPolicy(Q, r);
	Q->cost1[r] = FLT_MAX;
	maxorder->data[r] = INT_MAX;
	Q->color[r] = WHITE;
	pred->data[r] = NIL;
      }

      //---------------------------------------
      while(T[0] > 0){
	p = T[T[0]];
	T[0]--;
	
	px = p%pred->ncols;
	py = p/pred->ncols;
	for(i = 1; i < (sg->A)->n; i++){
	  qx = px + (sg->A)->dx[i];
	  qy = py + (sg->A)->dy[i];
	  //if(gft::Image32::IsValidPixel(pred, qx, qy)){
	  if(qx >= 0 && qx < pred->ncols &&
	     qy >= 0 && qy < pred->nrows){
	    q = qx + qy*pred->ncols;
	    if(p == pred->data[q]){
	      T[0]++;
	      T[T[0]] = q;

	      if (Q->color[q] == GRAY)
		gft::Heap32fi_lex::Delete_MinPolicy(Q, q);
	      Q->cost1[q] = FLT_MAX;
	      maxorder->data[q] = INT_MAX;
	      Q->color[q] = WHITE;
	      pred->data[q] = NIL;
	    }
	    else if(Q->cost1[q] != FLT_MAX){
	      if(Q->color[q] != GRAY){
		Q->color[q] = WHITE;
		gft::Heap32fi_lex::Insert_MinPolicy(Q, q);
	      }
	    }
	  }
	}
      }
      //---------------------------------------	
    }
    */


    
    void DOIFT_TreeRemoval_5(int *R,
			     sImageGraph *sg,
			     sHeap32fi_lex *Q,
			     sImage32 *pred,
			     sImage32 *maxorder,
			     int *T){
      gft::sBMap *Fcolor;
      int i,j,k,r,p,q;
      int px,py,qx,qy;
      int top = pred->n+1;
      
      T[0] = 0;
      for(j = 1; j <= R[0]; j++){
	r = R[j];
	T[0]++;
	T[T[0]] = r;

	if (Q->color[r] == GRAY)
	  gft::Heap32fi_lex::Delete_MinPolicy(Q, r);
	Q->cost1[r] = FLT_MAX;
	maxorder->data[r] = INT_MAX;
	Q->color[r] = WHITE;
	pred->data[r] = NIL;
      }

      //---------------------------------------
      while(T[0] > 0){
	p = T[T[0]];
	T[0]--;
	
	px = p%pred->ncols;
	py = p/pred->ncols;
	for(i = 1; i < (sg->A)->n; i++){
	  qx = px + (sg->A)->dx[i];
	  qy = py + (sg->A)->dy[i];
	  //if(gft::Image32::IsValidPixel(pred, qx, qy)){
	  if(qx >= 0 && qx < pred->ncols &&
	     qy >= 0 && qy < pred->nrows){
	    q = qx + qy*pred->ncols;
	    if(p == pred->data[q]){
	      T[0]++;
	      T[T[0]] = q;

	      if (Q->color[q] == GRAY)
		gft::Heap32fi_lex::Delete_MinPolicy(Q, q);
	      Q->cost1[q] = FLT_MAX;
	      maxorder->data[q] = INT_MAX;
	      Q->color[q] = WHITE;
	      pred->data[q] = NIL;
	    }
	    else if(Q->cost1[q] != FLT_MAX &&
		    Q->color[q] != GRAY){
	      top--;
	      T[top] = q;
	    }
	  }
	}
      }
      //---------------------------------------	
      //printf("per: %f\n",(float)(pred->n+1-top)/((float)(pred->n)));
      while(top < pred->n+1){
	q = T[top];
	top++;
	if(Q->cost1[q] != FLT_MAX && Q->color[q] != GRAY){
	  Q->color[q] = WHITE;
	  gft::Heap32fi_lex::Insert_MinPolicy(Q, q);	  
	}
      }
    }





    void DOIFT_TreeRemoval_5(int *R,
			     sGraph *graph,
			     sHeap32fi_lex *Q,
			     int *pred,
			     int *maxorder,
			     int *T){
      gft::sBMap *Fcolor;
      int i,j,k,r,p,q;
      int top = graph->nnodes+1;
      
      T[0] = 0;
      for(j = 1; j <= R[0]; j++){
	r = R[j];
	T[0]++;
	T[T[0]] = r;

	if (Q->color[r] == GRAY)
	  gft::Heap32fi_lex::Delete_MinPolicy(Q, r);
	Q->cost1[r] = FLT_MAX;
	maxorder[r] = INT_MAX;
	Q->color[r] = WHITE;
	pred[r] = NIL;
      }

      //---------------------------------------
      while(T[0] > 0){
	p = T[T[0]];
	T[0]--;
	
	for(i = 0; i < graph->nodes[p].outdegree; i++){
	  q = graph->nodes[p].adjList[i];

	  if(p == pred[q]){
	    T[0]++;
	    T[T[0]] = q;

	    if (Q->color[q] == GRAY)
	      gft::Heap32fi_lex::Delete_MinPolicy(Q, q);
	    Q->cost1[q] = FLT_MAX;
	    maxorder[q] = INT_MAX;
	    Q->color[q] = WHITE;
	    pred[q] = NIL;
	  }
	  else if(Q->cost1[q] != FLT_MAX &&
		  Q->color[q] != GRAY){
	    top--;
	    T[top] = q;
	  }
	}
      }
      //---------------------------------------	
      while(top < graph->nnodes+1){
	q = T[top];
	top++;
	if(Q->cost1[q] != FLT_MAX && Q->color[q] != GRAY){
	  Q->color[q] = WHITE;
	  gft::Heap32fi_lex::Insert_MinPolicy(Q, q);	  
	}
      }
    }





    //Errado, apagar depois
    void DOIFT_TreeRemoval__5(int *R,
			      sImageGraph *sg,
			      sHeap32fif_lex *Q,
			      sImage32 *pred,
			      sImage32 *maxorder,
			      int *T){
      gft::sBMap *Fcolor;
      int i,j,k,r,p,q;
      int px,py,qx,qy;
      int top = pred->n+1;
      
      T[0] = 0;
      for(j = 1; j <= R[0]; j++){
	r = R[j];
	T[0]++;
	T[T[0]] = r;

	if (Q->color[r] == GRAY)
	  gft::Heap32fif_lex::Delete_MinPolicy(Q, r);
	Q->cost1[r] = FLT_MAX;
	maxorder->data[r] = INT_MAX;
	Q->cost3[r] = FLT_MAX;
	Q->color[r] = WHITE;
	pred->data[r] = NIL;
      }

      //---------------------------------------
      while(T[0] > 0){
	p = T[T[0]];
	T[0]--;
	
	px = p%pred->ncols;
	py = p/pred->ncols;
	for(i = 1; i < (sg->A)->n; i++){
	  qx = px + (sg->A)->dx[i];
	  qy = py + (sg->A)->dy[i];
	  //if(gft::Image32::IsValidPixel(pred, qx, qy)){
	  if(qx >= 0 && qx < pred->ncols &&
	     qy >= 0 && qy < pred->nrows){
	    q = qx + qy*pred->ncols;
	    if(p == pred->data[q]){
	      T[0]++;
	      T[T[0]] = q;

	      if (Q->color[q] == GRAY)
		gft::Heap32fif_lex::Delete_MinPolicy(Q, q);
	      Q->cost1[q] = FLT_MAX;
	      maxorder->data[q] = INT_MAX;
	      Q->cost3[q] = FLT_MAX;
	      Q->color[q] = WHITE;
	      pred->data[q] = NIL;
	    }
	    else if(Q->cost1[q] != FLT_MAX &&
		    Q->color[q] != GRAY){
	      top--;
	      T[top] = q;
	    }
	  }
	}
      }
      //---------------------------------------	
      //printf("per: %f\n",(float)(pred->n+1-top)/((float)(pred->n)));
      while(top < pred->n+1){
	q = T[top];
	top++;
	if(Q->cost1[q] != FLT_MAX && Q->color[q] != GRAY){
	  Q->color[q] = WHITE;
	  gft::Heap32fif_lex::Insert_MinPolicy(Q, q);	  
	}
      }
      /*
      while(top < pred->n+1){
	q = T[top];
	top++;
	if(Q->cost1[q] != FLT_MAX){

	  while(Q->cost1[q] != Q->cost3[q])
	    q = pred->data[q];

	  if(Q->color[q] != GRAY){
	    Q->color[q] = WHITE;
	    maxorder->data[q]++;
	    gft::Heap32fif_lex::Insert_MinPolicy(Q, q);
	  }

	  //Reset subtrees of node q:
	  T[0] = 1;
	  T[1] = q;
	  while(T[0] > 0){
	    p = T[T[0]];
	    T[0]--;
	
	    px = p%pred->ncols;
	    py = p/pred->ncols;
	    for(i = 1; i < (sg->A)->n; i++){
	      qx = px + (sg->A)->dx[i];
	      qy = py + (sg->A)->dy[i];
	      //if(gft::Image32::IsValidPixel(pred, qx, qy)){
	      if(qx >= 0 && qx < pred->ncols &&
		 qy >= 0 && qy < pred->nrows){
		q = qx + qy*pred->ncols;
		if(p == pred->data[q]){
		  T[0]++;
		  T[T[0]] = q;
		  
		  if (Q->color[q] == GRAY)
		    gft::Heap32fif_lex::Delete_MinPolicy(Q, q);
		  Q->cost1[q] = FLT_MAX;
		  maxorder->data[q] = INT_MAX;
		  Q->cost3[q] = FLT_MAX;
		  Q->color[q] = WHITE;
		  pred->data[q] = NIL;
		}
	      }
	    }
	  }
	  }
	  }
      */
   }

    
    

    void DOIFT_2(sImageGraph *sg,
		 int *S, //new seeds only.
		 int *R, //new removal only.
		 sImage32 *label,
		 sHeap32fi_lex *Q,
		 sImage32 *pred,
		 sImage32 *maxorder,
		 int *iter){
      int *i_inv;
      int u_x,u_y,v_x,v_y;
      int p,q,i,j,tq;
      float w,Wm;
      int Removal[2];
      int maxord;
      int *T = NULL;
      
      i_inv = gft::AdjRel::InverseIndexes(sg->A);

      if(R != NULL && R[0] > 0){
	if(T == NULL)
	  T = (int *)malloc((pred->n+1)*sizeof(int));
	DOIFT_TreeRemoval_5(R, sg, Q, pred,
			    maxorder, T);
      }

      for(i = 1; i <= S[0]; i++){
	p = S[i];
	pred->data[p] = NIL;
	if (Q->color[p] != BLACK) {
	  gft::Heap32fi_lex::Update_MinPolicy(Q, p, 0.0, INT_MAX);
	} else {
	  Q->color[p] = WHITE;
	  Q->cost1[p] = 0.0;
	  maxorder->data[p] = INT_MAX;
	  gft::Heap32fi_lex::Insert_MinPolicy(Q, p);
	}
	Q->color[p] = GRAY;
      }

      //----------
      while(!Heap32fi_lex::IsEmpty(Q)) {
	Heap32fi_lex::Remove_MinPolicy(Q, &p);

	//q = pred->data[p];
	if(maxorder->data[p] == INT_MAX){
	   //&& (q == NIL || (Q->cost1[p] != Q->cost1[q]))){
	  (*iter)++;
	  maxorder->data[p] = *iter;
	}
	
	//-------------------------------------------------------
	u_x = p%label->ncols; //PixelX(label, p);
	u_y = p/label->ncols; //PixelY(label, p);
	
	for(i=1; i<sg->A->n; i++){
	  v_x = u_x + sg->A->dx[i];
	  v_y = u_y + sg->A->dy[i];
	  //if(Image32::IsValidPixel(label,v_x,v_y)){
	  if(v_x >= 0 && v_x < label->ncols &&
	     v_y >= 0 && v_y < label->nrows){	  
	    q = v_x + label->ncols*v_y;
	    
	    if(label->data[p]==0){
	      j = i_inv[i];
	      w = (sg->n_link[q])[j];
	    }
	    else
	      w = (sg->n_link[p])[i];
	    
	    Wm = MAX(w, Q->cost1[p]);
	    
	    if(Wm != w)
	      maxord = maxorder->data[p];
	    else
	      maxord = INT_MAX;
	    
	    if((Wm < Q->cost1[q] || (Wm == Q->cost1[q] && maxord < maxorder->data[q]))){
	      label->data[q] = label->data[p];
	      pred->data[q] = p;
	      
	      if (Q->color[q] != BLACK) {
		gft::Heap32fi_lex::Update_MinPolicy(Q, q, Wm, maxord);
	      } else {
		Q->color[q] = WHITE;
		Q->cost1[q] = Wm;
		maxorder->data[q] = maxord;
		gft::Heap32fi_lex::Insert_MinPolicy(Q, q);
	      }
	    }
	    else if(p == pred->data[q]){
	      if((label->data[q] != label->data[p])||
		 (Wm != Q->cost1[q]) ||
		 (maxorder->data[p] > maxorder->data[q])){
		if(T == NULL)
		  T = (int *)malloc((pred->n+1)*sizeof(int));

		Removal[0] = 1; Removal[1] = q;
		DOIFT_TreeRemoval_5(Removal, sg, Q, pred,
				    maxorder, T);
	      }
	    }
	  }
	}
	//-------------------------------------------------------
      }
      //----------
      free(i_inv);
      if(T != NULL)
	free(T);
    }




    void DOIFT__2(sImageGraph *sg,
		  int *S, //new seeds only.
		  int *R, //new removal only.
		  sImage32 *label,
		  sHeap32fi_lex *Q,
		  sImage32 *pred,
		  sImage32 *maxorder,
		  int *iter){
      int *i_inv;
      int u_x,u_y,v_x,v_y;
      int p,q,i,j,tq,cm;
      float w,Wm;
      int Removal[2];
      int maxord;
      int *T = NULL;
      
      i_inv = gft::AdjRel::InverseIndexes(sg->A);

      if(R != NULL && R[0] > 0){
	if(T == NULL)
	  T = (int *)malloc((pred->n+1)*sizeof(int));
	DOIFT_TreeRemoval_5(R, sg, Q, pred,
			    maxorder, T);
      }

      for(i = 1; i <= S[0]; i++){
	p = S[i];
	pred->data[p] = NIL;
	if(label->data[p]==0)
	  cm = 1; //1;
	else
	  cm = 2; //2;
	  
	if (Q->color[p] != BLACK) {
	  gft::Heap32fi_lex::Update_MinPolicy(Q, p, 0.0, cm); //INT_MAX);
	} else {
	  Q->color[p] = WHITE;
	  Q->cost1[p] = 0.0;
	  maxorder->data[p] = cm; //INT_MAX;
	  gft::Heap32fi_lex::Insert_MinPolicy(Q, p);
	}
	Q->color[p] = GRAY;
      }

      //----------
      while(!Heap32fi_lex::IsEmpty(Q)) {
	Heap32fi_lex::Remove_MinPolicy(Q, &p);

	/*
	if(maxorder->data[p] == INT_MAX){
	  (*iter)++;
	  maxorder->data[p] = *iter;
	}
	*/
	//-------------------------------------------------------
	u_x = p%label->ncols; //PixelX(label, p);
	u_y = p/label->ncols; //PixelY(label, p);
	
	for(i=1; i<sg->A->n; i++){
	  v_x = u_x + sg->A->dx[i];
	  v_y = u_y + sg->A->dy[i];
	  //if(Image32::IsValidPixel(label,v_x,v_y)){
	  if(v_x >= 0 && v_x < label->ncols &&
	     v_y >= 0 && v_y < label->nrows){	  
	    q = v_x + label->ncols*v_y;
	    
	    if(label->data[p]==0){
	      j = i_inv[i];
	      w = (sg->n_link[q])[j];
	      cm = 1; //1;
	    }
	    else{
	      w = (sg->n_link[p])[i];
	      cm = 2; //2;
	    }
	      
	    Wm = MAX(w, Q->cost1[p]);
	    
	    if(Wm != w)
	      maxord = maxorder->data[p];
	    else if( Q->cost1[p] == Wm )
	      maxord = maxorder->data[p] + 2;
	    else
	      maxord = cm;
	      
	    if((Wm < Q->cost1[q] || (Wm == Q->cost1[q] && maxord < maxorder->data[q]))){
	      label->data[q] = label->data[p];
	      pred->data[q] = p;
	      
	      if (Q->color[q] != BLACK) {
		gft::Heap32fi_lex::Update_MinPolicy(Q, q, Wm, maxord);
	      } else {
		Q->color[q] = WHITE;
		Q->cost1[q] = Wm;
		maxorder->data[q] = maxord;
		gft::Heap32fi_lex::Insert_MinPolicy(Q, q);
	      }
	    }
	    else if(p == pred->data[q]){
	      if((label->data[q] != label->data[p])||
		 (Wm != Q->cost1[q]) ||
		 (maxord > maxorder->data[q])){
		if(T == NULL)
		  T = (int *)malloc((pred->n+1)*sizeof(int));

		Removal[0] = 1; Removal[1] = q;
		DOIFT_TreeRemoval_5(Removal, sg, Q, pred,
				    maxorder, T);
		break;
	      }
	    }
	  }
	}
	//-------------------------------------------------------
      }
      //----------
      free(i_inv);
      if(T != NULL)
	free(T);
    }




    void DOIFT__2(sGraph *graph,
		  sGraph *transpose,
		  int *S,
		  int *R,
		  int *label,
		  sHeap32fi_lex *Q,
		  int *pred,
		  int *maxorder,
		  int *iter){
      int p,q,i,j,tq,cm;
      float w,Wm;
      int Removal[2];
      int maxord;
      int *T = NULL;
      sGraph *g;

      if(R != NULL && R[0] > 0){
	if(T == NULL)
	  T = (int *)malloc((graph->nnodes+1)*sizeof(int));
	DOIFT_TreeRemoval_5(R, graph, Q, pred,
			    maxorder, T);
      }

      for(i = 1; i <= S[0]; i++){
	p = S[i];
	pred[p] = NIL;
	if(label[p]==0)
	  cm = 1; //1;
	else
	  cm = 2; //2;
	  
	if (Q->color[p] != BLACK) {
	  gft::Heap32fi_lex::Update_MinPolicy(Q, p, 0.0, cm);
	} else {
	  Q->color[p] = WHITE;
	  Q->cost1[p] = 0.0;
	  maxorder[p] = cm;
	  gft::Heap32fi_lex::Insert_MinPolicy(Q, p);
	}
	Q->color[p] = GRAY;
      }

      //----------
      while(!Heap32fi_lex::IsEmpty(Q)) {
	Heap32fi_lex::Remove_MinPolicy(Q, &p);
	
	if(label[p]==0){ g = transpose; cm = 1; } //1;
	else{   	 g = graph;     cm = 2; } //2;
	//-------------------------------------------------------
	for(i = 0; i < g->nodes[p].outdegree; i++){
	  q = g->nodes[p].adjList[i];

	  w = g->nodes[p].Warcs[i];
	  
	  Wm = MAX(w, Q->cost1[p]);
	    
	  if(Wm != w)
	    maxord = maxorder[p];
	  else if( Q->cost1[p] == Wm )
	    maxord = maxorder[p] + 2;
	  else
	    maxord = cm;
	      
	  if((Wm < Q->cost1[q] || (Wm == Q->cost1[q] && maxord < maxorder[q]))){
	    label[q] = label[p];
	    pred[q] = p;
	      
	    if (Q->color[q] != BLACK) {
	      gft::Heap32fi_lex::Update_MinPolicy(Q, q, Wm, maxord);
	    } else {
	      Q->color[q] = WHITE;
	      Q->cost1[q] = Wm;
	      maxorder[q] = maxord;
	      gft::Heap32fi_lex::Insert_MinPolicy(Q, q);
	    }
	  }
	  else if(p == pred[q]){
	    if((label[q] != label[p])||
	       (Wm != Q->cost1[q]) ||
	       (maxord > maxorder[q])){
	      if(T == NULL)
		T = (int *)malloc((graph->nnodes+1)*sizeof(int));
	      
	      Removal[0] = 1; Removal[1] = q;
	      DOIFT_TreeRemoval_5(Removal, graph, Q, pred,
				  maxorder, T);
	      break;
	    }
	  }
	}
      
	//-------------------------------------------------------
      }
      //----------
      if(T != NULL)
	free(T);
    }




    void DOIFT__2p(sImageGraph *sg,
		   int *S, //new seeds only.
		   int *R, //new removal only.
		   sImage32 *label,
		   sHeap32fi_lex *Q,
		   sImage32 *pred,
		   sImage32 *maxorder,
		   int *iter){
      int *i_inv;
      int u_x,u_y,v_x,v_y;
      int p,q,i,j,tq,cm;
      float w,Wm;
      int Removal[2];
      int maxord;
      int *T = NULL;
      bool seenfirst_pq, seenfirst_qp, aux;
      
      i_inv = gft::AdjRel::InverseIndexes(sg->A);

      if(R != NULL && R[0] > 0){
	if(T == NULL)
	  T = (int *)malloc((pred->n+1)*sizeof(int));
	DOIFT_TreeRemoval_5(R, sg, Q, pred,
			    maxorder, T);
      }

      for(i = 1; i <= S[0]; i++){
	p = S[i];
	pred->data[p] = NIL;
	if(label->data[p]==0)
	  cm = 1; //1;
	else
	  cm = 2; //2;
	  
	if (Q->color[p] != BLACK) {
	  gft::Heap32fi_lex::Update_MinPolicy(Q, p, 0.0, cm);
	} else {
	  Q->color[p] = WHITE;
	  Q->cost1[p] = 0.0;
	  maxorder->data[p] = cm;
	  gft::Heap32fi_lex::Insert_MinPolicy(Q, p);
	}
	Q->color[p] = GRAY;
      }

      //----------
      while(!Heap32fi_lex::IsEmpty(Q)) {
	Heap32fi_lex::Remove_MinPolicy(Q, &p);
	//-------------------------------------------------------
	u_x = p%label->ncols; 
	u_y = p/label->ncols; 
	
	for(i=1; i<sg->A->n; i++){
	  v_x = u_x + sg->A->dx[i];
	  v_y = u_y + sg->A->dy[i];
	  if(v_x >= 0 && v_x < label->ncols &&
	     v_y >= 0 && v_y < label->nrows){	  
	    q = v_x + label->ncols*v_y;
	    
	    if(label->data[p]==0){
	      j = i_inv[i];
	      w = (sg->n_link[q])[j];
	      cm = 1; //1;
	    }
	    else{
	      w = (sg->n_link[p])[i];
	      cm = 2; //2;
	    }
	      
	    Wm = MAX(w, Q->cost1[p]);
	    
	    if(Wm != w)
	      maxord = maxorder->data[p];
	    else if( Q->cost1[p] == Wm )
	      maxord = maxorder->data[p] + 2;
	    else
	      maxord = cm;

	    seenfirst_pq = seenfirst_qp = false;
	    if(pred->data[q] != NIL){
	      aux = (w > Q->cost1[p] && Q->cost1[q] > Q->cost1[pred->data[q]] && w == Q->cost1[q]);
	      seenfirst_pq = aux && (Q->cost1[p] < Q->cost1[pred->data[q]] ||
				     (Q->cost1[p] == Q->cost1[pred->data[q]] &&
				      maxorder->data[p] < maxorder->data[pred->data[q]]));
	      seenfirst_qp = aux && (Q->cost1[p] > Q->cost1[pred->data[q]] ||
				     (Q->cost1[p] == Q->cost1[pred->data[q]] &&
				      maxorder->data[p] > maxorder->data[pred->data[q]]));
	    }
	    
	    if(Wm < Q->cost1[q] ||
	       (Wm == Q->cost1[q] && ((maxord < maxorder->data[q] && !seenfirst_qp) || seenfirst_pq))){
	      label->data[q] = label->data[p];
	      pred->data[q] = p;
	      if(Q->color[q] == GRAY)
		gft::Heap32fi_lex::Delete_MinPolicy(Q, q);
	      Q->color[q] = WHITE;
	      Q->cost1[q] = Wm;
	      maxorder->data[q] = maxord;
	      gft::Heap32fi_lex::Insert_MinPolicy(Q, q);
	    }
	    else if(p == pred->data[q]){
	      if((label->data[q] != label->data[p])||
		 (Wm != Q->cost1[q]) ||
		 (maxord > maxorder->data[q])){
		if(T == NULL)
		  T = (int *)malloc((pred->n+1)*sizeof(int));

		Removal[0] = 1; Removal[1] = q;
		DOIFT_TreeRemoval_5(Removal, sg, Q, pred,
				    maxorder, T);
		break;
	      }
	    }
	  }
	}
	//-------------------------------------------------------
      }
      //----------
      free(i_inv);
      if(T != NULL)
	free(T);
    }




    void DOIFT__2p(sGraph *graph,
		   sGraph *transpose,
		   int *S,
		   int *R,
		   int *label,
		   sHeap32fi_lex *Q,
		   int *pred,
		   int *maxorder,
		   int *iter){
      int p,q,i,j,tq,cm;
      float w,Wm;
      int Removal[2];
      int maxord;
      int *T = NULL;
      sGraph *g;
      bool seenfirst_pq, seenfirst_qp, aux;
      
      if(R != NULL && R[0] > 0){
	if(T == NULL)
	  T = (int *)malloc((graph->nnodes+1)*sizeof(int));
	DOIFT_TreeRemoval_5(R, graph, Q, pred,
			    maxorder, T);
      }

      for(i = 1; i <= S[0]; i++){
	p = S[i];
	pred[p] = NIL;
	if(label[p]==0)
	  cm = 1; //1;
	else
	  cm = 2; //2;
	  
	if (Q->color[p] != BLACK) {
	  gft::Heap32fi_lex::Update_MinPolicy(Q, p, 0.0, cm);
	} else {
	  Q->color[p] = WHITE;
	  Q->cost1[p] = 0.0;
	  maxorder[p] = cm;
	  gft::Heap32fi_lex::Insert_MinPolicy(Q, p);
	}
	Q->color[p] = GRAY;
      }

      //----------
      while(!Heap32fi_lex::IsEmpty(Q)) {
	Heap32fi_lex::Remove_MinPolicy(Q, &p);
	
	if(label[p]==0){ g = transpose; cm = 1; } //1;
	else{   	 g = graph;     cm = 2; } //2;
	//-------------------------------------------------------
	for(i = 0; i < g->nodes[p].outdegree; i++){
	  q = g->nodes[p].adjList[i];

	  w = g->nodes[p].Warcs[i];
	  
	  Wm = MAX(w, Q->cost1[p]);
	    
	  if(Wm != w)
	    maxord = maxorder[p];
	  else if( Q->cost1[p] == Wm )
	    maxord = maxorder[p] + 2;
	  else
	    maxord = cm;

	  seenfirst_pq = seenfirst_qp = false;
	  if(pred[q] != NIL){
	    aux = (w > Q->cost1[p] && Q->cost1[q] > Q->cost1[pred[q]] && w == Q->cost1[q]);
	    seenfirst_pq = aux && (Q->cost1[p] < Q->cost1[pred[q]] ||
				   (Q->cost1[p] == Q->cost1[pred[q]] &&
				    maxorder[p] < maxorder[pred[q]]));
	    seenfirst_qp = aux && (Q->cost1[p] > Q->cost1[pred[q]] ||
				   (Q->cost1[p] == Q->cost1[pred[q]] &&
				    maxorder[p] > maxorder[pred[q]]));
	  }
	  
	  if(Wm < Q->cost1[q] ||
	     (Wm == Q->cost1[q] && ((maxord < maxorder[q] && !seenfirst_qp) || seenfirst_pq))){
	    label[q] = label[p];
	    pred[q] = p;
	    if(Q->color[q] == GRAY)
	      gft::Heap32fi_lex::Delete_MinPolicy(Q, q);
	    Q->color[q] = WHITE;
	    Q->cost1[q] = Wm;
	    maxorder[q] = maxord;
	    gft::Heap32fi_lex::Insert_MinPolicy(Q, q);
	  }
	  else if(p == pred[q]){
	    if((label[q] != label[p])||
	       (Wm != Q->cost1[q]) ||
	       (maxord > maxorder[q])){
	      if(T == NULL)
		T = (int *)malloc((graph->nnodes+1)*sizeof(int));
	      
	      Removal[0] = 1; Removal[1] = q;
	      DOIFT_TreeRemoval_5(Removal, graph, Q, pred,
				  maxorder, T);
	      break;
	    }
	  }
	}
      
	//-------------------------------------------------------
      }
      //----------
      if(T != NULL)
	free(T);
    }




    void DOIFT__2p_hierarchy(sGraph *graph,
			     sGraph *transpose,
			     int *S,
			     int *R,
			     int *label,
			     sHeap32fi_lex *Q,
			     int *pred,
			     int *maxorder,
			     sHeap32 *E,
			     int *H,
			     int *obj_area){
      int p,q,b,i,j,tq,cm;
      float w,Wm;
      int Removal[2];
      int maxord, T_area;
      int *T = NULL;
      sGraph *g;
      bool seenfirst_pq, seenfirst_qp, aux;

      T_area = *obj_area;
      
      if(R != NULL && R[0] > 0){
	if(T == NULL)
	  T = (int *)malloc((graph->nnodes+1)*sizeof(int));
	DOIFT_TreeRemoval_5(R, graph, Q, pred,
			    maxorder, T);
      }

      for(i = 1; i <= S[0]; i++){
	p = S[i];
	pred[p] = NIL;
	if(label[p]==0)
	  cm = 1; //1;
	else
	  cm = 2; //2;
	  
	if (Q->color[p] != BLACK) {
	  gft::Heap32fi_lex::Update_MinPolicy(Q, p, 0.0, cm);
	} else {
	  Q->color[p] = WHITE;
	  Q->cost1[p] = 0.0;
	  maxorder[p] = cm;
	  gft::Heap32fi_lex::Insert_MinPolicy(Q, p);
	}
	Q->color[p] = GRAY;
      }

      //----------
      while(!Heap32fi_lex::IsEmpty(Q)) {
	Heap32fi_lex::Remove_MinPolicy(Q, &p);
	
	if(label[p]==0){
	  g = transpose; cm = 1;
	  if(H[p] == 0){
	    (*obj_area)--;
	    H[p] = T_area;
	  }
	  if(E->color[p] == GRAY)
	    gft::Heap32::Delete_MaxPolicy(E, p);
	} //1;
	else{   	 g = graph;     cm = 2; } //2;
	//-------------------------------------------------------
	for(i = 0; i < g->nodes[p].outdegree; i++){
	  q = g->nodes[p].adjList[i];

	  w = g->nodes[p].Warcs[i];
	  
	  Wm = MAX(w, Q->cost1[p]);
	    
	  if(Wm != w)
	    maxord = maxorder[p];
	  else if( Q->cost1[p] == Wm )
	    maxord = maxorder[p] + 2;
	  else
	    maxord = cm;

	  seenfirst_pq = seenfirst_qp = false;
	  if(pred[q] != NIL){
	    aux = (w > Q->cost1[p] && Q->cost1[q] > Q->cost1[pred[q]] && w == Q->cost1[q]);
	    seenfirst_pq = aux && (Q->cost1[p] < Q->cost1[pred[q]] ||
				   (Q->cost1[p] == Q->cost1[pred[q]] &&
				    maxorder[p] < maxorder[pred[q]]));
	    seenfirst_qp = aux && (Q->cost1[p] > Q->cost1[pred[q]] ||
				   (Q->cost1[p] == Q->cost1[pred[q]] &&
				    maxorder[p] > maxorder[pred[q]]));
	  }
	  
	  if(Wm < Q->cost1[q] ||
	     (Wm == Q->cost1[q] && ((maxord < maxorder[q] && !seenfirst_qp) || seenfirst_pq))){
	    label[q] = label[p];
	    pred[q] = p;
	    if(Q->color[q] == GRAY)
	      gft::Heap32fi_lex::Delete_MinPolicy(Q, q);
	    Q->color[q] = WHITE;
	    Q->cost1[q] = Wm;
	    maxorder[q] = maxord;
	    gft::Heap32fi_lex::Insert_MinPolicy(Q, q);
	  }
	  else if((p == pred[q]) &&
		  ((label[q] != label[p])||
		   (Wm != Q->cost1[q]) ||
		   (maxord > maxorder[q]))){

	    if(label[q] == label[p])
	      printf("Ops...\n");
	    
	    if(T == NULL)
	      T = (int *)malloc((graph->nnodes+1)*sizeof(int));
	      
	    Removal[0] = 1; Removal[1] = q;
	    DOIFT_TreeRemoval_5(Removal, graph, Q, pred,
				maxorder, T);
	    break;
	  }
	  else if(label[p] == 0 && label[q] == 1){
	    if(E->color[q] != GRAY)
	      gft::Heap32::Insert_MaxPolicy(E, q);
	  }
	  else if(label[p] == 1 && label[q] == 0 && H[q] > 0){
	    if(E->color[p] != GRAY)
	      gft::Heap32::Insert_MaxPolicy(E, p);
	  }
	}
      
	//-------------------------------------------------------
      }
      //----------
      if(T != NULL)
	free(T);
    }
    



    int *OIFT_area_hierarchy(sGraph *graph,
			     sGraph *transpose,
			     int *S0,
			     int *S1){
			     //gft::sImage32 *SP){
      /*
      char filename[512];
      gft::sImage32 *temp = NULL;
      int sp,it = 0;
      temp = gft::Image32::Create(SP);
      */
      //-------------------------------------
      int *H = NULL, *Rset = NULL, *S = NULL;
      gft::sHeap32fi_lex *Q = NULL;
      gft::sHeap32 *E = NULL;
      float *cost1=NULL;
      int *energy=NULL, *maxorder=NULL, *label=NULL, *pred=NULL;
      int p,p_max,i, obj_area;
      H = gft::AllocIntArray(graph->nnodes);
      S = gft::AllocIntArray(graph->nnodes+1);
      cost1  = gft::AllocFloatArray(graph->nnodes);
      energy   = gft::AllocIntArray(graph->nnodes);
      maxorder = gft::AllocIntArray(graph->nnodes);
      label    = gft::AllocIntArray(graph->nnodes);
      pred     = gft::AllocIntArray(graph->nnodes);
      for(p=0; p < graph->nnodes; p++){
	H[p] = 0;
	cost1[p] = FLT_MAX;
	maxorder[p] = INT_MAX;
	pred[p] = NIL;
	label[p] = NIL;
      }
      Q = gft::Heap32fi_lex::Create(graph->nnodes,
				    cost1, maxorder);
      
      for(i=1; i<=S1[0]; i++){
	p = S1[i];
	label[p] = 1;
      }
      gft::ift::IFT_fmax(graph, S1, label, energy);
      //------------------------------------------
      /*
      for(p = 0; p < temp->n; p++){
	sp = SP->data[p];
	temp->data[p] = energy[sp];
      }
      gft::Image32::Write(temp, (char*)"energy.pgm");
      */
      //------------------------------------------
      for(i=1; i<=S1[0]; i++){
	p = S1[i];
	energy[p] = -1;
      }
      E = gft::Heap32::Create(graph->nnodes, energy);
      
      for(p=0; p < graph->nnodes; p++){
	label[p] = NIL;
      }
      for(i=1; i<=S1[0]; i++){
	p = S1[i];
	label[p] = 1;
      }
      for(i=1; i<=S0[0]; i++){
	p = S0[i];
	label[p] = 0;
      }

      S[0] = 0;
      for(i=1; i<=S1[0]; i++){
	S[0]++;
	S[S[0]] = S1[i];
      }
      for(i=1; i<=S0[0]; i++){
	S[0]++;
	S[S[0]] = S0[i];
      }
      //------------------------------
      //printf("#Internal seeds: %d\n", S1[0]);
      //------------------------------
      obj_area = graph->nnodes;
      while(obj_area > S1[0]){
      
	gft::ift::DOIFT__2p_hierarchy(graph,
				      transpose,
				      S,
				      Rset,
				      label,
				      Q,
				      pred,
				      maxorder,
				      E,
				      H,
				      &obj_area);
	//------------------------------------------
	/*
	it++;
	printf("it: %d, obj_area: %d\n", it, obj_area);
	for(p = 0; p < temp->n; p++){
	  sp = SP->data[p];
	  temp->data[p] = (E->color[sp] == GRAY)*255;
	}
	sprintf(filename, "objboundary_%02d.pgm", it);
	gft::Image32::Write(temp, filename);
	for(p = 0; p < temp->n; p++){
	  sp = SP->data[p];
	  temp->data[p] = label[sp]*255;
	}
	sprintf(filename, "label_%02d.pgm", it);
	gft::Image32::Write(temp, filename);
	*/
	//------------------------------------------
	
	Heap32::Remove_MaxPolicy(E, &p_max);
	
	S[0] = 1;
	S[1] = p_max;
	label[p_max] = 0;

      }

      for(i=1; i<=S1[0]; i++){
	p = S1[i];
	H[p] = S1[0];
      }
      
      gft::Heap32::Destroy(&E);      
      gft::Heap32fi_lex::Destroy(&Q);
      if(S  != NULL) gft::FreeIntArray(&S);
      if(energy   != NULL) gft::FreeIntArray(&energy);
      if(label    != NULL) gft::FreeIntArray(&label);
      if(pred     != NULL) gft::FreeIntArray(&pred);
      if(maxorder != NULL) gft::FreeIntArray(&maxorder);
      if(cost1 != NULL) gft::FreeFloatArray(&cost1);
      //--------------------------------------------
      //gft::Image32::Destroy(&temp);
      //--------------------------------------------
      return H;
    }



    void DOIFT__3(sImageGraph *sg,
		  int *S, //new seeds only.
		  int *R, //new removal only.
		  sImage32 *label,
		  sHeap32fif_lex *Q,
		  sImage32 *pred,
		  sImage32 *maxorder,
		  int *iter){
      int *i_inv;
      int u_x,u_y,v_x,v_y;
      int p,q,i,j,tq,cm;
      float w,Wm;
      int Removal[2];
      int maxord;
      int *T = NULL;
      bool seenfirst_pq, seenfirst_qp, eq;

      printf("iter: %d\n",*iter);
      
      T = (int *)malloc((pred->n+1)*sizeof(int));
      i_inv = gft::AdjRel::InverseIndexes(sg->A);

      if(R != NULL && R[0] > 0){
	if(T == NULL)
	  T = (int *)malloc((pred->n+1)*sizeof(int));
	DOIFT_TreeRemoval__5(R, sg, Q, pred,
			     maxorder, T);
      }
      
      /*
      if(*iter >= 100){
	for(p=0; p<label->n; p++){
	  if(Q->color[p] == GRAY)
	    label->data[p] += 2;
	  if(Q->cost3[p] == FLT_MAX)
	    label->data[p] = 4;
	}
	return;
      }
      */
      
      for(i = 1; i <= S[0]; i++){
	p = S[i];
	pred->data[p] = NIL;
	if(label->data[p]==0)
	  cm = 2; //2;
	else
	  cm = 2; //4;
	  
	if (Q->color[p] != BLACK) {
	  gft::Heap32fif_lex::Update_MinPolicy(Q, p, 0.0, cm, -1.0); //INT_MAX);
	} else {
	  Q->color[p] = WHITE;
	  Q->cost1[p] = 0.0;
	  maxorder->data[p] = cm; //INT_MAX;
	  Q->cost3[p] = -1.0;
	  gft::Heap32fif_lex::Insert_MinPolicy(Q, p);
	}
	Q->color[p] = GRAY;
      }

      //----------
      while(!Heap32fif_lex::IsEmpty(Q)) {
	Heap32fif_lex::Remove_MinPolicy(Q, &p);
	//maxorder->data[p]--;
	//-------------------------------------------------------
	/*
	printf("p: %d, c1: %f, c2: %d, c3: %f\n",
	       p, Q->cost1[p], Q->cost2[p], Q->cost3[p]);
	*/
	//-------------------------------------------------------
	u_x = p%label->ncols; //PixelX(label, p);
	u_y = p/label->ncols; //PixelY(label, p);
	
	for(i=1; i<sg->A->n; i++){
	  v_x = u_x + sg->A->dx[i];
	  v_y = u_y + sg->A->dy[i];
	  //if(Image32::IsValidPixel(label,v_x,v_y)){
	  if(v_x >= 0 && v_x < label->ncols &&
	     v_y >= 0 && v_y < label->nrows){	  
	    q = v_x + label->ncols*v_y;
	    
	    if(label->data[p]==0){
	      j = i_inv[i];
	      w = (sg->n_link[q])[j];
	      cm = 2; //2;
	    }
	    else{
	      w = (sg->n_link[p])[i];
	      cm = 2; //4;
	    }
	      
	    Wm = MAX(w, Q->cost1[p]);


	    /*
	    if(*iter < 100)
	      eq = (p == 1 && q == 5);
	    else
	      eq = (p == 6 && q == 2) || (p == 1 && q == 5);
	    */
	    
	    /*
	    if(Wm != w)
	      maxord = maxorder->data[p]; //+1
	    else if( Q->cost1[p] == Wm )
	      maxord = maxorder->data[p]; //+1 + 4
	    else
	    */
	    maxord = cm;

	    
	    //eq = w == Wm && w == Q->cost1[q] && w == Q->cost3[q];
	    //seenfirst_pq = eq && pred->data[q] != NIL && Q->cost1[p] < Q->cost1[pred->data[q]];
	    //seenfirst_qp = eq && pred->data[q] != NIL && Q->cost1[p] > Q->cost1[pred->data[q]];
	    seenfirst_pq = seenfirst_qp = false;
	    /*
	    if(Wm < Q->cost1[q] ||
	       (Wm == Q->cost1[q] && seenfirst_pq) ||
	       (Wm == Q->cost1[q] && !seenfirst_qp && maxord < maxorder->data[q]) ||
	       (Wm == Q->cost1[q] && !seenfirst_qp && maxord == maxorder->data[q] && w < Q->cost3[q])){
	    */

	    if (Q->color[q] == BLACK) continue;
	    /*
	    if(eq || Wm < Q->cost1[q]){
	       //	       ||
	       //(Wm == Q->cost1[q] && maxord < maxorder->data[q])){
	       */
	    if(Wm < Q->cost1[q] ||
	       (Wm == Q->cost1[q] && maxord < maxorder->data[q]) ||
	       (Wm == Q->cost1[q] && maxord == maxorder->data[q] && w < Q->cost3[q])){
	      label->data[q] = label->data[p];
	      pred->data[q] = p;
	      
	      if (Q->color[q] != BLACK) {
		gft::Heap32fif_lex::Update_MinPolicy(Q, q, Wm, maxord, w);
	      } else {
		Q->color[q] = WHITE;
		Q->cost1[q] = Wm;
		maxorder->data[q] = maxord;
		Q->cost3[q] = w;
		gft::Heap32fif_lex::Insert_MinPolicy(Q, q);
	      }
	    }
	    /*
	    else if((w == Wm && Q->cost1[q] == Wm && Q->cost3[q] == Wm) &&
		    pred->data[q] != NIL && Q->cost1[p] < Q->cost1[pred->data[q]]){
		Removal[0] = 1; Removal[1] = q;
		DOIFT_TreeRemoval__5(Removal, sg, Q, pred,
				     maxorder, T);
		break;
	    }
	    */
	    else if(p == pred->data[q]){
	      if((label->data[q] != label->data[p])||
		 (Wm != Q->cost1[q])){
		 //(maxord > maxorder->data[q])){
		Removal[0] = 1; Removal[1] = q;
		DOIFT_TreeRemoval__5(Removal, sg, Q, pred,
				     maxorder, T);
		break;
	      }
	    }
	  }
	}
	//-------------------------------------------------------
      }
      //----------
      free(i_inv);
      free(T);

      *iter = 100;
    }

    
    

  } /*end ift namespace*/
} /*end gft namespace*/


