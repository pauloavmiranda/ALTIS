/*
ALTIS consists of three steps - lungs-and-trachea extraction, seed estimation for each object, and seed labeling and delineation of each lung and trachea.

This implementation assumes that the patient orientation in the CT scan is from inferior to superior along the axial slices (z-axis), from right to left along the sagittal slices (x-axis), and from posterior to anterior along the coronal slices (y-axis).
*/

//#define APPDEBUG  1

#include "gft.h"


void FixOrientation(gft::sScene32 *scn,
		    bool *flippedX,
		    bool *flippedY,
		    bool *flippedZ){
  int icod, jcod, kcod;
  int i, j , k;
  *flippedX = false;
  *flippedY = false;
  *flippedZ = false;
  if(scn->nii_hdr != NULL){
    if(scn->nii_hdr->sform_code > 0)
      nifti_mat44_to_orientation(scn->nii_hdr->sto_xyz, &icod, &jcod, &kcod);
    else if(scn->nii_hdr->qform_code > 0)
      nifti_mat44_to_orientation(scn->nii_hdr->qto_xyz, &icod, &jcod, &kcod);
    else
      return;
    
    if(icod == NIFTI_L2R){
      //Invert the x-axis:
      gft::Scene32::FlipX(scn);
      *flippedX = true;
    }
    if(jcod == NIFTI_A2P){
      //Invert the y-axis:
      gft::Scene32::FlipY(scn);
      *flippedY = true;
    }
    if(kcod == NIFTI_S2I){
      //Invert the z-axis:
      gft::Scene32::FlipZ(scn);
      *flippedZ = true;
    }
  }
}



void UndoFixOrientation(gft::sScene32 *scn,
			bool flippedX,
			bool flippedY,
			bool flippedZ){
  if(flippedX)
    gft::Scene32::FlipX(scn);
  if(flippedY)
    gft::Scene32::FlipY(scn);
  if(flippedZ)
    gft::Scene32::FlipZ(scn);
}



bool GetLungsSeeds_aux(gft::sScene8 *bin,
		       gft::sAdjRel3 *A,
		       gft::sScene32 *Slabel){
  gft::sScene32 *label=NULL;
  int Lmax, l, x;
  int *area=NULL;
  float *CoG_x = NULL;
  int imax,imax2,i,p,n=bin->n;
  int nlb,nlb2;
  int *S = NULL;
  
  if(gft::Scene8::GetMaxVal(bin)==0)
    return false;
  label = gft::Scene8::LabelBinComp(bin, A);
  Lmax = gft::Scene32::GetMaxVal(label);
  area  = (int *)calloc(Lmax+1, sizeof(int));
  CoG_x = (float *)calloc(Lmax+1, sizeof(float));
  if(area == NULL || CoG_x == NULL)
    return false;
  
  for (p=0; p < n; p++){
    l = label->data[p];
    if(l > 0){
      x = gft::Scene32::GetAddressX(label, p);
      area[l]++;
      CoG_x[l] += x;
    }
  }

  for(i=1; i <= Lmax; i++)
    if(area[i] != 0)
      CoG_x[i] /= area[i];
  
  //Largest component:
  imax = 0;
  for (i=1; i <= Lmax; i++) 
    if (area[i]>area[imax])
      imax = i;

  //Second largest component:
  imax2 = 0;
  for (i=1; i <= Lmax; i++) 
    if (area[i]>area[imax2] && i != imax)
      imax2 = i;  

  //printf("imax2: %d\n", imax2);
  if(imax2 != 0){
  
    if(CoG_x[imax] > CoG_x[imax2]){
      nlb = 2;
      nlb2 = 1;
    }
    else{
      nlb = 1;
      nlb2 = 2;
    }
    
    for(p=0; p < n; p++){
      l = label->data[p];
      if(l == imax)
	Slabel->data[p] = nlb;
      else if(l == imax2)
	Slabel->data[p] = nlb2;
    }

  }
  
  gft::Scene32::Destroy(&label);
  free(area);
  free(CoG_x);

  return !(imax2 == 0);
}



void GetLungsSeeds(gft::sScene8 *bin,
		   gft::sAdjRel3 *A,
		   gft::sScene32 *Slabel,
		   float gamma){
  gft::sScene8 *I_Slungs = NULL;
  gft::sSet *S = NULL;
  int ntimes = 0;
  bool r = false;
  do{
    S = NULL;
    I_Slungs = gft::Scene8::ErodeBin(bin, &S, gamma, A);
    gft::Set::Destroy(&S);
    
    r = GetLungsSeeds_aux(I_Slungs, A, Slabel);
    gft::Scene8::Destroy(&I_Slungs);
    ntimes++;
    gamma *= 0.9;
  }while(!r && ntimes <= 3);
  //printf("r: %d, ntimes: %d\n", r, ntimes);
}




gft::sScene32 *GeodesicDistanceTransform(gft::sScene8 *I4,
					 gft::sAdjRel3 *A,
					 int *S){
  gft::sPQueue32 *Q=NULL;
  int i,j,p,q,cst;
  gft::Voxel u,v;
  gft::sScene32 *C;
  int Dmax = 0;
  float *D;
  int *Dint;
  
  if(S == NULL)
    return NULL;

  D = gft::AdjRel3::GetDistanceArray(A);
  Dint = (int *)calloc(A->n, sizeof(int));
  for(i = 0; i < A->n; i++){
    Dint[i] = ROUND(10.0*D[i]);
    Dmax = MAX(Dint[i], Dmax);    
  }
  C = gft::Scene32::Create(I4->xsize, I4->ysize, I4->zsize);
  for(p = 0; p < C->n; p++){
    if(I4->data[p] > 0)
      C->data[p] = INT_MAX;
    else
      C->data[p] = 0;
  }
  
  Q = gft::PQueue32::Create(10*(Dmax+1)+2, C->n, C->data);
  
  for(i = 1; i <= S[0]; i++){
    p = S[i];
    C->data[p] = 0;
    gft::PQueue32::InsertElem(&Q, p);
  }
  
  while(!gft::PQueue32::IsEmpty(Q)) {
    p = gft::PQueue32::RemoveMinFIFO(Q);
    u.c.x = gft::Scene32::GetAddressX(C, p);
    u.c.y = gft::Scene32::GetAddressY(C, p);
    u.c.z = gft::Scene32::GetAddressZ(C, p);
    
    for(i = 1; i < A->n; i++){
      v.v = u.v + A->d[i].v;
      if(gft::Scene32::IsValidVoxel(C, v)){
	q = gft::Scene32::GetVoxelAddress(C,v);
	if(Q->L.elem[q].color != BLACK &&
	   I4->data[q] > 0){
	  cst = C->data[p] + Dint[i];
	  
	  if(cst < C->data[q]){
	    if(Q->L.elem[q].color == GRAY)
	      gft::PQueue32::RemoveElem(Q, q);
	    C->data[q] = cst;
	    gft::PQueue32::InsertElem(&Q, q);
	  }
	}
      }
    }
  }
  gft::PQueue32::Destroy(&Q);
  gft::FreeFloatArray(&D);
  free(Dint);
  return C;
}



void SelectComponentFromTop(gft::sScene8 *bin,
			    gft::sAdjRel3 *A,
			    gft::sScene32 *Slabel){
  gft::sScene32 *label=NULL;
  int Lmax, l, z, i_zmax;
  int *area=NULL;
  float *CoG_z = NULL;
  int i,p,n=bin->n;
  int *S = NULL;
  
  if(gft::Scene8::GetMaxVal(bin)==0)
    return;
  label = gft::Scene8::LabelBinComp(bin, A);
  Lmax = gft::Scene32::GetMaxVal(label);
  area  = (int *)calloc(Lmax+1, sizeof(int));
  CoG_z = (float *)calloc(Lmax+1, sizeof(float));
  
  for (p=0; p < n; p++){
    l = label->data[p];
    if(l > 0){
      z = gft::Scene32::GetAddressZ(label, p);
      area[l]++;
      CoG_z[l] += z;
    }
  }

  for(i=1; i <= Lmax; i++)
    if(area[i] != 0)
      CoG_z[i] /= area[i];
  
  //Top component:
  i_zmax = 0;
  for (i=1; i <= Lmax; i++) 
    if (CoG_z[i] > CoG_z[i_zmax])
      i_zmax = i;

  for(p=0; p < n; p++){
    l = label->data[p];
    if(l == i_zmax)
      Slabel->data[p] = 3;
  }

  gft::Scene32::Destroy(&label);
  free(area);
  free(CoG_z);
}



gft::sScene32 *ComputeGradientImage(gft::sScene32 *I2,
				    gft::sAdjRel3 *A){
  gft::sScene32 *I5;
  gft::Voxel u,v;
  int p,q,i,Ip;
  float *D = NULL;
  float d, sum1, sum2;
  I5 = gft::Scene32::Create(I2);
  D = gft::AdjRel3::GetDistanceArray(A);
  for(p = 0; p < I2->n; p++){
    u.c.x = gft::Scene32::GetAddressX(I2, p);
    u.c.y = gft::Scene32::GetAddressY(I2, p);
    u.c.z = gft::Scene32::GetAddressZ(I2, p);
    sum1 = 0.0;
    sum2 = 0.0;
    Ip = I2->data[p];
    for(i = 1; i < A->n; i++){
      v.v = u.v + A->d[i].v;
      if(gft::Scene32::IsValidVoxel(I2, v)){
	q = gft::Scene32::GetVoxelAddress(I2, v);
	d = expf(-(D[i]*D[i])/6.0);
	sum1 += abs(I2->data[q]-Ip)*d;
	sum2 += d;
      }
    }
    I5->data[p] = ROUND(sum1/sum2);
  }
  gft::FreeFloatArray(&D);
  return I5;
}



gft::sScene32 *LabelingByIFTPeak_1st(gft::sScene32 *I5,
				     gft::sAdjRel3 *A,
				     gft::sScene32 *label){
  gft::sPQueue32 *Q=NULL;
  gft::sScene32 *value, *pred;
  gft::Voxel u,v;
  int p,q,i,Imax,cst;
  int *S=NULL;
  pred  = gft::Scene32::Create(I5);
  value = gft::Scene32::Create(I5);
  gft::Scene32::Fill(pred, NIL);
  for(p = 0; p < label->n; p++){
    if(label->data[p]==NIL) value->data[p] = INT_MAX;
    else                    value->data[p] = 0;
  }
  Imax = gft::Scene32::GetMaxVal(I5);
  Q = gft::PQueue32::Create(Imax+2, value->n, value->data);
  S = gft::Scene32::GetMarkers(label, A);
  for(i = 1; i <= S[0]; i++)
    gft::PQueue32::FastInsertElem(Q, S[i]);
  free(S);
  
  while(!gft::PQueue32::IsEmpty(Q)) {
    p = gft::PQueue32::FastRemoveMinFIFO(Q);
    u.c.x = gft::Scene32::GetAddressX(label, p);
    u.c.y = gft::Scene32::GetAddressY(label, p);
    u.c.z = gft::Scene32::GetAddressZ(label, p);	
    
    for(i=1; i<A->n; i++){
      v.v = u.v + A->d[i].v;
      if(gft::Scene32::IsValidVoxel(label, v)){
	q = gft::Scene32::GetVoxelAddress(label,v);
	if(Q->L.elem[q].color != BLACK){
	  cst = MAX(value->data[p], I5->data[q]);
	  if(cst < value->data[q]){
	    if(Q->L.elem[q].color == GRAY)
	      gft::PQueue32::FastRemoveElem(Q, q);
	    value->data[q] = cst;
	    label->data[q] = label->data[p];
	    pred->data[q] = p;
	    gft::PQueue32::FastInsertElem(Q, q);
	  }
	}
      }
    }
  }
  gft::Scene32::Destroy(&value);
  gft::PQueue32::Destroy(&Q);
  return pred;
}




void LabelingByIFTPeak_2nd(gft::sScene32 *I1,
			   gft::sAdjRel3 *A,
			   gft::sScene32 *pred,
			   gft::sScene32 *label_1st,
			   gft::sScene32 *label){
  gft::sPQueue32 *Q=NULL;
  gft::sScene32 *value;
  gft::Voxel u,v;
  int p,q,t,i,k,Imax,cst;
  int n = 10; //to define the n-th ancestor.
  int *S=NULL;
  bool terminal_voxel;
  value = gft::Scene32::Create(I1);
  for(p = 0; p < label->n; p++){
    if(label->data[p]==NIL) value->data[p] = INT_MAX;
    else                    value->data[p] = 0;
  }
  Imax = gft::Scene32::GetMaxVal(I1);
  Q = gft::PQueue32::Create(Imax+2, value->n, value->data);
  S = gft::Scene32::GetMarkers(label, A);
  for(i = 1; i <= S[0]; i++)
    gft::PQueue32::FastInsertElem(Q, S[i]);
  free(S);

  //To use larger seed sets inside the objects:
  for(p = 0; p < label->n; p++){
    if(label_1st->data[p] > 0){
      u.c.x = gft::Scene32::GetAddressX(label, p);
      u.c.y = gft::Scene32::GetAddressY(label, p);
      u.c.z = gft::Scene32::GetAddressZ(label, p);	
      
      terminal_voxel = false;
      for(i=1; i<A->n; i++){
	v.v = u.v + A->d[i].v;
	if(gft::Scene32::IsValidVoxel(label, v)){
	  q = gft::Scene32::GetVoxelAddress(label,v);
	  if(label_1st->data[p] != label_1st->data[q]){
	    terminal_voxel = true;
	    break;
	  }
	}
      }

      if(terminal_voxel){
	t = p;
	for(k = 0; k < n; k++){
	  if(pred->data[t] == NIL)
	    break;
	  t = pred->data[t];
	}
	//Insert t in Q:
	if(Q->L.elem[t].color != GRAY){
	  value->data[t] = 0;
	  //the new seed sets are labeled by the previous IFT.
	  label->data[t] = label_1st->data[t];
	  gft::PQueue32::FastInsertElem(Q, t);
	}
      }
      
    }
  }
  
  while(!gft::PQueue32::IsEmpty(Q)) {
    p = gft::PQueue32::FastRemoveMinFIFO(Q);
    u.c.x = gft::Scene32::GetAddressX(label, p);
    u.c.y = gft::Scene32::GetAddressY(label, p);
    u.c.z = gft::Scene32::GetAddressZ(label, p);	
    
    for(i=1; i<A->n; i++){
      v.v = u.v + A->d[i].v;
      if(gft::Scene32::IsValidVoxel(label, v)){
	q = gft::Scene32::GetVoxelAddress(label,v);
	if(Q->L.elem[q].color != BLACK){
	  cst = MAX(value->data[p], abs(I1->data[q] - I1->data[p]));
	  if(cst < value->data[q]){
	    if(Q->L.elem[q].color == GRAY)
	      gft::PQueue32::FastRemoveElem(Q, q);
	    value->data[q] = cst;
	    label->data[q] = label->data[p];
	    gft::PQueue32::FastInsertElem(Q, q);
	  }
	}
      }
    }
  }
  gft::Scene32::Destroy(&value);
  gft::PQueue32::Destroy(&Q);
}



double Dice(gft::sScene32 *segm, int lb, gft::sScene32 *gtruth){
  gft::sScene32 *bin;
  double dice;
  int p;
  if(gtruth == NULL)
    return -1.0;
  bin = gft::Scene32::Create(segm);
  for(p = 0; p < bin->n; p++){
    if(segm->data[p] == lb)
      bin->data[p] = 1;
    else
      bin->data[p] = 0;
  }
  dice = gft::Scene32::DiceSimilarity(bin, gtruth);
  gft::Scene32::Destroy(&bin);
  return dice;
}



void SaveSeedsByLabel(char *filename,
		      gft::sScene32 *label,
		      int l){
  int p,lb,x,y,z,id,nseeds;
  FILE *fp;
  fp = fopen(filename, "w");
  if (fp == NULL) {
    printf("Error writing seeds.\n");
    exit(1);
  }
  nseeds = 0;
  for(p = 0; p < label->n; p++){
    lb = label->data[p];
    if(lb == NIL) continue;
    nseeds++;
  }
  fprintf(fp, " %d", nseeds);
  for(p = 0; p < label->n; p++){
    lb = label->data[p];
    if(lb == NIL) continue;
    x = gft::Scene32::GetAddressX(label, p);
    y = gft::Scene32::GetAddressY(label, p);
    z = gft::Scene32::GetAddressZ(label, p);
    id = 1;
    if(lb == l) lb = 1;
    else        lb = 0;
    fprintf(fp, " %d %d %d %d %d", x, y, z, id, lb);
  }
  fclose(fp);
}




void SelectLargestCompFitForLungs(gft::sScene8 *bin, gft::sAdjRel3 *A){
  gft::sScene32 *label=NULL;
  int Lmax;
  int *area=NULL;
  float *sx = NULL, *sy = NULL, *sz = NULL;
  float *sdx2 = NULL, *sdy2 = NULL, *sdz2 = NULL;
  int imax,i,l,p,n=bin->n;
  
  if(gft::Scene8::GetMaxVal(bin)==0)
    return;

  label = gft::Scene8::LabelBinComp(bin, A);
  Lmax = gft::Scene32::GetMaxVal(label);
  area = (int *)gft::AllocIntArray(Lmax+1);
  sx = (float *)gft::AllocFloatArray(Lmax+1);
  sy = (float *)gft::AllocFloatArray(Lmax+1);
  sz = (float *)gft::AllocFloatArray(Lmax+1);
  sdx2 = (float *)gft::AllocFloatArray(Lmax+1);
  sdy2 = (float *)gft::AllocFloatArray(Lmax+1);
  sdz2 = (float *)gft::AllocFloatArray(Lmax+1);
  
  for (p=0; p < n; p++){
    l = label->data[p];
    if (l > 0){
      area[l]++;
      sx[l] += gft::Scene8::GetAddressX(bin, p);
      sy[l] += gft::Scene8::GetAddressY(bin, p);
      sz[l] += gft::Scene8::GetAddressZ(bin, p);
    }
  }
  for (p=0; p < n; p++){
    l = label->data[p];
    if (l > 0 && area[l] > 0){
      sdx2[l] += SQUARE(gft::Scene8::GetAddressX(bin, p) - sx[l]/area[l]);
      sdy2[l] += SQUARE(gft::Scene8::GetAddressY(bin, p) - sy[l]/area[l]);
      sdz2[l] += SQUARE(gft::Scene8::GetAddressZ(bin, p) - sz[l]/area[l]);
    }
  }
  
  int size;
  float sdev_x,sdev_y,sdev_z;
  //--------Remover----------------  
  /*
  int *lord = (int *)gft::AllocIntArray(Lmax+1);
  for(l = 0; l <= Lmax; l++)
    lord[l] = l;
  //Sort by area:
  int x,j;
  for(i = 0; i < Lmax; i++){
    x = lord[i+1];
    j = i;
    while(j >= 0 && area[lord[j]] > area[x]){
      lord[j+1] = lord[j];
      j--;
    }
    lord[j+1] = x;
  }

  size = area[lord[Lmax]];
  sdev_x = sqrtf(sdx2[lord[Lmax]]/size);
  sdev_y = sqrtf(sdy2[lord[Lmax]]/size);
  sdev_z = sqrtf(sdz2[lord[Lmax]]/size);
  printf("1st: area=%d, sdev_x: %f, sdev_y: %f, sdev_z: %f\n",
	 size, sdev_x, sdev_y, sdev_z);
  if(Lmax >= 2){
    size = area[lord[Lmax-1]];
    sdev_x = sqrtf(sdx2[lord[Lmax-1]]/size);
    sdev_y = sqrtf(sdy2[lord[Lmax-1]]/size);
    sdev_z = sqrtf(sdz2[lord[Lmax-1]]/size);
    printf("2nd: area=%d, sdev_x: %f, sdev_y: %f, sdev_z: %f\n",
	   size, sdev_x, sdev_y, sdev_z);
  }
  if(Lmax >= 3){
    size = area[lord[Lmax-2]];
    sdev_x = sqrtf(sdx2[lord[Lmax-2]]/size);
    sdev_y = sqrtf(sdy2[lord[Lmax-2]]/size);
    sdev_z = sqrtf(sdz2[lord[Lmax-2]]/size);
    printf("3th: area=%d, sdev_x: %f, sdev_y: %f, sdev_z: %f\n",
	   size, sdev_x, sdev_y, sdev_z);
  }
  gft::FreeIntArray(&lord);
  */
  //-------------------------------
  imax = 0;
  for (i=1; i <= Lmax; i++){
    size = area[i];
    sdev_y = sqrtf(sdy2[i]/size);
    sdev_z = sqrtf(sdz2[i]/size);
    if (area[i]>area[imax] && sdev_z/sdev_y < 3.0) //to remove the stretcher
      imax = i;
  }
  if(imax == 0){ //choose by area only
    for (i=1; i <= Lmax; i++)
      if (area[i]>area[imax])
	imax = i;      
  }
  
  for (p=0; p < n; p++)  
    if (label->data[p]!=imax)
      bin->data[p]=0;
  
  gft::Scene32::Destroy(&label);
  gft::FreeIntArray(&area);
  gft::FreeFloatArray(&sx);
  gft::FreeFloatArray(&sy);
  gft::FreeFloatArray(&sz);
  gft::FreeFloatArray(&sdx2);
  gft::FreeFloatArray(&sdy2);
  gft::FreeFloatArray(&sdz2);
}



int main(int argc, char **argv) {
  gft::sScene32 *scn, *I1, *I2, *I5, *TDG_C, *pred;
  gft::sScene32 *tmp, *label, *label_2nd;
  gft::sScene32 *segm_1st, *segm_2nd, *seeds;
  gft::sScene32 *gtruth_left = NULL, *gtruth_right = NULL;
  gft::sScene8 *I3, *It, *I4, *I4_prime, *I_St;
  gft::sAdjRel3 *A;
  int i, p, Imin, T, Cmax;
  int output_type;
  clock_t end, start;
  double totaltime, dice;
  char buffer[1024];
  char value[128];
  gft::sAttributeList *al;
  bool flippedX, flippedY, flippedZ;

  if (argc < 3){
    fprintf(stdout, "usage:\n");
    fprintf(stdout, "altis <volume> <output_type> [T=value] [left=file] [right=file]\n");
    fprintf(stdout, "\toutput_type......... 0 - segmentation labels,\n");
    fprintf(stdout, "\t                     1 - seeds only.\n");
    fprintf(stdout, "Optional parameters:\n");
    fprintf(stdout, "\tT................... threshold integer value\n");
    fprintf(stdout, "\t                     (if not specified Otsu is used).\n");
    fprintf(stdout, "\tleft................ ground truth for left lung.\n");
    fprintf(stdout, "\tright............... ground truth for right lung.\n");
    exit(0);
  }

  //Original image:
  scn = gft::Scene32::Read(argv[1]);
  output_type = atoi(argv[2]);
  
  buffer[0] = '\0';
  for(i = 3; i < argc; i++){
    strcat(buffer, (char *)" ");
    strcat(buffer, argv[i]);
  }
  al = gft::AttributeList::Create(buffer);

  if( gft::AttributeList::GetAttributeValue(al, (char *)"left", value) )
    gtruth_left = gft::Scene32::Read(value);
  if( gft::AttributeList::GetAttributeValue(al, (char *)"right", value) )
    gtruth_right = gft::Scene32::Read(value);

  start = clock();
  //-----------------------------------------  
  //To remove negative values:
  I1 = gft::Scene32::Clone(scn);
  FixOrientation(I1, &flippedX, &flippedY, &flippedZ);  
  Imin = gft::Scene32::GetMinimumValue(I1);
  if (Imin < 0) {
    for (p = 0; p < I1->n; p++)
      I1->data[p] += (-Imin);
  }

  //All scans were linearly interpolated to
  //the same voxel size, 1.25 x 1.25 x 1.25 mm^3.
  tmp = I1;
  I1 = gft::Scene32::LinearInterp(tmp, 1.25, 1.25, 1.25);
  gft::Scene32::Destroy(&tmp);

  //------------------------------------------------
  
  label = gft::Scene32::Create(I1);
  gft::Scene32::Fill(label, NIL);
  
  I2 = gft::Scene32::CloseHolesSliceBySlice_z(I1);
  for(p = 0; p < I1->n; p++){
    I2->data[p] = I2->data[p] - I1->data[p];
  }


  if( gft::AttributeList::GetAttributeValue(al, (char *)"T", value) )
    T = atoi(value);
  else
    T = ROUND(1.2*gft::Scene32::Otsu(I2)); //k=1.2
  
  I3 = gft::Scene32::Threshold(I2, T, INT_MAX); 

#ifdef APPDEBUG
  if(gft::Scene8::GetMaxVal(I3) == 0)
    printf("Warning: empty selection when applying thresholding.\n");
  gft::Scene8::Write(I3, (char *)"./debug/I3_1st.scn");
#endif

  A = gft::AdjRel3::Spheric(sqrtf(3.0 + 0.000001));
  SelectLargestCompFitForLungs(I3, A);

  float d = (I1->dx + I1->dy + I1->dz)/3.0;
  float gamma = 2.5/d; //2.5 mm

  It = gft::Scene8::CloseBin(I3, gamma, A);
  I4 = gft::Scene8::CloseHoles(It);
  gft::Scene8::Destroy(&It);
  
  float gamma_i = 18.75/d; //18.75 mm
  float gamma_e = 2.5/d; //2.5 mm

  gft::sSet *S = NULL;
  
  I4_prime = gft::Scene8::DilateBin(I4, &S, gamma_e, A);
  gft::Set::Destroy(&S);

  //-----------------------
  GetLungsSeeds(I4, A, label, gamma_i);
  //-----------------------
  
  int *SA = NULL;
  SA = gft::Scene32::GetMarkers(label, A);
  TDG_C = GeodesicDistanceTransform(I4, A, SA);
  Cmax = gft::Scene32::GetMaxVal(TDG_C);
  free(SA);
  
  float alpha = 0.7; //0.4; //0.5; //0.6; //0.7; // 0.9
  I_St =  gft::Scene32::Threshold(TDG_C, ROUND(alpha*Cmax), INT_MAX);
  SelectComponentFromTop(I_St, A, label);
  gft::Scene8::Destroy(&I_St);
  
  for(p = 0; p < label->n; p++)
    if(I4_prime->data[p] == 0)
      label->data[p] = 0;


  //Output resolution will be the same as the input volume unless
  //a ground truth is provided in a different resolution.
  gft::sScene32 *ref = scn;
  if(gtruth_left != NULL && output_type == 0)
    ref = gtruth_left;
  else if(gtruth_right != NULL && output_type == 0)  
    ref = gtruth_right;
  //-------------------
  
  if(output_type == 0){
    I5 = ComputeGradientImage(I2, A);
  
    label_2nd = gft::Scene32::Clone(label);
    pred = LabelingByIFTPeak_1st(I5, A, label);
    LabelingByIFTPeak_2nd(I1, A, pred, label, label_2nd);
    gft::Scene32::Destroy(&pred);

    UndoFixOrientation(label,     flippedX, flippedY, flippedZ);
    UndoFixOrientation(label_2nd, flippedX, flippedY, flippedZ);
    segm_1st = gft::Scene32::ScaleLabel(label,     ref, gft::linear);//gft::none
    segm_2nd = gft::Scene32::ScaleLabel(label_2nd, ref, gft::linear);
    seeds = NULL;
  }
  else{ //output_type == 1
    I5 = NULL;
    label_2nd = NULL;
    segm_1st = NULL;
    segm_2nd = NULL;

    UndoFixOrientation(label, flippedX, flippedY, flippedZ);
    seeds = gft::Scene32::ScaleLabel(label, ref, gft::linear);
  }
  //-----------------------------------------  
  end = clock();
  totaltime = ((double)(end - start)) / CLOCKS_PER_SEC;
  printf("Time: %f sec\n", totaltime);

  
#ifdef APPDEBUG
  if(I5 != NULL)
    gft::Scene32::Write(I5, (char *)"./debug/I5.nii.gz");
  gft::Scene8::Write(I4_prime, (char *)"./debug/I4_prime.scn");
  gft::Scene8::Write(I4, (char *)"./debug/I4.scn");
  gft::Scene8::Write(I3, (char *)"./debug/I3.scn");
  gft::Scene32::Write(I2, (char *)"./debug/I2.scn");
  gft::Scene32::Write(I2, (char *)"./debug/I2.nii.gz");
  gft::Scene32::Write(I1, (char *)"./debug/I1.nii.gz");  
  gft::Scene32::Write(label, (char *)"./debug/label.nii.gz");
  if(label_2nd != NULL)
    gft::Scene32::Write(label_2nd, (char *)"./debug/label_2nd.nii.gz");
  gft::Scene32::Write(TDG_C, (char *)"./debug/TDG_C.nii.gz");
#endif

  if(seeds != NULL){
    gft::Scene32::Write(seeds, (char *)"./out/seeds.nii.gz");
    SaveSeedsByLabel((char *)"./out/seeds_L.txt", seeds, 2);
    SaveSeedsByLabel((char *)"./out/seeds_R.txt", seeds, 1);
  }
  
  //if(segm_1st != NULL)
  //  gft::Scene32::Write(segm_1st, (char *)"./out/segm_1st.nii.gz");
  if(segm_2nd != NULL)
    gft::Scene32::Write(segm_2nd, (char *)"./out/segm_altis.nii.gz");

  if(output_type == 0){
    if( gtruth_left != NULL ){
      printf("-------Left lung-------\n");
      //dice = Dice(segm_1st, 2, gtruth_left);
      //printf("1st: Dice = %lf\n", dice);
      
      dice = Dice(segm_2nd, 2, gtruth_left);
      printf("Dice = %lf\n", dice);    
    }
    if( gtruth_right != NULL ){
      printf("-------Right lung------\n");    
      //dice = Dice(segm_1st, 1, gtruth_right);
      //printf("1st: Dice = %lf\n", dice);
      
      dice = Dice(segm_2nd, 1, gtruth_right);
      printf("Dice = %lf\n", dice);
    }
  }

  gft::AttributeList::Destroy(&al);
  if(gtruth_left != NULL)
    gft::Scene32::Destroy(&gtruth_left);
  if(gtruth_right != NULL)
    gft::Scene32::Destroy(&gtruth_right);  
  if(I5 != NULL)
    gft::Scene32::Destroy(&I5);
  gft::Scene8::Destroy(&I4_prime);
  gft::Scene8::Destroy(&I4);
  gft::Scene8::Destroy(&I3);
  gft::Scene32::Destroy(&I2);
  gft::Scene32::Destroy(&I1);
  gft::Scene32::Destroy(&scn);
  gft::Scene32::Destroy(&label);
  if(label_2nd != NULL)
    gft::Scene32::Destroy(&label_2nd);
  if(segm_1st != NULL)
    gft::Scene32::Destroy(&segm_1st);
  if(segm_2nd != NULL)
    gft::Scene32::Destroy(&segm_2nd);
  gft::Scene32::Destroy(&TDG_C);
  if(seeds != NULL)
    gft::Scene32::Destroy(&seeds);
  gft::AdjRel3::Destroy(&A);
  return 0;
}
