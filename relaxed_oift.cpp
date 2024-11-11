
#include "gft.h"

#include <algorithm>
#include <vector>

gft::sScene32 *get_dilation_border(gft::sScene32 *scn, float radius_sphere) {
    gft::sScene32 *dil = NULL;
    int p, q, n, max, i;
    gft::Voxel u, v;
    gft::sAdjRel3 *A = gft::AdjRel3::Spheric(radius_sphere);
    n = scn->n;
    dil = gft::Scene32::Create(scn);
    gft::Scene32::Fill(dil, 0);
    for (p = 0; p < n; p++) {
        if (scn->data[p] == 1) {
            u.c.x = gft::Scene32::GetAddressX(scn, p);
            u.c.y = gft::Scene32::GetAddressY(scn, p);
            u.c.z = gft::Scene32::GetAddressZ(scn, p);
            for (i = 1; i < A->n; i++) {
                v.v = u.v + A->d[i].v;
                if (gft::Scene32::IsValidVoxel(scn, v)) {
                    q = gft::Scene32::GetVoxelAddress(scn, v);
                    if (scn->data[q] == 0) {
                        dil->data[q] = 1;
                    }
                }
            }
        }
    }
    gft::AdjRel3::Destroy(&A);
    return dil;
}

struct value_position {
    int value;
    int position;
    value_position(int value, int position) : value(value), position(position) {}
};

void condition_percentil(gft::sScene32 *scn, gft::sScene32 *border, float percentil) {

    int p, q, n;
    std::vector<value_position> border_values;
    n = scn->n;
    for (p = 0; p < n; p++) {
        if (border->data[p] == 1)
            border_values.emplace_back(scn->data[p], p);
    }
    std::stable_sort(border_values.begin(), border_values.end(), [](const value_position &a, const value_position &b) { return a.value < b.value; });
    float index_percentile = (percentil / 100.0) * (border_values.size() - 1);
    int value_position = static_cast<int>(index_percentile);
    int num_voxels_border = border_values.size();
    for (p = value_position + 1; p < num_voxels_border; p++) {
        q = border_values[p].position;
        border->data[q] = 0;
    }
}

void dilation_conditional(gft::sScene32 *scn, gft::sScene32 *label, float radius_sphere, float percentile) {
    gft::sScene32 *dilation_border;
    int p, n;

    dilation_border = get_dilation_border(label, radius_sphere);
    // gft::Scene32::Write(dilation_border, (char *)"label_dilation_border.nii.gz");
    condition_percentil(scn, dilation_border, percentile);
    // gft::Scene32::Write(dilation_border, (char *)"label_dilation_border_percentile.nii.gz");

    n = scn->n;
    for (p = 0; p < n; p++) {
        if (dilation_border->data[p] == 1)
            label->data[p] = 1;
    }

    gft::Scene32::Destroy(&dilation_border);
}






int *LoadSeeds(char *filename,
	       gft::sScene32 *scn,
	       gft::sScene32 *gtruth,
	       gft::sScene32 *label){
  int nseeds,j,i,p,x,y,z,id,lb;
  int n_inconsistent_seeds = 0;
  int *S = NULL;
  FILE *fp;
  fp = fopen(filename, "r");
  if (fp == NULL) {
    printf("Error reading seeds.\n");
    exit(1);
  }
  fscanf(fp, " %d", &nseeds);
  S = (int *)calloc((nseeds + 1), sizeof(int));
  S[0] = nseeds;
  j = 0;
  for (i = 0; i < nseeds; i++) {
    fscanf(fp, " %d %d %d %d %d", &x, &y, &z, &id, &lb);
    if (gft::Scene32::IsValidVoxel(scn, x, y, z)) {
      p = gft::Scene32::GetVoxelAddress(scn, x, y, z);
      j++;
      S[j] = p;
      label->data[p] = lb;
      if(gtruth != NULL){
	if((gtruth->data[p] == 0 && lb != 0) ||
	   (gtruth->data[p] != 0 && lb == 0)){
	  n_inconsistent_seeds++;
	}
      }
    }
  }
  S[0] = j;
  fclose(fp);
  
  if(gtruth != NULL)
    printf("inconsistent_seeds: %d (%f\%)\n",
	   n_inconsistent_seeds,
	   (100.0*n_inconsistent_seeds)/nseeds);
  return S;
}



void SaveSeeds(char *filename,
	       int *S,
	       gft::sScene32 *label){
  int i,p,lb,x,y,z,id;
  FILE *fp;
  fp = fopen(filename, "w");
  if (fp == NULL) {
    printf("Error writing seeds.\n");
    exit(1);
  }
  fprintf(fp, " %d", S[0]);
  for(i = 1; i <= S[0]; i++){
    p = S[i];
    lb = label->data[p];
    x = gft::Scene32::GetAddressX(label, p);
    y = gft::Scene32::GetAddressY(label, p);
    z = gft::Scene32::GetAddressZ(label, p);
    id = 1;
    if(gft::Scene32::IsValidVoxel(label, x, y, z))
      fprintf(fp, " %d %d %d %d %d", x, y, z, id, lb);
  }
  fclose(fp);
}



int main(int argc, char **argv) {
    gft::sScene32 *scn, *fscn, *label, *W, *Wx, *Wy, *Wz;
    gft::sScene32 *gtruth = NULL, *slabel = NULL;
    gft::sAdjRel3 *A;
    char fileseeds[512];
    char filename[512];
    clock_t end, start;
    double totaltime, dice;
    int *S;
    int p, i, j, Imin;
    int niter = 10;
    float pol = -0.1;
    int percentile = 90;

    if (argc < 7) {
        fprintf(stdout, "usage:\n");
        fprintf(stdout, "relaxed_oift <volume> <file_seeds> <pol> <niter> <percentile> <out> [ground truth]\n");
        fprintf(stdout, "\t pol.......... Boundary polarity. It can be in the range [-1.0, 1.0]\n");
        fprintf(stdout, "\t niter........ Number of iterations of the relaxation procedure.\n");
	fprintf(stdout, "\t percentile... Percentage of voxels for conditional acquisition in dilation in the range [0, 100]\n");
	fprintf(stdout, "\t out.......... Base name of the output file with the segmentation result.\n");
        exit(0);
    }

    A = gft::AdjRel3::Spheric(1.0);
    scn = gft::Scene32::Read(argv[1]);

    label = gft::Scene32::Create(scn);
    gft::Scene32::Fill(label, NIL);

    pol = atof(argv[3]);
    niter = atoi(argv[4]);
    percentile = atoi(argv[5]);
    //printf("pol: %f, niter: %d, percentile: %d\n", pol, niter, percentile);

    if (argc == 8){
      gtruth = gft::Scene32::Read(argv[7]);
    }

    strcpy(fileseeds, argv[2]);
    S = LoadSeeds(fileseeds, scn, NULL, label);
    
    Imin = gft::Scene32::GetMinimumValue(scn);
    if (Imin < 0) {
      for (p = 0; p < scn->n; p++)
	scn->data[p] += (-Imin);
    }

    
    start = clock();

    fscn = gft::Scene32::GaussianBlur(scn);
    gft::Scene32::Destroy(&scn);
    scn = fscn;
    fscn = gft::Scene32::GaussianBlur(scn);
    gft::Scene32::Destroy(&scn);
    
    gft::ift::OIFT(A, fscn, pol * 100.0, S, label);

    //-----------------
    // - radius sphere adjacency
    // - percentile, % of dark intensities in dilation to be added as object
    if(percentile > 0)
      dilation_conditional(fscn, label, 1 /*radius sphere adj */, percentile /*percentile*/);
    //-----------------
    if(niter > 0)
      gft::ift::ORelax_1(A, fscn, pol * 100.0, S, label, niter);
    gft::Scene32::Destroy(&fscn);

    end = clock();
    totaltime = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Time: %f sec\n", totaltime);

    sprintf(filename, (char *)"./out/%s.nii.gz", argv[6]);
    gft::Scene32::Write(label, filename);

    if (argc == 8){
      slabel = gft::Scene32::ScaleLabel(label, gtruth, gft::linear);
      dice = gft::Scene32::DiceSimilarity(slabel, gtruth);
      printf("Dice: %lf\n", dice);
    }
   
    free(S);
    if(gtruth != NULL)
      gft::Scene32::Destroy(&gtruth);
    if(slabel != NULL)
      gft::Scene32::Destroy(&slabel);
    gft::Scene32::Destroy(&label);
    gft::AdjRel3::Destroy(&A);
    return 0;
}
