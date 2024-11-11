
#include "gft.h"


int main(int argc, char **argv) {
  char file_image[612];
  char file_gtruth_L[612];
  char file_gtruth_R[612];  
  char command[2048];
  int i, T = 0;
  char buffer[1024];
  char value[128];
  char T_str[128];
  gft::sAttributeList *al;

  if (argc < 2) {
    fprintf(stdout, "usage:\n");
    fprintf(stdout, "altis_roift <volume> [T=value] [left=file] [right=file]\n");
    fprintf(stdout, "Optional parameters:\n");
    fprintf(stdout, "\tT................... threshold integer value\n");
    fprintf(stdout, "\t                     (if not specified Otsu is used).\n");
    fprintf(stdout, "\tleft................ ground truth for left lung.\n");
    fprintf(stdout, "\tright............... ground truth for right lung.\n");
    return 0;
  }

  strcpy(file_image, argv[1]);

  buffer[0] = '\0';
  for(i = 2; i < argc; i++){
    strcat(buffer, (char *)" ");
    strcat(buffer, argv[i]);
  }
  al = gft::AttributeList::Create(buffer);

  T_str[0] = '\0';
  if( gft::AttributeList::GetAttributeValue(al, (char *)"T", value) ){
    T = atoi(value);
    sprintf(T_str, "T=%d", T);
  }

  file_gtruth_L[0] = '\0';
  file_gtruth_R[0] = '\0';
  if( gft::AttributeList::GetAttributeValue(al, (char *)"left", value) )
    sprintf(file_gtruth_L, "%s", value);
  if( gft::AttributeList::GetAttributeValue(al, (char *)"right", value) )
    sprintf(file_gtruth_R, "%s", value);

  printf("-----Seeds by altis----\n");
  sprintf(command, "./altis %s 1 %s", file_image, T_str);
  system(command);

  printf("-------Left lung-------\n");
  sprintf(command, "./relaxed_oift %s ./out/seeds_L.txt -0.1 10 90 segm_left_lung %s",
	  file_image, file_gtruth_L);
  system(command);

  printf("-------Right lung------\n");
  sprintf(command, "./relaxed_oift %s ./out/seeds_R.txt -0.1 10 90 segm_right_lung %s",
	  file_image, file_gtruth_R);
  system(command);
  
  gft::AttributeList::Destroy(&al);
  
  return 0;
}

