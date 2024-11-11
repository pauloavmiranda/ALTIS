#ifndef _GFT_DICOM_H_
#define _GFT_DICOM_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <errno.h>
#include <fcntl.h>

#include "gft_image32.h"
#include "gft_scene32.h"

namespace gft{

  typedef struct _dicom_field {
    int  group;
    int  key;
    int  len;
    char data[64];
    struct _dicom_field *next;
  } DicomField;

  typedef struct _dicom_slice {
    sImage32 *img;
    float  zslice;
    int    seqnum;
    float  thickness;
    char   equip[256];
    char   modality[256];
  } DicomSlice;

  typedef struct _imglist {
    DicomSlice *I;
    struct _imglist *next;
  } ImageList;


  /* data structures */
  
  DicomField * DFNew();
  void         DFKill(DicomField *df);
  
  DicomSlice  *DicomSliceNew(int w, int h);
  void         DicomSliceKill(DicomSlice *i);
  
  void         DicomSliceDump16(FILE *f, DicomSlice *S);
  int          DicomSliceCmp(DicomSlice *a, DicomSlice *b);


  /* single file processing */
  
  DicomSlice  *ReadDicomSlice(char *filename, FILE *rf, int report);
  int          ReadDicomField(FILE *f, DicomField **df);
  DicomField  *DicomGetField(DicomField *H, int group, int key);
  int          ConvInt32(unsigned char *x);
  int          ConvInt16(unsigned char *x);
  
  /* directory sorting and volume assembling */

  //Read all ".dcm" files from a directory.
  gft::sScene32 *DicomDirectoryConvert(char *inputsample,
				       int report);

  //Read all ".dcm" files from a directory having the same basename and numbering in the range [first, last].
  gft::sScene32 *DicomDirectoryConvert(char *inputsample, //char *outputscn, 
				       int first,int last, int report);
  
  //void         DicomDirectoryInfo(char *inputsample, int extract);
  
  void         DicomNameSplit(char *isample, char *odirpath, int odsz, 
			      char *obasename);
  //int          DicomSeqInfo(int index, Image2D **vol,int first,int last,int extract);
  //void         DicomExtractSeq(int index,int b,int e);
  
  void         FindSampleFileIfDirectory(char *sample);
  
  /* utilities */
  int   slicenumber(const char *s);


  typedef struct _vr {
    char id[3];
    int  maxl;
    int  fixed;
    int  longlen;
  } VR;
  

} //end gft namespace

#endif
