// utils.cpp

#include "utils.h"
#include "stdafx.h"
typedef unsigned long ulong;

void **fspace_2d( int row, int col, int wide)
{
	void **  array;
    char * bArray;
  	register unsigned long i;
    typedef unsigned long ulong;
    ulong WCol,nbytes;

    WCol=(ulong)col*(ulong)wide;
    nbytes=(ulong)row*WCol;

  	array= (void **) calloc(row,sizeof(void*));
  	if (array==NULL)   	return(NULL);

    bArray=(char*) calloc(nbytes,sizeof(char));
    if(bArray==NULL)   return NULL;

  	for (i=0;i<(ulong)row;i++) array[i]=(void*)(bArray+i*WCol);

  	return( array);
/*	void ** pY;
	register unsigned long j;
	
	pY = (void **) calloc(Y, sizeof(void**)); ASSERT(pY);
	char *pXY = (char *) calloc(Y*X, sizeof(char)*wide); ASSERT(pXY);
	for(j=0; j<(ulong)Y; j++) pY[j] = pXY+j*X;
	return(pY);*/
}

void * fspace_3d( int X, int Y, int Z, int wide)
{
	void ***  pZ;
	void **  pYZ;
  	unsigned long j,k;

	pZ = (void ***) calloc(Z,sizeof(void**)); ASSERT(pZ);
	pYZ = (void **) calloc(Z*Y,sizeof(void*)); ASSERT(pYZ);
	char *pXYZ = (char *) calloc(Z*Y*X, sizeof(char)*wide); ASSERT(pXYZ);
	for(k=0; k<(ulong)Z; k++) pZ[k] = pYZ+k*Y;
	for(j=0; j<(ulong)(Z*Y); j++) pYZ[j] = (void*)(pXYZ+j*X*wide);
	return(pZ);
}

void ffree(void** array)
{
      free(array[0]);
      free(array);
}

void ffree(void*** array)
{
	  free(array[0][0]);
	  free(array[0]);
	  free(array);
}