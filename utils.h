// utils.h

#ifndef _UTILS_H
#define _UTILS_H

#define SIGN(a)    ((a)<0 ? -1 : 1)
#define ZSIGN(a)   ((a)<0 ? -1 : ((a)==0 ? 0 : 1) )
#define MIN(a,b)   ((a)<(b)? (a) : (b))
#define MAX(a,b)   ((a)>(b)? (a) : (b))
#define ABS(a)		((a)<0? -(a) : (a))
#define ZSign(x)	ZSIGN(x)

#define  ONE			255
#define  ZERO			0
#define  Error			-1
#define  FixedLevel    256
#define PI		3.1415926535897932384626433832795f
#define INF     1.0e10

typedef struct tagVECTOR
{
	float x;
	float y;
	tagVECTOR operator = (CPoint p)
	{x=float(p.x); y=float(p.y); return *this; }
}VECTOR;

void  **fspace_2d( int row, int col, int wide);
void ***fspace_3d( int X, int Y, int Z, int wide);
void ffree(void** array);
void ffree(void*** array);

#endif
