// mexfun.cpp : 定义控制台应用程序的入口点。
//

#include "mex.h"
#include "mat.h"

/*
 * mexfunobj.c : MEX-file code computing the objective function of the
 *             MATLAB
 *
 *              
 * Author: TangZedong
 * $Version: 0.0.1 $
 */
double sum(const mxArray *A);
double avg(const mxArray *A);
double var(const mxArray *A,size_t i,size_t j);    /*计算方差函数*/
double min(double a,double b);    /*取最小值*/
void readmat(mxArray **img);    /*读取图像*/
void
mexFunction( int nlhs, mxArray *plhs[],
	     int nrhs, const mxArray *prhs[] )
{
	//参数处理
  (void) nlhs;     /* unused parameters */
  (void) plhs;

  if(nrhs==0)
    mexErrMsgTxt( "Function 'mexfunobj' not defined for variables of class 'double'.\n");

  else if(!mxIsDouble(prhs[0])) {
    const char str[]="Function 'mexfunobj' not defined for variables of class";
    char errMsg[100];
    
    sprintf(errMsg,"%s '%s'\n",str,mxGetClassName(prhs[0]));
    mexErrMsgTxt(errMsg);
  }
  //函数主体
  else {
    /* causes MATLAB to execute the string as an expression 
     * or statement the biggest limitation is that it won't 
     * evaluate any left-hand arguments except 'ans'
     */
	/*读取img矩阵*/
    const int dbon=1;
	 mxArray *img;
	 readmat(&img);
     
     if(dbon==1)
        printf("Pass p1");
     
	 const mxArray *u=prhs[0];
	 size_t m=mxGetM(img);
	 size_t n=mxGetN(img);

          if(dbon==1)
        printf("Pass p2");
     
	 if(mxGetN(u)!=m*n)
	 {
		 mexErrMsgTxt("error:number of img's elements must be equal to u's");
		 return;
	 }
     
          if(dbon==1)
        printf("Pass p3");
     
	 mxArray *u2=mxCreateDoubleMatrix(1,mxGetN(u),mxREAL);
	 mxArray *c1=mxCreateDoubleMatrix(m,m,mxREAL);
	 mxArray *c2=mxCreateDoubleMatrix(m,n,mxREAL);
	 double *pu=mxGetPr(u2);
     double *pc1=mxGetPr(c1);
     double *pc2=mxGetPr(c2);
     
          if(dbon==1)
        printf("Pass p4");
     
	/*计算C矩阵*/
	 /*C边缘一像素边框和内边缘相等*/
     
     //计算第一类C矩阵
	 for(size_t i=1;i<m-1;i++)
	 {
		 for(size_t j=1;j<n-1;j++)
		 {
			 *(pc1+j*m+i)=var(img,i,j);
		 }
	 }
     
     //u~=1-u
	 for(int i=0;i<m;i++)
	 {
		 for(int j=0;j<n;j++)
		 {
			 *(pu+i+j*m)=1-*(pu+i+j*m);
		 }
	 }
     
     //计算第二类C矩阵
	 for(size_t i=1;i<m-1;i++)
	 {
		 for(size_t j=1;j<n-1;j++)
		 {
			 *(pc2+j*m+i)=var(img,i,j);
		 }
	 }
     
          if(dbon==1)
        printf("Pass p5");
	/*计算J函数值*/
	/*用户空间清理*/
	//mxDestoryArray(c1);
	//mxDestoryArray(c2);
     //free(c1);
     //free(c2);
     plhs[0]=c1;
     plhs[1]=c2;
  }
}


void readmat(mxArray **img)
{
	MATFile *fmat;
	const char *name = "diffimg";
	fmat=matOpen("test.mat","r");
	if(fmat==NULL){
		mexErrMsgTxt( "make sure you have test.mat in current dir");
		return;
	}
	else
	{
		*img=matGetVariable(fmat,name);
	}
	matClose(fmat);
}

double min(double a,double b)    /*取最小值*/
{
	if(a>b)
		return b;
	else
		return a;
}

double sum(const mxArray *A)
{
	size_t m=mxGetM(A);
	size_t n=mxGetN(A);
	double *p=mxGetPr(A);
	double sum=0;
    
	for(size_t j=0;j<n;j++)
	{
		for(size_t i=0;i<m;i++)
		{
			sum += *(p+i+j*m);
		}
	}
    
	return sum;
}

double avg(const mxArray *A)
{
	size_t m=mxGetM(A);
	size_t n=mxGetN(A);
	return sum(A)/(m*n);
}

double var(const mxArray *A,size_t i,size_t j)
{
	mxArray *loc = mxCreateDoubleMatrix(3,3,mxREAL);
	double *ploc = mxGetPr(loc);
	size_t ml=mxGetM(loc),nl=mxGetN(loc);
	double *pimg = mxGetPr(A);
	size_t mi=mxGetM(A),ni=mxGetN(A);

	for(size_t x=i-1;x<=i+1;x++)
	{
		for(size_t y=j-1;y<=j+1;y++)
		{
			*(ploc+(x-i+1)+(y-j+1)*ml)=*(pimg+x+y*mi);
		}
	}

	double a=avg(loc);

	for(size_t x=0;x<3;x++)
		for(size_t y=0;y<3;y++)
			*(ploc+x+y*ml)=(*(ploc+x+y*ml)-a)*(*(ploc+x+y*ml)-a);

	double s = sum(loc);
	mxDestroyArray(loc);

	return s/(a*a);
}