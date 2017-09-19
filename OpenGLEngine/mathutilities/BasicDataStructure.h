#ifndef __WXJ_DATASTRUCTURE_H__
#define __WXJ_DATASTRUCTURE_H__

//#include <math.h>
#include "..\app.h"

///Generic vecoter
class XVec{
private:
	int size;
	double *data;
public:
	XVec(int size){
		if(size <= 0) throw;
		data = new double[size];
		for(int i = 0; i < size; i++) data[i] = 0;
		this->size = size;
	}
	~XVec(){ delete data; }

	inline double operator()(int i) const{ return data[i]; }
	inline double& operator()(int i){ return data[i]; }
};

///2d Vector
class XVec2{
private:
	double elem[2];
public:
	XVec2(){
		elem[0] = 0;
		elem[1] = 0;
	}
	XVec2(double x, double y){
		elem[0] = x;
		elem[1] = y;
	}

	inline double operator()(int i) const { return (elem[i]); }
	inline double& operator()(int i){ return (elem[i]); }

	inline XVec2 operator + (XVec2 x){return XVec2(elem[0] + x(0), elem[1] + x(1));}
	inline XVec2 operator - (XVec2 x){return XVec2(elem[0] - x(0), elem[1] - x(1));}

	inline double Dot(XVec2 &x){return elem[0] * x(0) + elem[1] * x(1);}

	inline double Normalize(){
		double norm = sqrt(elem[0]*elem[0] + elem[1]*elem[1]);
		if(norm != 0){
			elem[0] = elem[0] / (double)norm;
			elem[1] = elem[1] / (double)norm;
		}
		return norm;
	}

	inline friend XVec2 operator * (XVec2 x, double);
	inline friend XVec2 operator * (double, XVec2 x);

	inline double GetLength(){//not tested
		return sqrt(elem[0]*elem[0] + elem[1]*elem[1]);
	}
};

inline XVec2 operator * (XVec2 x,double y){ return XVec2(x(0) * y, x(1) * y);}
inline XVec2 operator * (double y, XVec2 x){return XVec2(x(0) * y, x(1) * y);}

///3d vector
class XVec3{
private:
	double elem[3];
public:
	XVec3(){
		elem[0] = 0;
		elem[1] = 0;
		elem[2] = 0;
	}
	XVec3(double x, double y, double z){
		elem[0] = x;
		elem[1] = y;
		elem[2] = z;
	}

	void Zero(){
		elem[0] = 0;
		elem[1] = 0;
		elem[2] = 0;
	}

	inline double operator()(int i) const { return (elem[i]); }
	inline double& operator()(int i){ return (elem[i]); }

	inline XVec3 operator + (XVec3 x){ return XVec3(elem[0] + x(0), elem[1] + x(1), elem[2] + x(2)); }
	inline XVec3 operator - (XVec3 x){ return XVec3(elem[0] - x(0), elem[1] - x(1), elem[2] - x(2)); }

	inline double Dot(XVec3 &x){ return elem[0] * x(0) + elem[1] * x(1) + elem[2] * x(2); }
	inline XVec3 Cross(XVec3 &x){ return XVec3(elem[1]*x(2) - elem[2]*x(1), elem[2]*x(0) - elem[0]*x(2), elem[0]*x(1) - elem[1]*x(0)); }

	inline double Normalize(){
		double norm = sqrt(elem[0]*elem[0] + elem[1]*elem[1] + elem[2]*elem[2]);
		if(norm != 0){
			elem[0] = elem[0] / (double)norm;
			elem[1] = elem[1] / (double)norm;
			elem[2] = elem[2] / (double)norm;
		}
		return norm;
	}

	inline double GetLength(){//not tested
		return sqrt(elem[0]*elem[0] + elem[1]*elem[1] + elem[2]*elem[2]);
	}

	inline friend XVec3 operator * (XVec3 x, double);
	inline friend XVec3 operator * (double, XVec3 x);
};
inline XVec3 operator * (XVec3 x,double y){ return XVec3(x(0) * y, x(1) * y, x(2) * y); }
inline XVec3 operator * (double y, XVec3 x){ return XVec3(x(0) * y, x(1) * y, x(2) * y); }


///general sparse matrix
///inplemented in sparsematrix.h

///22d matrix
class XMat22{
private:
	double elem[4];
public:
	XMat22(){
		elem[0] = 0; elem[1] = 0;
		elem[2] = 0; elem[3] = 0;
	}
	XMat22(double x11, double x12, double x21, double x22){
		elem[0] = x11; elem[1] = x12;
		elem[2] = x21; elem[3] = x22;
	}

	inline double operator()(int row, int col) const { return (elem[row*2 + col]); }
	inline double& operator()(int row, int col){ return (elem[row*2 + col]); }

	inline XMat22 operator + (XMat22 x){ return XMat22(elem[0] + x(0,0), elem[1] + x(0,1),  elem[2] + x(1,0), elem[3] + x(1,1)); }
	inline XMat22 operator - (XMat22 x){ return XMat22(elem[0] - x(0,0), elem[1] - x(0,1),  elem[2] - x(1,0), elem[3] - x(1,1)); }

	inline XVec2 operator * (XVec2 x){ 
		return XVec2(elem[0] * x(0) + elem[1] * x(1),
					 elem[2] * x(0) + elem[3] * x(1)); 
	}

	inline XMat22 operator * (XMat22 x){
		return XMat22(
			elem[0] * x(0,0) + elem[1] * x(1,0), 
			elem[0] * x(0,1) + elem[1] * x(1,1),  
			elem[2] * x(0,0) + elem[3] * x(1,0), 
			elem[2] * x(0,1) + elem[3] * x(1,1));
	}

	inline XMat22 Transpose(){
		return XMat22(elem[0], elem[2],
					  elem[1], elem[3]);
	}

	inline friend XMat22 operator * (XMat22 x, double);
	inline friend XMat22 operator * (double, XMat22 x);
};
inline XMat22 operator * (XMat22 x,double y){ return XMat22(x(0,0) * y, x(0,1) * y, x(1,0) * y, x(1,1) * y); }
inline XMat22 operator * (double y, XMat22 x){ return XMat22(x(0,0) * y, x(0,1) * y, x(1,0) * y, x(1,1) * y); }


///33d matrix
class XMat33{
private:
	double elem[9];
public:
	XMat33(){
		for(int i = 0; i< 9;i++) elem[i] = 0;
	}
	XMat33(double x11, double x12, double x13, 
			  double x21, double x22, double x23, 
			  double x31, double x32, double x33){
		elem[0] = x11; elem[1] = x12; elem[2] = x13;
		elem[3] = x21; elem[4] = x22; elem[5] = x23;
		elem[6] = x31; elem[7] = x32; elem[8] = x33;
	}

	inline double operator()(int row, int col) const { return (elem[row*3 + col]); }
	inline double& operator()(int row, int col){ return (elem[row*3 + col]); }

	inline XMat33 operator + (XMat33 x){ 
		return XMat33(elem[0] + x(0,0), elem[1] + x(0,1),  elem[2] + x(0,2),
					  elem[3] + x(1,0), elem[4] + x(1,1),  elem[5] + x(1,2),
					  elem[6] + x(2,0), elem[7] + x(2,1),  elem[8] + x(2,2)); 
	}
	inline XMat33 operator - (XMat33 x){ 
		return XMat33(elem[0] - x(0,0), elem[1] - x(0,1),  elem[2] - x(0,2),
					  elem[3] - x(1,0), elem[4] - x(1,1),  elem[5] - x(1,2),
					  elem[6] - x(2,0), elem[7] - x(2,1),  elem[8] - x(2,2)); 
	}

	inline XVec3 operator * (XVec3 x){ 
		return XVec3(elem[0] * x(0) + elem[1] * x(1) + elem[2] * x(2),
					 elem[3] * x(0) + elem[4] * x(1) + elem[5] * x(2),
					 elem[6] * x(0) + elem[7] * x(1) + elem[8] * x(2)); 
	}

	inline XMat33 operator * (XMat33 x){
		XMat33 r;
		for(int i = 0; i < 3; i++)
			for(int j = 0; j < 3; j++)
				for(int k = 0; k < 3; k++)
					r(i,j) = r(i,j) + elem[i*3+k] * x(k,j);
		return r;
	}	

	inline XMat33 Transpose(){
		return XMat33(elem[0], elem[3], elem[6],
					  elem[1], elem[4], elem[7],
					  elem[2], elem[5], elem[8]);
	}

	inline void Identity(){
		elem[0] = 1; elem[1] = 0; elem[2] = 0;
		elem[3] = 0; elem[4] = 1; elem[5] = 0;
		elem[6] = 0; elem[7] = 0; elem[8] = 1;
	}

	inline friend XMat33 operator * (XMat33 x, double);
	inline friend XMat33 operator * (double, XMat33 x);
};
inline XMat33 operator * (XMat33 x,double y){ 
	return XMat33(x(0,0) * y, x(0,1) * y, x(0,2) * y,
				  x(1,0) * y, x(1,1) * y, x(1,2) * y,
				  x(2,0) * y, x(2,1) * y, x(2,2) * y); 
}
inline XMat33 operator * (double y, XMat33 x){ 
	return XMat33(x(0,0) * y, x(0,1) * y, x(0,2) * y,
				  x(1,0) * y, x(1,1) * y, x(1,2) * y,
				  x(2,0) * y, x(2,1) * y, x(2,2) * y); 
}



///44d matrix
//need?

#endif