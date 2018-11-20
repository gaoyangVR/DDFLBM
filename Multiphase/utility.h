#ifndef UTILITY_H
#define UTILITY_H
#include<vector_types.h>

// #define NX 24
// #define NY 24
// #define NZ 96
#define MAXITER 200
#define FLIP_ALPHA 0.95f
#define M_PI       3.14159265358979323846
const float DEGtoRAD = 3.1415926f / 180.0f;

#define TYPEFLUID 0
#define TYPEAIR 1
#define TYPEBOUNDARY 2
#define TYPEVACUUM 3
#define TYPEAIRSOLO 4
#define TYPESOLID 5
#define TYPECNT 6

typedef unsigned int uint;
#define CELL_UNDEF 0xffffffff
#define  NTHREADS 32
#define UNDEF_TEMPERATURE -10000.0f

struct FlipConstant{
	int gnum;
	int3 gvnum;
	float samplespace;
	float dt;
	float3 gravity;
	float3 gmin, gmax, cellsize;
	float m0;
	float airm0;
	float waterrho, solidrho;

	float pradius;
	float3	triHashSize, triHashRes;		//triHashSize是HASH网格的大小;  triHashRes是每一个维度上有几个HASH网格,程序执行过程中不再变化
	float3 t_min, t_max;
	int triHashCells;			//预留的hash数组大小，程序执行过程中不再变化
	//for SPH-like part
	float poly6kern, spikykern, lapkern;

	//***************************LBM*****************************
	float Pr, Ra, R;
	float T0;
	float p0;
	float tau_f, tau_h;
	float niu;
	float cv;
	float wf, wh;
	float total_E, T_heat;
	//****************************************************************
};

//创建一个方便转换1维与3维数组的数据结构
struct farray{
	float* data;
	int xn, yn, zn;
	farray();
	void setdim(int _xn, int _yn, int _zn){ xn = _xn, yn = _yn, zn = _zn; }

	__host__ __device__ inline float &operator ()(int i, int j, int k)
	{
		return data[i*yn*zn + j*zn + k];
	}
	__host__ __device__ inline float &operator ()(int i)
	{
		return data[i];
	}
	__host__ __device__ inline float &operator [](int i)
	{
		return data[i];
	}
};

//创建一个方便转换1维与3维数组的数据结构
struct charray{
	char* data;
	int xn, yn, zn;
	charray();//{ data = NULL; /*xn=NX; yn=NY; zn=NZ;*/}
	void setdim(int _xn, int _yn, int _zn){ xn = _xn, yn = _yn, zn = _zn; }

	__host__ __device__ inline char &operator ()(int i, int j, int k)
	{
		return data[i*yn*zn + j*zn + k];
	}
	__host__ __device__ inline char &operator ()(int i)
	{
		return data[i];
	}
	__host__ __device__ inline char &operator [](int i)
	{
		return data[i];
	}
};

__host__ __device__ inline void getijk(int &i, int &j, int &k, int &idx, int w, int h, int d)
{
	i = idx / d / h;
	j = idx / d%h;
	k = idx%d;
}

enum ERENDERMODE{
	RENDER_PARTICLE = 0,
	RENDER_MC,
	RENDER_GRID,
	RENDER_ALL,
	RENDER_CNT
};

enum SIMULATIONMODE{
	SIMULATION_WATER = 0,
	SIMULATION_SOLIDCOUPLING,
	SIMULATION_SMOKE,
	SIMULATION_BUBBLE,
	SIMULATION_HEATONLY,
	SIMULATION_CNT,
	SIMULATION_LBMWATER
};

enum SCENE{
	SCENE_FLUIDSPHERE = 0,
	SCENE_SMOKE,
	SCENE_BOILING,
	SCENE_BOILING_HIGHRES,
	SCENE_MULTIBUBBLE,
	SCENE_DAMBREAK,
	SCENE_MELTING,
	SCENE_MELTINGPOUR,		//melting simulation by pouring water.
	SCENE_FREEZING,
	SCENE_INTERACTION,			//interact with small bubbles, i.e., sub-grid bubbles.
	SCENE_INTERACTION_HIGHRES,			//interact with small bubbles, i.e., sub-grid bubbles.
	SCENE_MELTANDBOIL,		//interact with big bubble
	SCENE_MELTANDBOIL_HIGHRES,		//interact with big bubble
	SCENE_HEATTRANSFER,
	SCENE_CNT,
	SCENE_LBMWATER,
	SCENE_ALL
};

enum VELOCITYMODEL{
	FLIP = 0,
	CIP,
	HYBRID,
	VELOCITYMODEL_CNT
};

enum ECOLORMODE{
	COLOR_PRESS = 0,
	COLOR_UX,
	COLOR_UY,
	COLOR_UZ,
	COLOR_DIV,	//4
	COLOR_PHI,
	COLOR_MARK,	//6
	COLOR_LS,	//7
	COLOR_TP,	//8
	COLOR_CNT
};

enum TIMESTAT
{
	TIME_DYNAMICS,
	TIME_TRANSITION,
	TIME_DISPLAY,
	TIME_TOTAL,
	TIME_COUNT
};

typedef struct AABB
{
	float xMin, xMax;
	float yMin, yMax;
	float zMin, zMax;
} *pAabb;

//0~ total blue, >=6~total red.
__host__ __device__ inline float3 mapColorBlue2Red(float v);

struct matrix4
{
	float m[16];
};
struct matrix3x3
{
	float x00, x01, x02;
	float x10, x11, x12;
	float x20, x21, x22;
};

	inline matrix3x3 operator+(matrix3x3 A, matrix3x3 B)
	{
		B.x00 += A.x00; B.x01 += A.x01; B.x02 += A.x02;
		B.x10 += A.x10; B.x11 += A.x11; B.x12 += A.x12;
		B.x20 += A.x20; B.x21 += A.x21; B.x22 += A.x22;
		return B;
	}
inline matrix3x3 operator*(matrix3x3 B, float b)
{
	matrix3x3 A;
	B.x00 *= b; B.x01 *= b; B.x02 *= b;
	B.x10 *= b; B.x11 *= b; B.x12 *= b;
	B.x20 *= b; B.x21 *= b; B.x22 *= b;
	return B;
}
inline matrix3x3 operator/(matrix3x3 B, float b)
{
	matrix3x3 A;
	B.x00 /= b; B.x01 /= b; B.x02 /= b;
	B.x10 /= b; B.x11 /= b; B.x12 /= b;
	B.x20 /= b; B.x21 /= b; B.x22 /= b;
	
	return B;
}

struct  float9
{
	float x0, x1, x2, x3, x4, x5, x6, x7, x8;
	//float x, y, z, x2, y2, z2, xy, yz, zx;
	/////xu
	float x[9];//
};

struct matrix9x9
{
	float x[9][9];
// 	float x00, x01, x02, x03, x04, x05, x06, x07, x08;
// 	float x10, x11, x12, x13, x14, x15, x16, x17, x18;
// 	float x20, x21, x22, x23, x24, x25, x26, x27, x28;
// 	float x30, x31, x32, x33, x34, x35, x36, x37, x38;
// 	float x40, x41, x42, x43, x44, x45, x46, x47, x48;
// 	float x50, x51, x52, x53, x54, x55, x56, x57, x58;
// 	float x60, x61, x62, x63, x64, x65, x66, x67, x68;
// 	float x70, x71, x72, x73, x74, x75, x76, x77, x78;
// 	float x80, x81, x82, x83, x84, x85, x86, x87, x88;
};

/////xu
struct matrix3x9
{
	float x[3][9];
	//float x00, x01, x02, x03, x04, x05, x06, x07, x08;
	//float x10, x11, x12, x13, x14, x15, x16, x17, x18;
	//float x20, x21, x22, x23, x24, x25, x26, x27, x28;
};

inline matrix3x9 operator+(matrix3x9 A, matrix3x9 B)
{
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 9; j++)
			B.x[i][j] += A.x[i][j];
	}
	return B;
}

inline matrix3x9 operator*(matrix3x9 B, float b)
{
	matrix3x9 A;
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 9; j++)
			B.x[i][j] *= b;
	}
	return B;
}

inline matrix3x9 operator/(matrix3x9 B, float b)
{
	matrix3x9 A;
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 9; j++)
			B.x[i][j] /= b;
	}
	return B;
}

inline matrix9x9 float9Multfloat9_9x9(float9 x)
{
	matrix9x9 A;
	for (int i = 0; i < 9; i++){
		for (int j = 0; j < 9; j++)
			A.x[i][j] = x.x[i] * x.x[j];
	}
	return A;
}

struct cluster
{
	int num;
	float3 restCM;
	float wsum;
	int snum;
	float3 evedis;
	bool havepar;
	float3 cm;

	matrix3x3 invRestMat;
	matrix3x3 A;
	float det;
	matrix3x3 mat;
	matrix3x3 Apq;
	matrix3x3 Aqq;
	matrix3x3 R, U, D;
	float detA;

	matrix9x9 A_bar;
	float9 q_bar;
	matrix9x9 Aqq_bar;
	matrix3x9 Apq_bar;
	matrix9x9 invRestMat_bar;
	matrix3x9 mat_bar;
	float det_bar;
	matrix3x9 ApqMultAqq_bar;
	matrix3x9 R_bar;

	//float corr_cluster
};
//

float determinant(matrix3x3 m);
/////xu
float determinant(matrix9x9 arcs, int n);
matrix9x9 getAStart(matrix9x9 arcs);
matrix9x9 inverse(matrix9x9 src);
float3 mat39Multfloat9(matrix3x9 a, float9 b);
//
matrix3x3 inverse(matrix3x3 m);
matrix3x3 mat3Multmat3(matrix3x3 a, matrix3x3 b);
float3 mat3Multfloat3(matrix3x3 a, float3 b);
matrix3x3 polarDecompositionStable(matrix3x3 mat, float eps);
float3 matRow(matrix3x3 m, int i);
float3 matCol(matrix3x3 m, int i);
matrix3x3 polarDecomposition(matrix3x3 mat, matrix3x3 R, matrix3x3 U, matrix3x3 D);

//*********************LBM parameters*******************
//struct LBMConstant {
//	int gnum;
//	int3 gvnum;
//	float samplespace;
//	float dt;
//	float3 gravity;
//	float3 gmin, gmax, cellsize;
//	float m0;
//	float waterrho;//	float solidrho;
//	float Pr, Ra, R;
//	float T0;
//	float p0;
//	float tau_f, tau_h;
//	float niu;
//	float cv;
//	float wf, wh;
//	float total_E, T_heat;
//// 	float pradius;
//// 	float3	triHashSize, triHashRes;		//triHashSize是HASH网格的大小;  triHashRes是每一个维度上有几个HASH网格,程序执行过程中不再变化
//// 	float3 t_min, t_max;
//// 	int triHashCells;			//预留的hash数组大小，程序执行过程中不再变化
//};
//********************************init LBM parameters****************************
//delta_T = 0.1;
//Pr = 0.7;
//Ra = 10000;
//RHO = 0.2;
//U = 0.20;
//R = 8.3144;	//气体普适常量
//T0 = 0.04;
//p0 = RHO * R *T0;
//tau_f = 0.50;	//????存疑
//tau_h = tau_f / Pr;
//niu = tau_f * p0;
//cv = (3 + 3)*R / 2.0;
//frame = 0;
//wf = 2 * delta_T / (2 * tau_f + delta_T);
//wh = 2 * delta_T / (2 * tau_h + delta_T);
//total_E = 0.0;
//T_heat = 200;
//********************************************************************************
//******************************************************

#endif
