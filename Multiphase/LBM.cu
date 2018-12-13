/// \file lbm3D.c
/// \brief D3Q19 model implementation.
//
//	Copyright (c) 2014, Christian B. Mendl
//	All rights reserved.
//	http://christian.mendl.net
//
//	This program is free software; you can redistribute it and/or
//	modify it under the terms of the Simplified BSD License
//	http://www.opensource.org/licenses/bsd-license.php
//
//	Reference:
//	  Nils Thuerey, Physically based animation of free surface flows with the lattice Boltzmann method,
//	  PhD thesis, University of Erlangen-Nuremberg (2007)
//_______________________________________________________________________________________________________________________
//
#include<stdlib.h>
#include<stdio.h>
#include <time.h>
#include<memory.h>
#include"spray.h"
#include<GL/freeglut.h>
#include"spray.h"
#include"utility.h"

static const int vel3Di[19][3] = {
	{ 0, 0, 0 },		// zero direction
	{ -1, 0, 0 },		// 6 directions with velocity 1
	{ 1, 0, 0 },
	{ 0, -1, 0 },
	{ 0, 1, 0 },
	{ 0, 0, -1 },
	{ 0, 0, 1 },
	{ -1, -1, 0 },		// 12 directions with velocity sqrt(2)
	{ -1, 1, 0 },
	{ 1, -1, 0 },
	{ 1, 1, 0 },
	{ 0, -1, -1 },
	{ 0, -1, 1 },
	{ 0, 1, -1 },
	{ 0, 1, 1 },
	{ -1, 0, -1 },
	{ -1, 0, 1 },
	{ 1, 0, -1 },
	{ 1, 0, 1 },
};

/// \brief Velocity vectors for D3Q19 as floating-point values
static const float3 vel3Dv[19] = {
	{ 0, 0, 0 },		// zero direction
	{ -1, 0, 0 },		// 6 directions with velocity 1
	{ 1, 0, 0 },
	{ 0, -1, 0 },
	{ 0, 1, 0 },
	{ 0, 0, -1 },
	{ 0, 0, 1 },
	{ -1, -1, 0 },		// 12 directions with velocity sqrt(2)
	{ -1, 1, 0 },
	{ 1, -1, 0 },
	{ 1, 1, 0 },
	{ 0, -1, -1 },
	{ 0, -1, 1 },
	{ 0, 1, -1 },
	{ 0, 1, 1 },
	{ -1, 0, -1 },
	{ -1, 0, 1 },
	{ 1, 0, -1 },
	{ 1, 0, 1 },
};

/// \brief Index of inverse direction
static const int invVel3D[19] = { 0, 2, 1, 4, 3, 6, 5, 10, 9, 8, 7, 14, 13, 12, 11, 18, 17, 16, 15 };


//_______________________________________________________________________________________________________________________
//
/*

/// \brief D3Q19 weights
static const real weights3D[19] = { (real)1. / 3,
(real)1. / 18, (real)1. / 18, (real)1. / 18, (real)1. / 18, (real)1. / 18, (real)1. / 18,
(real)1. / 36, (real)1. / 36, (real)1. / 36, (real)1. / 36, (real)1. / 36, (real)1. / 36, (real)1. / 36, (real)1. / 36, (real)1. / 36, (real)1. / 36, (real)1. / 36, (real)1. / 36 };

device inline void DFD3Q19_DeriveQuantities(LBM::dfD3Q19_t *df)
{
	int i;

	// calculate average density
	df->rho = 0;
	for (i = 0; i < 19; i++)
	{
		df->rho += df->f[i];
	}
	assert(df->rho >= 0);

	// calculate average velocity u
	df->u.x = 0;
	df->u.y = 0;
	df->u.z = 0;
	if (df->rho > 0)
	{
		for (i = 0; i < 19; i++)
		{
			df->u.x += df->f[i] * vel3Dv[i].x;
			df->u.y += df->f[i] * vel3Dv[i].y;
			df->u.z += df->f[i] * vel3Dv[i].z;
		}
		real s = 1 / df->rho;
		df->u.x *= s;
		df->u.y *= s;
		df->u.z *= s;
	}

	// rescale in case maximum velocity is exceeded
	real n = Vec3_Norm(df->u);
	if (n > v_max)
	{
		df->u = Vec3_ScalarMultiply(v_max / n, df->u);
	}
}*/