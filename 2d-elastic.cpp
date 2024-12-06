#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <time.h>
#include <math.h>

#include <vector>
#include <algorithm>

#include "alglibmisc.h"
using namespace alglib;

//#define DEBUG_PRINT

constexpr int n_part = 271;
constexpr int n_equ = n_part*4; //n_part*4
constexpr int n_max_vec = 6;

constexpr double Kf = 1.0;
constexpr double mu = 0.5;

double radiusHS = 0.22;
double radiusLS = 0.21;
double radiusAVG = 0.215;

double radius = radiusHS;
double L0 = 0.6;
double F[n_part], Fex, Fey;
double prob_LH[n_part];
double prob_HL[n_part];

constexpr long double DH = 3040.0000;
constexpr long double DS = 10.3000;
constexpr long double DE = 2000.0000;
constexpr long double ka = 1450.0000;
constexpr long double tau = 10.00;
double T[n_part];
double T_HS = 500.0;
double T_LS = 150.0; 

double timp;
double t_max = 1.0e4;
double t_init = 0.0;
double step_t = 1.0;
double eps = 1.0e-3;
double diference = 1.0;
double dif_x, dif_y;

struct Medium {
	double x, y, z, r;
};
struct Medium Mediu[n_part];
int no_of_neigh[n_part];
int neighbours[n_part][n_max_vec];
double sol[n_equ], sol_old[n_equ];

//*************************************************************************

void alglib_function_neighbours(void);
int Funct_Dopri(double time, double* input, double* deriv);
int Dopri5(double x, double xend, double eps, double hmax, double has, double* sol);
void Probabilities();

//*************************************************************************

int main()
{
	FILE *fp = fopen("hexagonal_system.txt", "r");
	float dummy;
	for (int i = 0; i < n_part; i++)
	{
		fscanf(fp, "%f %lf %lf %lf %lf\n", &dummy, &Mediu[i].x, &Mediu[i].y, &Mediu[i].z, &Mediu[i].r);
	}
	fclose(fp);

	alglib_function_neighbours();

	//////////////////////////
	
	for (int i = 0; i < n_part; i++) //Conditii initiale
	{
		sol[4 * i + 0] = Mediu[i].x;
		sol[4 * i + 1] = 0.0;
		sol[4 * i + 2] = Mediu[i].y;
		sol[4 * i + 3] = 0.0;

		sol_old[4 * i + 0] = sol[4 * i + 0];
		sol_old[4 * i + 1] = sol[4 * i + 1];
		sol_old[4 * i + 2] = sol[4 * i + 2];
		sol_old[4 * i + 3] = sol[4 * i + 3];
	}

	timp = t_init;
	while ((diference >= 1.0e-3) || ((timp - t_init) < 100.0 * step_t))
	{
		for (int i = 0; i < n_part; i++) //Conditii initiale
		{
			sol_old[4 * i + 0] = sol[4 * i + 0];
			sol_old[4 * i + 1] = sol[4 * i + 1];
			sol_old[4 * i + 2] = sol[4 * i + 2];
			sol_old[4 * i + 3] = sol[4 * i + 3];
		}

		Dopri5(timp, timp + step_t, eps, step_t, step_t / 4.0, &sol[0]);

		for (int i = 0; i < n_part; i++)
		{
			Mediu[i].x = sol[4 * i + 0];
			Mediu[i].y = sol[4 * i + 2];
		}
		timp += step_t;

		diference = 0.0;
		for (int i = 0; i < n_part; i++)
		{

			dif_x = fabs(sol[4 * i + 0] - sol_old[4 * i + 0]);
			dif_y = fabs(sol[4 * i + 2] - sol_old[4 * i + 2]);
			if (diference < dif_x)
			{
				diference = dif_x;
			}
			if (diference < dif_y)
			{
				diference = dif_y;
			}
		}
		printf("timp = %lf    dif=%lf\n", timp, diference);
	}


	double random_value;
	double nHS = 0.0;
	int HS = 0;
	int count_temp_step = 0;

	//////////////////////////////////////////////////////////////////////////
	///  MHL - TO UP
	//////////////////////////////////////////////////////////////////////////
	FILE *file = fopen("T_f(nHS).txt", "w");
	if (!file) {
		perror("Error opening file");
		return(0);
	}

	timp = t_init;
	for (int i = 0; i < n_part; i++) //Conditii initiale
	{
		sol[4 * i + 0] = Mediu[i].x;
		sol[4 * i + 1] = 0.0;
		sol[4 * i + 2] = Mediu[i].y;
		sol[4 * i + 3] = 0.0;

		sol_old[4 * i + 0] = sol[4 * i + 0];
		sol_old[4 * i + 1] = sol[4 * i + 1];
		sol_old[4 * i + 2] = sol[4 * i + 2];
		sol_old[4 * i + 3] = sol[4 * i + 3];
	}

	Dopri5(timp, timp + step_t, eps, step_t, step_t / 4.0, &sol[0]);
	for (int i = 0; i < n_part; i++)
	{
		Mediu[i].x = sol[4 * i + 0];
		Mediu[i].y = sol[4 * i + 2];
	}
	for (int i = 0; i < n_part; i++) 
	{
		T[i] = T_LS;
	}

	count_temp_step = 0;
	// Temperature decrease
	for (double temperature = T_LS; temperature <= T_HS; temperature += 1.0e-3)
	{
		count_temp_step++;
		for (int i = 0; i < n_part; i++) 
		{
			T[i] = temperature;
		}

		Probabilities();

		for (int i = 0; i < n_part; i++)
		{
			random_value = (double)rand() / ((double)RAND_MAX + 1);

			if ((Mediu[i].r > radiusAVG) && (prob_HL[i] > random_value))
			{
				Mediu[i].r = radiusLS;
				HS--;
			}
			else
			{
				if ((Mediu[i].r < radiusAVG) && (prob_LH[i] > random_value))
				{
					Mediu[i].r = radiusHS;
					HS++;
				}
			}

		}
		nHS = (double)HS / n_part;

		Dopri5(timp, timp + step_t, eps, step_t, step_t / 4.0, &sol[0]);

		for (int i = 0; i < n_part; i++)
		{
			Mediu[i].x = sol[4 * i + 0];
			Mediu[i].y = sol[4 * i + 2];
			sol_old[4 * i + 0] = sol[4 * i + 0];
			sol_old[4 * i + 1] = sol[4 * i + 1];
			sol_old[4 * i + 2] = sol[4 * i + 2];
			sol_old[4 * i + 3] = sol[4 * i + 3];
		}

		if ((count_temp_step % 100) == 0)
		{
			fprintf(file, "%.4f\t%.4f\n", temperature, nHS);
			printf("%.4f\t %.4f\n", temperature, nHS);
		}
	}
	fclose(file);


	//////////////////////////////////////////////////////////////////////////
	///  MHL - TO DOWN
	//////////////////////////////////////////////////////////////////////////
	file = fopen("T_f(nHS).txt", "a");
	if (!file) {
		perror("Error opening file");
		return(0);
	}

	nHS = 1.0;
	HS = n_part;

	timp = t_init;
	for (int i = 0; i < n_part; i++) //Conditii initiale
	{
		sol[4 * i + 0] = Mediu[i].x;
		sol[4 * i + 1] = 0.0;
		sol[4 * i + 2] = Mediu[i].y;
		sol[4 * i + 3] = 0.0;

		sol_old[4 * i + 0] = sol[4 * i + 0];
		sol_old[4 * i + 1] = sol[4 * i + 1];
		sol_old[4 * i + 2] = sol[4 * i + 2];
		sol_old[4 * i + 3] = sol[4 * i + 3];
	}

	Dopri5(timp, timp + step_t, eps, step_t, step_t / 4.0, &sol[0]);
	for (int i = 0; i < n_part; i++)
	{
		Mediu[i].x = sol[4 * i + 0];
		Mediu[i].y = sol[4 * i + 2];
	}

	count_temp_step = 0;
	// Temperature increase
	for (double temperature = T_HS; temperature >= T_LS; temperature -= 1.0e-3)
	{
		count_temp_step++;
		for (int i = 0; i < n_part; i++)
		{
			T[i] = temperature;
		}

		Probabilities();

		for (int i = 0; i < n_part; i++)
		{
			random_value = (double)rand() / ((double)RAND_MAX + 1);

			if ((Mediu[i].r > radiusAVG) && (prob_HL[i] > random_value))
			{
				Mediu[i].r = radiusLS;
				HS--;
			}
			else
			{
				if ((Mediu[i].r < radiusAVG) && (prob_LH[i] > random_value))
				{
					Mediu[i].r = radiusHS;
					HS++;
				}
			}

		}
		nHS = (double)HS / n_part;

		Dopri5(timp, timp + step_t, eps, step_t, step_t / 4.0, &sol[0]);

		for (int i = 0; i < n_part; i++)
		{
			Mediu[i].x = sol[4 * i + 0];
			Mediu[i].y = sol[4 * i + 2];
			sol_old[4 * i + 0] = sol[4 * i + 0];
			sol_old[4 * i + 1] = sol[4 * i + 1];
			sol_old[4 * i + 2] = sol[4 * i + 2];
			sol_old[4 * i + 3] = sol[4 * i + 3];
		}

		if ((count_temp_step % 100) == 0)
		{
			fprintf(file, "%.4f\t%.4f\n", temperature, nHS);
			printf("%.4f\t %.4f\n", temperature, nHS);
		}
	}


	//new radii file saving 
	file = fopen("radii.txt", "w");
	for (int i = 0; i < n_part; i++) {
		
		fprintf(file, "%lf\t  %lf\t  %lf \n ", Mediu[i].x, Mediu[i].y, Mediu[i].r);
	}
	fclose(file);
	//*************************************************************************
}

//*************************************************************************

void Probabilities()
{
	//*************************************************************************************************************
	//PROBABILITIES
	int index_neigh = 0;
	double Fe = 0.0;

	for (int i = 0; i < n_part; i++)
	{
		Fe = 0.0;
		for (int j = 0; j < no_of_neigh[i]; j++) 
		{
			index_neigh = neighbours[i][j];
			Fe += Kf * (sqrt(pow((Mediu[index_neigh].x - Mediu[i].x), 2.0) + pow((Mediu[index_neigh].y - Mediu[i].y), 2.0)) - L0 - Mediu[i].r - Mediu[index_neigh].r);
		}
		prob_LH[i] = (exp(-(DH - T[i] * DS) / 2.0 / T[i]) * exp(-((DE + ka * Fe) / T[i])) / tau);
		prob_HL[i] = (exp( (DH - T[i] * DS) / 2.0 / T[i]) * exp(-((DE - ka * Fe) / T[i])) / tau);
	}
}

//*************************************************************************

void alglib_function_neighbours(void)
{
	int	neighbours_max = 0, neighbours_med = 0;
	int i, j, local_index;
	real_2d_array a;

	a.setlength(n_part, 3);
	for (i = 0; i < n_part; i++)
	{
		a(i, 0) = Mediu[i].x;
		a(i, 1) = Mediu[i].y;
		a(i, 2) = Mediu[i].z;
	}

	integer_1d_array tags;
	tags.setlength(n_part);
	for (int i = 0; i < n_part; i++)
	{
		tags(i) = i;
	}

	ae_int_t nx = 3;
	ae_int_t ny = 0;
	ae_int_t normtype = 2;

	kdtree kdt;
	//kdtreebuild(a, nx, ny, normtype, kdt);
	kdtreebuildtagged(a, tags, nx, ny, normtype, kdt);

	real_1d_array x;
	x.setlength(3);
	real_2d_array r = "[[]]";

	integer_1d_array indexes;

	for (i = 0; i < n_part; i++)
	{
		x(0) = Mediu[i].x;
		x(1) = Mediu[i].y;
		x(2) = Mediu[i].z;

		ae_int_t k;

		k = kdtreequeryrnn(kdt, x, 1.1 * (2.0 * radius + L0), false);

		no_of_neigh[i] = (int)k;

// 		if (neighbours[i] + 1 > n_max_vec - 1)
// 		{
// 			printf("\n\n PREA MULTI VECINI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! \n\n");
// 			break;
// 		}
// 		else
		{
			kdtreequeryresultstags(kdt, indexes);

			//SORT FOR TYPE
			std::vector<int> myvector(no_of_neigh[i]);
			for (j = 0; j < no_of_neigh[i]; j++)
				myvector[j] = (int)indexes(j);

			std::sort(myvector.begin(), myvector.end());

#ifdef DEBUG_PRINT
			printf("particula: %d  cu  %d  vecini -> ", i, no_of_neigh[i]);
#endif
			for (j = 0; j < no_of_neigh[i]; j++)
			{
				local_index = myvector[j];// (int)indexes(j);

				neighbours[i][j] = local_index;  //<<<---- asta e vecin
#ifdef DEBUG_PRINT
				printf(" %d ", neighbours[i][j]);
#endif
			}
#ifdef DEBUG_PRINT
			printf("\n");
#endif
		}

		neighbours_med += no_of_neigh[i];
		if (no_of_neigh[i] > neighbours_max)	neighbours_max = no_of_neigh[i];
	}
	printf("Numar maxim de vecini: %d    Numar mediu de vecini: %f \n", neighbours_max, (double)neighbours_med / n_part);
}

//*************************************************************************

int Funct_Dopri(double time, double *input, double *deriv)
{
	long i, j, k, kk;
	double radical, Fe, Fex, Fey;
	double alungirea;
	double distanta_normala = L0;

	for (i = 0; i < n_part; i++)
	{
		k = 4 * i;
		Fe = 0.0;
		Fex = 0.0;
		Fey = 0.0;

		for (j = 0; j < no_of_neigh[i]; j++)
		{
			kk = 4 * neighbours[i][j];
			radical = sqrt((input[kk + 0] - input[k + 0]) * (input[kk + 0] - input[k + 0]) + (input[kk + 2] - input[k + 2]) * (input[kk + 2] - input[k + 2]) );

			alungirea = (radical - Mediu[neighbours[i][j]].r - Mediu[i].r - distanta_normala);

			Fe = Kf * alungirea;

			Fex += Fe * (input[kk + 0] - input[k + 0]) / radical;
			Fey += Fe * (input[kk + 2] - input[k + 2]) / radical;
		}

		deriv[k + 0] = input[k + 1];
		deriv[k + 1] = (Fex - mu * input[k + 1])/* / m*/;
		deriv[k + 2] = input[k + 3];
		deriv[k + 3] = (Fey - mu * input[k + 3])/* / m*/;
	}
	return(0);
}

//*************************************************************************

int Dopri5(double x, double xend, double eps, double hmax, double has, double* sol)

{
	constexpr double uround = 1.0e-7;
	int i;
	int n = n_equ, nmax, reject, nrejct, naccpt, nfcn, nstep;
	double posneg, auxmax, xph, er, auxmax1, auxmax2, denom, hnew, auxmin, aux, fac;

	double verif;

	//#if defined(GPU_aux) && defined(GPU)

	//#else
	//Memory allocation 
	double* k1, * k2, * k3, * k4, * k5, * ygrec1;
	k1 = (double*)malloc(n_equ * sizeof(double));
	k2 = (double*)malloc(n_equ * sizeof(double));
	k3 = (double*)malloc(n_equ * sizeof(double));
	k4 = (double*)malloc(n_equ * sizeof(double));
	k5 = (double*)malloc(n_equ * sizeof(double));
	ygrec1 = (double*)malloc(n_equ * sizeof(double));
	//#endif

#ifndef GPU_aux
	for (i = 0; i < n; i++)
	{
		ygrec1[i] = k1[i] = k2[i] = k3[i] = k4[i] = k5[i] = 0.0;
	}
#else
	Dopri5_00_init0 << <Blk, Thrd >> > (&k1[0], &k2[0], &k3[0], &k4[0], &k5[0], &ygrec1[0]);
	cudaDeviceSynchronize();
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "Dopri5_00_init0 launch failed: %s\n", cudaGetErrorString(cudaStatus));
		return(cudaStatus);
	}
#endif

	nmax = 3000;
	//posneg = sign(1.0, xend-x);
	((xend - x) >= 0.0) ? (posneg = 1.0) : (posneg = -1.0);

	//***************************
	//   initializare
	//***************************

	hmax = fabs(hmax);
	if (1.0e-4 > fabs(has))  auxmax = 1.0e-4;
	else auxmax = fabs(has);
	if (auxmax > hmax)       has = hmax;
	else has = auxmax;
	if (posneg > 0.0)        has = fabs(has);
	else has = -fabs(has);
	//eps    = Amax1r(eps, 7.0*uround);
	if (eps >= 7.0 * uround) (eps = 7.0 * uround);
	reject = 0;
	naccpt = 0;
	nrejct = 0;
	nfcn = 1;
	nstep = 0;

#ifndef GPU
	Funct_Dopri(x, &sol[0], &k1[0]);
#else
	Funct_DopriCUDA << <BlkD, ThrdD >> > (&k1[0], &sol[0], &Medium_k[0], &Medium_raza[0], &Position_Coef_vecin[0], &neighbours[0], Ext_force_x, Ext_force_y);
	cudaDeviceSynchronize();
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "Funct_DopriCUDA launch failed: %s\n", cudaGetErrorString(cudaStatus));
		return(cudaStatus);
}
#endif

	//****************************
	//   pas de integrare de baza
	//****************************

	while ((x - xend) * posneg + uround <= 0.0)
	{
		verif = x + 1.0e-1 * has - x;
		if ((nstep > nmax) || (verif == 0.0))
		{
			//AfxMessageBox("OUT!!!");

			return (1);
		}

		if ((x + has - xend) * posneg > 0.0)
			has = xend - x;
		nstep++;

		//****************************
		//   primele 6 etape
		//****************************

#ifndef GPU_DOPRI

		/////////////////////////////////////////////////////////////////////////////
#ifndef GPU_aux
		for (i = 0; i < n; i++)
			ygrec1[i] = sol[i] + has * 2.0e-1 * k1[i];
#else
		Dopri5_01_add1 << <Blk, Thrd >> > (&ygrec1[0], has, &sol[0], &k1[0]);
		cudaDeviceSynchronize();
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess)
		{
			fprintf(stderr, "Dopri5_01_add1 launch failed: %s\n", cudaGetErrorString(cudaStatus));
			return(cudaStatus);
	}
#endif
		//////////////////////////////////////////////////////////////////////////////
#ifndef GPU
		Funct_Dopri(x + 2.0e-1 * has, &ygrec1[0], &k2[0]);
#else
		Funct_DopriCUDA << <BlkD, ThrdD >> > (&k2[0], &ygrec1[0], &Medium_k[0], &Medium_raza[0], &Position_Coef_vecin[0], &neighbours[0], Ext_force_x, Ext_force_y);
		cudaDeviceSynchronize();
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess)
		{
			fprintf(stderr, "Funct_DopriCUDA launch failed: %s\n", cudaGetErrorString(cudaStatus));
			return(cudaStatus);
		}
#endif
		//////////////////////////////////////////////////////////////////////////////
#ifndef GPU_aux
		for (i = 0; i < n; i++)
			ygrec1[i] = sol[i] + has * ((3.0e0 / 4.0e1) * k1[i] + (9.0e0 / 4.0e1) * k2[i]);
#else
		Dopri5_02_add2 << <Blk, Thrd >> > (&ygrec1[0], has, &sol[0], &k1[0], &k2[0]);
		cudaDeviceSynchronize();
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess)
		{
			fprintf(stderr, "Dopri5_02_add2 launch failed: %s\n", cudaGetErrorString(cudaStatus));
			return(cudaStatus);
		}
#endif
		//////////////////////////////////////////////////////////////////////////////
#ifndef GPU
		Funct_Dopri(x + 3.0e-1 * has, &ygrec1[0], &k3[0]);
#else
		Funct_DopriCUDA << <BlkD, ThrdD >> > (&k3[0], &ygrec1[0], &Medium_k[0], &Medium_raza[0], &Position_Coef_vecin[0], &neighbours[0], Ext_force_x, Ext_force_y);
		cudaDeviceSynchronize();
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess)
		{
			fprintf(stderr, "Funct_DopriCUDA launch failed: %s\n", cudaGetErrorString(cudaStatus));
			return(cudaStatus);
		}
#endif
		//////////////////////////////////////////////////////////////////////////////
#ifndef GPU_aux
		for (i = 0; i < n; i++)
			ygrec1[i] = sol[i] + has * ((4.4e1 / 4.5e1) * k1[i] - (5.6e1 / 1.5e1) * k2[i] + (3.2e1 / 9.0e0) * k3[i]);
#else
		Dopri5_03_add3 << <Blk, Thrd >> > (&ygrec1[0], has, &sol[0], &k1[0], &k2[0], &k3[0]);
		cudaDeviceSynchronize();
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess)
		{
			fprintf(stderr, "Dopri5_03_add3 launch failed: %s\n", cudaGetErrorString(cudaStatus));
			return(cudaStatus);
		}
#endif
		//////////////////////////////////////////////////////////////////////////////
#ifndef GPU
		Funct_Dopri(x + 8.0e-1 * has, &ygrec1[0], &k4[0]);
#else
		Funct_DopriCUDA << <BlkD, ThrdD >> > (&k4[0], &ygrec1[0], &Medium_k[0], &Medium_raza[0], &Position_Coef_vecin[0], &neighbours[0], Ext_force_x, Ext_force_y);
		cudaDeviceSynchronize();
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess)
		{
			fprintf(stderr, "Funct_DopriCUDA launch failed: %s\n", cudaGetErrorString(cudaStatus));
			return(cudaStatus);
		}
#endif
		//////////////////////////////////////////////////////////////////////////////
#ifndef GPU_aux
		for (i = 0; i < n; i++)
			ygrec1[i] = sol[i] + has * ((1.9372e4 / 6.561e3) * k1[i] - (2.5360e4 / 2.187e3) * k2[i] + (6.4448e4 / 6.561e3) * k3[i] - (2.1200e2 / 7.290e2) * k4[i]);
#else
		Dopri5_04_add4 << <Blk, Thrd >> > (&ygrec1[0], has, &sol[0], &k1[0], &k2[0], &k3[0], &k4[0]);
		cudaDeviceSynchronize();
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess)
		{
			fprintf(stderr, "Dopri5_04_add4 launch failed: %s\n", cudaGetErrorString(cudaStatus));
			return(cudaStatus);
		}
#endif
		//////////////////////////////////////////////////////////////////////////////
#ifndef GPU
		Funct_Dopri(x + (8.0e0 / 9.0e0) * has, &ygrec1[0], &k5[0]);
#else
		Funct_DopriCUDA << <BlkD, ThrdD >> > (&k5[0], &ygrec1[0], &Medium_k[0], &Medium_raza[0], &Position_Coef_vecin[0], &neighbours[0], Ext_force_x, Ext_force_y);
		cudaDeviceSynchronize();
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess)
		{
			fprintf(stderr, "Funct_DopriCUDA launch failed: %s\n", cudaGetErrorString(cudaStatus));
			return(cudaStatus);
		}
#endif
		//////////////////////////////////////////////////////////////////////////////
#ifndef GPU_aux
		for (i = 0; i < n; i++)
			ygrec1[i] = sol[i] + has * ((9.0170e3 / 3.1680e3) * k1[i] - (3.5500e2 / 3.3000e1) * k2[i] + (4.6732e4 / 5.2470e3) * k3[i] + (4.9000e1 / 1.7600e2) * k4[i] - (5.1030e3 / 1.8656e4) * k5[i]);
#else
		Dopri5_05_add5 << <Blk, Thrd >> > (&ygrec1[0], has, &sol[0], &k1[0], &k2[0], &k3[0], &k4[0], &k5[0]);
		cudaDeviceSynchronize();
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess)
		{
			fprintf(stderr, "Dopri5_05_add5 launch failed: %s\n", cudaGetErrorString(cudaStatus));
			return(cudaStatus);
		}
#endif
		//////////////////////////////////////////////////////////////////////////////
		xph = x + has;
		//////////////////////////////////////////////////////////////////////////////
#ifndef GPU
		Funct_Dopri(xph, &ygrec1[0], &k2[0]);
#else
		Funct_DopriCUDA << <BlkD, ThrdD >> > (&k2[0], &ygrec1[0], &Medium_k[0], &Medium_raza[0], &Position_Coef_vecin[0], &neighbours[0], Ext_force_x, Ext_force_y);
		cudaDeviceSynchronize();
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess)
		{
			fprintf(stderr, "Funct_DopriCUDA launch failed: %s\n", cudaGetErrorString(cudaStatus));
			return(cudaStatus);
		}
#endif
		//////////////////////////////////////////////////////////////////////////////
#ifndef GPU_aux
		for (i = 0; i < n; i++)
			ygrec1[i] = sol[i] + has * ((3.500e1 / 3.840e2) * k1[i] + (1.100e1 / 8.400e1) * k2[i] + (5.000e2 / 1.113e3) * k3[i] + (1.250e2 / 1.920e2) * k4[i] - (2.187e3 / 6.784e3) * k5[i]);
#else
		Dopri5_06_add6 << <Blk, Thrd >> > (&ygrec1[0], has, &sol[0], &k1[0], &k2[0], &k3[0], &k4[0], &k5[0]);
		cudaDeviceSynchronize();
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess)
		{
			fprintf(stderr, "Dopri5_02_add2 launch failed: %s\n", cudaGetErrorString(cudaStatus));
			return(cudaStatus);
		}
#endif
		//////////////////////////////////////////////////////////////////////////////
		//*******************************
		//   calculul sumei intermediare
		//   pentru economie de memorie
		//*******************************
#ifndef GPU_aux
		for (i = 0; i < n; i++)
			k2[i] = (7.1000e1 / 5.7600e4) * k1[i] + (2.2000e1 / 5.2500e2) * k2[i] - (7.1000e1 / 1.6695e4) * k3[i] + (7.1000e1 / 1.9200e3) * k4[i] - (1.7263e4 / 3.3920e5) * k5[i];
#else
		Dopri5_07_add7 << <Blk, Thrd >> > (&k2[0], &k1[0], &k3[0], &k4[0], &k5[0]);
		cudaDeviceSynchronize();
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess)
		{
			fprintf(stderr, "Dopri5_07_add7 launch failed: %s\n", cudaGetErrorString(cudaStatus));
			return(cudaStatus);
		}
#endif
		//////////////////////////////////////////////////////////////////////////////
		//****************************
		//   ultima etapa
		//****************************
		//////////////////////////////////////////////////////////////////////////////
#ifndef GPU
		Funct_Dopri(xph, &ygrec1[0], &k3[0]);
#else
		Funct_DopriCUDA << <BlkD, ThrdD >> > (&k3[0], &ygrec1[0], &Medium_k[0], &Medium_raza[0], &Position_Coef_vecin[0], &neighbours[0], Ext_force_x, Ext_force_y);
		cudaDeviceSynchronize();
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess)
		{
			fprintf(stderr, "Funct_DopriCUDA launch failed: %s\n", cudaGetErrorString(cudaStatus));
			return(cudaStatus);
		}
#endif
		//////////////////////////////////////////////////////////////////////////////
#ifndef GPU_aux
		for (i = 0; i < n; i++)
			k4[i] = (k2[i] - (1.0e0 / 4.0e1) * k3[i]) * has;
#else
		Dopri5_08_add8 << <Blk, Thrd >> > (&k4[0], has, &k2[0], &k3[0]);
		cudaDeviceSynchronize();
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess)
		{
			fprintf(stderr, "Dopri5_08_add8 launch failed: %s\n", cudaGetErrorString(cudaStatus));
			return(cudaStatus);
		}
#endif
		//////////////////////////////////////////////////////////////////////////////
#else
		Dopri5_XX_FullStep << <BlkD, ThrdD >> > (&ygrec1[0], has, &sol[0], &k1[0], &k2[0], &k3[0], &k4[0], &k5[0], &Medium_raza[0], &Position_Coef_vecin[0], &neighbours[0]);
		cudaDeviceSynchronize();
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess)
		{
			fprintf(stderr, "Dopri5_XX_FullStep launch failed: %s\n", cudaGetErrorString(cudaStatus));
			return(cudaStatus);
		}
		xph = x + has;
#endif
		nfcn += 6;

		//****************************
		//   estimarea erorii
		//****************************

		er = 0.0e0;
		for (i = 0; i < n; i++)
		{
			//auxmax1 = Amax1r(1.0e-5, fabs(ygrec1[i]));
			//auxmax2 = Amax1r(fabs(sol[i]), 2.0e0 * uround / eps);
			//denom   = Amax1r(auxmax1, auxmax2);
			//er     += (k4[i] / denom) * (k4[i] / denom);

			(1.0e-5 >= fabs(ygrec1[i])) ? (auxmax1 = 1.0e-5) : (auxmax1 = fabs(ygrec1[i]));
			(fabs(sol[i]) >= 2.0e0 * uround / eps) ? (auxmax2 = fabs(sol[i])) : (auxmax2 = 2.0e0 * uround / eps);
			(auxmax1 >= auxmax2) ? (denom = auxmax1) : (denom = auxmax2);
			er += (k4[i] / denom) * (k4[i] / denom);
		}

		er = sqrt(er / n);

		//***************************************
		//   calculul variabilei hnew
		//   se cere 0.2d0 <= hnew / has <= 10.0d0
		//***************************************

		aux = exp(log(er / eps) / 5.0e0);
		//auxmin = Amin1r(5.0e0, aux / 9.0e-1);
		//fac    = Amax1r(1.0e-1, auxmin);
		(5.0e0 <= aux / 9.0e-1) ? (auxmin = 5.0e0) : (auxmin = aux / 9.0e-1);
		(1.0e-1 >= auxmin) ? (fac = 1.0e-1) : (fac = auxmin);
		hnew = has / fac;
		if (er <= eps)
		{
			//****************************
			//   pasul este acceptat
			//****************************

			naccpt++;
			for (i = 0; i < n; i++)
			{
				k1[i] = k3[i];
				sol[i] = ygrec1[i];
			}

			x = xph;
			if (fabs(hnew) > hmax)
				hnew = posneg * hmax;
			//auxmin = Amin1r(fabs(hnew), fabs(has));
			(fabs(hnew) <= fabs(has)) ? (auxmin = fabs(hnew)) : (auxmin = fabs(has));
			if (reject == 1)
				hnew = posneg * auxmin;
			reject = 0;
		}

		//****************************
		//   pasul nu este acceptat
		//****************************

		if (er > eps)
		{
			reject = 1;
			if (naccpt >= 1)
				nrejct++;
		}

		has = hnew;
	}

	free(k1); free(k2); free(k3); free(k4); free(k5); free(ygrec1);

	return(0);
}

//*************************************************************************
