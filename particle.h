/*
 * particle.h
 *
 *  Created on: Nov 18, 2015
 *      Author: reza
 */

#ifndef PARTICLE_H_
#define PARTICLE_H_

#include <iostream>
#include <fstream>
#include <vector>
#include <time.h>
#include "vectorclass/vector3d.h"

using namespace std;

inline double kernel(double dij, double hRcp){
//	double alpha = 0.167*hRcp;				//In 1d = (1./6)*(1./*h)
	double alpha = 0.114*hRcp*hRcp;			//In 2d = (1./6)*(15/(7*pi*h*h))
//	double alpha = 0.080*hRcp*hRcp*hRcp;	//In 3d = (1./6)*(3/(2*pi*h*h*h))

	double R = dij*hRcp;

	double W;
	if(R>2.0)
		W = 0.;
	else if(R>1.0){
		double a = (2.-R);
		W = alpha*a*a*a;
	}
	else{
		double a = (2.-R);
		double b = (1.-R);
		W = alpha*(a*a*a - 4.*b*b*b);
	}
	return W;
}

inline Vec3d gradKernel(Vec3d Rij, double hRcp){
	//	double alpha = -0.500*hRcp*hRcp;				//In 1d = (-3./h)*(1./6)*(1./*h)
		double alpha = -0.341*hRcp*hRcp*hRcp;			//In 2d = (-3./h)*(1./6)*(15/(7*pi*h*h))
	//	double alpha = -0.239*hRcp*hRcp*hRcp*hRcp;		//In 3d = (-3./h)*(1./6)*(3/(2*pi*h*h*h))

	double dij = vector_length(Rij);
	double R = dij*hRcp;

	double W;
	if(R>2.0)
		W = 0.;
	else if(R>1.0){
		double a = (2.-R);
		W = alpha*a*a;
	}
	else{
		double a = (2.-R);
		double b = (1.-R);
		W = alpha*(a*a - 4.*b*b);
	}

	Vec3d nij = normalize_vector(Rij);
	return W*nij;
}

struct particle{
	Vec3d r0;				//initial position
	Vec3d u;				//displacement
	Vec3d v;				//velocity
	Vec3d vleap;			//velocity leap
	Vec3d a;				//acceleration
	vector<double> gu;		//gradient of displacement
	vector<double> epsilon;	//strain
	vector<double> sigma;	//stress
	double m;				//mass
	double rho;				//density
	double V;				//Initial volume
	double E;				//Young modulus
	double nu;				//Poisson ratio
	double U;				//Strain energy
	int type;				//type
	particle(): gu(9),epsilon(9),sigma(9){};
};

class particleSystem{
private:
	vector<particle> _particle;

public:
	int particleSize();
	void particleResize(int n);

	particle getParticle(int n);
	void setParticle(int n, particle p);
	Vec3d getPosInit(int n);
	void setPosInit(int n, Vec3d r0);
	Vec3d getDisp(int n);
	void setDisp(int n, Vec3d u);
	vector<double> getGradDisp(int n);
	void setGradDisp(int n, vector<double> gu);
	vector<double> getEpsilon(int n);
	void setEpsilon(int n, vector<double> epsilon);
	vector<double> getSigma(int n);
	void setSigma(int n, vector<double> sigma);
	Vec3d getCurrPos(int n);
	Vec3d getVel(int n);
	void setVel(int n, Vec3d v);
	Vec3d getVelLeap(int n);
	void setVelLeap(int n, Vec3d vleap);
	Vec3d getAcc(int n);
	void setAcc(int n, Vec3d acc);
	double getMass(int n);
	void setMass(int n, double m);
	double getDensity(int n);
	void setDensity(int n, double rho);
	void addDensity(int n, double rho);
	double getVolume(int n);
	void setVolume(int n, double V);
	double getYoung(int n);
	void setYoung(int n, double E);
	double getPoisson(int n);
	void setPoisson(int n, double nu);
	double getStrainEnergy(int n);
	void setStrainEnergy(int n, double U);
	int getType(int n);
	void setType(int n, int t);
};

class sph{
private:
	particleSystem _pSys;
	double t;
	double tMax;
	double dt;
	double h;
	double hRcp;
	double g;

	time_t timeStart;

	int fileCounter;
	int dumpSkip;
	int loopCounter;

public:
	sph();
	double getT();
	double getTMax();
	int getDumpSkip();
	int getLoopCounter();

	void calcDensity();
	void calcGradDisp();
	void calcStrain();
	void calcForce();
	void timeIntegration();
	void dumpData();
};

#endif /* PARTICLE_H_ */
