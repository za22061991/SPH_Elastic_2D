/*
 * particle.cpp
 *
 *  Created on: Nov 18, 2015
 *      Author: reza
 */

#include "particle.h"

int particleSystem::particleSize(){
	return _particle.size();
}
void particleSystem::particleResize(int n){
	_particle.resize(n);
}

particle particleSystem::getParticle(int n){
	return _particle[n];
}
void particleSystem::setParticle(int n, particle p){
	_particle[n] = p;
}

Vec3d particleSystem::getPosInit(int n){
	return _particle[n].r0;
}
void particleSystem::setPosInit(int n, Vec3d r0){
	_particle[n].r0 = r0;
}

Vec3d particleSystem::getDisp(int n){
	return _particle[n].u;
}
void particleSystem::setDisp(int n, Vec3d u){
	_particle[n].u = u;
}

vector<double> particleSystem::getGradDisp(int n){
	return _particle[n].gu;
}
void particleSystem::setGradDisp(int n, vector<double> gu){
	_particle[n].gu = gu;
}

vector<double> particleSystem::getEpsilon(int n){
	return _particle[n].epsilon;
}
void particleSystem::setEpsilon(int n, vector<double> epsilon){
	_particle[n].epsilon = epsilon;
}

vector<double> particleSystem::getSigma(int n){
	return _particle[n].sigma;
}
void particleSystem::setSigma(int n, vector<double> sigma){
	_particle[n].sigma = sigma;
}

Vec3d particleSystem::getCurrPos(int n){
	return (_particle[n].r0 + _particle[n].u);
}

Vec3d particleSystem::getVel(int n){
	return _particle[n].v;
}
void particleSystem::setVel(int n, Vec3d v){
	_particle[n].v = v;
}

Vec3d particleSystem::getVelLeap(int n){
	return _particle[n].vleap;
}
void particleSystem::setVelLeap(int n, Vec3d vl){
	_particle[n].vleap = vl;
}

Vec3d particleSystem::getAcc(int n){
	return _particle[n].a;
}
void particleSystem::setAcc(int n, Vec3d a){
	_particle[n].a = a;
}

double particleSystem::getMass(int n){
	return _particle[n].m;
}
void particleSystem::setMass(int n, double m){
	_particle[n].m = m;
}

double particleSystem::getDensity(int n){
	return _particle[n].rho;
}
void particleSystem::setDensity(int n, double rho){
	_particle[n].rho = rho;
}
void particleSystem::addDensity(int n, double rho){
	_particle[n].rho += rho;
}

double particleSystem::getVolume(int n){
	return _particle[n].V;
}
void particleSystem::setVolume(int n, double V){
	_particle[n].V = V;
}

double particleSystem::getYoung(int n){
	return _particle[n].E;
}
void particleSystem::setYoung(int n, double E){
	_particle[n].E = E;
}

double particleSystem::getPoisson(int n){
	return _particle[n].nu;
}
void particleSystem::setPoisson(int n, double nu){
	_particle[n].nu = nu;
}

double particleSystem::getStrainEnergy(int n){
	return _particle[n].U;
}
void particleSystem::setStrainEnergy(int n, double U){
	_particle[n].U = U;
}

int particleSystem::getType(int n){
	return _particle[n].type;
}
void particleSystem::setType(int n, int t){
	_particle[n].type = t;
}

/*
//Initialization for roller
sph::sph(){
	//Initializing parameters
	double initRad = 0.05;		//Particle "radius"
	double r1 = 0.3;			//Inner radius of roller
	double r2 = 1.;				//Outer radius of roller

	g = -9.8;
	t = 0.;
	tMax = 1.;
	dt = 0.005;
	h = 4.*initRad;
	hRcp = 1./h;
	dumpSkip = 0.01/dt;
	fileCounter = 0;
	loopCounter = 0;

	timeStart = clock();

	//Initializing particles system
	double initFlag = 1;
	double radAnnulus = r1;
	int annulusCounter = 0;
	while(initFlag){
		int nInit = _pSys.particleSize();
		int n = (0.5*M_PI*radAnnulus/initRad);
		n *= 2;
		double dTeta = 2.*M_PI/n;
		_pSys.particleResize(nInit+n);

		for(int i=0; i<n; i++){
			particle pTemp;
			pTemp.r0 = Vec3d(radAnnulus*cos(i*dTeta), radAnnulus*sin(i*dTeta), 0.);
			pTemp.u = Vec3d(0., 0., 0.);
			pTemp.v = Vec3d(0., 0., 0.);
			pTemp.vleap = Vec3d(0., 0., 0.);
			pTemp.a = Vec3d(0., 0., 0.);
			pTemp.m = 0.001;
			pTemp.E = 100.;
			pTemp.nu = 0.4;
			if(annulusCounter == 0)
				pTemp.type = 1;
			else if(i==0)
				pTemp.type = 2;
			else
				pTemp.type = 0;
			_pSys.setParticle(nInit+i, pTemp);
		}
		radAnnulus += 2.*(3./M_PI)*initRad;
		annulusCounter++;
		if(radAnnulus>r2)
			initFlag = 0;
	}

	int n = _pSys.particleSize();

	FILE *fout;
	fout=fopen("structure.dat", "w+");
	for(int i=0; i<n; i++){
		Vec3d posTemp = _pSys.getPosInit(i);
		fprintf(fout, "%f\t%f\n", posTemp[0], posTemp[1]);
	}
	fclose(fout);

	dumpData();

	calcDensity();
}
*/

sph::sph(){
	//Initializing parameters
	double initRad = 0.05;		//Particle "radius"
	int nx = 20;
	int ny = 5;

	g = -9.8;
	t = 0.;
	tMax = 1.;
	dt = 0.0001;
	h = 4.*initRad;
	hRcp = 1./h;
	dumpSkip = 0.01/dt;
	fileCounter = 0;
	loopCounter = 0;

	timeStart = clock();

	//Initializing particles system
	_pSys.particleResize(nx*ny);
	for(int i=0; i<ny; i++){
		for(int j=0; j<nx; j++){
			particle pTemp;
			pTemp.r0 = Vec3d((j*2.+1.)*initRad, (i*2.+1.)*initRad, 0.);
			pTemp.u = Vec3d(0., 0., 0.);
			pTemp.v = Vec3d(0., 0., 0.);
			pTemp.vleap = Vec3d(0., 0., 0.);
			pTemp.a = Vec3d(0., 0., 0.);
			pTemp.m = 0.01;
			pTemp.E = 100.;
			pTemp.nu = 0.49;
			if(i == 0)
				pTemp.type = 1;
			else
				pTemp.type = 0;
			_pSys.setParticle(i*nx+j, pTemp);
		}
	}

	int n = _pSys.particleSize();

	FILE *fout;
	fout=fopen("structure.dat", "w+");
	for(int i=0; i<n; i++){
		Vec3d posTemp = _pSys.getPosInit(i);
		fprintf(fout, "%f\t%f\n", posTemp[0], posTemp[1]);
	}
	fclose(fout);

	dumpData();

	calcDensity();
}

double sph::getT(){
	return t;
}

double sph::getTMax(){
	return tMax;
}

int sph::getDumpSkip(){
	return dumpSkip;
}

int sph::getLoopCounter(){
	return loopCounter;
}

void sph::calcDensity(){
	int NNodes = _pSys.particleSize();

	for(int i=0; i<NNodes; i++){
		double rho0 = _pSys.getMass(i)*kernel(0., hRcp);
		_pSys.setDensity(i, rho0);
	}

	for(int i=0; i<NNodes; i++){
		Vec3d ri = _pSys.getPosInit(i);
		double mi = _pSys.getMass(i);
		double rhoi = _pSys.getDensity(i);
		for(int j=i+1; j<NNodes; j++){
			Vec3d rj = _pSys.getPosInit(j);
			double mj = _pSys.getMass(j);

			double dij = vector_length(ri-rj);
			double Wij = kernel(dij, hRcp);
			rhoi += mj*Wij;
			_pSys.addDensity(j, mi*Wij);
		}
		_pSys.setDensity(i, rhoi);
		_pSys.setVolume(i, mi/rhoi);
	}
}

void sph::calcGradDisp(){
	int NNodes = _pSys.particleSize();

	for(int i=0; i<NNodes; i++){
		vector<double> U(9, 0.);
		_pSys.setGradDisp(i, U);
	}

	for(int i=0; i<NNodes; i++){
		Vec3d ri = _pSys.getPosInit(i);
		Vec3d ui = _pSys.getDisp(i);
		vector<double> gui(_pSys.getGradDisp(i));
		double voli = _pSys.getVolume(i);
		for(int j=i+1; j<NNodes; j++){
			Vec3d rj = _pSys.getPosInit(j);
			Vec3d uj = _pSys.getDisp(j);
			vector<double> guj(_pSys.getGradDisp(j));
			double volj = _pSys.getVolume(j);

			Vec3d Rij = ri-rj;
			Vec3d uji = uj-ui;
			Vec3d gWij = gradKernel(Rij, hRcp);

			for(int k=0; k<3; k++){
				for(int l=0; l<3; l++){
					double gukl = uji[k]*gWij[l];
					gui[3*k+l] += volj*gukl;
					guj[3*k+l] += voli*gukl;
				}
			}
			_pSys.setGradDisp(j, guj);
		}
		_pSys.setGradDisp(i, gui);
	}
}

void sph::calcStrain(){
	int NNodes = _pSys.particleSize();

	for(int i=0; i<NNodes; i++){
		vector<double> gui(_pSys.getGradDisp(i));
		double Ei = _pSys.getYoung(i);
		double nui = _pSys.getPoisson(i);
		double voli = _pSys.getVolume(i);

		//Make gradient of displacement into Jacobian
		for(int j=0; j<3; j++){
			gui[3*j+j] += 1.;
		}

		//Calculate strain
		vector<double> epsilon(9,0.);
		for(int j=0; j<3; j++){
			for(int k=0; k<3; k++){
				for(int l=0; l<3; l++){
					epsilon[3*j+k] += gui[3*l+j]*gui[3*l+k];
				}
				if(j==k)
					epsilon[3*j+k] -= 1.;
				epsilon[3*j+k] *= 0.5;
			}
		}
		_pSys.setEpsilon(i, epsilon);

		//Calculate stress
		double CConst = Ei/((1.+nui)*(1.-2.*nui));
		vector<double> sigma(9);
		sigma[0] = CConst*((1.-nui)*epsilon[0] + nui*(epsilon[4]+epsilon[8]));
		sigma[4] = CConst*((1.-nui)*epsilon[4] + nui*(epsilon[0]+epsilon[8]));
		sigma[8] = CConst*((1.-nui)*epsilon[8] + nui*(epsilon[0]+epsilon[4]));
		for(int j=1; j<8; j++)
			if(j!=4)
				sigma[j] = CConst*(0.5-nui)*epsilon[j];
		_pSys.setSigma(i, sigma);

		//Calculate strain energy
		double Ui = 0.;
		for(int j=0; j<9; j++)
			Ui += sigma[j]*epsilon[j];
		Ui *= 0.5*voli;
		_pSys.setStrainEnergy(i, Ui);
	}
}

void sph::calcForce(){
	int NNodes = _pSys.particleSize();

	for(int i=0; i<NNodes; i++){
		_pSys.setAcc(i, Vec3d(0.,g,0.));
	}

	for(int i=0; i<NNodes; i++){
		Vec3d ri = _pSys.getPosInit(i);
		Vec3d ai = _pSys.getAcc(i);
		vector<double> gui(_pSys.getGradDisp(i));
		vector<double> sigmai(_pSys.getSigma(i));
		vector<double> epsiloni(_pSys.getEpsilon(i));
		double mi = _pSys.getMass(i);
		double mircp = 1./mi;
		double voli = _pSys.getVolume(i);
		for(int j=i+1; j<NNodes; j++){
			Vec3d rj = _pSys.getPosInit(j);
			Vec3d aj = _pSys.getAcc(j);
			vector<double> guj(_pSys.getGradDisp(j));
			vector<double> sigmaj(_pSys.getSigma(j));
			vector<double> epsilonj(_pSys.getEpsilon(j));
			double mj = _pSys.getMass(j);
			double mjrcp = 1./mj;
			double volj = _pSys.getVolume(j);

			Vec3d Rij = ri-rj;
			Vec3d dij = -voli*volj*gradKernel(Rij, hRcp);

			//Gradient U add with I
			vector<double> guit(gui), gujt(guj);
			for(int k=0; k<3; k++){
				guit[3*k+k] += 1.;
				gujt[3*k+k] += 1.;
			}

			//Multiplication of (gradU+I)*sigma
			vector<double> Mi(9,0.), Mj(9,0.);
			for(int k=0; k<3; k++){
				for(int l=0; l<3; l++){
					for(int m=0; m<3; m++){
						Mi[3*k+l] += guit[3*k+m]*sigmai[3*m+l];
						Mj[3*k+l] += gujt[3*k+m]*sigmaj[3*m+l];
					}
				}
			}

			//Calculate acceleration
			vector<double> aai(3,0.);
			vector<double> aaj(3,0.);
			for(int k=0; k<3; k++){
				aai[k] = 0.;
				aaj[k] = 0.;
				for(int l=0; l<3; l++){
					aai[k] -= Mj[3*k+l]*dij[l];
					aaj[k] += Mi[3*k+l]*dij[l];
				}
				aai[k] *= mircp;
				aaj[k] *= mjrcp;
			}
			ai += Vec3d(aai[0], aai[1], aai[2]);
			aj += Vec3d(aaj[0], aaj[1], aaj[2]);
			_pSys.setAcc(j, aj);
		}
		_pSys.setAcc(i, ai);
	}
}

void sph::timeIntegration(){
	int NNodes = _pSys.particleSize();

	if(t == 0.){
//		cout << "Initializing for leap velocity..." << endl;
		for(int i=0; i<NNodes; i++){
			_pSys.setVelLeap(i, _pSys.getVel(i)-0.5*_pSys.getVel(i)*dt);
		}
	}
	else{
		for(int i=0; i<NNodes; i++){
			int ti = _pSys.getType(i);
			if(ti != 1){
				Vec3d ai = _pSys.getAcc(i);
				Vec3d vli = _pSys.getVelLeap(i);
				Vec3d ui = _pSys.getDisp(i);
				vli += ai*dt;
				ui += vli*dt;
				_pSys.setVelLeap(i, vli);
				_pSys.setDisp(i, ui);
				_pSys.setVel(i, vli+0.5*ai*dt);
			}
//			//For rolling movement
//			else{
//				Vec3d ri = _pSys.getCurrPos(i);
//				Vec3d ui = _pSys.getDisp(i);
//				Vec3d wi(0., 0., 1.*sin(2.*M_PI*t));
////				Vec3d wi(0., 0., 2.*(1.-exp(-2.*M_PI*t)));
//
//				Vec3d vlli = _pSys.getVelLeap(i);
//				Vec3d vli = cross_product(wi, ri);
//				Vec3d ai = (vli - vlli) / dt;
//				ui += vli*dt;
//
//				_pSys.setVelLeap(i, vli);
//				_pSys.setVel(i, vli+0.5*ai*dt);
//				_pSys.setDisp(i, ui);
//			}
			//For hung particle
			else{
				_pSys.setVelLeap(i, Vec3d(0.,0.,0.));
				_pSys.setDisp(i, Vec3d(0.,0.,0.));
				_pSys.setVel(i, Vec3d(0.,0.,0.));
			}
		}
	}

	t += dt;
	loopCounter++;
}

void sph::dumpData(){
	int NNodes = _pSys.particleSize();

	FILE*fout;
	char fileName[30];
	sprintf(fileName, "dat/dat%05i.dat", fileCounter);
	fout = fopen(fileName, "w+");
	for(int i=0; i<NNodes; i++){
		Vec3d ri = _pSys.getCurrPos(i);
		fprintf(fout, "%f\t%f\n", ri[0], ri[1]);
	}
	fclose(fout);
	fileCounter++;
	loopCounter = 0;

	cout << "TSimulation = " << t << "\tTRun = " <<
			double(clock()-timeStart)/double(CLOCKS_PER_SEC) << " s." << endl;
}

