//============================================================================
// Name        : SPH_Elastic_2D.cpp
// Author      : Reza
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include "particle.h"

using namespace std;

int main(){
	sph _sph;

	while(_sph.getT() < _sph.getTMax()){
//		_sph.calcDensity();
		_sph.calcGradDisp();
		_sph.calcStrain();
		_sph.calcForce();
		_sph.timeIntegration();

		if(_sph.getLoopCounter() == _sph.getDumpSkip())
			_sph.dumpData();
	}

	return 0;
}
