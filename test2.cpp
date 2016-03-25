#include "Potential_HarmonicAngle.h"
#include "Potential_StillingerWeber.h"
#include <cmath>
#include <iostream>
#include <iomanip>

using namespace std;

Potential_HarmonicAngle potential;
//Potential_StillingerWeber potential;
array3d rij{-0.9, 0.4, -0.1};
array3d rik{0.9, 0.3, 0.1};
constexpr double finiteDist = 1.e-10;


void test_1stDerivative ()
{
	auto force = potential.Force(rij, rik);
	cout << force[0][0] << " " <<force[0][1] << " " << force[0][2] << "\t"
	<< force[1][0] << " " <<force[1][1] << " " << force[1][2] << endl;

	force = potential.Force(rij, rik, FiniteDifference_t());
	cout << force[0][0] << " " <<force[0][1] << " " << force[0][2] << "\t"
	<< force[1][0] << " " <<force[1][1] << " " << force[1][2] << endl;

	double energy = potential.Energy(rij, rik);
	double e;
	rij[0] += finiteDist;
	e = potential.Energy(rij, rik);
	cout << (e-energy)/finiteDist << " ";
	rij[0] -= finiteDist;
	rij[1] += finiteDist;
	e = potential.Energy(rij, rik);
	cout << (e-energy)/finiteDist << " ";
	rij[1] -= finiteDist;
	rij[2] += finiteDist;
	e = potential.Energy(rij, rik);
	cout << (e-energy)/finiteDist << "\t";
	rij[2] -= finiteDist;
	rik[0] += finiteDist;
	e = potential.Energy(rij, rik);
	cout << (e-energy)/finiteDist << " ";
	rik[0] -= finiteDist;
	rik[1] += finiteDist;
	e = potential.Energy(rij, rik);
	cout << (e-energy)/finiteDist << " ";
	rik[1] -= finiteDist;
	rik[2] += finiteDist;
	e = potential.Energy(rij, rik);
	cout << (e-energy)/finiteDist << "\n";
	rik[2] -= finiteDist;
}

void test_2ndDerivative ()
{
	auto hessian = potential.Hessian(rij, rik);
	cout << hessian[0][0][0] << " " << hessian[0][0][1] << " " << hessian[0][0][2] << "\t"
	<< hessian[1][0][0] << " " << hessian[1][0][1] << " " << hessian[1][0][2] << "\t"
	<< hessian[2][0][0] << " " << hessian[2][0][1] << " " << hessian[2][0][2] << "\n"
	<< hessian[0][1][0] << " " << hessian[0][1][1] << " " << hessian[0][1][2] << "\t"
	<< hessian[1][1][0] << " " << hessian[1][1][1] << " " << hessian[1][1][2] << "\t"
	<< hessian[2][1][0] << " " << hessian[2][1][1] << " " << hessian[2][1][2] << "\n"
	<< hessian[0][2][0] << " " << hessian[0][2][1] << " " << hessian[0][2][2] << "\t"
	<< hessian[1][2][0] << " " << hessian[1][2][1] << " " << hessian[1][2][2] << "\t"
	<< hessian[2][2][0] << " " << hessian[2][2][1] << " " << hessian[2][2][2] << "\n"<< endl;

	hessian = potential.Hessian(rij, rik, FiniteDifference_t());
	cout << hessian[0][0][0] << " " << hessian[0][0][1] << " " << hessian[0][0][2] << "\t"
	<< hessian[1][0][0] << " " << hessian[1][0][1] << " " << hessian[1][0][2] << "\t"
	<< hessian[2][0][0] << " " << hessian[2][0][1] << " " << hessian[2][0][2] << "\n"
	<< hessian[0][1][0] << " " << hessian[0][1][1] << " " << hessian[0][1][2] << "\t"
	<< hessian[1][1][0] << " " << hessian[1][1][1] << " " << hessian[1][1][2] << "\t"
	<< hessian[2][1][0] << " " << hessian[2][1][1] << " " << hessian[2][1][2] << "\n"
	<< hessian[0][2][0] << " " << hessian[0][2][1] << " " << hessian[0][2][2] << "\t"
	<< hessian[1][2][0] << " " << hessian[1][2][1] << " " << hessian[1][2][2] << "\t"
	<< hessian[2][2][0] << " " << hessian[2][2][1] << " " << hessian[2][2][2] << "\n"<< endl;

	// ij, jk
	auto force = potential.Force(rij, rik);
	auto forceI = Scale( Addition(force[0], force[1]), -1. );
	rij[0] -= finiteDist;
	rik[0] -= finiteDist;
	auto dForce = potential.Force(rij, rik);
	auto dForceI = Scale( Addition(dForce[0], dForce[1]), -1.);
	cout <<	(dForceI[0]-forceI[0])/finiteDist << " "
	<<	(dForceI[1]-forceI[1])/finiteDist << " "
	<<	(dForceI[2]-forceI[2])/finiteDist << "\t"
	<<	(dForce[0][0]-force[0][0])/finiteDist << " "
	<<	(dForce[0][1]-force[0][1])/finiteDist << " "
	<<	(dForce[0][2]-force[0][2])/finiteDist << "\t"
	<<	(dForce[1][0]-force[1][0])/finiteDist << " "
	<<	(dForce[1][1]-force[1][1])/finiteDist << " "
	<<	(dForce[1][2]-force[1][2])/finiteDist << "\n";
	rij[0] += finiteDist;
	rik[0] += finiteDist;
	rij[1] -= finiteDist;
	rik[1] -= finiteDist;
	dForce = potential.Force(rij, rik);
	dForceI = Scale( Addition(dForce[0], dForce[1]), -1.);
	cout <<	(dForceI[0]-forceI[0])/finiteDist << " "
	<<	(dForceI[1]-forceI[1])/finiteDist << " "
	<<	(dForceI[2]-forceI[2])/finiteDist << "\t"
	<<	(dForce[0][0]-force[0][0])/finiteDist << " "
	<<	(dForce[0][1]-force[0][1])/finiteDist << " "
	<<	(dForce[0][2]-force[0][2])/finiteDist << "\t"
	<<	(dForce[1][0]-force[1][0])/finiteDist << " "
	<<	(dForce[1][1]-force[1][1])/finiteDist << " "
	<<	(dForce[1][2]-force[1][2])/finiteDist << "\n";
	rij[1] += finiteDist;
	rik[1] += finiteDist;
	rij[2] -= finiteDist;
	rik[2] -= finiteDist;
	dForce = potential.Force(rij, rik);
	dForceI = Scale( Addition(dForce[0], dForce[1]), -1.);
	cout <<	(dForceI[0]-forceI[0])/finiteDist << " "
	<<	(dForceI[1]-forceI[1])/finiteDist << " "
	<<	(dForceI[2]-forceI[2])/finiteDist << "\t"
	<<	(dForce[0][0]-force[0][0])/finiteDist << " "
	<<	(dForce[0][1]-force[0][1])/finiteDist << " "
	<<	(dForce[0][2]-force[0][2])/finiteDist << "\t"
	<<	(dForce[1][0]-force[1][0])/finiteDist << " "
	<<	(dForce[1][1]-force[1][1])/finiteDist << " "
	<<	(dForce[1][2]-force[1][2])/finiteDist << "\n";
	rij[2] += finiteDist;
	rik[2] += finiteDist;

	rij[0] += finiteDist;
	dForce = potential.Force(rij, rik);
	dForceI = Scale( Addition(dForce[0], dForce[1]), -1.);
	cout <<	(dForceI[0]-forceI[0])/finiteDist << " "
	<<	(dForceI[1]-forceI[1])/finiteDist << " "
	<<	(dForceI[2]-forceI[2])/finiteDist << "\t"
	<<	(dForce[0][0]-force[0][0])/finiteDist << " "
	<<	(dForce[0][1]-force[0][1])/finiteDist << " "
	<<	(dForce[0][2]-force[0][2])/finiteDist << "\t"
	<<	(dForce[1][0]-force[1][0])/finiteDist << " "
	<<	(dForce[1][1]-force[1][1])/finiteDist << " "
	<<	(dForce[1][2]-force[1][2])/finiteDist << "\n";
	rij[0] -= finiteDist;
	rij[1] += finiteDist;
	dForce = potential.Force(rij, rik);
	dForceI = Scale( Addition(dForce[0], dForce[1]), -1.);
	cout <<	(dForceI[0]-forceI[0])/finiteDist << " "
	<<	(dForceI[1]-forceI[1])/finiteDist << " "
	<<	(dForceI[2]-forceI[2])/finiteDist << "\t"
	<<	(dForce[0][0]-force[0][0])/finiteDist << " "
	<<	(dForce[0][1]-force[0][1])/finiteDist << " "
	<<	(dForce[0][2]-force[0][2])/finiteDist << "\t"
	<<	(dForce[1][0]-force[1][0])/finiteDist << " "
	<<	(dForce[1][1]-force[1][1])/finiteDist << " "
	<<	(dForce[1][2]-force[1][2])/finiteDist << "\n";
	rij[1] -= finiteDist;
	rij[2] += finiteDist;
	dForce = potential.Force(rij, rik);
	dForceI = Scale( Addition(dForce[0], dForce[1]), -1.);
	cout <<	(dForceI[0]-forceI[0])/finiteDist << " "
	<<	(dForceI[1]-forceI[1])/finiteDist << " "
	<<	(dForceI[2]-forceI[2])/finiteDist << "\t"
	<<	(dForce[0][0]-force[0][0])/finiteDist << " "
	<<	(dForce[0][1]-force[0][1])/finiteDist << " "
	<<	(dForce[0][2]-force[0][2])/finiteDist << "\t"
	<<	(dForce[1][0]-force[1][0])/finiteDist << " "
	<<	(dForce[1][1]-force[1][1])/finiteDist << " "
	<<	(dForce[1][2]-force[1][2])/finiteDist << "\n";
	rij[2] -= finiteDist;

	rik[0] += finiteDist;
	dForce = potential.Force(rij, rik);
	dForceI = Scale( Addition(dForce[0], dForce[1]), -1.);
	cout <<	(dForceI[0]-forceI[0])/finiteDist << " "
	<<	(dForceI[1]-forceI[1])/finiteDist << " "
	<<	(dForceI[2]-forceI[2])/finiteDist << "\t"
	<<	(dForce[0][0]-force[0][0])/finiteDist << " "
	<<	(dForce[0][1]-force[0][1])/finiteDist << " "
	<<	(dForce[0][2]-force[0][2])/finiteDist << "\t"
	<<	(dForce[1][0]-force[1][0])/finiteDist << " "
	<<	(dForce[1][1]-force[1][1])/finiteDist << " "
	<<	(dForce[1][2]-force[1][2])/finiteDist << "\n";
	rik[0] -= finiteDist;
	rik[1] += finiteDist;
	dForce = potential.Force(rij, rik);
	dForceI = Scale( Addition(dForce[0], dForce[1]), -1.);
	cout <<	(dForceI[0]-forceI[0])/finiteDist << " "
	<<	(dForceI[1]-forceI[1])/finiteDist << " "
	<<	(dForceI[2]-forceI[2])/finiteDist << "\t"
	<<	(dForce[0][0]-force[0][0])/finiteDist << " "
	<<	(dForce[0][1]-force[0][1])/finiteDist << " "
	<<	(dForce[0][2]-force[0][2])/finiteDist << "\t"
	<<	(dForce[1][0]-force[1][0])/finiteDist << " "
	<<	(dForce[1][1]-force[1][1])/finiteDist << " "
	<<	(dForce[1][2]-force[1][2])/finiteDist << "\n";
	rik[1] -= finiteDist;
	rik[2] += finiteDist;
	dForce = potential.Force(rij, rik);
	dForceI = Scale( Addition(dForce[0], dForce[1]), -1.);
	cout <<	(dForceI[0]-forceI[0])/finiteDist << " "
	<<	(dForceI[1]-forceI[1])/finiteDist << " "
	<<	(dForceI[2]-forceI[2])/finiteDist << "\t"
	<<	(dForce[0][0]-force[0][0])/finiteDist << " "
	<<	(dForce[0][1]-force[0][1])/finiteDist << " "
	<<	(dForce[0][2]-force[0][2])/finiteDist << "\t"
	<<	(dForce[1][0]-force[1][0])/finiteDist << " "
	<<	(dForce[1][1]-force[1][1])/finiteDist << " "
	<<	(dForce[1][2]-force[1][2])/finiteDist << endl;
	rik[2] -= finiteDist;
}

int main ()
{
cout<<setprecision(10);
	test_1stDerivative();
cout << endl;
	test_2ndDerivative();
}

