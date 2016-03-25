#include "Potential_StillingerWeber.h"
#include "Potential_HarmonicBond.h"
#include "Potential_LennardJones.h"
#include "Potential_Buckingham.h"
#include <cmath>
#include <iomanip>
#include <iostream>

using namespace std;

//Potential_StillingerWeber potential;
//Potential_HarmonicBond potential;
//Potential_LennardJones potential;
Potential_Buckingham potential;
array3d rij{-0.9, 0.4, -0.1};
constexpr double finiteDist = 1.e-10;


void test_1stDerivative ()
{
	auto force = potential.Force(rij);
	cout << force[0] << " " << force[1] << " " << force[2] << endl;

	force = potential.Force(rij, FiniteDifference_t());
	cout << force[0] << " " << force[1] << " " << force[2] << endl;

	double energy = potential.Energy(rij);
	double e;
	rij[0] += finiteDist;
	e = potential.Energy(rij);
	cout << (e-energy)/finiteDist << " ";
	rij[0] -= finiteDist;
	rij[1] += finiteDist;
	e = potential.Energy(rij);
	cout << (e-energy)/finiteDist << " ";
	rij[1] -= finiteDist;
	rij[2] += finiteDist;
	e = potential.Energy(rij);
	cout << (e-energy)/finiteDist << endl;
	rij[2] -= finiteDist;
}

void test_2ndDerivative ()
{
	auto hessian = potential.Hessian(rij);
	cout << hessian[0][0] << " " << hessian[0][1] << " " << hessian[0][2] << "\n"
	<< hessian[1][0] << " " << hessian[1][1] << " " << hessian[1][2] << "\n"
	<< hessian[2][0] << " " << hessian[2][1] << " " << hessian[2][2] << endl<<endl;

	hessian = potential.Hessian(rij, FiniteDifference_t());
	cout << hessian[0][0] << " " << hessian[0][1] << " " << hessian[0][2] << "\n"
	<< hessian[1][0] << " " << hessian[1][1] << " " << hessian[1][2] << "\n"
	<< hessian[2][0] << " " << hessian[2][1] << " " << hessian[2][2] << endl<<endl;

	// compare to drjdri
	auto force = potential.Force(rij);
	rij[0] -= finiteDist;
	auto dForce = potential.Force(rij);
	cout << (dForce[0]-force[0])/finiteDist << " "
		<< (dForce[1]-force[1])/finiteDist << " "
		<< (dForce[2]-force[2])/finiteDist << endl;
	rij[0] += finiteDist;
	rij[1] -= finiteDist;
	dForce = potential.Force(rij);
	cout << (dForce[0]-force[0])/finiteDist << " "
		<< (dForce[1]-force[1])/finiteDist << " "
		<< (dForce[2]-force[2])/finiteDist << endl;
	rij[1] += finiteDist;
	rij[2] -= finiteDist;
	dForce = potential.Force(rij);
	cout << (dForce[0]-force[0])/finiteDist << " "
		<< (dForce[1]-force[1])/finiteDist << " "
		<< (dForce[2]-force[2])/finiteDist << endl;
	rij[2] += finiteDist;
}

int main ()
{
cout << setprecision(10);
	test_1stDerivative();
cout << endl;
	test_2ndDerivative();
}

