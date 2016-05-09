#include "Potential_HarmonicDihedral.h"
#include <iomanip>
#include <iostream>

using namespace std;


Potential_HarmonicDihedral potential;
array3d rij{0.4, -0.4, 0.9};
array3d rjk{0.4, 0.38, -0.4};
array3d rkl{0.1, 0.7, 0.6};

constexpr double FiniteDist = 1.e-10;
constexpr double InverseFiniteDist = 1.e10;

void test_1stDerivative ()
{
	auto energy = potential.Energy(rij, rjk, rkl);
	cout << energy<<endl;

	auto force = potential.Force(rij, rjk, rkl);
	cout << "potentialPack:" << endl;
	cout << "fi: " << force[0][0] << " " << force[0][1] << " " << force[0][2] << endl;
	cout << "fj: " << -force[0][0]-force[1][0] << " " << -force[0][1]-force[1][1] << " " << -force[0][2]-force[1][2] << endl;
	cout << "fk: " << force[1][0]-force[2][0] << " " << force[1][1]-force[2][1] << " " << force[1][2]-force[2][2] << endl;
	cout << "fl: " << force[2][0] << " " << force[2][1] << " " << force[2][2] << endl;
	cout << endl;

	force = potential.Force(rij, rjk, rkl, FiniteDifference_t());
	cout << "potentialPack (finite difference):" << endl;
	cout << "fi: " << force[0][0] << " " << force[0][1] << " " << force[0][2] << endl;
	cout << "fj: " << -force[0][0]-force[1][0] << " " << -force[0][1]-force[1][1] << " " << -force[0][2]-force[1][2] << endl;
	cout << "fk: " << force[1][0]-force[2][0] << " " << force[1][1]-force[2][1] << " " << force[1][2]-force[2][2] << endl;
	cout << "fl: " << force[2][0] << " " << force[2][1] << " " << force[2][2] << endl;
	cout << endl;

	cout << "finite difference:" << endl;
	cout << "fi: ";
	rij[0] -= FiniteDist;
	cout << (potential.Energy(rij, rjk, rkl)-energy) * InverseFiniteDist << " ";
	rij[0] += FiniteDist;
	rij[1] -= FiniteDist;
	cout << (potential.Energy(rij, rjk, rkl)-energy) * InverseFiniteDist << " ";
	rij[1] += FiniteDist;
	rij[2] -= FiniteDist;
	cout << (potential.Energy(rij, rjk, rkl)-energy) * InverseFiniteDist << endl;
	rij[2] += FiniteDist;

	cout << "fj: ";
	rij[0] += FiniteDist;
	rjk[0] -= FiniteDist;
	cout << (potential.Energy(rij, rjk, rkl)-energy) * InverseFiniteDist << " ";
	rij[0] -= FiniteDist;
	rjk[0] += FiniteDist;
	rij[1] += FiniteDist;
	rjk[1] -= FiniteDist;
	cout << (potential.Energy(rij, rjk, rkl)-energy) * InverseFiniteDist << " ";
	rij[1] -= FiniteDist;
	rjk[1] += FiniteDist;
	rij[2] += FiniteDist;
	rjk[2] -= FiniteDist;
	cout << (potential.Energy(rij, rjk, rkl)-energy) * InverseFiniteDist << endl;
	rij[2] -= FiniteDist;
	rjk[2] += FiniteDist;

	cout << "fk: ";
	rjk[0] += FiniteDist;
	rkl[0] -= FiniteDist;
	cout << (potential.Energy(rij, rjk, rkl)-energy) * InverseFiniteDist << " ";
	rjk[0] -= FiniteDist;
	rkl[0] += FiniteDist;
	rjk[1] += FiniteDist;
	rkl[1] -= FiniteDist;
	cout << (potential.Energy(rij, rjk, rkl)-energy) * InverseFiniteDist << " ";
	rjk[1] -= FiniteDist;
	rkl[1] += FiniteDist;
	rjk[2] += FiniteDist;
	rkl[2] -= FiniteDist;
	cout << (potential.Energy(rij, rjk, rkl)-energy) * InverseFiniteDist << endl;
	rjk[2] -= FiniteDist;
	rkl[2] += FiniteDist;

	cout << "fl: ";
	rkl[0] += FiniteDist;
	cout << (potential.Energy(rij, rjk, rkl)-energy) * InverseFiniteDist << " ";
	rkl[0] -= FiniteDist;
	rkl[1] += FiniteDist;
	cout << (potential.Energy(rij, rjk, rkl)-energy) * InverseFiniteDist << " ";
	rkl[1] -= FiniteDist;
	rkl[2] += FiniteDist;
	cout << (potential.Energy(rij, rjk, rkl)-energy) * InverseFiniteDist << endl;
	rkl[2] -= FiniteDist;
}


int main ()
{
	test_1stDerivative();
}

