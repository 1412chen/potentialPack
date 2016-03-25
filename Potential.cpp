#include "Potential.h"
#include <utility>
#include <cmath>
#include <iostream>

using namespace std;

constexpr double FiniteDistance = 1.e-10;
constexpr double InverseFiniteDistance = 1./FiniteDistance;

double
Potential::Energy (
	const array3d& rij
) const noexcept
{
	return ObjectiveFunction( Norm(rij) );
}


double
Potential::Energy (
	const array3d& rij,
	const array3d& rik
) const noexcept
{
	auto rij_ = Norm(rij);
	auto rik_ = Norm(rik);
	auto cosTheta = CosTheta (
		UnitVector(rij, 1./rij_),
		UnitVector(rik, 1./rik_)
	);
	return ObjectiveFunction( cosTheta, rij_, rik_ );	
}


double
Potential::Energy (
	const array3d& rij,
	const array3d& rjk,
	const array3d& rkl
) const noexcept
{
	auto rij_ = Norm(rij);
	auto rjk_ = Norm(rjk);
	auto rkl_ = Norm(rkl);
	auto nij = UnitVector(rij, 1./rij_);
	auto njk = UnitVector(rjk, 1./rjk_);
	auto nkl = UnitVector(rkl, 1./rkl_);
	auto cosThetaIJK = CosTheta(nij, njk);
	auto cosThetaJKL = CosTheta(njk, nkl);
	auto cosPhi = CosPhi(nij, njk, nkl);
	return ObjectiveFunction(cosPhi, rij_, rjk_, rkl_, cosThetaIJK, cosThetaJKL);
}

//---------------------------------------------------------------------------//

array3d
Potential::Force (
	const array3d& rij
) const noexcept
{
	auto _1stDeriFunc = [this](double rij) {
		return this->_1stDerivative(rij);
	};
	return ForceImp(rij, _1stDeriFunc);
}


array3d
Potential::Force (
	const array3d& rij,
	FiniteDifference_t
) const noexcept
{
	auto _1stDeriFunc = [this](double rij_) {
		return this->Potential::_1stDerivative(rij_);
	};
	return ForceImp(rij, _1stDeriFunc);
}


std::array<array3d, 2>
Potential::Force (
	const array3d& rij,
	const array3d& rik
) const noexcept
{
	auto _1stDeriFunc = [this](double cosTheta, double rij_, double rik_) {
		return this->_1stDerivative(cosTheta, rij_, rik_);
	};
	return ForceImp(rij, rik, _1stDeriFunc);
}


std::array<array3d, 2>
Potential::Force (
	const array3d& rij,
	const array3d& rik,
	FiniteDifference_t
) const noexcept
{
	auto _1stDeriFunc = [this](double cosTheta, double rij_, double rik_) {
		return this->Potential::_1stDerivative(cosTheta, rij_, rik_);
	};
	return ForceImp(rij, rik, _1stDeriFunc);
}


/*std::array<array3d, 3>
Potential::Force (
	const array3d& rij,
	const array3d& rjk,
	const array3d& rkl
) const noexcept
{
	auto _1stDeriFunc = [this](double cosPhi, double rij_, double rjk_, double rkl_, double cosThetaIJK, double cosThetaJKL) {
		return this->_1stDerivative(cosPhi, rij_, rjk_, rkl_, cosThetaIJK, cosThetaJKL);
	};
	return ForceImp(rij, rjk, rkl, _1stDeriFunc);
}


std::array<array3d, 3>
Potential::Force (
	const array3d& rij,
	const array3d& rjk,
	const array3d& rkl,
	FiniteDifference_t
) const noexcept
{
	auto _1stDeriFunc = [this](double cosPhi, double rij_, double rjk_, double rkl_, double cosThetaIJK, double cosThetaJKL) {
		return this->Potential::_1stDerivative(cosPhi, rij_, rjk_, rkl_, cosThetaIJK, cosThetaJKL);
	};
	return ForceImp(rij, rjk, rkl, _1stDeriFunc);
}*/


matrix3d
Potential::Hessian (
	const array3d& rij
) const noexcept
{
	auto _1stDeriFunc = [this](double rij_) {
		return this->_1stDerivative(rij_);
	};
	auto _2ndDeriFunc = [this](double rij_) {
		return this->_2ndDerivative(rij_);
	};
	return HessianImp(rij, _1stDeriFunc, _2ndDeriFunc);
}


matrix3d
Potential::Hessian (
	const array3d& rij,
	FiniteDifference_t
) const noexcept
{
	auto _1stDeriFunc = [this](double rij_) {
		return this->Potential::_1stDerivative(rij_);
	};
	auto _2ndDeriFunc = [this](double rij_) {
		return this->Potential::_2ndDerivative(rij_);
	};
	return HessianImp(rij, _1stDeriFunc, _2ndDeriFunc);
}


std::array<matrix3d, 3>
Potential::Hessian (
	const array3d& rij,
	const array3d& rik
) const noexcept
{
	auto _1stDeriFunc = [this](double cosTheta, double rij_, double rik_) {
		return this->_1stDerivative(cosTheta, rij_, rik_);
	};
	auto _2ndDeriFunc = [this](double cosTheta, double rij_, double rik_) {
		return this->_2ndDerivative(cosTheta, rij_, rik_);
	};
	return HessianImp(rij, rik, _1stDeriFunc, _2ndDeriFunc);
}


std::array<matrix3d, 3>
Potential::Hessian (
	const array3d& rij,
	const array3d& rik,
	FiniteDifference_t
) const noexcept
{
	auto _1stDeriFunc = [this](double cosTheta, double rij_, double rik_) {
		return this->Potential::_1stDerivative(cosTheta, rij_, rik_);
	};
	auto _2ndDeriFunc = [this](double cosTheta, double rij_, double rik_) {
		return this->Potential::_2ndDerivative(cosTheta, rij_, rik_);
	};
	return HessianImp(rij, rik, _1stDeriFunc, _2ndDeriFunc);
}


/*std::array<matrix3d, 6>
Potential::Hessian (
	const array3d& rij,
	const array3d& rjk,
	const array3d& rkl
) const noexcept
{
	auto _1stDeriFunc = [this](double cosPhi, double rij_, double rjk_, double rkl_, double cosThetaIJK, double cosThetaJKL) {
		return this->_1stDerivative(cosPhi, rij_, rjk_, rkl_, cosThetaIJK, cosThetaJKL);
	};
	auto _2ndDeriFunc = [this](double cosPhi, double rij_, double rjk_, double rkl_, double cosThetaIJK, double cosThetaJKL) {
		return this->_2ndDerivative(cosPhi, rij_, rjk_, rkl_, cosThetaIJK, cosThetaJKL);
	};
	return HessianImp(rij, rjk, rkl, _1stDeriFunc, _2ndDeriFunc);
}


std::array<matrix3d, 6>
Potential::Hessian (
	const array3d& rij,
	const array3d& rjk,
	const array3d& rkl,
	FiniteDifference_t
) const noexcept
{
	auto _1stDeriFunc = [this](double cosPhi, double rij_, double rjk_, double rkl_, double cosThetaIJK, double cosThetaJKL) {
		return this->Potential::_1stDerivative(cosPhi, rij_, rjk_, rkl_, cosThetaIJK, cosThetaJKL);
	};
	auto _2ndDeriFunc = [this](double cosPhi, double rij_, double rjk_, double rkl_, double cosThetaIJK, double cosThetaJKL) {
		return this->Potential::_2ndDerivative(cosPhi, rij_, rjk_, rkl_, cosThetaIJK, cosThetaJKL);
	};
	return HessianImp(rij, rjk, rkl, _1stDeriFunc, _2ndDeriFunc);
}*/

//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//

array3d
Potential::ForceImp (
	const array3d& rij,
	std::function< Pair_1stDerivative_t(double) > _1stDeriFunc
) const noexcept
{
	auto rij_ = Norm(rij);
	auto _1stDeri = _1stDeriFunc(rij_);
	return Scale( UnitVector(rij, 1./rij_), _1stDeri );
}


std::array<array3d, 2>
Potential::ForceImp (
	const array3d& rij,
	const array3d& rik,
	std::function< Angle_1stDerivative_t(double, double, double) > _1stDeriFunc
) const noexcept
{
	using namespace Angle_Namespace;

	auto rij_ = Norm(rij);
	auto inverseRij = 1./rij_;
	auto rik_ = Norm(rik);
	auto inverseRik = 1./rik_;
	auto nij = UnitVector(rij, inverseRij);
	auto nik = UnitVector(rik, inverseRik);
	auto cosTheta = CosTheta(nij, nik);
	auto _1stDeri = _1stDeriFunc(cosTheta, rij_, rik_);
	auto termIJ_1 = _1stDeri[DCosTheta] * inverseRij;
	auto termIJ_2 = _1stDeri[DRij] - termIJ_1*cosTheta;
	auto termIK_1 = _1stDeri[DCosTheta] * inverseRik;
	auto termIK_2 = _1stDeri[DRik] - termIK_1*cosTheta;
	auto forceIJ = Addition (
		Scale(nik, termIJ_1),
		Scale(nij, termIJ_2)
	);
	auto forceIK = Addition (
		Scale(nij, termIK_1),
		Scale(nik, termIK_2)
	);
	return {{ move(forceIJ), move(forceIK) }};
}


/*std::array<array3d, 3>
Potential::ForceImp (
	const array3d& rij,
	const array3d& rjk,
	const array3d& rkl,
	std::function< Dihedral_1stDerivative_t(double, double, double, double, double, double) _1stDeriFunc
) const noexcept
{
	using namespace Dihedral_Namespace;

	;
}*/

//---------------------------------------------------------------------------//

matrix3d
Potential::HessianImp (
	const array3d& rij,
	std::function< Pair_1stDerivative_t(double) > _1stDeriFunc,
	std::function< Pair_2ndDerivative_t(double) > _2ndDeriFunc
) const noexcept
{
	auto rij_ = Norm(rij);
	auto inverseRij = 1./rij_;
	auto nij = UnitVector(rij, inverseRij);
	auto _1stDeri = _1stDeriFunc(rij_);
	auto _2ndDeri = _2ndDeriFunc(rij_);
	auto nijnij = OuterDot(nij, nij);
	return Addition (
		Scale(nijnij, -_2ndDeri),
		Scale(Minus(nijnij, Identity()), _1stDeri*inverseRij)
	);
}


std::array<matrix3d, 3>
Potential::HessianImp (
	const array3d& rij,
	const array3d& rik,
	std::function< Angle_1stDerivative_t(double, double, double) > _1stDeriFunc,
	std::function< Angle_2ndDerivative_t(double, double, double) > _2ndDeriFunc
) const noexcept
{
	using namespace Angle_Namespace;

	auto rij_ = Norm(rij);
	auto inverseRij = 1. / rij_;
	auto rik_ = Norm(rik);
	auto inverseRik = 1. / rik_;
	auto nij = UnitVector(rij, inverseRij);
	auto nik = UnitVector(rik, inverseRik);
	auto cosTheta = CosTheta(nij, nik);
	auto nij_nij = OuterDot(nij, nij);
	auto nik_nik = OuterDot(nik, nik);

	auto cosTheta_drj = Scale( Addition(nik, Scale(nij, -cosTheta)), inverseRij );
	auto cosTheta_drk = Scale( Addition(nij, Scale(nik, -cosTheta)), inverseRik );
	auto cosTheta_dri = Scale( Addition( cosTheta_drj, cosTheta_drk ), -1. );

	auto cosTheta_dri_drj = Scale( Addition (
		OuterDot( nij, cosTheta_drj ),
		OuterDot( cosTheta_dri, -nij ),
		Scale( Minus(nik_nik, Identity()), inverseRik ),
		Scale( Minus(nij_nij, Identity()), -cosTheta*inverseRij )
	), inverseRij );
	auto cosTheta_drj_drk = Scale( Addition (
		Scale( Minus(Identity(), nij_nij), inverseRij ),
		OuterDot( cosTheta_drj, -nik )
	), inverseRik );
	auto cosTheta_dri_drk = Scale( Addition (
		OuterDot( nik, cosTheta_drk ),
		OuterDot( cosTheta_dri, -nik ),
		Scale( Minus(nij_nij, Identity()), inverseRij ),
		Scale( Minus(nik_nik, Identity()), -cosTheta*inverseRik )
	), inverseRik );

	auto rij_dri_drj = Scale( Minus(nij_nij, Identity()), inverseRij );
	auto rik_dri_drk = Scale( Minus(nik_nik, Identity()), inverseRik );

	//-------------------------------------//

	auto _1stDeri = _1stDeriFunc(cosTheta, rij_, rik_);
	auto _2ndDeri = _2ndDeriFunc(cosTheta, rij_, rik_);
cout<<_1stDeri[0]<<" "<<_1stDeri[1]<<" "<<_1stDeri[2]<<"\n"
<<_2ndDeri[0]<<" "<<_2ndDeri[1]<<" "<<_2ndDeri[2]<<" "<<_2ndDeri[3]<<" "<<_2ndDeri[4]<<" "<<_2ndDeri[5]<<endl;

	//-------------------------------------//

	auto termIJ_1L = Addition (
		Scale(nij, -_2ndDeri[DRij_DRij]),
		Scale(nik, -_2ndDeri[DRij_DRik]),
		Scale(cosTheta_dri, _2ndDeri[DCosTheta_DRij])
	);
	auto& termIJ_1R = nij;
	auto termIJ_2 = Scale(rij_dri_drj, _1stDeri[DRij]);
	auto termIJ_3L = Addition (
		Scale(nij, -_2ndDeri[DCosTheta_DRij]),
		Scale(nik, -_2ndDeri[DCosTheta_DRik]),
		Scale(cosTheta_dri, _2ndDeri[DCosTheta_DCosTheta])
	);
	auto& termIJ_3R = cosTheta_drj;
	auto termIJ_4 = Scale(cosTheta_dri_drj, _1stDeri[DCosTheta]);
	auto hessianIJ = Addition (
		OuterDot(termIJ_1L, termIJ_1R),
		termIJ_2,
		OuterDot(termIJ_3L, termIJ_3R),
		termIJ_4
	);

	//-------------------------------------//

	auto termIK_1L = Addition (
		Scale(nij, -_2ndDeri[DRij_DRik]),
		Scale(nik, -_2ndDeri[DRij_DRij]),
		Scale(cosTheta_dri, _2ndDeri[DCosTheta_DRik])
	);
	auto& termIK_1R = nik;
	auto termIK_2 = Scale(rik_dri_drk, _1stDeri[DRik]);
	auto& termIK_3L = termIJ_3L;/*Addition (
		Scale(cosTheta_dri, _2ndDeri[DCosTheta_DCosTheta]),
		Scale(nij, _2ndDeri[DCosTheta_DRij]),
		Scale(nik, _2ndDeri[DCosTheta_DRik])
	);*/
	auto& termIK_3R = cosTheta_drk;
	auto termIK_4 = Scale( cosTheta_dri_drk, _1stDeri[DCosTheta] );
	auto hessianIK = Addition (
		OuterDot(termIK_1L, termIK_1R),
		termIK_2,
		OuterDot(termIK_3L, termIK_3R),
		termIK_4
	);

	//-------------------------------------//

	auto termJK_1L = Addition (
		Scale(nij, _2ndDeri[DRij_DRik]),
		Scale(cosTheta_drj, _2ndDeri[DCosTheta_DRik])
	);
	auto& termJK_1R = nik;
	auto termJK_2L = Addition (
		Scale(cosTheta_drj, _2ndDeri[DCosTheta_DCosTheta]),
		Scale(nij, _2ndDeri[DCosTheta_DRij])
	);
	auto& termJK_2R = cosTheta_drk;
	auto termJK_3 = Scale(cosTheta_drj_drk, _1stDeri[DCosTheta]);
	auto hessianJK = Addition (
		OuterDot(termJK_1L, termJK_1R),
		OuterDot(termJK_2L, termJK_2R),
		termJK_3
	);

	return { move(hessianIJ), move(hessianIK), move(hessianJK) };
}


/*std::array<matrix3d, 6>
Potential::HessianImp (
	const array3d& rij,
	const array3d& rjk,
	const array3d& rkl,
	std::function< Dihedral_1stDerivative_t(double, double, double, double, double, double) > _1stDeriFunc,
	std::function< Dihedral_2ndDerivative_t(double, double, double, double, double, double) > _2ndDeriFunc
) const noexcept
{
	using namespace Dihedral_Namespace;

	;
}*/

//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//

double
Potential::Norm (
	const array3d& rij
) const noexcept
{
	return sqrt( Square(rij) );
}


array3d
Potential::UnitVector (
	const array3d& rij
) const noexcept
{
	auto rij_ = Norm(rij);
	return UnitVector(rij, 1./rij_);
}


array3d
Potential::UnitVector (
	const array3d& rij,
	double inverseRij
) const noexcept
{
	return Scale(rij, inverseRij);
}


double
Potential::CosTheta (
	const array3d& nij,
	const array3d& nik
) const noexcept
{
	return Dot(nij, nik);
}


double
Potential::CosPhi (
	const array3d& nij,
	const array3d& njk,
	const array3d& nkl
) const noexcept
{
	auto normalIJK = Cross(nij, njk);
	auto normalJKL = Cross(njk, nkl);
	return CosTheta(normalIJK, normalJKL);
}

//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//

Potential::Pair_1stDerivative_t
Potential::_1stDerivative (
	double rij
) const noexcept
{
	auto currentObjFun = ObjectiveFunction(rij);
	auto incrementObjFun = ObjectiveFunction(rij+FiniteDistance);
	return (incrementObjFun-currentObjFun) * InverseFiniteDistance;
}


Potential::Angle_1stDerivative_t
Potential::_1stDerivative (
	double cosTheta,
	double rij,
	double rik
) const noexcept
{
	auto currentObjFun = ObjectiveFunction(cosTheta, rij, rik);
	auto incrementObjFun_CosTheta = ObjectiveFunction(cosTheta+FiniteDistance, rij, rik);
	auto incrementObjFun_Rij = ObjectiveFunction(cosTheta, rij+FiniteDistance, rik);
	auto incrementObjFun_Rik = ObjectiveFunction(cosTheta, rij, rik+FiniteDistance);
	return {
		(incrementObjFun_CosTheta-currentObjFun) * InverseFiniteDistance,
		(incrementObjFun_Rij-currentObjFun) * InverseFiniteDistance,
		(incrementObjFun_Rik-currentObjFun) * InverseFiniteDistance
	};
}


Potential::Dihedral_1stDerivative_t
Potential::_1stDerivative (
	double cosPhi,
	double rij,
	double rjk,
	double rkl,
	double cosThetaIJK,
	double cosThetaJKL
) const noexcept
{
	auto currentObjFun = ObjectiveFunction(cosPhi, rij, rjk, rkl, cosThetaIJK, cosThetaJKL);
	auto incrementObjFun_CosPhi = ObjectiveFunction(cosPhi+FiniteDistance, rij, rjk, rkl, cosThetaIJK, cosThetaJKL);
	auto incrementObjFun_Rij = ObjectiveFunction(cosPhi, rij+FiniteDistance, rjk, rkl, cosThetaIJK, cosThetaJKL);
	auto incrementObjFun_Rjk = ObjectiveFunction(cosPhi, rij, rjk+FiniteDistance, rkl, cosThetaIJK, cosThetaJKL);
	auto incrementObjFun_Rkl = ObjectiveFunction(cosPhi, rij, rjk, rkl+FiniteDistance, cosThetaIJK, cosThetaJKL);
	auto incrementObjFun_CosThetaIJK = ObjectiveFunction(cosPhi, rij, rjk, rkl, cosThetaIJK+FiniteDistance, cosThetaJKL);
	auto incrementObjFun_CosThetaJKL = ObjectiveFunction(cosPhi, rij, rjk, rkl, cosThetaIJK, cosThetaJKL+FiniteDistance);
	return {
		(incrementObjFun_CosPhi-currentObjFun) * InverseFiniteDistance,
		(incrementObjFun_Rij-currentObjFun) * InverseFiniteDistance,
		(incrementObjFun_Rjk-currentObjFun) * InverseFiniteDistance,
		(incrementObjFun_Rkl-currentObjFun) * InverseFiniteDistance,
		(incrementObjFun_CosThetaIJK-currentObjFun) * InverseFiniteDistance,
		(incrementObjFun_CosThetaJKL-currentObjFun) * InverseFiniteDistance,
	};
}


Potential::Pair_2ndDerivative_t
Potential::_2ndDerivative (
	double rij
) const noexcept
{
	auto current1stDeri = _1stDerivative(rij);
	auto increment1stDeri = _1stDerivative(rij+FiniteDistance);
	return (increment1stDeri-current1stDeri) * InverseFiniteDistance;
}


Potential::Angle_2ndDerivative_t
Potential::_2ndDerivative (
	double cosTheta,
	double rij,
	double rik
) const noexcept
{
	auto current1stDeri = _1stDerivative(cosTheta, rij, rik);
	auto increment1stDeri_CosTheta = _1stDerivative(cosTheta+FiniteDistance, rij, rik);
	auto increment1stDeri_Rij = _1stDerivative(cosTheta, rij+FiniteDistance, rik);
	auto increment1stDeri_Rik = _1stDerivative(cosTheta, rij, rik+FiniteDistance);
	return {
		(increment1stDeri_CosTheta[0]-current1stDeri[0]) * InverseFiniteDistance,
		(increment1stDeri_CosTheta[1]-current1stDeri[1]) * InverseFiniteDistance,
		(increment1stDeri_CosTheta[2]-current1stDeri[2]) * InverseFiniteDistance,
		(increment1stDeri_Rij[1]-current1stDeri[1]) * InverseFiniteDistance,
		(increment1stDeri_Rij[2]-current1stDeri[2]) * InverseFiniteDistance,
		(increment1stDeri_Rik[2]-current1stDeri[2]) * InverseFiniteDistance
	};
}


Potential::Dihedral_2ndDerivative_t
Potential::_2ndDerivative (
	double cosPhi,
	double rij,
	double rjk,
	double rkl,
	double cosThetaIJK,
	double cosThetaJKL
) const noexcept
{
	auto current1stDeri = _1stDerivative(cosPhi, rij, rjk, rkl, cosThetaIJK, cosThetaJKL);
	auto increment1stDeri_CosPhi = _1stDerivative(cosPhi+FiniteDistance, rij, rjk, rkl, cosThetaIJK, cosThetaJKL);
	auto increment1stDeri_Rij = _1stDerivative(cosPhi, rij+FiniteDistance, rjk, rkl, cosThetaIJK, cosThetaJKL);
	auto increment1stDeri_Rjk = _1stDerivative(cosPhi, rij, rjk+FiniteDistance, rkl, cosThetaIJK, cosThetaJKL);
	auto increment1stDeri_Rkl = _1stDerivative(cosPhi, rij, rjk, rkl+FiniteDistance, cosThetaIJK, cosThetaJKL);
	auto increment1stDeri_CosThetaIJK = _1stDerivative(cosPhi, rij, rjk, rkl, cosThetaIJK+FiniteDistance, cosThetaJKL);
	auto increment1stDeri_CosThetaJKL = _1stDerivative(cosPhi, rij, rjk, rkl, cosThetaIJK, cosThetaJKL+FiniteDistance);
	return {
		(increment1stDeri_CosPhi[0]-current1stDeri[0]) * InverseFiniteDistance,
		(increment1stDeri_CosPhi[1]-current1stDeri[1]) * InverseFiniteDistance,
		(increment1stDeri_CosPhi[2]-current1stDeri[2]) * InverseFiniteDistance,
		(increment1stDeri_CosPhi[3]-current1stDeri[3]) * InverseFiniteDistance,
		(increment1stDeri_CosPhi[4]-current1stDeri[4]) * InverseFiniteDistance,
		(increment1stDeri_CosPhi[5]-current1stDeri[5]) * InverseFiniteDistance,
		(increment1stDeri_Rij[1]-current1stDeri[1]) * InverseFiniteDistance,
		(increment1stDeri_Rij[2]-current1stDeri[2]) * InverseFiniteDistance,
		(increment1stDeri_Rij[3]-current1stDeri[3]) * InverseFiniteDistance,
		(increment1stDeri_Rij[4]-current1stDeri[4]) * InverseFiniteDistance,
		(increment1stDeri_Rij[5]-current1stDeri[5]) * InverseFiniteDistance,
		(increment1stDeri_Rjk[2]-current1stDeri[2]) * InverseFiniteDistance,
		(increment1stDeri_Rjk[3]-current1stDeri[3]) * InverseFiniteDistance,
		(increment1stDeri_Rjk[4]-current1stDeri[4]) * InverseFiniteDistance,
		(increment1stDeri_Rjk[5]-current1stDeri[5]) * InverseFiniteDistance,
		(increment1stDeri_Rkl[3]-current1stDeri[3]) * InverseFiniteDistance,
		(increment1stDeri_Rkl[4]-current1stDeri[4]) * InverseFiniteDistance,
		(increment1stDeri_Rkl[5]-current1stDeri[5]) * InverseFiniteDistance,
		(increment1stDeri_CosThetaIJK[4]-current1stDeri[4]) * InverseFiniteDistance,
		(increment1stDeri_CosThetaIJK[5]-current1stDeri[5]) * InverseFiniteDistance,
		(increment1stDeri_CosThetaJKL[5]-current1stDeri[5]) * InverseFiniteDistance,
	};
}

