#include "PotentialDihedral.h"
#include <iostream>
using namespace std;

constexpr double FiniteDistance = 1.e-13;
constexpr double InverseFiniteDistance = 1./FiniteDistance;


double
PotentialDihedral::Energy (
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
	auto cosPhi = CosPhi(rij, rjk, rkl);
	return ObjectiveFunction(cosPhi, rij_, rjk_, rkl_, cosThetaIJK, cosThetaJKL);
}

//---------------------------------------------------------------------------//

std::array<array3d, 3>
PotentialDihedral::Force (
	const array3d& rij,
	const array3d& rjk,
	const array3d& rkl
) const noexcept
{
	auto _1stDeriFunc = [this](double cosPhi, double rij_, double rjk_, double rkl_, double cosThetaIJK, double cosThetaJKL) noexcept {
		return this->_1stDerivative(cosPhi, rij_, rjk_, rkl_, cosThetaIJK, cosThetaJKL);
	};
	return ForceImp(rij, rjk, rkl, _1stDeriFunc);
}


std::array<array3d, 3>
PotentialDihedral::Force (
	const array3d& rij,
	const array3d& rjk,
	const array3d& rkl,
	FiniteDifference_t
) const noexcept
{
	auto _1stDeriFunc = [this](double cosPhi, double rij_, double rjk_, double rkl_, double cosThetaIJK, double cosThetaJKL) noexcept {
		return this->PotentialDihedral::_1stDerivative(cosPhi, rij_, rjk_, rkl_, cosThetaIJK, cosThetaJKL);
	};
	return ForceImp(rij, rjk, rkl, _1stDeriFunc);
}

//---------------------------------------------------------------------------//

std::array<matrix3d, 6>
PotentialDihedral::Hessian (
	const array3d& rij,
	const array3d& rjk,
	const array3d& rkl
) const noexcept
{
	auto _1stDeriFunc = [this](double cosPhi, double rij_, double rjk_, double rkl_, double cosThetaIJK, double cosThetaJKL) noexcept {
		return this->_1stDerivative(cosPhi, rij_, rjk_, rkl_, cosThetaIJK, cosThetaJKL);
	};
	auto _2ndDeriFunc = [this](double cosPhi, double rij_, double rjk_, double rkl_, double cosThetaIJK, double cosThetaJKL) noexcept {
		return this->_2ndDerivative(cosPhi, rij_, rjk_, rkl_, cosThetaIJK, cosThetaJKL);
	};
	return HessianImp(rij, rjk, rkl, _1stDeriFunc, _2ndDeriFunc);
}


std::array<matrix3d, 6>
PotentialDihedral::Hessian (
	const array3d& rij,
	const array3d& rjk,
	const array3d& rkl,
	FiniteDifference_t
) const noexcept
{
	auto _1stDeriFunc = [this](double cosPhi, double rij_, double rjk_, double rkl_, double cosThetaIJK, double cosThetaJKL) noexcept {
		return this->PotentialDihedral::_1stDerivative(cosPhi, rij_, rjk_, rkl_, cosThetaIJK, cosThetaJKL);
	};
	auto _2ndDeriFunc = [this](double cosPhi, double rij_, double rjk_, double rkl_, double cosThetaIJK, double cosThetaJKL) noexcept {
		return this->PotentialDihedral::_2ndDerivative(cosPhi, rij_, rjk_, rkl_, cosThetaIJK, cosThetaJKL);
	};
	return HessianImp(rij, rjk, rkl, _1stDeriFunc, _2ndDeriFunc);
}

//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//

std::array<array3d, 3>
PotentialDihedral::ForceImp (
	const array3d& rij,
	const array3d& rjk,
	const array3d& rkl,
	_1stD_Func _1stDeriFunc
) const noexcept
{
	using namespace Dihedral_Namespace;

	auto rij_ = Norm(rij);
	auto inverseRij_ = 1./rij_;
	auto nij = UnitVector(rij, inverseRij_);
	auto rjk_ = Norm(rjk);
	auto inverseRjk_ = 1./rjk_;
	auto njk = UnitVector(rjk, inverseRjk_);
	auto rkl_ = Norm(rkl);
	auto inverseRkl_ = 1./rkl_;
	auto nkl = UnitVector(rkl, inverseRkl_);
	auto cosThetaIJK = CosTheta(nij, njk);
	auto cosThetaJKL = CosTheta(njk, nkl);
	auto nijk = Cross(rij, rjk);
	auto njkl = Cross(rjk, rkl);
	auto inverseNijk_ = 1./Norm(nijk);
	auto inverseNjkl_ = 1./Norm(njkl);
	auto cosPhi = CosTheta(UnitVector(nijk, inverseNijk_), UnitVector(njkl, inverseNjkl_));
	auto _1stDeri = _1stDeriFunc(cosPhi, rij_, rjk_, rkl_, cosThetaIJK, cosThetaJKL);

	auto inverseNijk_Sq = inverseNijk_*inverseNijk_;
	auto inverseNjkl_Sq = inverseNjkl_*inverseNjkl_;
	auto inverseNijk_Njkl_ = inverseNijk_*inverseNjkl_;
	auto rijCnijk = Cross(rij, nijk);
	auto rijCnjkl = Cross(rij, njkl);
	auto rjkCnijk = Cross(rjk, nijk);
	auto rjkCnjkl = Cross(rjk, njkl);
	auto rklCnijk = Cross(rkl, nijk);
	auto rklCnjkl = Cross(rkl, njkl);

	auto forceIJ = Addition (
		Scale(nij, -_1stDeri[DRij] +
			_1stDeri[DCosThetaIJK]*cosThetaIJK*inverseRij_),
		Scale(njk, -_1stDeri[DCosThetaIJK]*inverseRij_),
		Scale(rjkCnijk, _1stDeri[DCosPhi]*cosPhi*inverseNijk_Sq),
		Scale(rjkCnjkl, -_1stDeri[DCosPhi]*inverseNijk_Njkl_)
	);
	auto forceKJ = Addition (
		Scale(njk, _1stDeri[DRjk] -
			_1stDeri[DCosThetaIJK]*cosThetaJKL * inverseRjk_),
		Scale(nkl, _1stDeri[DCosThetaJKL]*inverseRjk_),
		Scale(rijCnijk, _1stDeri[DCosPhi]*cosPhi*inverseNijk_Sq),
		Scale(rijCnjkl, -_1stDeri[DCosPhi]*inverseNijk_Njkl_),
		Scale(rklCnijk, _1stDeri[DCosPhi]*inverseNijk_Njkl_),
		Scale(rklCnjkl, -_1stDeri[DCosPhi]*cosPhi*inverseNjkl_Sq)
	);
	auto forceLK = Addition (
		Scale(njk, -_1stDeri[DCosThetaJKL]*inverseRkl_),
		Scale(nkl, _1stDeri[DRkl] -
			_1stDeri[DCosThetaJKL]*cosThetaJKL*inverseRkl_),
		Scale(rjkCnijk, -_1stDeri[DCosPhi]*inverseNijk_Njkl_),
		Scale(rjkCnjkl, _1stDeri[DCosPhi]*cosPhi*inverseNjkl_Sq)
	);
	return {{ move(forceIJ), move(forceKJ), move(forceLK) }};
}

//---------------------------------------------------------------------------//

std::array<matrix3d, 6>
PotentialDihedral::HessianImp (
	const array3d& rij,
	const array3d& rjk,
	const array3d& rkl,
	_1stD_Func _1stDeriFunc,
	_2ndD_Func _2ndDeriFunc
) const noexcept
{
	using namespace Dihedral_Namespace;

	auto rij_ = Norm(rij);
	auto inverseRij_ = 1./rij_;
	auto nij = UnitVector(rij, inverseRij_);
	auto rjk_ = Norm(rjk);
	auto inverseRjk_ = 1./rjk_;
	auto njk = UnitVector(rjk, inverseRjk_);
	auto rkl_ = Norm(rkl);
	auto inverseRkl_ = 1./rkl_;
	auto nkl = UnitVector(rkl, inverseRkl_);
	auto cosThetaIJK = CosTheta(nij, njk);
	auto cosThetaJKL = CosTheta(njk, nkl);
	auto nijk = Cross(rij, rjk);
	auto njkl = Cross(rjk, rkl);
	auto inverseNijk_ = 1./Norm(nijk);
	auto inverseNjkl_ = 1./Norm(njkl);
	auto cosPhi = CosTheta(UnitVector(nijk, inverseNijk_), UnitVector(njkl, inverseNjkl_));
	auto _1stDeri = _1stDeriFunc(cosPhi, rij_, rjk_, rkl_, cosThetaIJK, cosThetaJKL);
	auto _2ndDeri = _2ndDeriFunc(cosPhi, rij_, rjk_, rkl_, cosThetaIJK, cosThetaJKL);

/*	auto inverseNijk_Njkl_ = inverseNijk_*inverseNjkl_;
	auto inverseNijk_Sq = inverseNijk_*inverseNijk_;
	auto inverseNjkl_Sq = inverseNjkl_*inverseNjkl_;
	auto cosPhi_dri = Addition (
		Scale(Cross(rjk, njkl), -inverseNijk_Njkl_),
		Scale(Cross(rjk, nijk), cosPhi*inverseNijk_Sq)
	);
	auto cosPhi_drkj = Addition (
		Addition (
			Scale(Cross(rij, njkl), -inverseNijk_Njkl_),
			Scale(Cross(rkl, nijk), inverseNijk_Njkl_)
		),
		Scale(Cross(rij, nijk), cosPhi*inverseNijk_Sq),
		Scale(Cross(rkl, njkl), -cosPhi*inverseNjkl_Sq)
	);
	auto cosPhi_drj = Scale(Addition(cosPhi_dri, cosPhi_drkj), -1.);
	auto cosPhi_dri_drj = ;
	auto cosThetaIJK_dri = Addition (
		Scale(nij, cosThetaIJK*inverseRij_),
		Scale(njk, -inverseRij_)
	);
	auto cosThetaIJK_dri_drj = ;
	matrix3d hessianIJ = Addition (
		OuterDot( Addition (
			Scale(cosPhi_dri, _2ndDeri[DCosPhi_DCosPhi]),
			Scale(nij, -_2ndDeri[DCosPhi_DRij]),
			Scale(cosThetaIJK_dri, _2ndDeri[DCosPhi_DCosThetaIJK])
		), cosPhi_drj),
		OuterDot( Addition (
			Scale(cosPhi_dri, _2ndDeri[DCosPhi_DRij]),
			Scale(nij, _2ndDeri[DRij_DRij]),
			Scale(cosThetaIJK_dri, _2ndDeri[DRij_DCosThetaIJK])
		), nij),
		OuterDot( Addition (
			Scale(cosPhi_dri, _2ndDeri[DCosPhi_Drjk]),
			Scale(nij, _2ndDeri[DRij_DRjk]),
			Scale(cosThetaIJK_dri, _2ndDeri[DRjk_DCosThetaIJK])
		), Scale(njk, -1.)),
		OuterDot( Addition (
			Scale(cosPhi_dri, _2ndDeri[DCosPhi_DCosPhi]),
			Scale(nij, _2ndDeri[DRij_DCosThetaIJK]),
			Scale(cosThetaIJK_dri, _2ndDeri[DCosThetaIJK_DCosThetaIJK])
		), cosThetaIJK_drj),
		OuterDot( Addition (
			Scale(cosPhi_dri, _2ndDeri[DCosPhi_DCosThetaJKL]),
			Scale(nij, _2ndDeri[DRij_DCosThetaJKL]),
			Scale(cosThetaIJK_dri, _2ndDeri[DCosThetaIJK_DCosThetaJKL])
		), cosThetaJKL_drj),
		Scale(cosPhi_dri_drj, _1stDeri[DCosPhi]),
		Scale(Identity(), -_1stDeri[DRij]),
		Scale(cosThetaIJK_dri_drj, _1stDeri[DCosThetaIJK]),
	);*/
matrix3d hessianIJ{};
	matrix3d hessianIK{};
	matrix3d hessianIL{};
	matrix3d hessianJK{};
	matrix3d hessianJL{};
	matrix3d hessianKL{};
	return {{
		move(hessianIJ), move(hessianIK), move(hessianIL),
		move(hessianJK), move(hessianJL),
		move(hessianKL)
	}};
}

//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//

PotentialDihedral::_1stDerivative_t
PotentialDihedral::_1stDerivative (
	double cosPhi,
	double rij_,
	double rjk_,
	double rkl_,
	double cosThetaIJK,
	double cosThetaJKL
) const noexcept
{
	auto currentObjFun = ObjectiveFunction(cosPhi, rij_, rjk_, rkl_, cosThetaIJK, cosThetaJKL);
	auto incrementObjFun_CosPhi = ObjectiveFunction(cosPhi+FiniteDistance, rij_, rjk_, rkl_, cosThetaIJK, cosThetaJKL);
	auto incrementObjFun_Rij = ObjectiveFunction(cosPhi, rij_+FiniteDistance, rjk_, rkl_, cosThetaIJK, cosThetaJKL);
	auto incrementObjFun_Rjk = ObjectiveFunction(cosPhi, rij_, rjk_+FiniteDistance, rkl_, cosThetaIJK, cosThetaJKL);
	auto incrementObjFun_Rkl = ObjectiveFunction(cosPhi, rij_, rjk_, rkl_+FiniteDistance, cosThetaIJK, cosThetaJKL);
	auto incrementObjFun_CosThetaIJK = ObjectiveFunction(cosPhi, rij_, rjk_, rkl_, cosThetaIJK+FiniteDistance, cosThetaJKL);
	auto incrementObjFun_CosThetaJKL = ObjectiveFunction(cosPhi, rij_, rjk_, rkl_, cosThetaIJK, cosThetaJKL+FiniteDistance);
	return {
		(incrementObjFun_CosPhi-currentObjFun) * InverseFiniteDistance,
		(incrementObjFun_Rij-currentObjFun) * InverseFiniteDistance,
		(incrementObjFun_Rjk-currentObjFun) * InverseFiniteDistance,
		(incrementObjFun_Rkl-currentObjFun) * InverseFiniteDistance,
		(incrementObjFun_CosThetaIJK-currentObjFun) * InverseFiniteDistance,
		(incrementObjFun_CosThetaJKL-currentObjFun) * InverseFiniteDistance,
	};
}

PotentialDihedral::_2ndDerivative_t
PotentialDihedral::_2ndDerivative (
	double cosPhi,
	double rij_,
	double rjk_,
	double rkl_,
	double cosThetaIJK,
	double cosThetaJKL
) const noexcept
{
	auto current1stDeri = _1stDerivative(cosPhi, rij_, rjk_, rkl_, cosThetaIJK, cosThetaJKL);
	auto increment1stDeri_CosPhi = _1stDerivative(cosPhi+FiniteDistance, rij_, rjk_, rkl_, cosThetaIJK, cosThetaJKL);
	auto increment1stDeri_Rij = _1stDerivative(cosPhi, rij_+FiniteDistance, rjk_, rkl_, cosThetaIJK, cosThetaJKL);
	auto increment1stDeri_Rjk = _1stDerivative(cosPhi, rij_, rjk_+FiniteDistance, rkl_, cosThetaIJK, cosThetaJKL);
	auto increment1stDeri_Rkl = _1stDerivative(cosPhi, rij_, rjk_, rkl_+FiniteDistance, cosThetaIJK, cosThetaJKL);
	auto increment1stDeri_CosThetaIJK = _1stDerivative(cosPhi, rij_, rjk_, rkl_, cosThetaIJK+FiniteDistance, cosThetaJKL);
	auto increment1stDeri_CosThetaJKL = _1stDerivative(cosPhi, rij_, rjk_, rkl_, cosThetaIJK, cosThetaJKL+FiniteDistance);
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


double
PotentialDihedral::CosNPhi (
	unsigned n,
	double cosPhi
) const noexcept
{
	double cosSqPhi;
	switch( n )
	{
	case 0:
		return 1.;
	case 1:
		return cosPhi;
	case 2:
		return 2.*cosPhi*cosPhi-1.;
	case 3:
		return (4.*cosPhi*cosPhi-3.)*cosPhi;
	case 4:
		cosSqPhi = cosPhi * cosPhi;
		return 8.*cosSqPhi*(cosSqPhi-1) + 1.;
	case 5:
		cosSqPhi = cosPhi * cosPhi;
		return cosPhi*(cosSqPhi*(16.*cosSqPhi-20.)+5.);
	case 6:
		cosSqPhi = cosPhi * cosPhi;
		return cosSqPhi*(cosSqPhi*(32.*cosSqPhi-48.)+18.)-1.;
	}
	return cos( n*acos(cosPhi) );
}


double
PotentialDihedral::_1stDCosNPhi (
	unsigned n,
	double cosPhi
) const noexcept
{
	double cosPhiSq;
	switch( n )
	{
	case 0:
		return 0.;
	case 1:
		return 1.;
	case 2:
		return 4.*cosPhi;
	case 3:
		return 12.*cosPhi*cosPhi - 3.;
	case 4:
		return cosPhi*(32.*cosPhi*cosPhi-16.);
	case 5:
		cosPhiSq = cosPhi * cosPhi;
		return 80.*cosPhiSq*cosPhiSq - 60.*cosPhiSq + 5.;
	case 6:
		cosPhiSq = cosPhi * cosPhi;
		return cosPhi*(cosPhiSq*(192.*cosPhiSq-192.)+36.);
	}
	auto phi = acos(cosPhi);
	return n*sin(n*phi)/sin(phi);
}


double
PotentialDihedral::_2ndDCosNPhi (
	unsigned n,
	double cosPhi
) const noexcept
{
	double cosPhiSq;
	switch( n )
	{
	case 0:
	case 1:
		return 0.;
	case 2:
		return 4.;
	case 3:
		return 24.*cosPhi;
	case 4:
		return 96.*cosPhi*cosPhi-16.;
	case 5:
		cosPhiSq = cosPhi * cosPhi;
		return cosPhi*(320.*cosPhiSq - 120.);
	case 6:
		cosPhiSq = cosPhi * cosPhi;
		return cosPhiSq*(cosPhiSq*960.-576.)+36.;
	}
	auto phi = acos(cosPhi);
	auto sinPhi = sin(phi);
	auto sinNPhi = sin(n*phi);
	auto sinPhiSq = sinPhi*sinPhi;
	return -n*n*CosNPhi(n, cosPhi)/sinPhiSq + n*sinNPhi*cosPhi/(sinPhi*sinPhiSq);
}

