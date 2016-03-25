#ifndef POTENTIAL_H_INCLUDED
#define POTENTIAL_H_INCLUDED

/* definition:
	rij = rj - ri
	fij = force at i due to j
	Hij = potential second derivative about ri and rj

	pair: i-j
	angle: j-i-k
	dihedral: i-j-k-l
*/

#include "array3d.h"
//#include "PotentialType.h"
#include <functional>

struct FiniteDifference_t {};

namespace Angle_Namespace {
	enum {DCosTheta, DRij, DRik};
	enum {DCosTheta_DCosTheta, DCosTheta_DRij, DCosTheta_DRik, DRij_DRij, DRij_DRik, DRik_DRik};
	using _1stDerivative_t = std::array<double, 3>;
	using _2ndDerivative_t = std::array<double, 6>;
};
namespace Dihedral_Namespace {
	enum {DCosPhi, DRij, DRjk, DRkl, DCosThetaIJK, DCosThetaJKL};
	enum {	DCosPhi_DCosPhi, DCosPhi_DRij, DCosPhi_DRjk, DCosPhi_DRkl, DCosPhi_DCosThetaIJK, DCosPhi_DCosThetaJKL,
		DRij_DRij, DRij_DRjk, DRij_DRkl, DRij_DCosThetaIJK, DRij_DCosThetaJKL,
		DRjk_DRjk, DRjk_DRkl, DRjk_DCosThetaIJK, DRjk_DCosThetaJKL,
		DRkl_DRkl, DRkl_DCosThetaIJK, DRkl_DCosThetaJKL,
		DCosThetaIJK_DCosThetaIJK, DCosThetaIJK_DCosThetaJKL,
		DCosThetaJKL_DCosThetaJKL
	};
	using _1stDerivative_t = std::array<double, 6>;
	using _2ndDerivative_t = std::array<double, 21>;
};


class Potential {
public:
	using Pair_1stDerivative_t = double;
	using Pair_2ndDerivative_t = double;
	using Angle_1stDerivative_t = Angle_Namespace::_1stDerivative_t;
	using Angle_2ndDerivative_t = Angle_Namespace::_2ndDerivative_t;
	using Dihedral_1stDerivative_t = Dihedral_Namespace::_1stDerivative_t;
	using Dihedral_2ndDerivative_t = Dihedral_Namespace::_2ndDerivative_t;


	double
	Energy (
		const array3d& rij
	) const noexcept;


	double
	Energy (
		const array3d& rij,
		const array3d& rik
	) const noexcept;


	double
	Energy (
		const array3d& rij,
		const array3d& rjk,
		const array3d& rkl
	) const noexcept;

	//-------------------------------------//

	array3d // fij
	Force (
		const array3d& rij
	) const noexcept;


	array3d
	Force (
		const array3d& rij,
		FiniteDifference_t
	) const noexcept;


	std::array<array3d, 2> // fij, fik
	Force (
		const array3d& rij,
		const array3d& rik
	) const noexcept;


	std::array<array3d, 2>
	Force (
		const array3d& rij,
		const array3d& rik,
		FiniteDifference_t
	) const noexcept;


	std::array<array3d, 3> // fi, fjk, fl
	Force (
		const array3d& rij,
		const array3d& rjk,
		const array3d& rkl
	) const noexcept;


	std::array<array3d, 3>
	Force (
		const array3d& rij,
		const array3d& rjk,
		const array3d& rkl,
		FiniteDifference_t
	) const noexcept;

	//-------------------------------------//

	matrix3d // Hij
	Hessian (
		const array3d& rij
	) const noexcept;


	matrix3d
	Hessian (
		const array3d& rij,
		FiniteDifference_t
	) const noexcept;


	std::array<matrix3d, 3> // Hij, Hik, Hjk
	Hessian (
		const array3d& rij,
		const array3d& rik
	) const noexcept;


	std::array<matrix3d, 3>
	Hessian (
		const array3d& rij,
		const array3d& rik,
		FiniteDifference_t
	) const noexcept;


	std::array<matrix3d, 6> // Hij, Hik, Hil, Hjk, Hjl, Hkl
	Hessian (
		const array3d& rij,
		const array3d& rjk,
		const array3d& rkl
	) const noexcept;


	std::array<matrix3d, 6>
	Hessian (
		const array3d& rij,
		const array3d& rjk,
		const array3d& rkl,
		FiniteDifference_t
	) const noexcept;

protected:
	double
	Norm (
		const array3d& rij
	) const noexcept;


	array3d
	UnitVector (
		const array3d& rij
	) const noexcept;


	array3d
	UnitVector (
		const array3d& rij,
		double inverseRij
	) const noexcept;


	double
	CosTheta (
		const array3d& nij,
		const array3d& nik
	) const noexcept;


	double
	CosPhi (
		const array3d& nij,
		const array3d& njk,
		const array3d& nkl
	) const noexcept;

private:
	array3d
	ForceImp (
		const array3d& rij,
		std::function< Pair_1stDerivative_t(double) > _1stDeriFunc
	) const noexcept;


	std::array<array3d, 2>
	ForceImp (
		const array3d& rij,
		const array3d& rik,
		std::function< Angle_1stDerivative_t(double, double, double) > _1stDeriFunc
	) const noexcept;


	std::array<array3d, 3>
	ForceImp (
		const array3d& rij,
		const array3d& rjk,
		const array3d& rkl,
		std::function< Dihedral_1stDerivative_t(double, double, double, double, double, double) > _1stDeriFunc
	) const noexcept;


	matrix3d
	HessianImp (
		const array3d& rij,
		std::function< Pair_1stDerivative_t(double) > _1stDeriFunc,
		std::function< Pair_2ndDerivative_t(double) > _2ndDeriFunc
	) const noexcept;


	std::array<matrix3d, 3>
	HessianImp (
		const array3d& rij,
		const array3d& rik,
		std::function< Angle_1stDerivative_t(double, double, double) > _1stDeriFunc,
		std::function< Angle_2ndDerivative_t(double, double, double) > _2ndDeriFunc
	) const noexcept;


	std::array<matrix3d, 6>
	HessianImp (
		const array3d& rij,
		const array3d& rjk,
		const array3d& rkl,
		std::function< Dihedral_1stDerivative_t(double, double, double, double, double, double) > _1stDeriFunc,
		std::function< Dihedral_2ndDerivative_t(double, double, double, double, double, double) > _2ndDeriFunc
	) const noexcept;

	//-------------------------------------//
	//-------------------------------------//
	//-------------------------------------//

	virtual double
	ObjectiveFunction (
		double /*rij*/
	) const noexcept
	{	return 0.; }


	virtual double
	ObjectiveFunction (
		double /*cosTheta*/,
		double /*rij*/,
		double /*rik*/
	) const noexcept
	{	return 0.; }


	virtual double
	ObjectiveFunction (
		double /*cosPhi*/,
		double /*rij*/,
		double /*rjk*/,
		double /*rkl*/,
		double /*cosThetaIJK*/,
		double /*cosThetaJKL*/
	) const noexcept
	{	return 0.; }


	virtual Pair_1stDerivative_t // drij
	_1stDerivative (
		double rij
	) const noexcept;


	virtual Angle_1stDerivative_t // dcosTheta, drij, drik
	_1stDerivative (
		double cosTheta,
		double rij,
		double rik
	) const noexcept;


	virtual Dihedral_1stDerivative_t // dcosPhi, drij, drik, drkl, dcosThetaIJK, dcosThetaJKL
	_1stDerivative (
		double cosPhi,
		double rij,
		double rjk,
		double rkl,
		double cosThetaIJK,
		double cosThetaJKL
	) const noexcept;


	virtual Pair_2ndDerivative_t
	_2ndDerivative (
		double rij
	) const noexcept;


	virtual Angle_2ndDerivative_t	// dcosTheta_dcosTheta, dcosTheta_drij, dcosTheta_drik
	_2ndDerivative (		// drij_drij, drij_drik
		double cosTheta,	// drik_drik
		double rij,
		double rik
	) const noexcept;


	virtual Dihedral_2ndDerivative_t
	_2ndDerivative (
		double cosPhi,
		double rij,
		double rjk,
		double rkl,
		double cosThetaIJK,
		double cosThetaJKL
	) const noexcept;
};

#endif // POTENTIAL_H_INCLUDED
