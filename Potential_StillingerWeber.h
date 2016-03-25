#include "Potential.h"

class Potential_StillingerWeber : public Potential {
public:
	Potential_StillingerWeber (
		double epsilon = 2.1683,
		double sigma = 2.0951,
		double a = 1.8,
		double lambda = 21.,
		double gamma = 1.2,
		double cosTheta0 = -1./3.,
		double A = 7.049556277,
		double B = 0.6022245584,
		double p = 4.,
		double q = 0.
	) noexcept;

protected:
	double
	ObjectiveFunction (
		double rij
	) const noexcept override;


	double
	ObjectiveFunction (
		double cosTheta,
		double rij,
		double rik
	) const noexcept override;

	//-------------------------------------//

	double
	_1stDerivative (
		double rij
	) const noexcept override;


	Angle_Namespace::_1stDerivative_t
	_1stDerivative (
		double cosTheta,
		double rij,
		double rik
	) const noexcept override;

	//-------------------------------------//

	double
	_2ndDerivative (
		double rij
	) const noexcept override;


	Angle_Namespace::_2ndDerivative_t
	_2ndDerivative (
		double cosTheta,
		double rij,
		double rik
	) const noexcept override;

	//-------------------------------------//

	double
	FactorPair (
		double rij,
		double minus = 0.
	) const noexcept;

private:
	double m_epsilon;
	double m_sigma;
	double m_a;
	double m_lambda;
	double m_gamma;
	double m_cosTheta0;
	double m_A;
	double m_B;
	double m_p;
	double m_q;
};

