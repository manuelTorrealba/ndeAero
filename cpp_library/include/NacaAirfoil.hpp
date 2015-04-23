/**
 * File:   NacaAirfoil.hpp
 * Author: kenobi
 *
 * Created on Apr 17, 2015, 11:02 AM
 */

#ifndef INCLUDE_NACA_AIRFOIL_HPP
#define INCLUDE_NACA_AIRFOIL_HPP

#include "EquationSolver.hpp"
#include "SmartPtr.hpp"
#include "ThinAirfoil.hpp"
#include "Vector.hpp"

namespace nde {


	/**
	  * Naca Airfoil structure
	  */
	class NacaAirfoil {
	public:
      NacaAirfoil(double chord, const Vector<unsigned int>& naca_codenum);
		NacaAirfoil(const NacaAirfoil& naca_airfoil);

		virtual ~NacaAirfoil();
		
		double getChord() const {
			return _chord;
		}

      double top(double x) const {
      	return camber(x) + thickness(x);
      }

      double bottom(double x) const {
      	return camber(x) - thickness(x);
		}

	protected:
		double camber(double x) const;
		double thickness(double x) const;
		double dcamberdx(double x) const;

	private:
   	double _chord;
		Vector<unsigned int> _naca_codenum;
		ThinAirfoil* _naca_thin_airfoil; // this pointer is initialized in ctor.
											// and will be destroyed in dtor.

		unsigned int _num_digits;
		double _max_thickness;

		void helperCtor();

	};


	/**
	  * Naca airfoil coordinates for 4 and 5-digits
	  */
	class NacaFourDigits: public ThinAirfoil {
	public:
		NacaFourDigits(Vector<unsigned int> digits);

		virtual double camber(double t) const;
		virtual	double dCamberDx(double t) const;

	private:
		Vector<unsigned int> _digits;
		double _max_camber_val;
		double _max_camber_x; // maximum camber position

	};


	class NacaFiveDigits: public ThinAirfoil {
	public:
		NacaFiveDigits(Vector<unsigned int> digits);

		virtual double camber(double t) const;
		virtual	double dCamberDx(double t) const;

	private:
		Vector<unsigned int> _digits;
		double _max_camber_val;
		double _max_camber_x; // maximum camber position

		void helperCtor(double& max_camber_val, double& max_camber_x) const;

	};


	class NacaFiveDigitsReflx: public ThinAirfoil, public EquationSolver {
	public:
		NacaFiveDigitsReflx(Vector<unsigned int> digits);

		virtual double camber(double t) const;
		virtual	double dCamberDx(double t) const;

	protected:
		virtual double solverF(double x);

	private:
		Vector<unsigned int> _digits;
		double _max_camber_val;
		double _max_camber_x; // maximum camber position
		unsigned int _solver_objective; // = 1, _max_camber_x
												  // = 2, _max_camber_val

	};

	class NacaFiveDigitsMaxCamber: public EquationSolver {
	public:
			NacaFiveDigitsMaxCamber(double xf);
			double calcMaxCamberX();

	protected:
		virtual double solverF(double x);

	private:
		double _xf;
	};

} /* end of namespace nde */

#endif

