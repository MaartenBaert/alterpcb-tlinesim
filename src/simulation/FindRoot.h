/*
Copyright (C) 2016  The AlterPCB team
Contact: Maarten Baert <maarten-baert@hotmail.com>

This file is part of AlterPCB.

AlterPCB is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

AlterPCB is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this AlterPCB.  If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once

#include "Basics.h"
#include "MiscMath.h"

// TODO: remove
//#include <iostream>
//#include <iomanip>

// Find the root of func(x) for x in the interval [a, b] using Brent's method. This method uses a combination of
// the bisection method, secant method and inverse quadratic interpolation method .
// The sign of func(a) and func(b) must be opposite.
template<typename Func>
real_t FindRootBracketed(Func &&func, real_t a, real_t b, real_t fa, real_t fb, real_t xtol, real_t ftol) {
	if(!std::isfinite(a))
		throw std::runtime_error("Input a is non-finite.");
	if(!std::isfinite(b))
		throw std::runtime_error("Input b is non-finite.");
	if(!std::isfinite(fa))
		throw std::runtime_error("Input fa is non-finite.");
	if(!std::isfinite(fb))
		throw std::runtime_error("Input fb is non-finite.");
	if(!FinitePositive(xtol))
		throw std::runtime_error("Input xtol must be positive.");
	if(!FinitePositive(ftol))
		throw std::runtime_error("Input ftol must be positive.");
	if((fa > 0.0) == (fb > 0.0))
		throw std::runtime_error("Inputs fa and fb should have opposite signs.");
	if(fabs(fa) < fabs(fb)) {
		std::swap(a, b);
		std::swap(fa, fb);
	}
	//std::cerr.precision(12);
	real_t c = a, fc = fa, d = 0.0, delta = 4.0 * std::numeric_limits<real_t>::epsilon() * fmax(fabs(a), fabs(b));
	bool mflag = true;
	while(fabs(a - b) > xtol && fabs(fb) > ftol) {
		real_t ab = a * 0.75 + b * 0.25; // overflow-safe
		real_t bcd = (mflag)? fabs(b - c) : fabs(c - d);
		real_t s;
		if(fa != fc && fb != fc) {
			// use inverse quadratic interpolation method
			s = (a * fb * fc * (fb - fc) - b * fa * fc * (fa - fc) + c * fa * fb * (fa - fb)) / ((fa - fb) * (fb - fc) * (fa - fc));
		} else {
			// use secant method
			s = b - fb * (b - a) / (fb - fa);
		}
		//std::cerr << "mflag " << (s <= fmin(ab, b)) << (s >= fmax(ab, b)) << (fabs(s - b) >= bcd * 0.5) << (bcd < delta) << "   ";
		mflag = (!std::isfinite(s) || s <= fmin(ab, b) || s >= fmax(ab, b) || fabs(s - b) >= bcd * 0.5 || bcd < delta);
		if(mflag) {
			// use bisection method
			s = a * 0.5 + b * 0.5; // overflow-safe
			//std::cerr << "bisection   " << std::setw(15) << a << " " << std::setw(15) << b << " " << std::setw(15) << c << " " << std::setw(15) << d << " " << std::setw(15) << s << std::endl;
		} else {
			if(fa != fc && fb != fc) {
				//std::cerr << "quadratic   " << std::setw(15) << a << " " << std::setw(15) << b << " " << std::setw(15) << c << " " << std::setw(15) << d << " " << std::setw(15) << s << std::endl;
			} else {
				//std::cerr << "secant      " << std::setw(15) << a << " " << std::setw(15) << b << " " << std::setw(15) << c << " " << std::setw(15) << d << " " << std::setw(15) << s << std::endl;
			}
		}
		d = c;
		c = b;
		fc = fb;
		real_t fs = func(s);
		if(!std::isfinite(fs)) {
			throw std::runtime_error("Function result is non-finite.");
		}
		if((fa > 0.0) != (fs > 0.0)) {
			b = s;
			fb = fs;
		} else {
			a = s;
			fa = fs;
		}
		if(fabs(fa) < fabs(fb)) {
			std::swap(a, b);
			std::swap(fa, fb);
		}
	}
	return b;
}

template<typename Func>
real_t FindRootBracketed(Func &&func, real_t a, real_t b, real_t xtol, real_t ftol) {
	return FindRootBracketed(std::forward<Func>(func), a, b, func(a), func(b), xtol, ftol);
}

template<typename Func>
real_t FindRootRelative(Func &&func, real_t x, real_t xreltol, real_t ftol, real_t search_range) {
	if(!FinitePositive(x))
		throw std::runtime_error("Input x should be positive.");
	if(!FinitePositive(xreltol))
		throw std::runtime_error("Input xreltol should be positive.");
	if(!FinitePositive(ftol))
		throw std::runtime_error("Input ftol should be positive.");
	if(!FinitePositive(search_range))
		throw std::runtime_error("Input search_range should be positive.");
	real_t fx = func(x);
	if(!std::isfinite(fx)) {
		throw std::runtime_error("Function result is non-finite.");
	}
	if(fabs(fx) <= ftol)
		return x;
	real_t a = x, fa = fx, b = x, fb = fx;
	//std::cerr.precision(12);
	for( ; ; ) {

		// try lower value
		x = a * 0.5;
		fx = func(x);
		if(!std::isfinite(fx)) {
			throw std::runtime_error("Function result is non-finite.");
		}
		//std::cerr << "lower       " << std::setw(15) << a << " " << std::setw(15) << x << std::endl;
		if((fx > 0.0) != (fa > 0.0)) {
			return FindRootBracketed(std::forward<Func>(func), a, x, fa, fx, fmax(fabs(a), fabs(x)) * xreltol, ftol);
		}
		a = x;
		fa = fx;

		// try higher value
		x = b * 2.0;
		fx = func(x);
		if(!std::isfinite(fx)) {
			throw std::runtime_error("Function result is non-finite.");
		}
		//std::cerr << "higher      " << std::setw(15) << b << " " << std::setw(15) << x << std::endl;
		if((fx > 0.0) != (fb > 0.0)) {
			return FindRootBracketed(std::forward<Func>(func), b, x, fb, fx, fmax(fabs(b), fabs(x)) * xreltol, ftol);
		}
		b = x;
		fb = fx;

		// check search range
		if(b > a * search_range) {
			throw std::runtime_error("No sign change found in search range.");
		}

	}
}
