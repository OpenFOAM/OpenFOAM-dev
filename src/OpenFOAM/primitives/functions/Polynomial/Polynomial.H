/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::Polynomial

Description
    Polynomial templated on size (order):

    \verbatim
        poly = sum(coeffs[i]*x^i) + logCoeff*log(x)
    \endverbatim

    where <tt> 0 <= i <= N </tt>

    - integer powers, starting at zero
    - \c value(x) to evaluate the poly for a given value
    - \c derivative(x) returns derivative at value
    - \c integral(x1, x2) returns integral between two scalar values
    - \c integral() to return a new, integral coeff polynomial
      - increases the size (order)
    - \c integralMinus1() to return a new, integral coeff polynomial where
      the base poly starts at order -1

SourceFiles
    Polynomial.C

\*---------------------------------------------------------------------------*/

#ifndef Polynomial_H
#define Polynomial_H

#include "word.H"
#include "scalar.H"
#include "Ostream.H"
#include "VectorSpace.H"
#include <type_traits>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
template<int PolySize>
class Polynomial;

// Forward declaration of friend functions
template<int PolySize>
Ostream& operator<<
(
    Ostream&,
    const Polynomial<PolySize>&
);


/*---------------------------------------------------------------------------*\
                        Class Polynomial Declaration
\*---------------------------------------------------------------------------*/

template<int PolySize>
class Polynomial
:
    public VectorSpace<Polynomial<PolySize>, scalar, PolySize>
{
    static_assert(PolySize > 0, "Size must be positive (non-zero)");

    // Private data

        //- Include the log term? - only activated using integralMinus1()
        bool logActive_;

        //- Log coefficient - only activated using integralMinus1()
        scalar logCoeff_;


public:

    typedef Polynomial<PolySize> polyType;

    typedef Polynomial<PolySize+1> intPolyType;


    // Constructors

        //- Construct null, with all coefficients = 0.0
        Polynomial();

        //- Copy constructor
        Polynomial(const Polynomial&);

        //- Construct from C-array of coefficients
        explicit Polynomial(const scalar coeffs[PolySize]);

        //- Construct from a list of coefficients
        explicit Polynomial(const UList<scalar>& coeffs);

        //- Construct from Istream
        Polynomial(Istream&);

        //- Construct from name and Istream
        Polynomial(const word& name, Istream&);


    // Member Functions

        // Access

            //- Return true if the log term is active
            bool logActive() const;

            //- Return the log coefficient
            scalar logCoeff() const;


        // Evaluation

            //- Return polynomial value
            scalar value(const scalar x) const;

            //- Return derivative of the polynomial at the given x
            scalar derivative(const scalar x) const;

            //- Return integral between two values
            scalar integral(const scalar x1, const scalar x2) const;

            //- Return integral coefficients.
            //  Argument becomes zero'th element (constant of integration)
            intPolyType integral(const scalar intConstant = 0.0) const;

            //- Return integral coefficients when lowest order is -1.
            //  Argument becomes zero'th element (constant of integration)
            polyType integralMinus1(const scalar intConstant = 0.0) const;


    //- Ostream Operator
    friend Ostream& operator<< <PolySize>
    (
        Ostream&,
        const Polynomial&
    );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "Polynomial.C"
    #include "PolynomialIO.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
