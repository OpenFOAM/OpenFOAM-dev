/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2020 OpenFOAM Foundation
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
    Foam::rPolynomial

Description
    Reciprocal polynomial equation of state for liquids and solids

    \f[
        1/\rho = C_0 + C_1 T + C_2 T^2 - C_3 p - C_4 p T
    \f]

    This polynomial for the reciprocal of the density provides a much better fit
    than the equivalent polynomial for the density and has the advantage that it
    support coefficient mixing to support liquid and solid mixtures in an
    efficient manner.

Usage
    \table
        Property     | Description
        C            | Density polynomial coefficients
    \endtable

    Example of the specification of the equation of state for pure water:
    \verbatim
    equationOfState
    {
        C (0.001278 -2.1055e-06 3.9689e-09 4.3772e-13 -2.0225e-16);
    }
    \endverbatim
    Note: This fit is based on the small amount of data which is freely
    available for the range 20-65degC and 1-100bar.

SourceFiles
    rPolynomialI.H
    rPolynomial.C

\*---------------------------------------------------------------------------*/

#ifndef rPolynomial_H
#define rPolynomial_H

#include "autoPtr.H"
#include "VectorSpace.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie> class rPolynomial;

template<class Specie>
inline rPolynomial<Specie> operator+
(
    const rPolynomial<Specie>&,
    const rPolynomial<Specie>&
);

template<class Specie>
inline rPolynomial<Specie> operator*
(
    const scalar,
    const rPolynomial<Specie>&
);

template<class Specie>
inline rPolynomial<Specie> operator==
(
    const rPolynomial<Specie>&,
    const rPolynomial<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const rPolynomial<Specie>&
);


/*---------------------------------------------------------------------------*\
                         Class rPolynomial Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class rPolynomial
:
    public Specie
{
    // Private Data

        class coeffList
        :
            public VectorSpace<coeffList, scalar, 5>
        {
            public:

            // Constructors

                //- Construct null
                inline coeffList()
                {}

                //- Construct from Istream
                inline coeffList(Istream& is)
                :
                    VectorSpace<coeffList, scalar, 5>(is)
                {}
        };


        //- Density coefficients
        coeffList C_;


public:

    // Constructors

        //- Construct from components
        inline rPolynomial
        (
            const Specie& sp,
            const coeffList& coeffs
        );

        //- Construct from dictionary
        rPolynomial(const dictionary& dict);

        //- Construct as named copy
        inline rPolynomial(const word& name, const rPolynomial&);

        //- Construct and return a clone
        inline autoPtr<rPolynomial> clone() const;

        // Selector from dictionary
        inline static autoPtr<rPolynomial> New(const dictionary& dict);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "rPolynomial<" + word(Specie::typeName_()) + '>';
        }


        // Fundamental properties

            //- Is the equation of state is incompressible i.e. rho != f(p)
            static const bool incompressible = false;

            //- Is the equation of state is isochoric i.e. rho = const
            static const bool isochoric = false;

            //- Return density [kg/m^3]
            inline scalar rho(scalar p, scalar T) const;

            //- Return enthalpy contribution [J/kg]
            inline scalar H(const scalar p, const scalar T) const;

            //- Return Cp contribution [J/(kg K]
            inline scalar Cp(scalar p, scalar T) const;

            //- Return internal energy contribution [J/kg]
            inline scalar E(const scalar p, const scalar T) const;

            //- Return Cv contribution [J/(kg K]
            inline scalar Cv(scalar p, scalar T) const;

            //- Return entropy contribution to the integral of Cp/T [J/kg/K]
            inline scalar Sp(const scalar p, const scalar T) const;

            //- Return entropy contribution to the integral of Cv/T [J/kg/K]
            inline scalar Sv(const scalar p, const scalar T) const;

            //- Return compressibility [s^2/m^2]
            inline scalar psi(scalar p, scalar T) const;

            //- Return compression factor []
            inline scalar Z(scalar p, scalar T) const;

            //- Return (Cp - Cv) [J/(kg K]
            inline scalar CpMCv(scalar p, scalar T) const;


        // IO

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member Operators

        inline void operator+=(const rPolynomial&);
        inline void operator*=(const scalar);


    // Friend operators

        friend rPolynomial operator+ <Specie>
        (
            const rPolynomial&,
            const rPolynomial&
        );

        friend rPolynomial operator* <Specie>
        (
            const scalar s,
            const rPolynomial&
        );

        friend rPolynomial operator== <Specie>
        (
            const rPolynomial&,
            const rPolynomial&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie>
        (
            Ostream&,
            const rPolynomial&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "rPolynomialI.H"

#ifdef NoRepository
    #include "rPolynomial.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
