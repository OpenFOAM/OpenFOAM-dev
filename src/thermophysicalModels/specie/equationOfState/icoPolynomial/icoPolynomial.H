/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::icoPolynomial

Description
    Incompressible, polynomial form of equation of state, using a polynomial
    function for density.

SourceFiles
    icoPolynomialI.H
    icoPolynomial.C

\*---------------------------------------------------------------------------*/

#ifndef icoPolynomial_H
#define icoPolynomial_H

#include "autoPtr.H"
#include "Polynomial.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie, int PolySize>
class icoPolynomial;

template<class Specie, int PolySize>
icoPolynomial<Specie, PolySize> operator+
(
    const icoPolynomial<Specie, PolySize>&,
    const icoPolynomial<Specie, PolySize>&
);

template<class Specie, int PolySize>
icoPolynomial<Specie, PolySize> operator-
(
    const icoPolynomial<Specie, PolySize>&,
    const icoPolynomial<Specie, PolySize>&
);

template<class Specie, int PolySize>
icoPolynomial<Specie, PolySize> operator*
(
    const scalar,
    const icoPolynomial<Specie, PolySize>&
);

template<class Specie, int PolySize>
icoPolynomial<Specie, PolySize> operator==
(
    const icoPolynomial<Specie, PolySize>&,
    const icoPolynomial<Specie, PolySize>&
);

template<class Specie, int PolySize>
Ostream& operator<<
(
    Ostream&,
    const icoPolynomial<Specie, PolySize>&
);


/*---------------------------------------------------------------------------*\
                        Class icoPolynomial Declaration
\*---------------------------------------------------------------------------*/

template<class Specie, int PolySize=8>
class icoPolynomial
:
    public Specie
{
    // Private data

        //- Density polynomial coefficients
        //  Note: input in [kg/m3], but internally uses [kg/m3/kmol]
        Polynomial<PolySize> rhoCoeffs_;


public:

    // Constructors

        //- Construct from components
        inline icoPolynomial
        (
            const Specie& sp,
            const Polynomial<PolySize>& rhoPoly
        );

        //- Construct from Istream
        icoPolynomial(Istream&);

        //- Construct from dictionary
        icoPolynomial(const dictionary& dict);

        //- Construct as copy
        inline icoPolynomial(const icoPolynomial&);

        //- Construct as named copy
        inline icoPolynomial(const word& name, const icoPolynomial&);

        //- Construct and return a clone
        inline autoPtr<icoPolynomial> clone() const;

        // Selector from Istream
        inline static autoPtr<icoPolynomial> New(Istream& is);

        // Selector from dictionary
        inline static autoPtr<icoPolynomial> New(const dictionary& dict);


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "icoPolynomial<" + word(Specie::typeName_()) + '>';
        }


        // Fundamental properties

            //- Is the equation of state is incompressible i.e. rho != f(p)
            static const bool incompressible = true;

            //- Is the equation of state is isochoric i.e. rho = const
            static const bool isochoric = false;

            //- Return density [kg/m^3]
            inline scalar rho(scalar p, scalar T) const;

            //- Return compressibility rho/p [s^2/m^2]
            inline scalar psi(scalar p, scalar T) const;

            //- Return compression factor []
            inline scalar Z(scalar p, scalar T) const;

            //- Return (cp - cv) [J/(kmol K]
            inline scalar cpMcv(scalar p, scalar T) const;


        // IO

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member operators

        inline icoPolynomial& operator=(const icoPolynomial&);
        inline void operator+=(const icoPolynomial&);
        inline void operator-=(const icoPolynomial&);

        inline void operator*=(const scalar);


    // Friend operators

        friend icoPolynomial operator+ <Specie, PolySize>
        (
            const icoPolynomial&,
            const icoPolynomial&
        );

        friend icoPolynomial operator- <Specie, PolySize>
        (
            const icoPolynomial&,
            const icoPolynomial&
        );

        friend icoPolynomial operator* <Specie, PolySize>
        (
            const scalar s,
            const icoPolynomial&
        );

        friend icoPolynomial operator== <Specie, PolySize>
        (
            const icoPolynomial&,
            const icoPolynomial&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie, PolySize>
        (
            Ostream&,
            const icoPolynomial&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeIcoPolynomial(PolySize)                                          \
                                                                             \
defineTemplateTypeNameAndDebugWithName                                       \
(                                                                            \
    icoPolynomial<Specie, PolySize>,                                         \
    "icoPolynomial<"#PolySize">",                                            \
    0                                                                        \
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "icoPolynomialI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "icoPolynomial.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
