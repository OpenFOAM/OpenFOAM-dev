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
    Foam::polynomialTransport

Description
    Transport package using polynomial functions for mu and kappa

SourceFiles
    polynomialTransportI.H
    polynomialTransport.C

\*---------------------------------------------------------------------------*/

#ifndef polynomialTransport_H
#define polynomialTransport_H

#include "Polynomial.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo, int PolySize> class polynomialTransport;

template<class Thermo, int PolySize>
inline polynomialTransport<Thermo, PolySize> operator+
(
    const polynomialTransport<Thermo, PolySize>&,
    const polynomialTransport<Thermo, PolySize>&
);

template<class Thermo, int PolySize>
inline polynomialTransport<Thermo, PolySize> operator-
(
    const polynomialTransport<Thermo, PolySize>&,
    const polynomialTransport<Thermo, PolySize>&
);

template<class Thermo, int PolySize>
inline polynomialTransport<Thermo, PolySize> operator*
(
    const scalar,
    const polynomialTransport<Thermo, PolySize>&
);

template<class Thermo, int PolySize>
inline polynomialTransport<Thermo, PolySize> operator==
(
    const polynomialTransport<Thermo, PolySize>&,
    const polynomialTransport<Thermo, PolySize>&
);

template<class Thermo, int PolySize>
Ostream& operator<<
(
    Ostream&,
    const polynomialTransport<Thermo, PolySize>&
);


/*---------------------------------------------------------------------------*\
                     Class polynomialTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo, int PolySize=8>
class polynomialTransport
:
    public Thermo
{
    // Private data

        //- Dynamic viscosity polynomial coefficients
        //  Note: input in [Pa.s], but internally uses [Pa.s/kmol]
        Polynomial<PolySize> muCoeffs_;

        //- Thermal conductivity polynomial coefficients
        //  Note: input in [W/m/K], but internally uses [W/m/K/kmol]
        Polynomial<PolySize> kappaCoeffs_;


    // Private Member Functions

        //- Construct from components
        inline polynomialTransport
        (
            const Thermo& t,
            const Polynomial<PolySize>& muPoly,
            const Polynomial<PolySize>& kappaPoly
        );


public:

    // Constructors

        //- Construct copy
        inline polynomialTransport(const polynomialTransport&);

        //- Construct as named copy
        inline polynomialTransport(const word&, const polynomialTransport&);

        //- Construct from Istream
        polynomialTransport(Istream& is);

        //- Construct from dictionary
        polynomialTransport(const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<polynomialTransport> clone() const;

        // Selector from Istream
        inline static autoPtr<polynomialTransport> New(Istream& is);

        // Selector from dictionary
        inline static autoPtr<polynomialTransport> New(const dictionary& dict);


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "polynomial<" + Thermo::typeName() + '>';
        }

        //- Dynamic viscosity [kg/ms]
        inline scalar mu(const scalar p, const scalar T) const;

        //- Thermal conductivity [W/mK]
        inline scalar kappa(const scalar p, const scalar T) const;

        //- Thermal diffusivity of enthalpy [kg/ms]
        inline scalar alphah(const scalar p, const scalar T) const;

        // Species diffusivity
        //inline scalar D(const scalar p, const scalar T) const;

        //- Write to Ostream
        void write(Ostream& os) const;


    // Member operators

        inline polynomialTransport& operator=(const polynomialTransport&);
        inline void operator+=(const polynomialTransport&);
        inline void operator-=(const polynomialTransport&);
        inline void operator*=(const scalar);


    // Friend operators

        friend polynomialTransport operator+ <Thermo, PolySize>
        (
            const polynomialTransport&,
            const polynomialTransport&
        );

        friend polynomialTransport operator- <Thermo, PolySize>
        (
            const polynomialTransport&,
            const polynomialTransport&
        );

        friend polynomialTransport operator* <Thermo, PolySize>
        (
            const scalar,
            const polynomialTransport&
        );

        friend polynomialTransport operator== <Thermo, PolySize>
        (
            const polynomialTransport&,
            const polynomialTransport&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo, PolySize>
        (
            Ostream&,
            const polynomialTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "polynomialTransportI.H"

#ifdef NoRepository
#   include "polynomialTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
