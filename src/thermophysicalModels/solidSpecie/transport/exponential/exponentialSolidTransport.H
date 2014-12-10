/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::exponentialSolidTransport

Description
    Exponential properties for solid heat transport
    Templated into a given thermodynamics package.

SourceFiles
    exponentialSolidTransportI.H
    exponentialSolidTransport.C

\*---------------------------------------------------------------------------*/

#ifndef exponentialSolidTransport_H
#define exponentialSolidTransport_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Thermo> class exponentialSolidTransport;

template<class Thermo>
inline exponentialSolidTransport<Thermo> operator*
(
    const scalar,
    const exponentialSolidTransport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const exponentialSolidTransport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                  Class exponentialSolidTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class exponentialSolidTransport
:
    public Thermo
{
    // Private data

        //- Constant thermal coefficient.
        scalar kappa0_;

        //- Exponent coefficient
        scalar n0_;

        //- Reference temperature
        scalar Tref_;


    // Private Member Functions

        //- Construct from components
        inline exponentialSolidTransport
        (
            const Thermo& t,
            const scalar kappa0,
            const scalar n0,
            const scalar Tref
        );


public:

    // Constructors

        //- Construct as named copy
        inline exponentialSolidTransport
        (
            const word&,
            const exponentialSolidTransport&
        );

        //- Construct from dictionary
        exponentialSolidTransport(const dictionary&);

        // Selector from dictionary
        inline static autoPtr<exponentialSolidTransport> New
        (
            const dictionary& dict
        );


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "exponential<" + Thermo::typeName() + '>';
        }

        //- Is the thermal conductivity isotropic
        static const bool isotropic = true;

        //- Thermal conductivity [W/mK]
        inline scalar kappa(const scalar p, const scalar T) const;

        //- Thermal conductivity [W/mK]
        inline vector Kappa(const scalar p, const scalar T) const;

         //- Dynamic viscosity [kg/ms]
        inline scalar mu(const scalar p, const scalar T) const;

        //- Thermal diffusivity of enthalpy [kg/ms]
        inline scalar alphah(const scalar p, const scalar T) const;

        //- Write to Ostream
        void write(Ostream& os) const;


    // Member operators

        inline exponentialSolidTransport& operator=
        (
            const exponentialSolidTransport&
        );
        inline void operator+=(const exponentialSolidTransport&);
        inline void operator-=(const exponentialSolidTransport&);


    // Friend operators

        friend exponentialSolidTransport operator* <Thermo>
        (
            const scalar,
            const exponentialSolidTransport&
        );

    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const exponentialSolidTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "exponentialSolidTransportI.H"

#ifdef NoRepository
#   include "exponentialSolidTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
