/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::constTransport

Description
    Constant properties Transport package.
    Templated into a given thermodynamics package (needed for thermal
    conductivity).

SourceFiles
    constTransportI.H
    constTransport.C

\*---------------------------------------------------------------------------*/

#ifndef constTransport_H
#define constTransport_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo> class constTransport;

template<class Thermo>
inline constTransport<Thermo> operator+
(
    const constTransport<Thermo>&,
    const constTransport<Thermo>&
);

template<class Thermo>
inline constTransport<Thermo> operator*
(
    const scalar,
    const constTransport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const constTransport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                       Class constTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class constTransport
:
    public Thermo
{
    // Private Data

        //- Constant dynamic viscosity [Pa.s]
        scalar mu_;

        //- Reciprocal Prandtl Number []
        scalar rPr_;


    // Private Member Functions

        //- Construct from components
        inline constTransport
        (
            const Thermo& t,
            const scalar mu,
            const scalar Pr
        );


public:

    // Constructors

        //- Construct as named copy
        inline constTransport(const word&, const constTransport&);

        //- Construct from dictionary
        constTransport(const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<constTransport> clone() const;

        // Selector from dictionary
        inline static autoPtr<constTransport> New(const dictionary& dict);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "const<" + Thermo::typeName() + '>';
        }

        //- Dynamic viscosity [kg/m/s]
        inline scalar mu(const scalar p, const scalar T) const;

        //- Thermal conductivity [W/m/K]
        inline scalar kappa(const scalar p, const scalar T) const;

        // Species diffusivity
        // inline scalar D(const scalar p, const scalar T) const;

        //- Write to Ostream
        void write(Ostream& os) const;


    // Member Operators

        inline void operator+=(const constTransport&);
        inline void operator*=(const scalar);


    // Friend operators

        friend constTransport operator+ <Thermo>
        (
            const constTransport&,
            const constTransport&
        );

        friend constTransport operator* <Thermo>
        (
            const scalar,
            const constTransport&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const constTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "constTransportI.H"

#ifdef NoRepository
    #include "constTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
