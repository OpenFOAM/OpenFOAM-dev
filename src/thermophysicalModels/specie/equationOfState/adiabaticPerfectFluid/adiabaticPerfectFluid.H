/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
    Foam::adiabaticPerfectFluid

Description
    AdiabaticPerfect gas equation of state.

SourceFiles
    adiabaticPerfectFluidI.H
    adiabaticPerfectFluid.C

\*---------------------------------------------------------------------------*/

#ifndef adiabaticPerfectFluid_H
#define adiabaticPerfectFluid_H

#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie> class adiabaticPerfectFluid;

template<class Specie>
inline adiabaticPerfectFluid<Specie> operator+
(
    const adiabaticPerfectFluid<Specie>&,
    const adiabaticPerfectFluid<Specie>&
);

template<class Specie>
inline adiabaticPerfectFluid<Specie> operator-
(
    const adiabaticPerfectFluid<Specie>&,
    const adiabaticPerfectFluid<Specie>&
);

template<class Specie>
inline adiabaticPerfectFluid<Specie> operator*
(
    const scalar,
    const adiabaticPerfectFluid<Specie>&
);

template<class Specie>
inline adiabaticPerfectFluid<Specie> operator==
(
    const adiabaticPerfectFluid<Specie>&,
    const adiabaticPerfectFluid<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const adiabaticPerfectFluid<Specie>&
);


/*---------------------------------------------------------------------------*\
                           Class adiabaticPerfectFluid Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class adiabaticPerfectFluid
:
    public Specie
{
    // Private data

        //- Reference pressure
        scalar p0_;

        //- Reference density
        scalar rho0_;

        //- The isentropic exponent
        scalar gamma_;

        //- Pressure offset for a stiffened gas
        scalar B_;


public:

    // Constructors

        //- Construct from components
        inline adiabaticPerfectFluid
        (
            const Specie& sp,
            const scalar p0,
            const scalar rho0,
            const scalar gamma,
            const scalar B
        );

        //- Construct from Istream
        adiabaticPerfectFluid(Istream&);

        //- Construct from dictionary
        adiabaticPerfectFluid(const dictionary& dict);

        //- Construct as named copy
        inline adiabaticPerfectFluid
        (
            const word& name,
            const adiabaticPerfectFluid&
        );

        //- Construct and return a clone
        inline autoPtr<adiabaticPerfectFluid> clone() const;

        // Selector from Istream
        inline static autoPtr<adiabaticPerfectFluid> New(Istream& is);

        // Selector from dictionary
        inline static autoPtr<adiabaticPerfectFluid> New
        (
            const dictionary& dict
        );


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "adiabaticPerfectFluid<" + word(Specie::typeName_()) + '>';
        }


        // Fundamental properties

            //- Is the equation of state is incompressible i.e. rho != f(p)
            static const bool incompressible = false;

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

        inline void operator+=(const adiabaticPerfectFluid&);
        inline void operator-=(const adiabaticPerfectFluid&);

        inline void operator*=(const scalar);


    // Friend operators

        friend adiabaticPerfectFluid operator+ <Specie>
        (
            const adiabaticPerfectFluid&,
            const adiabaticPerfectFluid&
        );

        friend adiabaticPerfectFluid operator- <Specie>
        (
            const adiabaticPerfectFluid&,
            const adiabaticPerfectFluid&
        );

        friend adiabaticPerfectFluid operator* <Specie>
        (
            const scalar s,
            const adiabaticPerfectFluid&
        );

        friend adiabaticPerfectFluid operator== <Specie>
        (
            const adiabaticPerfectFluid&,
            const adiabaticPerfectFluid&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie>
        (
            Ostream&,
            const adiabaticPerfectFluid&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "adiabaticPerfectFluidI.H"

#ifdef NoRepository
#   include "adiabaticPerfectFluid.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
