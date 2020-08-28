/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2020 OpenFOAM Foundation
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
    Foam::Boussinesq

Description
    Incompressible gas equation of state using the Boussinesq approximation for
    the density as a function of temperature only:

    \verbatim
        rho = rho0*(1 - beta*(T - T0))
    \endverbatim

    Coefficient mixing is very inaccurate and not supported,
    so this equation of state is not applicable to mixtures.

SourceFiles
    BoussinesqI.H
    Boussinesq.C

\*---------------------------------------------------------------------------*/

#ifndef Boussinesq_H
#define Boussinesq_H

#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie> class Boussinesq;

template<class Specie>
inline Boussinesq<Specie> operator+
(
    const Boussinesq<Specie>&,
    const Boussinesq<Specie>&
);

template<class Specie>
inline Boussinesq<Specie> operator*
(
    const scalar,
    const Boussinesq<Specie>&
);

template<class Specie>
inline Boussinesq<Specie> operator==
(
    const Boussinesq<Specie>&,
    const Boussinesq<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const Boussinesq<Specie>&
);


/*---------------------------------------------------------------------------*\
                         Class Boussinesq Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class Boussinesq
:
    public Specie
{
    // Private Data

        //- Reference density
        scalar rho0_;

        //- Reference temperature
        scalar T0_;

        //- Thermal expansion coefficient
        scalar beta_;


public:

    // Constructors

        //- Construct from components
        inline Boussinesq
        (
            const Specie& sp,
            const scalar rho0,
            const scalar T0,
            const scalar beta
        );

        //- Construct from dictionary
        Boussinesq(const dictionary& dict);

        //- Construct as named copy
        inline Boussinesq
        (
            const word& name,
            const Boussinesq&
        );

        //- Construct and return a clone
        inline autoPtr<Boussinesq> clone() const;

        // Selector from dictionary
        inline static autoPtr<Boussinesq> New
        (
            const dictionary& dict
        );


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return
                "Boussinesq<"
              + word(Specie::typeName_()) + '>';
        }


        // Fundamental properties

            //- Is the equation of state is incompressible i.e. rho != f(p)
            static const bool incompressible = true;

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

        inline void operator+=(const Boussinesq&);
        inline void operator*=(const scalar);


    // Friend operators

        friend Boussinesq operator+ <Specie>
        (
            const Boussinesq&,
            const Boussinesq&
        );

        friend Boussinesq operator* <Specie>
        (
            const scalar s,
            const Boussinesq&
        );

        friend Boussinesq operator== <Specie>
        (
            const Boussinesq&,
            const Boussinesq&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie>
        (
            Ostream&,
            const Boussinesq&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "BoussinesqI.H"

#ifdef NoRepository
    #include "Boussinesq.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
