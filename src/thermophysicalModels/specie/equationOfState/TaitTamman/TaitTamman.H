/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2023 OpenFOAM Foundation
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
    Foam::TaitTamman

Description
    Equation of state for liquid polymers:

    \verbatim
        v  = v0*(1 - C*log(1 + p/B))
        v0 = b1 + b2*dT,
        B  = b3 * exp(-b4*dT)
        dT = T - T0
    \endverbatim

    C = 0.0894 is a universal constant

Usage
    \table
        Property     | Description
        T0           | Reference temperature [K]
        b1           | Reference volume at T=T0 [m^3/kg]
        b2           | Coefficient for volume change with temperature [m^3/kg-K]
        b3           | Coefficient for pressure dependence [Pa]
        b4           | Coefficient for pressure dependence [1/K]
    \endtable

    Example specification of the Tait-Tamman equation of state
    for LNG:
    \verbatim
    equationOfState
    {
        T0          414.5;
        b1          12.74e-4;
        b2          10.26e-7;
        b3          9.263e7;
        b4          4.941e-3;
    }
    \endverbatim

SourceFiles
    TaitTammanI.H
    TaitTamman.C

\*---------------------------------------------------------------------------*/

#ifndef TaitTamman_H
#define TaitTamman_H

#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie> class TaitTamman;

template<class Specie>
inline TaitTamman<Specie> operator+
(
    const TaitTamman<Specie>&,
    const TaitTamman<Specie>&
);

template<class Specie>
inline TaitTamman<Specie> operator*
(
    const scalar,
    const TaitTamman<Specie>&
);

template<class Specie>
inline TaitTamman<Specie> operator==
(
    const TaitTamman<Specie>&,
    const TaitTamman<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const TaitTamman<Specie>&
);


/*---------------------------------------------------------------------------*\
                    Class TaitTamman Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class TaitTamman
:
    public Specie
{
    // Private Data

    //- Reference temperature
    scalar T0_;

    //- Reference volume at T=T0 [m^3/kg]
    scalar b1_;

    //- Coefficient for volume change with temperature [m^3/kg-K]
    scalar b2_;

    //- Coefficient for pressure dependence [Pa]
    scalar b3_;

    //- Coefficient for pressure dependence [1/K]
    scalar b4_;

    static scalar C_;

public:

    // Constructors

        //- Construct from components
        inline TaitTamman
        (
            const Specie& sp,
            const scalar T0,
            const scalar b1,
            const scalar b2,
            const scalar b3,
            const scalar b4
        );

        //- Construct from name and dictionary
        TaitTamman(const word& name, const dictionary& dict);

        //- Construct as named copy
        inline TaitTamman
        (
            const word& name,
            const TaitTamman&
        );

        //- Construct and return a clone
        inline autoPtr<TaitTamman> clone() const;


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "TaitTamman<" + word(Specie::typeName_()) + '>';
        }


        // Fundamental properties

            //- Is the equation of state is incompressible i.e. rho != f(p)
            static const bool incompressible = false;

            //- Is the equation of state is isochoric i.e. rho = const
            static const bool isochoric = false;

            //- Return density [kg/m^3]
            inline scalar rho(scalar p, scalar T) const;

            //- Return enthalpy contribution [J/kg]
            inline scalar h(const scalar p, const scalar T) const;

            //- Return Cp contribution [J/(kg K]
            inline scalar Cp(scalar p, scalar T) const;

            //- Return internal energy contribution [J/kg]
            inline scalar e(const scalar p, const scalar T) const;

            //- Return Cv contribution [J/(kg K]
            inline scalar Cv(scalar p, scalar T) const;

            // //- Return entropy contribution to the integral of Cp/T [J/kg/K]
            // inline scalar sp(const scalar p, const scalar T) const;

            // //- Return entropy contribution to the integral of Cv/T [J/kg/K]
            // inline scalar sv(const scalar p, const scalar T) const;

            //- Return compressibility [s^2/m^2]
            inline scalar psi(scalar p, scalar T) const;

            // //- Return compression factor []
            // inline scalar Z(scalar p, scalar T) const;

            //- Return (Cp - Cv) [J/(kg K]
            inline scalar CpMCv(scalar p, scalar T) const;

            // //- Return volumetric coefficient of thermal expansion [1/T]
            // inline scalar alphav(const scalar p, const scalar T) const;


        // IO

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member Operators

        inline void operator+=(const TaitTamman&);
        inline void operator*=(const scalar);


    // Friend operators

        friend TaitTamman operator+ <Specie>
        (
            const TaitTamman&,
            const TaitTamman&
        );

        friend TaitTamman operator* <Specie>
        (
            const scalar s,
            const TaitTamman&
        );

        friend TaitTamman operator== <Specie>
        (
            const TaitTamman&,
            const TaitTamman&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie>
        (
            Ostream&,
            const TaitTamman&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "TaitTammanI.H"

#ifdef NoRepository
    #include "TaitTamman.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
