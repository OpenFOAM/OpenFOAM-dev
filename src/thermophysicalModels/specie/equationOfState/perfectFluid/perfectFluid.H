/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2023 OpenFOAM Foundation
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
    Foam::perfectFluid

Description
    Simple extension of the perfect gas equation of state to liquids by
    the addition of a constant density off-set:

    \verbatim
        rho = rho0 + p/(R*T)
    \endverbatim

    Coefficient mixing is very inaccurate and not supported, so this equation of
    state is not applicable to mixtures.

    This equation of state is rather inaccurate and has been superseded by
    rPolynomial which is much more accurate and supports mixtures.

Usage
    \table
        Property     | Description
        rho0         | Reference density
        R            | Fluid constant
    \endtable

    Example specification of the perfectFluid equation of state:
    \verbatim
    equationOfState
    {
        rho0        1000;
        R           100;
    }
    \endverbatim

SourceFiles
    perfectFluidI.H
    perfectFluid.C

\*---------------------------------------------------------------------------*/

#ifndef perfectFluid_H
#define perfectFluid_H

#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie> class perfectFluid;

template<class Specie>
inline perfectFluid<Specie> operator+
(
    const perfectFluid<Specie>&,
    const perfectFluid<Specie>&
);

template<class Specie>
inline perfectFluid<Specie> operator*
(
    const scalar,
    const perfectFluid<Specie>&
);

template<class Specie>
inline perfectFluid<Specie> operator==
(
    const perfectFluid<Specie>&,
    const perfectFluid<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const perfectFluid<Specie>&
);


/*---------------------------------------------------------------------------*\
                        Class perfectFluid Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class perfectFluid
:
    public Specie
{
    // Private Data

        //- Fluid constant
        scalar R_;

        //- The reference density
        scalar rho0_;


public:

    // Constructors

        //- Construct from components
        inline perfectFluid
        (
            const Specie& sp,
            const scalar R,
            const scalar rho0
        );

        //- Construct from name and dictionary
        perfectFluid(const word& name, const dictionary& dict);

        //- Construct as named copy
        inline perfectFluid(const word& name, const perfectFluid&);

        //- Construct and return a clone
        inline autoPtr<perfectFluid> clone() const;


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "perfectFluid<" + word(Specie::typeName_()) + '>';
        }


        // Fundamental properties

            //- Is the equation of state is incompressible i.e. rho != f(p)
            static const bool incompressible = false;

            //- Is the equation of state is isochoric i.e. rho = const
            static const bool isochoric = false;

            //- Return fluid constant [J/kg/K]
            inline scalar R() const;

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

            //- Return entropy contribution to the integral of Cp/T [J/kg/K]
            inline scalar sp(const scalar p, const scalar T) const;

            //- Return entropy contribution to the integral of Cv/T [J/kg/K]
            inline scalar sv(const scalar p, const scalar T) const;

            //- Return compressibility [s^2/m^2]
            inline scalar psi(scalar p, scalar T) const;

            //- Return compression factor []
            inline scalar Z(scalar p, scalar T) const;

            //- Return (Cp - Cv) [J/(kg K]
            inline scalar CpMCv(scalar p, scalar T) const;

            //- Return volumetric coefficient of thermal expansion [1/T]
            inline scalar alphav(const scalar p, const scalar T) const;


        // IO

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member Operators

        inline void operator+=(const perfectFluid&);
        inline void operator*=(const scalar);


    // Friend operators

        friend perfectFluid operator+ <Specie>
        (
            const perfectFluid&,
            const perfectFluid&
        );

        friend perfectFluid operator* <Specie>
        (
            const scalar s,
            const perfectFluid&
        );

        friend perfectFluid operator== <Specie>
        (
            const perfectFluid&,
            const perfectFluid&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie>
        (
            Ostream&,
            const perfectFluid&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "perfectFluidI.H"

#ifdef NoRepository
    #include "perfectFluid.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
