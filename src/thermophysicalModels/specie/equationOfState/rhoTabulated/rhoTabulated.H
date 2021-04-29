/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2021 OpenFOAM Foundation
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
    Foam::rhoTabulated

Description
    Incompressible of equation of state using uniform tabulated
    density vs pressure and temperature.

    Thermodynamic contributions derived from tabulated density data are not
    sufficiently accurate and hence the H, Cp, E, Cv, Sp, Sv and CpMCv are not
    implemented.  Thus this equation of state should be used in conjunction with
    hTabulated or eTabulated thermodynamics in which both the energy and heat
    capacity are tabulated w.r.t. pressure and temperature avoiding the need
    for equation of state contributions.

Usage
    \table
        Property     | Description
        rho          | Density vs pressure and temperature table
    \endtable

    Example of the specification of the equation of state:
    \verbatim
    equationOfState
    {
        rho
        {
            pLow    1e4;
            pHigh   5e5;

            Tlow    200;
            Thigh   1500;

            values
            <m> <n>
            (
                (..........)
                .
                .
                .
                (..........)
            );
        }
    }
    \endverbatim

SourceFiles
    rhoTabulatedI.H
    rhoTabulated.C

See also
    Foam::Function1s::NonUniformTable

\*---------------------------------------------------------------------------*/

#ifndef rhoTabulated_H
#define rhoTabulated_H

#include "UniformTable2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie>
class rhoTabulated;

template<class Specie>
inline rhoTabulated<Specie> operator+
(
    const rhoTabulated<Specie>&,
    const rhoTabulated<Specie>&
);

template<class Specie>
inline rhoTabulated<Specie> operator*
(
    const scalar,
    const rhoTabulated<Specie>&
);

template<class Specie>
inline rhoTabulated<Specie> operator==
(
    const rhoTabulated<Specie>&,
    const rhoTabulated<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const rhoTabulated<Specie>&
);


/*---------------------------------------------------------------------------*\
                        Class rhoTabulated Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class rhoTabulated
:
    public Specie
{
    typedef Function2s::UniformTable<scalar> table2D;


    // Private Data

        //- Density table [kg/m^3]
        table2D rho_;


public:

    // Constructors

        //- Construct from components
        inline rhoTabulated
        (
            const Specie& sp,
            const table2D& rho
        );

        //- Construct from dictionary
        rhoTabulated(const dictionary& dict);

        //- Construct as named copy
        inline rhoTabulated(const word& name, const rhoTabulated&);

        //- Construct and return a clone
        inline autoPtr<rhoTabulated> clone() const;

        // Selector from dictionary
        inline static autoPtr<rhoTabulated> New(const dictionary& dict);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "rhoTabulated<" + word(Specie::typeName_()) + '>';
        }


        // Fundamental properties

            //- Is the equation of state is incompressible i.e. rho != f(p)
            static const bool incompressible = false;

            //- Is the equation of state is isochoric i.e. rho = const
            static const bool isochoric = false;

            //- Return density [kg/m^3]
            inline scalar rho(scalar p, scalar T) const;

            //- Return enthalpy contribution [J/kg]
            //  Not implemented
            inline scalar H(const scalar p, const scalar T) const;

            //- Return Cp contribution [J/(kg K]
            //  Not implemented
            inline scalar Cp(scalar p, scalar T) const;

            //- Return internal energy contribution [J/kg]
            //  Not implemented
            inline scalar E(const scalar p, const scalar T) const;

            //- Return Cv contribution [J/(kg K]
            //  Not implemented
            inline scalar Cv(scalar p, scalar T) const;

            //- Return entropy contribution to the integral of Cp/T [J/kg/K]
            //  Not implemented
            inline scalar Sp(const scalar p, const scalar T) const;

            //- Return entropy contribution to the integral of Cv/T [J/kg/K]
            //  Not implemented
            inline scalar Sv(const scalar p, const scalar T) const;

            //- Return compressibility [s^2/m^2]
            inline scalar psi(scalar p, scalar T) const;

            //- Return compression factor []
            inline scalar Z(scalar p, scalar T) const;

            //- Return (Cp - Cv) [J/(kg K]
            //  Not implemented
            inline scalar CpMCv(scalar p, scalar T) const;


        // IO

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member Operators

        inline void operator+=(const rhoTabulated&);
        inline void operator*=(const scalar);


    // Friend operators

        friend rhoTabulated operator+ <Specie>
        (
            const rhoTabulated&,
            const rhoTabulated&
        );

        friend rhoTabulated operator* <Specie>
        (
            const scalar s,
            const rhoTabulated&
        );

        friend rhoTabulated operator== <Specie>
        (
            const rhoTabulated&,
            const rhoTabulated&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie>
        (
            Ostream&,
            const rhoTabulated&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "rhoTabulatedI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "rhoTabulated.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
