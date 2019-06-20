/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2019 OpenFOAM Foundation
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
    Foam::rhoConst

Description
    RhoConst (rho = const) of state.

SourceFiles
    rhoConstI.H
    rhoConst.C

\*---------------------------------------------------------------------------*/

#ifndef rhoConst_H
#define rhoConst_H

#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie> class rhoConst;

template<class Specie>
inline rhoConst<Specie> operator+
(
    const rhoConst<Specie>&,
    const rhoConst<Specie>&
);

template<class Specie>
inline rhoConst<Specie> operator*
(
    const scalar,
    const rhoConst<Specie>&
);

template<class Specie>
inline rhoConst<Specie> operator==
(
    const rhoConst<Specie>&,
    const rhoConst<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const rhoConst<Specie>&
);


/*---------------------------------------------------------------------------*\
                        Class rhoConst Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class rhoConst
:
    public Specie
{
    // Private Data

        //- Density
        scalar rho_;


public:

    // Constructors

        //- Construct from components
        inline rhoConst(const Specie& sp, const scalar rho);

        //- Construct from dictionary
        rhoConst(const dictionary& dict);

        //- Construct as named copy
        inline rhoConst(const word& name, const rhoConst&);

        //- Construct and return a clone
        inline autoPtr<rhoConst> clone() const;


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "rhoConst<" + word(Specie::typeName_()) + '>';
        }


        // Fundamental properties

            //- Is the equation of state is incompressible i.e. rho != f(p)
            static const bool incompressible = true;

            //- Is the equation of state is isochoric i.e. rho = const
            static const bool isochoric = true;

            //- Return density [kg/m^3]
            inline scalar rho(scalar p, scalar T) const;

            //- Return enthalpy departure [J/kg]
            inline scalar H(const scalar p, const scalar T) const;

            //- Return Cp departure [J/(kg K]
            inline scalar Cp(scalar p, scalar T) const;

            //- Return internal energy departure [J/kg]
            inline scalar E(const scalar p, const scalar T) const;

            //- Return Cv departure [J/(kg K]
            inline scalar Cv(scalar p, scalar T) const;

            //- Return entropy [J/kg/K]
            inline scalar S(const scalar p, const scalar T) const;

            //- Return compressibility rho/p [s^2/m^2]
            inline scalar psi(scalar p, scalar T) const;

            //- Return compression factor []
            inline scalar Z(scalar p, scalar T) const;

            //- Return (Cp - Cv) [J/(kg K]
            inline scalar CpMCv(scalar p, scalar T) const;


        // IO

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member Operators

        inline void operator+=(const rhoConst&);
        inline void operator*=(const scalar);


    // Friend operators

        friend rhoConst operator+ <Specie>
        (
            const rhoConst&,
            const rhoConst&
        );

        friend rhoConst operator* <Specie>
        (
            const scalar s,
            const rhoConst&
        );

        friend rhoConst operator== <Specie>
        (
            const rhoConst&,
            const rhoConst&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie>
        (
            Ostream&,
            const rhoConst&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "rhoConstI.H"

#ifdef NoRepository
    #include "rhoConst.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
