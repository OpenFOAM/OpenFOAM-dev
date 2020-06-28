/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2020 OpenFOAM Foundation
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
    Foam::linear

Description
    Linear equation of state with constant compressibility

    \verbatim
        rho = rho0 + psi*p
    \endverbatim

    Coefficient mixing is very inaccurate and not supported,
    so this equation of state is not applicable to mixtures.

SourceFiles
    linearI.H
    linear.C

\*---------------------------------------------------------------------------*/

#ifndef linear_H
#define linear_H

#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie> class linear;

template<class Specie>
inline linear<Specie> operator+
(
    const linear<Specie>&,
    const linear<Specie>&
);

template<class Specie>
inline linear<Specie> operator*
(
    const scalar,
    const linear<Specie>&
);

template<class Specie>
inline linear<Specie> operator==
(
    const linear<Specie>&,
    const linear<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const linear<Specie>&
);


/*---------------------------------------------------------------------------*\
                           Class linear Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class linear
:
    public Specie
{
    // Private Data

        //- Compressibility
        scalar psi_;

        //- The reference density
        scalar rho0_;


public:

    // Constructors

        //- Construct from components
        inline linear
        (
            const Specie& sp,
            const scalar psi,
            const scalar rho0
        );

        //- Construct from dictionary
        linear(const dictionary& dict);

        //- Construct as named copy
        inline linear(const word& name, const linear&);

        //- Construct and return a clone
        inline autoPtr<linear> clone() const;

        // Selector from dictionary
        inline static autoPtr<linear> New(const dictionary& dict);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "linear<" + word(Specie::typeName_()) + '>';
        }


        // Fundamental properties

            //- Is the equation of state is incompressible i.e. rho != f(p)
            static const bool incompressible = false;

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

        inline void operator+=(const linear&);
        inline void operator*=(const scalar);


    // Friend operators

        friend linear operator+ <Specie>
        (
            const linear&,
            const linear&
        );

        friend linear operator* <Specie>
        (
            const scalar s,
            const linear&
        );

        friend linear operator== <Specie>
        (
            const linear&,
            const linear&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie>
        (
            Ostream&,
            const linear&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "linearI.H"

#ifdef NoRepository
    #include "linear.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
