/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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
    Foam::NSRDS

Description
    Incompressible equation of state using NSRDS function number 105 to
    calculate density as a function of the temperature.

    Reference:
    \verbatim
    Daubert, T. E., & Danner, R. P. (1985).
    Data compilation tables of properties of pure compounds.
    \endverbatim

Usage
    \table
        Property     | Description
        a            | Coefficient
        b            | Coefficient
        c            | Coefficient
        d            | Coefficient
    \endtable

    Example specification of the NSRDS equation of state for water, with
    coefficients as given by Daubert & Danner:
    \verbatim
    equationOfState
    {
        a               98.343885;
        b               0.30542;
        c               647.13;
        d               0.081;
    }
    \endverbatim

SourceFiles
    NSRDSI.H
    NSRDS.C

See also
    Foam::Function1s::NonUniformTable

\*---------------------------------------------------------------------------*/

#ifndef NSRDS_H
#define NSRDS_H

#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie>
class NSRDS;

template<class Specie>
NSRDS<Specie> operator+(const NSRDS<Specie>&, const NSRDS<Specie>&);

template<class Specie>
NSRDS<Specie> operator*(const scalar, const NSRDS<Specie>&);

template<class Specie>
NSRDS<Specie> operator==(const NSRDS<Specie>&, const NSRDS<Specie>&);

template<class Specie>
Ostream& operator<<(Ostream&, const NSRDS<Specie>&);


/*---------------------------------------------------------------------------*\
                        Class NSRDS Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class NSRDS
:
    public Specie
{
    // Private Data

        //- Coefficients
        scalar a_, b_, c_, d_;


    // Private Constructors

        //- Construct from name and parent dictionary and dictionary
        NSRDS
        (
            const word& name,
            const dictionary& dict,
            const dictionary& subDict
        );


public:

    // Constructors

        //- Construct from name and dictionary
        NSRDS(const word& name, const dictionary& dict);

        //- Construct as named copy
        inline NSRDS
        (
            const word& name,
            const NSRDS&
        );

        //- Construct and return a clone
        inline autoPtr<NSRDS> clone() const;


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "NSRDS<" + word(Specie::typeName_()) + '>';
        }


        // Fundamental properties

            //- Is the equation of state is incompressible i.e. rho != f(p)
            static const bool incompressible = true;

            //- Is the equation of state is isochoric i.e. rho = const
            static const bool isochoric = false;

            //- Return density [kg/m^3]
            inline scalar rho(const scalar p, const scalar T) const;

            //- Return enthalpy contribution [J/kg]
            inline scalar h(const scalar p, const scalar T) const;

            //- Return Cp contribution [J/(kg K]
            inline scalar Cp(const scalar p, const scalar T) const;

            //- Return internal energy contribution [J/kg]
            inline scalar e(const scalar p, const scalar T) const;

            //- Return Cv contribution [J/(kg K]
            inline scalar Cv(const scalar p, const scalar T) const;

            //- Return entropy contribution to the integral of Cp/T [J/kg/K]
            inline scalar sp(const scalar p, const scalar T) const;

            //- Return entropy contribution to the integral of Cv/T [J/kg/K]
            inline scalar sv(const scalar p, const scalar T) const;

            //- Return compressibility [s^2/m^2]
            inline scalar psi(const scalar p, const scalar T) const;

            //- Return compression factor []
            inline scalar Z(const scalar p, const scalar T) const;

            //- Return (Cp - Cv) [J/(kg K]
            inline scalar CpMCv(const scalar p, const scalar T) const;

            //- Return volumetric coefficient of thermal expansion [1/T]
            inline scalar alphav(const scalar p, const scalar T) const;


        // IO

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member Operators

        inline void operator+=(const NSRDS&);
        inline void operator*=(const scalar);


    // Friend operators

        friend NSRDS operator+ <Specie>(const NSRDS&, const NSRDS&);
        friend NSRDS operator* <Specie>(const scalar s, const NSRDS&);
        friend NSRDS operator== <Specie>(const NSRDS&, const NSRDS&);


    // Ostream Operator

        friend Ostream& operator<< <Specie>(Ostream&, const NSRDS&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "NSRDSI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "NSRDS.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
