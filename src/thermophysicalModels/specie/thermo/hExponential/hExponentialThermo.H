/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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
    Foam::hExponentialThermo

Description
    Exponential properties thermodynamics package
    templated into the equationOfState.

SourceFiles
    hExponentialThermoI.H
    hExponentialThermo.C

\*---------------------------------------------------------------------------*/

#ifndef hExponentialThermo_H
#define hExponentialThermo_H

#include "scalar.H"
#include "dictionary.H"
#include "specie.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class EquationOfState> class hExponentialThermo;

template<class EquationOfState>
inline hExponentialThermo<EquationOfState> operator+
(
    const hExponentialThermo<EquationOfState>&,
    const hExponentialThermo<EquationOfState>&
);

template<class EquationOfState>
inline hExponentialThermo<EquationOfState> operator-
(
    const hExponentialThermo<EquationOfState>&,
    const hExponentialThermo<EquationOfState>&
);

template<class EquationOfState>
inline hExponentialThermo<EquationOfState> operator*
(
    const scalar,
    const hExponentialThermo<EquationOfState>&
);


template<class EquationOfState>
inline hExponentialThermo<EquationOfState> operator==
(
    const hExponentialThermo<EquationOfState>&,
    const hExponentialThermo<EquationOfState>&
);


template<class EquationOfState>
Ostream& operator<<
(
    Ostream&,
    const hExponentialThermo<EquationOfState>&
);


/*---------------------------------------------------------------------------*\
                         Class hExponentialThermo Declaration
\*---------------------------------------------------------------------------*/

template<class EquationOfState>
class hExponentialThermo
:
    public EquationOfState
{

    // Private data

        scalar c0_;
        scalar n0_;
        scalar Tref_;
        scalar Hf_;


    // Private Member Functions

        //- Check given temperature is within the range of the fitted coeffs
        inline void checkT(const scalar T) const;

        //- Integrate Cp expression
        inline scalar integrateCp(const scalar T) const;

        //- Construct from components
        inline hExponentialThermo
        (
            const EquationOfState& st,
            const scalar c0,
            const scalar n0,
            const scalar Tref,
            const scalar Hf
        );

public:

    // Constructors

        //- Construct from Istream
        hExponentialThermo(Istream&);

        //- Construct from dictionary
        hExponentialThermo(const dictionary&);

        //- Construct as a named copy
        inline hExponentialThermo
        (
            const word&,
            const hExponentialThermo&
        );

         //- Construct and return a clone
        inline autoPtr<hExponentialThermo> clone() const;

        //- Selector from Istream
        inline static autoPtr<hExponentialThermo> New(Istream& is);

        //- Selector from dictionary
        inline static autoPtr<hExponentialThermo> New(const dictionary& dict);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "hExponential<" + EquationOfState::typeName() + '>';
        }

        //- Limit the temperature to be in the range Tlow_ to Thigh_
        inline scalar limit(const scalar T) const;


        // Fundamental properties

            //- Heat capacity at constant pressure [J/(kg K)]
            inline scalar cp(const scalar p, const scalar T) const;

            //- Absolute enthalpy [J/kmol]
            inline scalar ha(const scalar p, const scalar T) const;

            //- Sensible enthalpy [J/kg]
            inline scalar hs(const scalar p, const scalar T) const;

            //- chemical enthalpy [J/kg]
            inline scalar hc() const;

            //- Entropy [J/(kmol K)]
            inline scalar s(const scalar p, const scalar T) const;


    // Member operators

        inline void operator+=(const hExponentialThermo&);
        inline void operator-=(const hExponentialThermo&);


    // Friend operators

        friend hExponentialThermo operator+ <EquationOfState>
        (
            const hExponentialThermo&,
            const hExponentialThermo&
        );

        friend hExponentialThermo operator- <EquationOfState>
        (
            const hExponentialThermo&,
            const hExponentialThermo&
        );

        friend hExponentialThermo operator* <EquationOfState>
        (
            const scalar,
            const hExponentialThermo&
        );


        friend hExponentialThermo operator== <EquationOfState>
        (
            const hExponentialThermo&,
            const hExponentialThermo&
        );


    // Ostream Operator

        friend Ostream& operator<< <EquationOfState>
        (
            Ostream&,
            const hExponentialThermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#ifdef NoRepository
#   include "hExponentialThermoI.H"
#   include "hExponentialThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
