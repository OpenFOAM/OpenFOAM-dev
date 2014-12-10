/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::hConstThermo

Description
    Constant properties thermodynamics package
    templated into the equationOfState.

SourceFiles
    hConstThermoI.H
    hConstThermo.C

\*---------------------------------------------------------------------------*/

#ifndef hConstThermo_H
#define hConstThermo_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class EquationOfState> class hConstThermo;

template<class EquationOfState>
inline hConstThermo<EquationOfState> operator+
(
    const hConstThermo<EquationOfState>&,
    const hConstThermo<EquationOfState>&
);

template<class EquationOfState>
inline hConstThermo<EquationOfState> operator-
(
    const hConstThermo<EquationOfState>&,
    const hConstThermo<EquationOfState>&
);

template<class EquationOfState>
inline hConstThermo<EquationOfState> operator*
(
    const scalar,
    const hConstThermo<EquationOfState>&
);

template<class EquationOfState>
inline hConstThermo<EquationOfState> operator==
(
    const hConstThermo<EquationOfState>&,
    const hConstThermo<EquationOfState>&
);

template<class EquationOfState>
Ostream& operator<<
(
    Ostream&,
    const hConstThermo<EquationOfState>&
);


/*---------------------------------------------------------------------------*\
                           Class hConstThermo Declaration
\*---------------------------------------------------------------------------*/

template<class EquationOfState>
class hConstThermo
:
    public EquationOfState
{
    // Private data

        scalar Cp_;
        scalar Hf_;


    // Private Member Functions

        //- Construct from components
        inline hConstThermo
        (
            const EquationOfState& st,
            const scalar cp,
            const scalar hf
        );


public:

    // Constructors

        //- Construct from Istream
        hConstThermo(Istream&);

        //- Construct from dictionary
        hConstThermo(const dictionary& dict);

        //- Construct as named copy
        inline hConstThermo(const word&, const hConstThermo&);

        //- Construct and return a clone
        inline autoPtr<hConstThermo> clone() const;

        //- Selector from Istream
        inline static autoPtr<hConstThermo> New(Istream& is);

        //- Selector from dictionary
        inline static autoPtr<hConstThermo> New(const dictionary& dict);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "hConst<" + EquationOfState::typeName() + '>';
        }

        //- Limit the temperature to be in the range Tlow_ to Thigh_
        inline scalar limit(const scalar T) const;


        // Fundamental properties

            //- Heat capacity at constant pressure [J/(kmol K)]
            inline scalar cp(const scalar p, const scalar T) const;

            //- Absolute Enthalpy [J/kmol]
            inline scalar ha(const scalar p, const scalar T) const;

            //- Sensible enthalpy [J/kmol]
            inline scalar hs(const scalar p, const scalar T) const;

            //- Chemical enthalpy [J/kmol]
            inline scalar hc() const;

            //- Entropy [J/(kmol K)]
            inline scalar s(const scalar p, const scalar T) const;


        // I-O

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member operators

        inline void operator+=(const hConstThermo&);
        inline void operator-=(const hConstThermo&);


    // Friend operators

        friend hConstThermo operator+ <EquationOfState>
        (
            const hConstThermo&,
            const hConstThermo&
        );

        friend hConstThermo operator- <EquationOfState>
        (
            const hConstThermo&,
            const hConstThermo&
        );

        friend hConstThermo operator* <EquationOfState>
        (
            const scalar,
            const hConstThermo&
        );

        friend hConstThermo operator== <EquationOfState>
        (
            const hConstThermo&,
            const hConstThermo&
        );


    // IOstream Operators

        friend Ostream& operator<< <EquationOfState>
        (
            Ostream&,
            const hConstThermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "hConstThermoI.H"

#ifdef NoRepository
#   include "hConstThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
