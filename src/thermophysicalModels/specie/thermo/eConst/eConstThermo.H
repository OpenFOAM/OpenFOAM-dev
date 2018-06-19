/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::eConstThermo

Description
    Constant properties thermodynamics package templated on an equation of
    state

SourceFiles
    eConstThermoI.H
    eConstThermo.C

\*---------------------------------------------------------------------------*/

#ifndef eConstThermo_H
#define eConstThermo_H

#include "thermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class EquationOfState> class eConstThermo;

template<class EquationOfState>
inline eConstThermo<EquationOfState> operator+
(
    const eConstThermo<EquationOfState>&,
    const eConstThermo<EquationOfState>&
);

template<class EquationOfState>
inline eConstThermo<EquationOfState> operator*
(
    const scalar,
    const eConstThermo<EquationOfState>&
);

template<class EquationOfState>
inline eConstThermo<EquationOfState> operator==
(
    const eConstThermo<EquationOfState>&,
    const eConstThermo<EquationOfState>&
);

template<class EquationOfState>
Ostream& operator<<
(
    Ostream&,
    const eConstThermo<EquationOfState>&
);


/*---------------------------------------------------------------------------*\
                        Class eConstThermo Declaration
\*---------------------------------------------------------------------------*/

template<class EquationOfState>
class eConstThermo
:
    public EquationOfState
{
    // Private data

        //- Heat capacity at constant volume
        //  Note: input in [J/(kg K)], but internally uses [J/(kmol K)]
        scalar Cv_;

        //- Heat of formation
        //  Note: input in [J/kg], but internally uses [J/kmol]
        scalar Hf_;


    // Private Member Functions

        //- Construct from components
        inline eConstThermo
        (
            const EquationOfState& st,
            const scalar cv,
            const scalar hf
        );


public:

    // Constructors

        //- Construct from dictionary
        eConstThermo(const dictionary& dict);

        //- Construct as named copy
        inline eConstThermo(const word&, const eConstThermo&);

        //- Construct and return a clone
        inline autoPtr<eConstThermo> clone() const;

        // Selector from dictionary
        inline static autoPtr<eConstThermo> New(const dictionary& dict);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "eConst<" + EquationOfState::typeName() + '>';
        }

        //- Limit the temperature to be in the range Tlow_ to Thigh_
        inline scalar limit(const scalar T) const;


        // Fundamental properties

            //- Heat capacity at constant pressure [J/(kg K)]
            inline scalar Cp(const scalar p, const scalar T) const;

            //- Absolute Enthalpy [J/kg]
            inline scalar Ha(const scalar p, const scalar T) const;

            //- Sensible enthalpy [J/kg]
            inline scalar Hs(const scalar p, const scalar T) const;

            //- Chemical enthalpy [J/kg]
            inline scalar Hc() const;

            //- Entropy [J/(kg K)]
            inline scalar S(const scalar p, const scalar T) const;


        // Derivative term used for Jacobian

            //- Derivative of Gibbs free energy w.r.t. temperature
            inline scalar dGdT(const scalar p, const scalar T) const;

            //- Temperature derivative of heat capacity at constant pressure
            inline scalar dCpdT(const scalar p, const scalar T) const;



        // I-O

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member operators

        inline void operator+=(const eConstThermo&);


    // Friend operators

        friend eConstThermo operator+ <EquationOfState>
        (
            const eConstThermo&,
            const eConstThermo&
        );

        friend eConstThermo operator* <EquationOfState>
        (
            const scalar,
            const eConstThermo&
        );

        friend eConstThermo operator== <EquationOfState>
        (
            const eConstThermo&,
            const eConstThermo&
        );


    // IOstream Operators

        friend Ostream& operator<< <EquationOfState>
        (
            Ostream&,
            const eConstThermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "eConstThermoI.H"

#ifdef NoRepository
    #include "eConstThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
