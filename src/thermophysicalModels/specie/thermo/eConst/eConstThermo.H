/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
    Internal energy based thermodynamics package using a constant heat capacity
    at constant volume:

    \verbatim
        e = Cv*(T - Tref) + Esref
    \endverbatim

Usage
    \table
        Property     | Description
        Cv           | Constant Heat capacity at constant volume [J/kg/K]
        Tref         | Reference temperature [K] (defaults to Tstd)
        Esref        | Reference sensible internal energy [J/kg] (defaults to 0)
        Hf           | Heat of formation [J/kg]
    \endtable

    Example specification of eConstThermo for air:
    \verbatim
    thermodynamics
    {
        Cv          724;
        Hf          0;
    }
    \endverbatim

SourceFiles
    eConstThermoI.H
    eConstThermo.C

\*---------------------------------------------------------------------------*/

#ifndef eConstThermo_H
#define eConstThermo_H

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
    // Private Data

        //- Heat capacity at constant volume [J/kg/K]
        scalar Cv_;

        //- Heat of formation [J/kg]
        scalar Hf_;

        //- Reference temperature around which to linearise [K]
        scalar Tref_;

        //- Reference sensible enthalpy around which to linearise [J/kg]
        scalar Esref_;


    // Private Member Functions

        //- Construct from components
        inline eConstThermo
        (
            const EquationOfState& st,
            const scalar Cv,
            const scalar Hf,
            const scalar Tref,
            const scalar Esref
        );


public:

    // Constructors

        //- Construct from dictionary
        eConstThermo(const dictionary& dict);

        //- Construct as named copy
        inline eConstThermo(const word&, const eConstThermo&);

        //- Construct and return a clone
        inline autoPtr<eConstThermo> clone() const;

        //- Selector from dictionary
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

            //- Heat capacity at constant volume [J/kg/K]
            inline scalar Cv(const scalar p, const scalar T) const;

            //- Sensible internal energy [J/kg]
            inline scalar Es(const scalar p, const scalar T) const;

            //- Absolute internal energy [J/kg]
            inline scalar Ea(const scalar p, const scalar T) const;

            //- Enthalpy of formation [J/kg]
            inline scalar Hf() const;

            //- Entropy [J/kg/K]
            inline scalar S(const scalar p, const scalar T) const;

            //- Gibbs free energy of the mixture in the standard state [J/kg]
            inline scalar Gstd(const scalar T) const;

            #include "EtoHthermo.H"


        // Derivative term used for Jacobian

            //- Temperature derivative of heat capacity at constant pressure
            inline scalar dCpdT(const scalar p, const scalar T) const;


        // I-O

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member Operators

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
