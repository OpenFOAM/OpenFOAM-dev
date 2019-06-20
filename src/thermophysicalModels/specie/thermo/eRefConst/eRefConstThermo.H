/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2019 OpenFOAM Foundation
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
    Foam::eRefConstThermo

Description
    Constant properties thermodynamics package
    templated into the EquationOfState.

SourceFiles
    eRefConstThermoI.H
    eRefConstThermo.C

\*---------------------------------------------------------------------------*/

#ifndef eRefConstThermo_H
#define eRefConstThermo_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class EquationOfState> class eRefConstThermo;

template<class EquationOfState>
inline eRefConstThermo<EquationOfState> operator+
(
    const eRefConstThermo<EquationOfState>&,
    const eRefConstThermo<EquationOfState>&
);

template<class EquationOfState>
inline eRefConstThermo<EquationOfState> operator*
(
    const scalar,
    const eRefConstThermo<EquationOfState>&
);

template<class EquationOfState>
inline eRefConstThermo<EquationOfState> operator==
(
    const eRefConstThermo<EquationOfState>&,
    const eRefConstThermo<EquationOfState>&
);

template<class EquationOfState>
Ostream& operator<<
(
    Ostream&,
    const eRefConstThermo<EquationOfState>&
);


/*---------------------------------------------------------------------------*\
                           Class eRefConstThermo Declaration
\*---------------------------------------------------------------------------*/

template<class EquationOfState>
class eRefConstThermo
:
    public EquationOfState
{
    // Private Data

        scalar Cv_;
        scalar Hf_;
        scalar Tref_;
        scalar Eref_;


    // Private Member Functions

        //- Construct from components
        inline eRefConstThermo
        (
            const EquationOfState& st,
            const scalar cv,
            const scalar hf,
            const scalar tref,
            const scalar eref
        );


public:

    // Constructors

        //- Construct from dictionary
        eRefConstThermo(const dictionary& dict);

        //- Construct as named copy
        inline eRefConstThermo(const word&, const eRefConstThermo&);

        //- Construct and return a clone
        inline autoPtr<eRefConstThermo> clone() const;

        //- Selector from dictionary
        inline static autoPtr<eRefConstThermo> New(const dictionary& dict);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "eRefConst<" + EquationOfState::typeName() + '>';
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

            //- Chemical enthalpy [J/kg]
            inline scalar Hc() const;

            //- Entropy [J/kg/K]
            inline scalar S(const scalar p, const scalar T) const;

            #include "EtoHthermo.H"


        // Derivative term used for Jacobian

            //- Derivative of Gibbs free energy w.r.t. temperature
            inline scalar dGdT(const scalar p, const scalar T) const;

            //- Temperature derivative of heat capacity at constant pressure
            inline scalar dCpdT(const scalar p, const scalar T) const;


        // I-O

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member Operators

        inline void operator+=(const eRefConstThermo&);


    // Friend operators

        friend eRefConstThermo operator+ <EquationOfState>
        (
            const eRefConstThermo&,
            const eRefConstThermo&
        );

        friend eRefConstThermo operator* <EquationOfState>
        (
            const scalar,
            const eRefConstThermo&
        );

        friend eRefConstThermo operator== <EquationOfState>
        (
            const eRefConstThermo&,
            const eRefConstThermo&
        );


    // IOstream Operators

        friend Ostream& operator<< <EquationOfState>
        (
            Ostream&,
            const eRefConstThermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "eRefConstThermoI.H"

#ifdef NoRepository
    #include "eRefConstThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
