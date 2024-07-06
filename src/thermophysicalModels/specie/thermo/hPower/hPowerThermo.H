/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2024 OpenFOAM Foundation
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
    Foam::hPowerThermo

Description
    Enthalpy based thermodynamics package using a power function of temperature
    for the constant heat capacity at constant volume which is particularly
    suitable for solids at low temperatures:

    \verbatim
        Cp = c0*pow(T/Tref, n0)
    \endverbatim

Usage
    \table
        Property     | Description
        c0           | Reference heat capacity at constant pressure [J/kg/K]
        n0           | Exponent of the power function
        Tref         | Reference temperature [K]
        hf           | Heat of formation [J/kg]
    \endtable

    Example specification of hPowerThermo:
    \verbatim
    thermodynamics
    {
        c0          230;
        Tref        470;
        n0          3;
        hf          0;
    }
    \endverbatim

SourceFiles
    hPowerThermoI.H
    hPowerThermo.C

\*---------------------------------------------------------------------------*/

#ifndef hPowerThermo_H
#define hPowerThermo_H

#include "scalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class EquationOfState> class hPowerThermo;

template<class EquationOfState>
inline hPowerThermo<EquationOfState> operator+
(
    const hPowerThermo<EquationOfState>&,
    const hPowerThermo<EquationOfState>&
);

template<class EquationOfState>
inline hPowerThermo<EquationOfState> operator*
(
    const scalar,
    const hPowerThermo<EquationOfState>&
);


template<class EquationOfState>
inline hPowerThermo<EquationOfState> operator==
(
    const hPowerThermo<EquationOfState>&,
    const hPowerThermo<EquationOfState>&
);


template<class EquationOfState>
Ostream& operator<<
(
    Ostream&,
    const hPowerThermo<EquationOfState>&
);


/*---------------------------------------------------------------------------*\
                        Class hPowerThermo Declaration
\*---------------------------------------------------------------------------*/

template<class EquationOfState>
class hPowerThermo
:
    public EquationOfState
{
    // Private Data

        scalar c0_;
        scalar n0_;
        scalar Tref_;
        scalar hf_;


    // Private Member Functions

        //- Check given temperature is within the range of the fitted coeffs
        inline void checkT(const scalar T) const;


public:

    // Constructors

        //- Construct from components
        inline hPowerThermo
        (
            const EquationOfState& st,
            const scalar c0,
            const scalar n0,
            const scalar Tref,
            const scalar hf
        );

        //- Construct from name and dictionary
        hPowerThermo(const word& name, const dictionary& dict);

        //- Construct as a named copy
        inline hPowerThermo
        (
            const word&,
            const hPowerThermo&
        );

         //- Construct and return a clone
        inline autoPtr<hPowerThermo> clone() const;


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "hPower<" + EquationOfState::typeName() + '>';
        }

        //- Limit the temperature to be in the range Tlow_ to Thigh_
        inline scalar limit(const scalar T) const;


        // Fundamental properties

            //- Heat capacity at constant pressure [J/kg/K]
            inline scalar Cp(const scalar p, const scalar T) const;

            //- Absolute enthalpy [J/kg]
            inline scalar ha(const scalar p, const scalar T) const;

            //- Sensible enthalpy [J/kg]
            inline scalar hs(const scalar p, const scalar T) const;

            //- Enthalpy of formation [J/kg]
            inline scalar hf() const;

            //- Entropy [J/kg/K]
            inline scalar s(const scalar p, const scalar T) const;

            //- Gibbs free energy of the mixture in the standard state [J/kg]
            inline scalar gStd(const scalar T) const;

            #include "HtoEthermo.H"


        // Derivative term used for Jacobian

            //- Temperature derivative of heat capacity at constant pressure
            inline scalar dCpdT(const scalar p, const scalar T) const;


    // Member Operators

        inline void operator+=(const hPowerThermo&);


    // Friend operators

        friend hPowerThermo operator+ <EquationOfState>
        (
            const hPowerThermo&,
            const hPowerThermo&
        );

        friend hPowerThermo operator* <EquationOfState>
        (
            const scalar,
            const hPowerThermo&
        );


        friend hPowerThermo operator== <EquationOfState>
        (
            const hPowerThermo&,
            const hPowerThermo&
        );


    // Ostream Operator

        friend Ostream& operator<< <EquationOfState>
        (
            Ostream&,
            const hPowerThermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "hPowerThermoI.H"
    #include "hPowerThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
