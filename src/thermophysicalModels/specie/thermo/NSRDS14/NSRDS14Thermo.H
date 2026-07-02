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
    Foam::NSRDS14Thermo

Description
    Enthalpy based thermodynamics package using NSRDS function number 114 to
    calculate specific heat capacity at constant pressure as a function of the
    temperature.

    Note: We don't have a reference for this function or the coefficients. They
    are not in the book from which all the other NSRDS data is taken (Daubert &
    Danner (1985)). This largely exists in order to maintain the behaviour of
    existing tutorial cases (i.e., the aachenBomb). Use with caution.

Usage
    \table
        Property     | Description
        hf           | Heat of formation
        sf           | Standard entropy
        Tc           | Coefficient
        a            | Coefficient
        b            | Coefficient
        c            | Coefficient
        d            | Coefficient
    \endtable

    Example specification of NSRDS14 thermodynamics for n-heptane:
    \verbatim
    thermodynamics
    {
        hf              -2237496.08846327;
        sf              0;
        Tc              540.2;
        a               6.11976102401216;
        b               3137.69909384855;
        c               182.274175063868;
        d               -254.530511150515;
    }
    \endverbatim

SourceFiles
    NSRDS14ThermoI.H
    NSRDS14Thermo.C

See also
    Foam::FixedPolynomial

\*---------------------------------------------------------------------------*/

#ifndef NSRDS14Thermo_H
#define NSRDS14Thermo_H

#include "scalar.H"
#include "FixedLaurentPolynomial.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class EquationOfState>
class NSRDS14Thermo;

template<class EquationOfState>
Ostream& operator<<(Ostream&, const NSRDS14Thermo<EquationOfState>&);


/*---------------------------------------------------------------------------*\
                      Class NSRDS14Thermo Declaration
\*---------------------------------------------------------------------------*/

template<class EquationOfState>
class NSRDS14Thermo
:
    public EquationOfState
{
    // Private Data

        //- Heat of formation
        scalar hf_;

        //- Standard entropy
        scalar sf_;

        //- Coefficients
        scalar Tc_, a_, b_, c_, d_;

        //- Derived polynomial coefficients
        FixedLaurentPolynomial<scalar, -1, 6> p_;

        //- Reference sensible enthalpy
        scalar hsRef_;


    // Private Constructors

        //- Construct from name and parent dictionary and dictionary
        NSRDS14Thermo
        (
            const word& name,
            const dictionary& dict,
            const dictionary& subDict
        );


public:

    // Constructors

        //- Construct from name and dictionary
        NSRDS14Thermo(const word& name, const dictionary& dict);

        //- Construct as a named copy
        inline NSRDS14Thermo(const word&, const NSRDS14Thermo&);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "NSRDS14<" + EquationOfState::typeName() + '>';
        }

        //- Limit the temperature to be in the range Tlow_ to Thigh_
        inline scalar limit(const scalar) const;


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


        // I-O

            //- Write to Ostream
            void write(Ostream& os) const;


    // Ostream Operator

        friend Ostream& operator<< <EquationOfState>
        (
            Ostream&,
            const NSRDS14Thermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "NSRDS14ThermoI.H"

#ifdef NoRepository
    #include "NSRDS14Thermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
