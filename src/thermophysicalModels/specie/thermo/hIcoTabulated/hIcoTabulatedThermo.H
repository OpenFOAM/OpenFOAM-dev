/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2023 OpenFOAM Foundation
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
    Foam::hIcoTabulatedThermo

Description
    Enthalpy based thermodynamics package using non-uniform tabulated data for
    heat capacity vs temperature.

Usage

    \table
        Property     | Description
        hf           | Heat of formation
        sf           | Standard entropy
        Cp           | Specific heat at constant pressure vs temperature table
    \endtable

    Example of the specification of the thermodynamic properties:
    \verbatim
    thermodynamics
    {
        hf              0;
        sf              0;
        Cp
        {
            values
            (
                (200     1005)
                (350     1010)
                (400     1020)
            );
        }
    }
    \endverbatim

SourceFiles
    hIcoTabulatedThermoI.H
    hIcoTabulatedThermo.C

See also
    Foam::Function1s::UniformTable

\*---------------------------------------------------------------------------*/

#ifndef hIcoTabulatedThermo_H
#define hIcoTabulatedThermo_H

#include "integratedNonUniformTable1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class EquationOfState>
class hIcoTabulatedThermo;

template<class EquationOfState>
Ostream& operator<<
(
    Ostream&,
    const hIcoTabulatedThermo<EquationOfState>&
);


/*---------------------------------------------------------------------------*\
                      Class hIcoTabulatedThermo Declaration
\*---------------------------------------------------------------------------*/

template<class EquationOfState>
class hIcoTabulatedThermo
:
    public EquationOfState
{
    // Private Typedefs

        //- Table type
        typedef Function1s::integratedNonUniformTable integratedNonUniformTable;


    // Private Data

        //- Heat of formation
        scalar hf_;

        //- Standard entropy
        scalar sf_;

        //- Specific heat at constant pressure table [J/kg/K]
        integratedNonUniformTable Cp_;


public:

    // Constructors

        //- Construct from name and dictionary
        hIcoTabulatedThermo(const word& name, const dictionary& dict);

        //- Construct as a named copy
        inline hIcoTabulatedThermo(const word&, const hIcoTabulatedThermo&);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "hIcoTabulated<" + EquationOfState::typeName() + '>';
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
            const hIcoTabulatedThermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "hIcoTabulatedThermoI.H"

#ifdef NoRepository
    #include "hIcoTabulatedThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
