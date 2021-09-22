/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2021 OpenFOAM Foundation
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
    Foam::eIcoTabulatedThermo

Description
    Internal energy based thermodynamics package using non-uniform tabulated
    data for heat capacity vs temperature.

Usage
    \table
        Property     | Description
        Hf           | Heat of formation
        Sf           | Standard entropy
        Cv           | Specific heat at constant volume vs temperature table
    \endtable

    Example of the specification of the thermodynamic properties:
    \verbatim
    thermodynamics
    {
        Hf              0;
        Sf              0;
        Cv
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
    eIcoTabulatedThermoI.H
    eIcoTabulatedThermo.C

See also
    Foam::Function1s::UniformTable

\*---------------------------------------------------------------------------*/

#ifndef eIcoTabulatedThermo_H
#define eIcoTabulatedThermo_H

#include "integratedNonUniformTable1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class EquationOfState>
class eIcoTabulatedThermo;

template<class EquationOfState>
Ostream& operator<<
(
    Ostream&,
    const eIcoTabulatedThermo<EquationOfState>&
);


/*---------------------------------------------------------------------------*\
                      Class eIcoTabulatedThermo Declaration
\*---------------------------------------------------------------------------*/

template<class EquationOfState>
class eIcoTabulatedThermo
:
    public EquationOfState
{
    typedef Function1s::integratedNonUniformTable integratedNonUniformTable;


    // Private Data

        //- Heat of formation
        scalar Hf_;

        //- Standard entropy
        scalar Sf_;

        //- Specific heat at constant volume table [J/kg/K]
        integratedNonUniformTable Cv_;


public:

    // Constructors

        //- Construct from dictionary
        eIcoTabulatedThermo(const dictionary& dict);

        //- Construct as a named copy
        inline eIcoTabulatedThermo(const word&, const eIcoTabulatedThermo&);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "eIcoTabulated<" + EquationOfState::typeName() + '>';
        }

        //- Limit the temperature to be in the range Tlow_ to Thigh_
        inline scalar limit(const scalar) const;


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


    // Ostream Operator

        friend Ostream& operator<< <EquationOfState>
        (
            Ostream&,
            const eIcoTabulatedThermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "eIcoTabulatedThermoI.H"

#ifdef NoRepository
    #include "eIcoTabulatedThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
