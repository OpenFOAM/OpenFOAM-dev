/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2020 OpenFOAM Foundation
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
    Foam::thermophysicalProperties

Description
    Base-class for thermophysical properties of solids, liquids and gases
    providing an interface compatible with the templated thermodynamics
    packages.

SourceFiles
    thermophysicalPropertiesI.H
    thermophysicalProperties.C

\*---------------------------------------------------------------------------*/

#ifndef thermophysicalProperties_H
#define thermophysicalProperties_H

#include "dictionary.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
class thermophysicalProperties;
Ostream& operator<<(Ostream& os, const thermophysicalProperties& l);

/*---------------------------------------------------------------------------*\
                  Class thermophysicalProperties Declaration
\*---------------------------------------------------------------------------*/

class thermophysicalProperties
{
    // Private Data

        //- Molecular weight [kg/kmol]
        scalar W_;


public:

    TypeName("thermophysicalProperties");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            thermophysicalProperties,
            ,
            (),
            ()
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            thermophysicalProperties,
            dictionary,
            (const dictionary& dict),
            (dict)
        );


    // Constructors

        //- Construct from molecular weight
        thermophysicalProperties(scalar W);

        //- Construct from dictionary
        thermophysicalProperties(const dictionary& dict);


    // Selectors

        //- Return a pointer to a new thermophysicalProperties created from name
        static autoPtr<thermophysicalProperties> New(const word& name);

        //- Return a pointer to a new thermophysicalProperties
        //  created from dictionary
        static autoPtr<thermophysicalProperties> New(const dictionary& dict);


    //- Destructor
    virtual ~thermophysicalProperties()
    {}


    // Member Functions

        // Physical constants which define the specie

            //- Molecular weight [kg/kmol]
            inline scalar W() const;

            //- Limit the temperature to be in the range Tlow_ to Thigh_
            inline scalar limit(const scalar T) const;


        // Fundamental equation of state properties

            //- Liquid density [kg/m^3]
            virtual scalar rho(scalar p, scalar T) const = 0;

            //- Liquid compressibility [s^2/m^2]
            //  Note: currently it is assumed the liquid is incompressible
            virtual scalar psi(scalar p, scalar T) const = 0;

            //- Return (Cp - Cv) [J/(kg K]
            //  Note: currently it is assumed the liquid is incompressible
            //  so CpMCv 0
            virtual scalar CpMCv(scalar p, scalar T) const = 0;


        // Fundamental thermodynamic properties

            //- Heat capacity at constant pressure [J/kg/K]
            virtual scalar Cp(const scalar p, const scalar T) const = 0;

            //- Sensible enthalpy [J/kg]
            virtual scalar Hs(const scalar p, const scalar T) const = 0;

            //- Enthalpy of formation [J/kg]
            virtual scalar Hf() const = 0;

            //- Absolute enthalpy [J/kg]
            virtual scalar Ha(const scalar p, const scalar T) const = 0;

            // Entropy [J/kg/K]
            virtual scalar S(const scalar p, const scalar T) const = 0;


        // Physical properties

            //- Liquid viscosity [Pa s]
            virtual scalar mu(scalar p, scalar T) const = 0;

            //- Liquid thermal conductivity  [W/m/K]
            virtual scalar kappa(scalar p, scalar T) const = 0;


    // I-O

        //- Read and set the properties present it the given dictionary
        void readIfPresent(const dictionary& dict);

        //- Write the function coefficients
        virtual void write(Ostream& os) const = 0;

        //- Ostream Operator
        friend Ostream& operator<<
        (
            Ostream& os,
            const thermophysicalProperties& l
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "thermophysicalPropertiesI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
