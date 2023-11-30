/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    Foam::solidProperties

Description
    The thermophysical properties of a solid

SourceFiles
    solidProperties.C

\*---------------------------------------------------------------------------*/

#ifndef solidProperties_H
#define solidProperties_H

#include "typeInfo.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
class solidProperties;
Ostream& operator<<(Ostream&, const solidProperties&);

/*---------------------------------------------------------------------------*\
                       Class solidProperties Declaration
\*---------------------------------------------------------------------------*/

class solidProperties
{
    // Private Data

        //- Density [kg/m^3]
        scalar rho_;

        //- Specific heat capacity [J/kg/K]
        scalar Cp_;

        //- Thermal conductivity [W/m/K]
        scalar kappa_;

        //- Heat of formation [J/kg]
        scalar hf_;

        //- Emissivity
        scalar emissivity_;


public:

    //- Runtime type information
    TypeName("solid");


    // Declare run-time constructor selection tables

    declareRunTimeSelectionTable
    (
        autoPtr,
        solidProperties,
        ,
        (),
        ()
    );

    declareRunTimeSelectionTable
    (
        autoPtr,
        solidProperties,
        dictionary,
        (const dictionary& dict),
        (dict)
    );


    // Constructors

        //- Construct from components
        solidProperties
        (
            scalar rho,
            scalar Cp,
            scalar kappa,
            scalar hf,
            scalar emissivity
        );

        //- Construct from dictionary
        solidProperties(const dictionary& dict);

        //- Construct and return clone
        virtual autoPtr<solidProperties> clone() const
        {
            return autoPtr<solidProperties>(new solidProperties(*this));
        }


    // Selectors

        //- Return a pointer to a new solidProperties created from name
        static autoPtr<solidProperties> New(const word& name);

        //- Return a pointer to a new solidProperties created from dictionary
        static autoPtr<solidProperties> New(const dictionary& dict);


    //- Destructor
    virtual ~solidProperties()
    {}


    // Member Functions

        // Physical constants which define the solidProperties

            //- Density [kg/m^3]
            inline scalar rho() const;

            //- Specific heat capacity [J/kg/K]
            inline scalar Cp() const;

            //- Heat of formation [J/kg]
            inline scalar hf() const;

            //- Sensible enthalpy [J/kg]
            inline scalar hs(const scalar T) const;

            //- Absolute enthalpy [J/kg]
            inline scalar ha(scalar T) const;

            //- Sensible internal energy [J/kg]
            inline scalar es(const scalar T) const;

            //- Absolute internal energy [J/kg]
            inline scalar ea(scalar T) const;

            //- Thermal conductivity [W/m/K]
            inline scalar kappa() const;

            //- Emissivity []
            inline scalar emissivity() const;


    // I-O

        //- Read and set the properties present it the given dictionary
        void readIfPresent(const dictionary& dict);

        //- Write the solidProperties properties
        virtual void write(Ostream& os) const;

        //- Ostream Operator
        friend Ostream& operator<<(Ostream& os, const solidProperties& s);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "solidPropertiesI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
