/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2022 OpenFOAM Foundation
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
    Foam::interfaceSaturationTemperatureModel

Description
    Wrapper around saturationTemperatureModel to facilitate convenient
    construction on interfaces

SourceFiles
    interfaceSaturationTemperatureModel.C
    interfaceSaturationTemperatureModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef interfaceSaturationTemperatureModel_H
#define interfaceSaturationTemperatureModel_H

#include "volFields.H"
#include "dictionary.H"
#include "phaseInterface.H"
#include "saturationTemperatureModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class interfaceSaturationTemperatureModel Declaration
\*---------------------------------------------------------------------------*/

class interfaceSaturationTemperatureModel
:
    public regIOobject
{
    // Private Data

        //- The saturation temperature model
        autoPtr<saturationTemperatureModel> saturationModel_;

        //- The interface
        const phaseInterface interface_;


public:

    //- Runtime type information
    TypeName("saturationTemperatureModel");


    //- Declare runtime construction
    declareRunTimeSelectionTable
    (
        autoPtr,
        interfaceSaturationTemperatureModel,
        dictionary,
        (
            const dictionary& dict,
            const phaseInterface& interface
        ),
        (dict, interface)
    );


    // Constructors

        //- Construct from a dictionary and an interface
        interfaceSaturationTemperatureModel
        (
            const dictionary& dict,
            const phaseInterface& interface
        );

    // Selectors

        static autoPtr<interfaceSaturationTemperatureModel> New
        (
            const dictionary& dict,
            const phaseInterface& interface,
            const bool outer=true
        );

    //- Destructor
    virtual ~interfaceSaturationTemperatureModel();


    // Member Functions

        //- Access the interface
        const phaseInterface& interface() const;

        //- Saturation temperature
        tmp<volScalarField> Tsat(const volScalarField& p) const;

        // Dummy write for regIOobject
        bool writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
