/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2023 OpenFOAM Foundation
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
    Foam::wallLubricationModel

Description
    Model for the wall lubrication force between two phases

SourceFiles
    wallLubricationModel.C
    wallLubricationModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef wallLubricationModel_H
#define wallLubricationModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "wallDependentModel.H"
#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "BlendedInterfacialModel.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class wallLubricationModel Declaration
\*---------------------------------------------------------------------------*/

class wallLubricationModel
:
    public wallDependentModel
{
protected:

    // Protected member functions

        //- Zero-gradient wall-lubrication force at walls
        tmp<volVectorField> zeroGradWalls(tmp<volVectorField>) const;


public:

    //- Runtime type information
    TypeName("wallLubricationModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            wallLubricationModel,
            dictionary,
            (
                const dictionary& dict,
                const phaseInterface& interface
            ),
            (dict, interface)
        );


    // Static Data Members

        //- Coefficient dimensions
        static const dimensionSet dimF;

        //- Does this model require correcting on fixed flux boundaries?
        static const bool correctFixedFluxBCs = true;


    // Constructors

        //- Construct from a dictionary and an interface
        wallLubricationModel
        (
            const dictionary& dict,
            const phaseInterface& interface
        );


    //- Destructor
    virtual ~wallLubricationModel();


    // Selectors

        static autoPtr<wallLubricationModel> New
        (
            const dictionary& dict,
            const phaseInterface& interface,
            const bool outer=true
        );


    // Member Functions

        //- Return wall lubrication force
        virtual tmp<volVectorField> F() const = 0;

        //- Return face wall lubrication force
        virtual tmp<surfaceScalarField> Ff() const = 0;
};


/*---------------------------------------------------------------------------*\
                   Class blendedWallLubricationModel Declaration
\*---------------------------------------------------------------------------*/

class blendedWallLubricationModel
:
    public BlendedInterfacialModel<wallLubricationModel>
{
public:

    // Constructors

        //- Inherit base class constructors
        using
            BlendedInterfacialModel<wallLubricationModel>::
            BlendedInterfacialModel;


    // Selectors

        static autoPtr<blendedWallLubricationModel> New
        (
            const dictionary& dict,
            const phaseInterface& interface
        );


    // Member Functions

        //- Return wall lubrication force
        tmp<volVectorField> F() const;

        //- Return face wall lubrication force
        tmp<surfaceScalarField> Ff() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
