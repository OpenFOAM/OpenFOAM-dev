/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2020 OpenFOAM Foundation
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

namespace Foam
{

class phasePair;

/*---------------------------------------------------------------------------*\
                    Class wallLubricationModel Declaration
\*---------------------------------------------------------------------------*/

class wallLubricationModel
:
    public wallDependentModel
{
protected:

    // Protected data

        //- Phase pair
        const phasePair& pair_;


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
                const phasePair& pair
            ),
            (dict, pair)
        );


    // Static Data Members

        //- Coefficient dimensions
        static const dimensionSet dimF;


    // Constructors

        //- Construct from components
        wallLubricationModel
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~wallLubricationModel();


    // Selectors

        static autoPtr<wallLubricationModel> New
        (
            const dictionary& dict,
            const phasePair& pair
        );


    // Member Functions

        //- Return phase-intensive wall lubrication force
        virtual tmp<volVectorField> Fi() const = 0;

        //- Return wall lubrication force
        virtual tmp<volVectorField> F() const;

        //- Return face wall lubrication force
        virtual tmp<surfaceScalarField> Ff() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
