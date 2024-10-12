/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2024 OpenFOAM Foundation
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
    Foam::wallHeatTransferCoeffModel

Description
    Abstract base class for run time selection of heat transfer coefficient
    models.

SourceFiles
    wallHeatTransferCoeffModel.C

\*---------------------------------------------------------------------------*/

#ifndef wallHeatTransferCoeffModel_H
#define wallHeatTransferCoeffModel_H

#include "incompressibleMomentumTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class wallHeatTransferCoeffModel Declaration
\*---------------------------------------------------------------------------*/

class wallHeatTransferCoeffModel
{

public:

    //- Runtime type information
    TypeName("wallHeatTransferCoeffModel");


    // Declare run-time selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            wallHeatTransferCoeffModel,
            word,
            (
                const word& name,
                const fvMesh& mesh,
                const dictionary& dict
            ),
            (name, mesh, dict)
        );


    // Constructors

        //- Construct from name, mesh and dict
        wallHeatTransferCoeffModel
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary&
        )
        {}

        //- Disallow default bitwise copy construction
        wallHeatTransferCoeffModel(const wallHeatTransferCoeffModel&) = delete;

        //- Clone
        autoPtr<wallHeatTransferCoeffModel> clone() const
        {
            NotImplemented;
            return autoPtr<wallHeatTransferCoeffModel>(nullptr);
        }


    // Selectors

        //- Return a reference to the selected subset
        static autoPtr<wallHeatTransferCoeffModel> New
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary&
        );


    //- Destructor
    virtual ~wallHeatTransferCoeffModel()
    {}


    // Member Functions

        //- Read the wallHeatTransferCoeffModel data
        virtual bool read(const dictionary&) = 0;

        //- Interface for heat transfer coefficient
        virtual tmp<volScalarField> htcByRhoCp
        (
            const momentumTransportModel& mmtm,
            const labelHashSet& patches
        ) const = 0;

        //- Calculate wall shear stress
        tmp<surfaceVectorField> tau
        (
            const momentumTransportModel& mmtm,
            const fvMesh& mesh
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const wallHeatTransferCoeffModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
