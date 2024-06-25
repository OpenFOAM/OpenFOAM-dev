/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2024 OpenFOAM Foundation
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
    Foam::fv::heatTransferCoefficientModel

Description
    Base class for heat transfer coefficient modelling used in heat transfer
    fvModels

SourceFiles
    heatTransferCoefficientModel.C

\*---------------------------------------------------------------------------*/

#ifndef heatTransferCoefficientModel_H
#define heatTransferCoefficientModel_H

#include "volFields.H"
#include "interRegionModel.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                          Class heatTransferCoefficientModel Declaration
\*---------------------------------------------------------------------------*/

class heatTransferCoefficientModel
{
protected:

    // Protected Member Data

        //- Reference to the mesh
        const fvMesh& mesh_;


public:

    //- Runtime type information
    TypeName("heatTransferCoefficientModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            heatTransferCoefficientModel,
            mesh,
            (
                const dictionary& dict,
                const fvMesh& mesh
            ),
            (dict, mesh)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            heatTransferCoefficientModel,
            model,
            (
                const dictionary& dict,
                const interRegionModel& model
            ),
            (dict, model)
        );


    // Constructors

        //- Construct from dictionary and mesh
        heatTransferCoefficientModel
        (
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Construct from dictionary and model
        heatTransferCoefficientModel
        (
            const word& modelType,
            const dictionary& dict,
            const interRegionModel& model
        );

        //- Disallow default bitwise copy construction
        heatTransferCoefficientModel
        (
            const heatTransferCoefficientModel&
        ) = delete;


    // Selectors

        //- Select from dictionary and mesh
        static autoPtr<heatTransferCoefficientModel> New
        (
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Select from dictionary and model
        static autoPtr<heatTransferCoefficientModel> New
        (
            const dictionary& dict,
            const interRegionModel& model
        );


    //- Destructor
    virtual ~heatTransferCoefficientModel();


    // Member Functions

        //- Get the heat transfer coefficient
        virtual tmp<volScalarField> htc() const = 0;

        //- Correct the heat transfer coefficient
        virtual void correct() = 0;

        //- Read dictionary
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
