/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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
    Foam::fv::heatTransferModel

Description
    Base class for heat transfer coefficient modelling used in heat transfer
    fvModels. Area per unit volume [1/m] (AoV) must be provided as a value in
    the coefficients dictionary or as a field in constant.

SourceFiles
    heatTransferModel.C

\*---------------------------------------------------------------------------*/

#ifndef heatTransferModel_H
#define heatTransferModel_H

#include "volFields.H"
#include "interRegionModel.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                          Class heatTransferModel Declaration
\*---------------------------------------------------------------------------*/

class heatTransferModel
{
    // Private Member Data

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Dictionary containing source coefficients
        dictionary coeffs_;

        //- Area per unit volume [1/m]
        dimensionedScalar AoV_;

        //- Area per unit volume [1/m]
        autoPtr<volScalarField> AoVPtr_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs();


public:

    //- Runtime type information
    TypeName("heatTransferModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            heatTransferModel,
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
            heatTransferModel,
            model,
            (
                const dictionary& dict,
                const interRegionModel& model
            ),
            (dict, model)
        );


    // Constructors

        //- Construct from dictionary and mesh
        heatTransferModel
        (
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Construct from dictionary and model
        heatTransferModel
        (
            const word& modelType,
            const dictionary& dict,
            const interRegionModel& model
        );

        //- Disallow default bitwise copy construction
        heatTransferModel(const heatTransferModel&) = delete;


    // Selectors

        //- Select from dictionary and mesh
        static autoPtr<heatTransferModel> New
        (
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Select from dictionary and model
        static autoPtr<heatTransferModel> New
        (
            const dictionary& dict,
            const interRegionModel& model
        );


    //- Destructor
    virtual ~heatTransferModel();


    // Member Functions

        //- Return reference to the mesh
        inline const fvMesh& mesh() const
        {
            return mesh_;
        }

        //- Return coeffs dictionary
        inline const dictionary& coeffs() const
        {
            return coeffs_;
        }

        //- Get the area per unit volume
        tmp<volScalarField> AoV() const;

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
