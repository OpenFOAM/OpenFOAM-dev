/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2023 OpenFOAM Foundation
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
    Foam::fv::heatTransferAoV

Description
    Class to handle area per unit volume [1/m] (AoV) for heat transfer
    fvModels, which must be provided as a value in the coefficients dictionary
    or as a field in constant.

SourceFiles
    heatTransferAoV.C

\*---------------------------------------------------------------------------*/

#ifndef heatTransferAoV_H
#define heatTransferAoV_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                          Class heatTransferAoV Declaration
\*---------------------------------------------------------------------------*/

class heatTransferAoV
{
    // Private Member Data

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Area per unit volume [1/m]
        dimensionedScalar AoV_;

        //- Area per unit volume [1/m]
        autoPtr<volScalarField> AoVPtr_;


    // Private Member Functions

        //- Read the area per unit volume
        void readCoeffs(const dictionary& dict);


public:

    // Constructors

        //- Construct from dictionary and mesh
        heatTransferAoV
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~heatTransferAoV();


    // Member Functions

        //- Get the area per unit volume
        tmp<volScalarField> AoV() const;

        //- Read dictionary
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
