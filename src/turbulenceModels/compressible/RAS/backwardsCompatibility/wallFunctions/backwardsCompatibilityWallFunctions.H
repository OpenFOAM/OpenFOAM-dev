/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::compressible

Description
    Auto creation of fields to provide backwards compatibility with
    runtime selectable wall functions

SourceFiles
    backwardsCompatibilityWallFunctions.C
    backwardsCompatibilityWallFunctionsTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef backwardsCompatibilityWallFunctions_H
#define backwardsCompatibilityWallFunctions_H

#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
    //- mut
    tmp<volScalarField> autoCreateMut
    (
        const word& fieldName,
        const fvMesh& mesh
    );

    //- mut for Low-Reynolds number models
    tmp<volScalarField> autoCreateLowReMut
    (
        const word& fieldName,
        const fvMesh& mesh
    );

    //- alphat
    tmp<volScalarField> autoCreateAlphat
    (
        const word& fieldName,
        const fvMesh& mesh
    );

    //- epsilon
    tmp<volScalarField> autoCreateEpsilon
    (
        const word& fieldName,
        const fvMesh& mesh
    );

    //- omega
    tmp<volScalarField> autoCreateOmega
    (
        const word& fieldName,
        const fvMesh& mesh
    );

    //- k
    tmp<volScalarField> autoCreateK
    (
        const word& fieldName,
        const fvMesh& mesh
    );

    //- Q
    tmp<volScalarField> autoCreateQ
    (
        const word& fieldName,
        const fvMesh& mesh
    );

    //- R
    tmp<volSymmTensorField> autoCreateR
    (
        const word& fieldName,
        const fvMesh& mesh
    );

    //- Helper function to create the new field
    template<class Type, class PatchType>
    tmp<GeometricField<Type, fvPatchField, volMesh> >
    autoCreateWallFunctionField
    (
        const word& fieldName,
        const fvMesh& mesh
    );


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "backwardsCompatibilityWallFunctionsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
