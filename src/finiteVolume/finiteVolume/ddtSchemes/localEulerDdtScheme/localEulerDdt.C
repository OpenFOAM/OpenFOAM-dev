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

\*---------------------------------------------------------------------------*/

#include "localEulerDdtScheme.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::word Foam::fv::localEulerDdt::rDeltaTName("rDeltaT");
Foam::word Foam::fv::localEulerDdt::rDeltaTfName("rDeltaTf");
Foam::word Foam::fv::localEulerDdt::rSubDeltaTName("rSubDeltaT");

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool Foam::fv::localEulerDdt::enabled(const fvMesh& mesh)
{
    return
        word(mesh.schemes().ddt("default"))
     == fv::localEulerDdtScheme<scalar>::typeName;
}


const Foam::volScalarField& Foam::fv::localEulerDdt::localRDeltaT
(
    const fvMesh& mesh
)
{
    return mesh.objectRegistry::lookupObject<volScalarField>
    (
        mesh.time().subCycling() ? rSubDeltaTName : rDeltaTName
    );
}


const Foam::surfaceScalarField& Foam::fv::localEulerDdt::localRDeltaTf
(
    const fvMesh& mesh
)
{
    return mesh.objectRegistry::lookupObject<surfaceScalarField>
    (
        rDeltaTfName
    );
}


Foam::tmp<Foam::volScalarField> Foam::fv::localEulerDdt::localRSubDeltaT
(
    const fvMesh& mesh,
    const label nAlphaSubCycles
)
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            rSubDeltaTName,
            nAlphaSubCycles
           *mesh.objectRegistry::lookupObject<volScalarField>
            (
                rDeltaTName
            )
        )
    );
}


// ************************************************************************* //
