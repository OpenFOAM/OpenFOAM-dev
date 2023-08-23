/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

#include "fvMeshStitcherTools.H"
#include "surfaceFields.H"
#include "calculatedFvsPatchField.H"
#include "nonConformalFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField::Boundary>
Foam::fvMeshStitcherTools::origNcMagSfb(const fvMesh& mesh)
{
    const fvBoundaryMesh& fvbm = mesh.boundary();

    const surfaceScalarField::Boundary& magSfb =
        fvbm.mesh().magSf().boundaryField();

    tmp<surfaceScalarField::Boundary> tresult
    (
        new surfaceScalarField::Boundary
        (
            fvbm,
            surfaceScalarField::Internal::null(),
            calculatedFvsPatchField<scalar>::typeName
        )
    );

    surfaceScalarField::Boundary& result = tresult.ref();

    result == 0;

    forAll(fvbm, ncPatchi)
    {
        const fvPatch& fvp = fvbm[ncPatchi];

        if (!isA<nonConformalFvPatch>(fvp)) continue;

        const nonConformalFvPatch& ncFvp =
            refCast<const nonConformalFvPatch>(fvp);

        const label origPatchi = ncFvp.origPatchID();
        const fvPatch& origFvp = ncFvp.origPatch();

        result[origPatchi] +=
            fieldRMapSum
            (
                magSfb[ncPatchi],
                origFvp.size(),
                ncFvp.polyFaces(),
                origFvp.start()
            );
    }

    return tresult;
}


// ************************************************************************* //
