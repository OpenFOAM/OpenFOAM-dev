/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2023 OpenFOAM Foundation
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

#include "MPLIC.H"
#include "MPLICcell.H"
#include "volPointInterpolation.H"
#include "syncTools.H"
#include "slicedSurfaceFields.H"
#include "upwind.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(MPLIC, 0);

    surfaceInterpolationScheme<scalar>::addMeshFluxConstructorToTable<MPLIC>
        addMPLICScalarMeshFluxConstructorToTable_;
}


// * * * * * * * * * * * * * * * Private Functions  * * * * * * * * * * * * * //

void Foam::MPLIC::setCellAlphaf
(
    const label celli,
    const scalarField& phi,
    scalarField& alphaf,
    boolList& correctedFaces,
    const DynamicList<scalar>& cellAlphaf,
    const fvMesh& mesh
) const
{
    // Face owners reference
    const labelList& own = mesh.faceOwner();

    // The cell face labels
    const labelList& cFaces = mesh.cells()[celli];

    // Fill the alphaf with surface interpolation in direction of the flow
    forAll(cFaces, i)
    {
        const label facei = cFaces[i];
        const scalar phiSigni = sign(phi[facei]);

        if
        (
            (own[facei] == celli && phiSigni == 1)
         || (own[facei] != celli && phiSigni == -1)
        )
        {
            alphaf[facei] = cellAlphaf[i];
            correctedFaces[facei] = true;
        }
    }
}


Foam::tmp<Foam::surfaceScalarField> Foam::MPLIC::surfaceAlpha
(
    const volScalarField& alpha,
    const surfaceScalarField& phi,
    scalarField& initAlphaf,
    const bool unweighted,
    const scalar tol,
    const bool isMPLIC
) const
{
    // Finite volume mesh reference
    const fvMesh& mesh = alpha.mesh();

    // Reference to primitive mesh
    const primitiveMesh& primMesh = mesh;

    // Velocity field reference
    const volVectorField& U
    (
        mesh.lookupObject<const volVectorField>
        (
            IOobject::groupName("U", phi.group())
        )
    );

    // Interpolate alpha from volume to the points of the mesh
    const scalarField alphap
    (
        volPointInterpolation::New(mesh).interpolate(alpha)
    );

    // Interpolate U from cell centres to the points of the mesh
    vectorField Up;

    if (!unweighted)
    {
        Up = volPointInterpolation::New(mesh).interpolate(U);
    }

    // Flatten down phi flux field
    const scalarField splicedPhi
    (
        slicedSurfaceScalarField
        (
            IOobject
            (
                "splicedPhi",
                mesh.time().name(),
                mesh
            ),
            phi,
            false
        ).splice()
    );

    scalarField alphaf(mesh.nFaces(), 0);

    // Mark which faces are corrected by MPLIC
    boolList correctedFaces(mesh.nFaces(), false);

    // Construct class for cell cut
    MPLICcell cutCell(unweighted, isMPLIC);

    // Loop through all the cells
    forAll(mesh.cells(), celli)
    {
        if (alpha[celli] < (1 - tol) && alpha[celli] > tol)
        {
            // Store cell information
            const MPLICcellStorage cellInfo
            (
                primMesh,
                alphap,
                Up,
                alpha[celli],
                U[celli],
                celli
            );

            // Volume ratio matching algorithm
            if (cutCell.matchAlpha(cellInfo))
            {
                // Fill cutCell.alphaf() with face values from this cell
                setCellAlphaf
                (
                    celli,
                    splicedPhi,
                    alphaf,
                    correctedFaces,
                    cutCell.alphaf(),
                    mesh
                );
            }
        }
    }

    // Synchronise across the processor and cyclic patches
    syncTools::syncFaceList(mesh, alphaf, plusEqOp<scalar>());
    syncTools::syncFaceList(mesh, correctedFaces, plusEqOp<bool>());

    // Correct selected faces
    forAll(correctedFaces, facei)
    {
        if (correctedFaces[facei])
        {
            initAlphaf[facei] = alphaf[facei];
        }
    }

    if (!mesh.conformal())
    {
        FatalErrorInFunction
            << "The " << type() << " scheme is not compatible with "
            << "non-conformal meshes" << exit(FatalError);
    }

    // Convert the alphaPhi spliced field into a surfaceScalarField
    tmp<surfaceScalarField> tsplicedAlpha
    (
        surfaceScalarField::New
        (
            "alphaf",
            slicedSurfaceScalarField
            (
                IOobject
                (
                    "alphaf",
                    mesh.time().name(),
                    mesh
                ),
                mesh,
                dimless,
                initAlphaf,
                false
            ),
            fvsPatchField<scalar>::calculatedType()
        )
    );
    surfaceScalarField& splicedAlpha = tsplicedAlpha.ref();

    forAll(mesh.boundary(), patchi)
    {
        if (alpha.boundaryField()[patchi].fixesValue())
        {
            splicedAlpha.boundaryFieldRef()[patchi] =
                alpha.boundaryField()[patchi];
        }
    }

    return tsplicedAlpha;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField> Foam::MPLIC::interpolate
(
    const VolField<scalar>& vf
) const
{
    tmp<surfaceScalarField> tvff(upwind<scalar>(mesh(), phi_).interpolate(vf));

    scalarField splicedTvff
    (
        slicedSurfaceScalarField
        (
            IOobject
            (
                "splicedTvff",
                mesh().time().name(),
                mesh()
            ),
            tvff,
            false
        ).splice()
    );

    return surfaceAlpha(vf, phi_, splicedTvff, true, 1e-6);
}


// ************************************************************************* //
