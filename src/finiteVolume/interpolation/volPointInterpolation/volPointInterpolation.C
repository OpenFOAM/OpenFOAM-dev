/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "volPointInterpolation.H"
#include "fvMesh.H"
#include "volFields.H"
#include "pointFields.H"
#include "demandDrivenData.H"
#include "pointConstraints.H"
#include "surfaceFields.H"
#include "syncTools.H"
#include "cpuTime.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(volPointInterpolation, 0);
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::volPointInterpolation::volPointInterpolation(const fvMesh& vm)
:
    DemandDrivenMeshObject
    <
        fvMesh,
        DeletableMeshObject,
        volPointInterpolation
    >(vm),
    pointWeights_(mesh().points().size()),
    boundary_
    (
        SubList<face>
        (
            mesh().faces(),
            mesh().nFaces() - mesh().nInternalFaces(),
            mesh().nInternalFaces()
        ),
        mesh().points()
    ),
    boundaryPointWeights_(boundary_.meshPoints().size()),
    boundaryPointNbrWeights_(boundary_.meshPoints().size())
{
    if (debug)
    {
        Pout<< "volPointInterpolation::volPointInterpolation(const fvMesh&) : "
            << "constructing weighting factors"
            << endl;
    }

    const pointField& points = mesh().points();
    const labelListList& pointCells = mesh().pointCells();
    const polyBoundaryMesh& pbm = mesh().boundaryMesh();
    const fvBoundaryMesh& fvbm = mesh().boundary();

    // Cache calls to patch coupled flags
    boolList isCoupledPolyPatch(pbm.size(), false);
    boolList isCoupledFvPatch(fvbm.size(), false);
    forAll(isCoupledFvPatch, patchi)
    {
        isCoupledPolyPatch[patchi] = pbm[patchi].coupled();
        isCoupledFvPatch[patchi] = fvbm[patchi].coupled();
    }

    // Determine the factor to which a point has its values set by adjacent
    // boundary faces, rather than connected cells
    scalarList pointBoundaryFactor(mesh().nPoints(), Zero);
    forAll(boundary_.meshPoints(), bPointi)
    {
        const label pointi = boundary_.meshPoints()[bPointi];

        const labelList& pFaces = boundary_.pointFaces()[bPointi];

        forAll(pFaces, pPointFacei)
        {
            // Poly indices
            const label patchi = pbm.patchIndices()[pFaces[pPointFacei]];
            const label patchFacei =
                pbm.patchFaceIndices()[pFaces[pPointFacei]];

            // FV indices
            const labelUList patches =
                mesh().polyBFacePatches()[pFaces[pPointFacei]];
            const labelUList patchFaces =
                mesh().polyBFacePatchFaces()[pFaces[pPointFacei]];

            scalar nonCoupledMagSf = 0;
            forAll(patches, i)
            {
                if
                (
                    !isCoupledPolyPatch[patches[i]]
                 && !isCoupledFvPatch[patches[i]]
                )
                {
                    nonCoupledMagSf += fvbm[patches[i]].magSf()[patchFaces[i]];
                }
            }

            pointBoundaryFactor[pointi] =
                max
                (
                    pointBoundaryFactor[pointi],
                    nonCoupledMagSf/pbm[patchi].magFaceAreas()[patchFacei]
                );
        }
    }
    syncTools::syncPointList
    (
        mesh(),
        pointBoundaryFactor,
        maxEqOp<scalar>(),
        scalar(0)
    );

    // Calculate inverse distances between cell centres and points
    // and store in the weighting factor array
    forAll(points, pointi)
    {
        if (pointBoundaryFactor[pointi] > 1 - rootSmall) continue;

        pointWeights_[pointi].setSize(pointCells[pointi].size());

        const scalar f = pointBoundaryFactor[pointi];

        forAll(pointCells[pointi], pointCelli)
        {
            const label celli = pointCells[pointi][pointCelli];

            pointWeights_[pointi][pointCelli] =
                (1 - f)/mag(points[pointi] - mesh().C()[celli]);
        }
    }

    // Get the cell centres on the other side of coupled boundaries
    typename volVectorField::Boundary CBnf
    (
        mesh().boundary(),
        volVectorField::Internal::null(),
        calculatedFvPatchField<vector>::typeName
    );
    forAll(fvbm, patchi)
    {
        if
        (
            !isCoupledPolyPatch[patchi]
         && isCoupledFvPatch[patchi]
        )
        {
            CBnf[patchi] = fvbm[patchi].Cn() + fvbm[patchi].delta();
        }
    }

    // Calculate inverse distanced between boundary face centres and points and
    // store in the boundary weighting factor arrays
    forAll(boundary_.meshPoints(), bPointi)
    {
        const label pointi = boundary_.meshPoints()[bPointi];

        const labelList& pFaces = boundary_.pointFaces()[bPointi];

        boundaryPointWeights_[bPointi].resize(pFaces.size());
        boundaryPointNbrWeights_[bPointi].resize(pFaces.size());

        forAll(pFaces, bPointFacei)
        {
            // Poly indices
            const label patchi = pbm.patchIndices()[pFaces[bPointFacei]];
            const label patchFacei =
                pbm.patchFaceIndices()[pFaces[bPointFacei]];

            // FV indices
            const labelUList patches =
                mesh().polyBFacePatches()[pFaces[bPointFacei]];
            const labelUList patchFaces =
                mesh().polyBFacePatchFaces()[pFaces[bPointFacei]];

            boundaryPointWeights_[bPointi][bPointFacei].resize
            (
                patches.size(),
                0
            );
            boundaryPointNbrWeights_[bPointi][bPointFacei].resize
            (
                patches.size(),
                0
            );

            forAll(patches, i)
            {
                const scalar a =
                    fvbm[patches[i]].magSf()[patchFaces[i]]
                   /pbm[patchi].magFaceAreas()[patchFacei];

                const scalar f = pointBoundaryFactor[pointi];

                // If FV coupled only, add a weight to the neighbouring cell.
                // This is necessary because point synchronisation will not sum
                // the contributions across this interface as would be the case
                // with a poly coupled interface.
                if
                (
                    !isCoupledPolyPatch[patches[i]]
                 && isCoupledFvPatch[patches[i]]
                )
                {
                    const point C = CBnf[patches[i]][patchFaces[i]];
                    boundaryPointNbrWeights_[bPointi][bPointFacei][i] =
                        (1 - f)*a/mag(points[pointi] - C);
                }

                // If not coupled, add a weight to the boundary value
                if
                (
                    !isCoupledPolyPatch[patches[i]]
                 && !isCoupledFvPatch[patches[i]]
                )
                {
                    const point Cf = fvbm[patches[i]].Cf()[patchFaces[i]];
                    boundaryPointWeights_[bPointi][bPointFacei][i] =
                        f*a/mag(points[pointi] - Cf);
                }
            }
        }
    }

    // Construct a sum of weights
    pointScalarField sumWeights
    (
        IOobject
        (
            "volPointSumWeights",
            mesh().polyMesh::instance(),
            mesh()
        ),
        pointMesh::New(mesh()),
        dimensionedScalar(dimless, 0)
    );

    // Add the internal weights
    forAll(pointWeights_, pointi)
    {
        forAll(pointWeights_[pointi], i)
        {
            sumWeights[pointi] += pointWeights_[pointi][i];
        }
    }

    // Add the boundary weights
    forAll(boundary_.meshPoints(), bPointi)
    {
        const label pointi = boundary_.meshPoints()[bPointi];
        forAll(boundaryPointWeights_[bPointi], i)
        {
            forAll(boundaryPointWeights_[bPointi][i], j)
            {
                sumWeights[pointi] += boundaryPointWeights_[bPointi][i][j];
                sumWeights[pointi] += boundaryPointNbrWeights_[bPointi][i][j];
            }
        }
    }

    // Synchronise over conformal couplings
    syncTools::syncPointList
    (
        mesh(),
        sumWeights,
        plusEqOp<scalar>(),
        scalar(0)
    );

    // Normalise internal weights
    forAll(pointWeights_, pointi)
    {
        forAll(pointWeights_[pointi], i)
        {
            pointWeights_[pointi][i] /= sumWeights[pointi];
        }
    }

    // Normalise boundary weights
    forAll(boundary_.meshPoints(), bPointi)
    {
        const label pointi = boundary_.meshPoints()[bPointi];
        forAll(boundaryPointWeights_[bPointi], i)
        {
            forAll(boundaryPointWeights_[bPointi][i], j)
            {
                boundaryPointWeights_[bPointi][i][j] /= sumWeights[pointi];
                boundaryPointNbrWeights_[bPointi][i][j] /= sumWeights[pointi];
            }
        }
    }
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

Foam::volPointInterpolation::~volPointInterpolation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::volPointInterpolation::interpolateDisplacement
(
    const volVectorField& vf,
    pointVectorField& pf
) const
{
    interpolateUnconstrained(vf, pf);

    // Apply displacement constraints
    pointConstraints::New(pf.mesh()).constrainDisplacement(pf);
}


// ************************************************************************* //
