/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "inverseFaceDistanceDiffusivity.H"
#include "surfaceFields.H"
#include "FvFaceCellWave.H"
#include "fvWallPoint.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(inverseFaceDistanceDiffusivity, 0);

    addToRunTimeSelectionTable
    (
        motionDiffusivity,
        inverseFaceDistanceDiffusivity,
        Istream
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::inverseFaceDistanceDiffusivity::inverseFaceDistanceDiffusivity
(
    const fvMesh& mesh,
    Istream& mdData
)
:
    motionDiffusivity(mesh),
    patchNames_(mdData)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::inverseFaceDistanceDiffusivity::~inverseFaceDistanceDiffusivity()
{}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField>
Foam::inverseFaceDistanceDiffusivity::operator()() const
{
    const labelHashSet patchIDs(mesh().boundaryMesh().patchSet(patchNames_));

    const surfaceVectorField::Boundary& CfBf = mesh().Cf().boundaryField();

    // Create changed face information
    List<labelPair> changedPatchAndFaces;
    List<fvWallPoint> changedFacesInfo;
    {
        label changedFacei = 0;
        forAllConstIter(labelHashSet, patchIDs, iter)
        {
            const label patchi = iter.key();
            changedFacei += mesh().boundary()[patchi].size();
        }

        changedPatchAndFaces.resize(changedFacei);
        changedFacesInfo.resize(changedFacei);

        changedFacei = 0;
        forAllConstIter(labelHashSet, patchIDs, iter)
        {
            const label patchi = iter.key();

            forAll(mesh().boundary()[patchi], patchFacei)
            {
                changedPatchAndFaces[changedFacei] =
                    labelPair(patchi, patchFacei);
                changedFacesInfo[changedFacei] =
                    fvWallPoint(CfBf[patchi][patchFacei], 0);

                changedFacei ++;
            }
        }
    }

    // Initialise wave storage
    List<fvWallPoint> internalFaceInfo(mesh().nInternalFaces());
    List<List<fvWallPoint>> patchFaceInfo
    (
        FvFaceCellWave<fvWallPoint>::template
        sizesListList<List<List<fvWallPoint>>>
        (
            FvFaceCellWave<fvWallPoint>::template
            listListSizes(mesh().boundary()),
            fvWallPoint()
        )
    );
    List<fvWallPoint> cellInfo(mesh().nCells());

    // Wave through the mesh
    FvFaceCellWave<fvWallPoint> wave
    (
        mesh(),
        changedPatchAndFaces,
        changedFacesInfo,
        internalFaceInfo,
        patchFaceInfo,
        cellInfo,
        mesh().globalData().nTotalCells() + 1 // max iterations
    );

    // Create the diffusivity field
    tmp<surfaceScalarField> tfaceDiffusivity
    (
        new surfaceScalarField
        (
            IOobject
            (
                "faceDiffusivity",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar(dimless, 1.0)
        )
    );
    surfaceScalarField& faceDiffusivity = tfaceDiffusivity.ref();

    // Convert waved distance data into diffusivities
    forAll(internalFaceInfo, facei)
    {
        const scalar dist = internalFaceInfo[facei].dist(wave.data());
        faceDiffusivity[facei] = 1/dist;
    }
    forAll(patchFaceInfo, patchi)
    {
        // Use cell distance on faces that are part of the patch set. This
        // avoids divide-by-zero issues.
        const bool useCellDist = patchIDs.found(patchi);

        const labelUList& patchCells = mesh().boundary()[patchi].faceCells();

        forAll(patchFaceInfo[patchi], patchFacei)
        {
            const scalar dist =
                useCellDist
              ? cellInfo[patchCells[patchFacei]].dist(wave.data())
              : patchFaceInfo[patchi][patchFacei].dist(wave.data());

            faceDiffusivity.boundaryFieldRef()[patchi][patchFacei] = 1/dist;
        }
    }

    return tfaceDiffusivity;
}


// ************************************************************************* //
