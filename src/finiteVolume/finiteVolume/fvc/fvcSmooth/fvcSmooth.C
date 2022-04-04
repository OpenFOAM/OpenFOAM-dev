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

#include "fvcSmooth.H"
#include "volFields.H"
#include "FvFaceCellWave.H"
#include "smoothData.H"
#include "sweepData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::fvc::smooth
(
    volScalarField& field,
    const scalar coeff
)
{
    const fvMesh& mesh = field.mesh();
    scalar maxRatio = 1 + coeff;

    DynamicList<labelPair> changedPatchAndFaces(mesh.nFaces()/100 + 100);
    DynamicList<smoothData> changedFacesInfo(changedPatchAndFaces.size());

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    forAll(owner, facei)
    {
        const label own = owner[facei];
        const label nbr = neighbour[facei];

        // Check if owner value much larger than neighbour value or vice versa
        if (field[own] > maxRatio*field[nbr])
        {
            changedPatchAndFaces.append(labelPair(-1, facei));
            changedFacesInfo.append(smoothData(field[own]));
        }
        else if (field[nbr] > maxRatio*field[own])
        {
            changedPatchAndFaces.append(labelPair(-1, facei));
            changedFacesInfo.append(smoothData(field[nbr]));
        }
    }

    // Insert all faces of coupled patches - FvFaceCellWave will correct them
    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& patch = mesh.boundary()[patchi];

        if (patch.coupled())
        {
            forAll(patch, patchFacei)
            {
                const label own = patch.faceCells()[patchFacei];

                changedPatchAndFaces.append(labelPair(patchi, patchFacei));
                changedFacesInfo.append(smoothData(field[own]));
            }
        }
    }

    changedPatchAndFaces.shrink();
    changedFacesInfo.shrink();

    // Set initial field on cells
    List<smoothData> cellData(mesh.nCells());
    forAll(field, celli)
    {
        cellData[celli] = field[celli];
    }

    // Set initial field on faces
    List<smoothData> internalFaceData(mesh.nInternalFaces());
    List<List<smoothData>> patchFaceData
    (
        FvFaceCellWave<smoothData, smoothData::trackData>::template
        sizesListList<List<List<smoothData>>>
        (
            FvFaceCellWave<smoothData, smoothData::trackData>::template
            listListSizes(mesh.boundary()),
            smoothData()
        )
    );

    // Create track data
    smoothData::trackData td;
    td.maxRatio = maxRatio;

    // Propagate information over whole domain
    FvFaceCellWave<smoothData, smoothData::trackData> smoothData
    (
        mesh,
        changedPatchAndFaces,
        changedFacesInfo,
        internalFaceData,
        patchFaceData,
        cellData,
        mesh.globalData().nTotalCells() + 1, // max iterations
        td
    );

    forAll(field, celli)
    {
        field[celli] = cellData[celli].value();
    }

    field.correctBoundaryConditions();
}


void Foam::fvc::spread
(
    volScalarField& field,
    const volScalarField& alpha,
    const label nLayers,
    const scalar alphaDiff,
    const scalar alphaMax,
    const scalar alphaMin
)
{
    const fvMesh& mesh = field.mesh();

    DynamicList<labelPair> changedPatchAndFaces(mesh.nFaces()/100 + 100);
    DynamicList<smoothData> changedFacesInfo(changedPatchAndFaces.size());

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    forAll(owner, facei)
    {
        const label own = owner[facei];
        const label nbr = neighbour[facei];

        // Check if owner value much larger than neighbour value or vice versa
        if (mag(alpha[own] - alpha[nbr]) > alphaDiff)
        {
            changedPatchAndFaces.append(labelPair(-1, facei));
            changedFacesInfo.append(smoothData(max(field[own], field[nbr])));
        }
    }

    // Insert all faces of coupled patches - FvFaceCellWave will correct them
    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& patch = mesh.boundary()[patchi];

        if (patch.coupled())
        {
            forAll(patch, patchFacei)
            {
                const label own = patch.faceCells()[patchFacei];

                const scalarField alphapn
                (
                    alpha.boundaryField()[patchi].patchNeighbourField()
                );

                if (mag(alpha[own] - alphapn[patchFacei]) > alphaDiff)
                {
                    changedPatchAndFaces.append(labelPair(patchi, patchFacei));
                    changedFacesInfo.append(smoothData(field[own]));
                }
            }
        }
    }

    changedPatchAndFaces.shrink();
    changedFacesInfo.shrink();

    // Set initial field on cells
    List<smoothData> cellData(mesh.nCells());
    forAll(field, celli)
    {
        cellData[celli] = field[celli];
    }

    // Set initial field on faces
    List<smoothData> internalFaceData(mesh.nInternalFaces());
    List<List<smoothData>> patchFaceData
    (
        FvFaceCellWave<smoothData, smoothData::trackData>::template
        sizesListList<List<List<smoothData>>>
        (
            FvFaceCellWave<smoothData, smoothData::trackData>::template
            listListSizes(mesh.boundary()),
            smoothData()
        )
    );

    // Create track data
    smoothData::trackData td;
    td.maxRatio = 1.0;

    // Propagate information over whole domain
    FvFaceCellWave<smoothData, smoothData::trackData> smoothData
    (
        mesh,
        internalFaceData,
        patchFaceData,
        cellData,
        td
    );

    smoothData.setFaceInfo(changedPatchAndFaces, changedFacesInfo);

    smoothData.iterate(nLayers);

    forAll(field, celli)
    {
        field[celli] = cellData[celli].value();
    }

    field.correctBoundaryConditions();
}


void Foam::fvc::sweep
(
    volScalarField& field,
    const volScalarField& alpha,
    const label nLayers,
    const scalar alphaDiff
)
{
    const fvMesh& mesh = field.mesh();

    DynamicList<labelPair> changedPatchAndFaces(mesh.nFaces()/100 + 100);
    DynamicList<sweepData> changedFacesInfo(changedPatchAndFaces.size());

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const surfaceVectorField& Cf = mesh.Cf();
    const surfaceVectorField::Boundary& CfBf = mesh.Cf().boundaryField();

    forAll(owner, facei)
    {
        const label own = owner[facei];
        const label nbr = neighbour[facei];

        // Check if owner value much larger than neighbour value or vice versa
        if (mag(alpha[own] - alpha[nbr]) > alphaDiff)
        {
            changedPatchAndFaces.append(labelPair(-1, facei));
            changedFacesInfo.append
            (
                sweepData(max(field[own], field[nbr]), Cf[facei])
            );
        }
    }

    // Insert all faces of coupled patches - FvFaceCellWave will correct them
    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& patch = mesh.boundary()[patchi];

        if (patch.coupled())
        {
            forAll(patch, patchFacei)
            {
                const label own = patch.faceCells()[patchFacei];

                const scalarField alphapn
                (
                    alpha.boundaryField()[patchi].patchNeighbourField()
                );

                if (mag(alpha[own] - alphapn[patchFacei]) > alphaDiff)
                {
                    changedPatchAndFaces.append(labelPair(patchi, patchFacei));
                    changedFacesInfo.append
                    (
                        sweepData(field[own], CfBf[patchi][patchFacei])
                    );
                }
            }
        }
    }

    changedPatchAndFaces.shrink();
    changedFacesInfo.shrink();

    // Set initial field on cells
    List<sweepData> cellData(mesh.nCells());

    // Set initial field on faces
    List<sweepData> internalFaceData(mesh.nInternalFaces());
    List<List<sweepData>> patchFaceData
    (
        FvFaceCellWave<sweepData>::template
        sizesListList<List<List<sweepData>>>
        (
            FvFaceCellWave<sweepData>::template
            listListSizes(mesh.boundary()),
            sweepData()
        )
    );

    // Propagate information over whole domain
    FvFaceCellWave<sweepData> sweepData
    (
        mesh,
        internalFaceData,
        patchFaceData,
        cellData
    );

    sweepData.setFaceInfo(changedPatchAndFaces, changedFacesInfo);

    sweepData.iterate(nLayers);

    forAll(field, celli)
    {
        if (cellData[celli].valid(sweepData.data()))
        {
            field[celli] = max(field[celli], cellData[celli].value());
        }
    }

    field.correctBoundaryConditions();
}


// ************************************************************************* //
