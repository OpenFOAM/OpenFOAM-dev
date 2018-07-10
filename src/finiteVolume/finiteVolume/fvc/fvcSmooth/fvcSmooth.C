/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
#include "FaceCellWave.H"
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

    DynamicList<label> changedFaces(mesh.nFaces()/100 + 100);
    DynamicList<smoothData> changedFacesInfo(changedFaces.size());

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    forAll(owner, facei)
    {
        const label own = owner[facei];
        const label nbr = neighbour[facei];

        // Check if owner value much larger than neighbour value or vice versa
        if (field[own] > maxRatio*field[nbr])
        {
            changedFaces.append(facei);
            changedFacesInfo.append(smoothData(field[own]));
        }
        else if (field[nbr] > maxRatio*field[own])
        {
            changedFaces.append(facei);
            changedFacesInfo.append(smoothData(field[nbr]));
        }
    }

    // Insert all faces of coupled patches - FaceCellWave will correct them
    forAll(mesh.boundaryMesh(), patchi)
    {
        const polyPatch& patch = mesh.boundaryMesh()[patchi];

        if (patch.coupled())
        {
            forAll(patch, patchFacei)
            {
                label facei = patch.start() + patchFacei;
                label own = mesh.faceOwner()[facei];

                changedFaces.append(facei);
                changedFacesInfo.append(smoothData(field[own]));
            }
        }
    }

    changedFaces.shrink();
    changedFacesInfo.shrink();

    // Set initial field on cells
    List<smoothData> cellData(mesh.nCells());

    forAll(field, celli)
    {
        cellData[celli] = field[celli];
    }

    // Set initial field on faces
    List<smoothData> faceData(mesh.nFaces());

    smoothData::trackData td;
    td.maxRatio = maxRatio;

    // Propagate information over whole domain
    FaceCellWave<smoothData, smoothData::trackData> smoothData
    (
        mesh,
        changedFaces,
        changedFacesInfo,
        faceData,
        cellData,
        mesh.globalData().nTotalCells(),   // max iterations
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

    DynamicList<label> changedFaces(mesh.nFaces()/100 + 100);
    DynamicList<smoothData> changedFacesInfo(changedFaces.size());

    // Set initial field on cells
    List<smoothData> cellData(mesh.nCells());

    forAll(field, celli)
    {
        cellData[celli] = field[celli];
    }

    // Set initial field on faces
    List<smoothData> faceData(mesh.nFaces());

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    forAll(owner, facei)
    {
        const label own = owner[facei];
        const label nbr = neighbour[facei];

        if (mag(alpha[own] - alpha[nbr]) > alphaDiff)
        {
            changedFaces.append(facei);
            changedFacesInfo.append
            (
                smoothData(max(field[own], field[nbr]))
            );
        }
    }

    // Insert all faces of coupled patches - FaceCellWave will correct them
    forAll(mesh.boundaryMesh(), patchi)
    {
        const polyPatch& patch = mesh.boundaryMesh()[patchi];

        if (patch.coupled())
        {
            forAll(patch, patchFacei)
            {
                label facei = patch.start() + patchFacei;
                label own = mesh.faceOwner()[facei];

                const scalarField alphapn
                (
                    alpha.boundaryField()[patchi].patchNeighbourField()
                );

                if (mag(alpha[own] - alphapn[patchFacei]) > alphaDiff)
                {
                    changedFaces.append(facei);
                    changedFacesInfo.append(smoothData(field[own]));
                }
            }
        }
    }

    changedFaces.shrink();
    changedFacesInfo.shrink();

    smoothData::trackData td;
    td.maxRatio = 1.0;

    // Propagate information over whole domain
    FaceCellWave<smoothData, smoothData::trackData> smoothData
    (
        mesh,
        faceData,
        cellData,
        td
    );

    smoothData.setFaceInfo(changedFaces, changedFacesInfo);

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

    DynamicList<label> changedFaces(mesh.nFaces()/100 + 100);
    DynamicList<sweepData> changedFacesInfo(changedFaces.size());

    // Set initial field on cells
    List<sweepData> cellData(mesh.nCells());

    // Set initial field on faces
    List<sweepData> faceData(mesh.nFaces());

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const vectorField& Cf = mesh.faceCentres();

    forAll(owner, facei)
    {
        const label own = owner[facei];
        const label nbr = neighbour[facei];

        if (mag(alpha[own] - alpha[nbr]) > alphaDiff)
        {
            changedFaces.append(facei);
            changedFacesInfo.append
            (
                sweepData(max(field[own], field[nbr]), Cf[facei])
            );
        }
    }

    // Insert all faces of coupled patches - FaceCellWave will correct them
    forAll(mesh.boundaryMesh(), patchi)
    {
        const polyPatch& patch = mesh.boundaryMesh()[patchi];

        if (patch.coupled())
        {
            forAll(patch, patchFacei)
            {
                label facei = patch.start() + patchFacei;
                label own = mesh.faceOwner()[facei];

                const scalarField alphapn
                (
                    alpha.boundaryField()[patchi].patchNeighbourField()
                );

                if (mag(alpha[own] - alphapn[patchFacei]) > alphaDiff)
                {
                    changedFaces.append(facei);
                    changedFacesInfo.append
                    (
                        sweepData(field[own], Cf[facei])
                    );
                }
            }
        }
    }

    changedFaces.shrink();
    changedFacesInfo.shrink();

    // Propagate information over whole domain
    FaceCellWave<sweepData> sweepData
    (
        mesh,
        faceData,
        cellData
    );

    sweepData.setFaceInfo(changedFaces, changedFacesInfo);

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
