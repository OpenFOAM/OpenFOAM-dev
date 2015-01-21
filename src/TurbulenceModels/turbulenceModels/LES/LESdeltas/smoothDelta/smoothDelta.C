/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "smoothDelta.H"
#include "addToRunTimeSelectionTable.H"
#include "FaceCellWave.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{
    defineTypeNameAndDebug(smoothDelta, 0);
    addToRunTimeSelectionTable(LESdelta, smoothDelta, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::LESModels::smoothDelta::setChangedFaces
(
    const polyMesh& mesh,
    const volScalarField& delta,
    DynamicList<label>& changedFaces,
    DynamicList<deltaData>& changedFacesInfo
)
{
    for (label faceI = 0; faceI < mesh.nInternalFaces(); faceI++)
    {
        scalar ownDelta = delta[mesh.faceOwner()[faceI]];

        scalar neiDelta = delta[mesh.faceNeighbour()[faceI]];

        // Check if owner delta much larger than neighbour delta or vice versa

        if (ownDelta > maxDeltaRatio_ * neiDelta)
        {
            changedFaces.append(faceI);
            changedFacesInfo.append(deltaData(ownDelta));
        }
        else if (neiDelta > maxDeltaRatio_ * ownDelta)
        {
            changedFaces.append(faceI);
            changedFacesInfo.append(deltaData(neiDelta));
        }
    }

    // Insert all faces of coupled patches no matter what. Let FaceCellWave
    // sort it out.
    forAll(mesh.boundaryMesh(), patchI)
    {
        const polyPatch& patch = mesh.boundaryMesh()[patchI];

        if (patch.coupled())
        {
            forAll(patch, patchFaceI)
            {
                label meshFaceI = patch.start() + patchFaceI;

                scalar ownDelta = delta[mesh.faceOwner()[meshFaceI]];

                changedFaces.append(meshFaceI);
                changedFacesInfo.append(deltaData(ownDelta));
            }
        }
    }

    changedFaces.shrink();
    changedFacesInfo.shrink();
}


void Foam::LESModels::smoothDelta::calcDelta()
{
    const fvMesh& mesh = turbulenceModel_.mesh();

    const volScalarField& geometricDelta = geometricDelta_();

    // Fill changed faces with info
    DynamicList<label> changedFaces(mesh.nFaces()/100 + 100);
    DynamicList<deltaData> changedFacesInfo(changedFaces.size());

    setChangedFaces(mesh, geometricDelta, changedFaces, changedFacesInfo);

    // Set initial field on cells.
    List<deltaData> cellDeltaData(mesh.nCells());

    forAll(geometricDelta, cellI)
    {
        cellDeltaData[cellI] = geometricDelta[cellI];
    }

    // Set initial field on faces.
    List<deltaData> faceDeltaData(mesh.nFaces());


    // Propagate information over whole domain.
    FaceCellWave<deltaData, scalar> deltaCalc
    (
        mesh,
        changedFaces,
        changedFacesInfo,
        faceDeltaData,
        cellDeltaData,
        mesh.globalData().nTotalCells()+1,  // max iterations
        maxDeltaRatio_
    );

    forAll(delta_, cellI)
    {
        delta_[cellI] = cellDeltaData[cellI].delta();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LESModels::smoothDelta::smoothDelta
(
    const word& name,
    const turbulenceModel& turbulence,
    const dictionary& dict
)
:
    LESdelta(name, turbulence),
    geometricDelta_
    (
        LESdelta::New
        (
            "geometricDelta",
            turbulence,
            dict.subDict(type() + "Coeffs")
        )
    ),
    maxDeltaRatio_
    (
        readScalar(dict.subDict(type() + "Coeffs").lookup("maxDeltaRatio"))
    )
{
    calcDelta();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::LESModels::smoothDelta::read(const dictionary& dict)
{
    const dictionary& coeffsDict(dict.subDict(type() + "Coeffs"));

    geometricDelta_().read(coeffsDict);
    coeffsDict.lookup("maxDeltaRatio") >> maxDeltaRatio_;
    calcDelta();
}


void Foam::LESModels::smoothDelta::correct()
{
    geometricDelta_().correct();

    if (turbulenceModel_.mesh().changing())
    {
        calcDelta();
    }
}


// ************************************************************************* //
