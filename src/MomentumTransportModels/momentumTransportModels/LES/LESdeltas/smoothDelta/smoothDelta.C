/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    for (label facei = 0; facei < mesh.nInternalFaces(); facei++)
    {
        scalar ownDelta = delta[mesh.faceOwner()[facei]];

        scalar neiDelta = delta[mesh.faceNeighbour()[facei]];

        // Check if owner delta much larger than neighbour delta or vice versa

        if (ownDelta > maxDeltaRatio_ * neiDelta)
        {
            changedFaces.append(facei);
            changedFacesInfo.append(deltaData(ownDelta));
        }
        else if (neiDelta > maxDeltaRatio_ * ownDelta)
        {
            changedFaces.append(facei);
            changedFacesInfo.append(deltaData(neiDelta));
        }
    }

    // Insert all faces of coupled patches no matter what. Let FaceCellWave
    // sort it out.
    forAll(mesh.boundaryMesh(), patchi)
    {
        const polyPatch& patch = mesh.boundaryMesh()[patchi];

        if (patch.coupled())
        {
            forAll(patch, patchFacei)
            {
                label meshFacei = patch.start() + patchFacei;

                scalar ownDelta = delta[mesh.faceOwner()[meshFacei]];

                changedFaces.append(meshFacei);
                changedFacesInfo.append(deltaData(ownDelta));
            }
        }
    }

    changedFaces.shrink();
    changedFacesInfo.shrink();
}


void Foam::LESModels::smoothDelta::calcDelta()
{
    const fvMesh& mesh = momentumTransportModel_.mesh();

    const volScalarField& geometricDelta = geometricDelta_();

    // Fill changed faces with info
    DynamicList<label> changedFaces(mesh.nFaces()/100 + 100);
    DynamicList<deltaData> changedFacesInfo(changedFaces.size());

    setChangedFaces(mesh, geometricDelta, changedFaces, changedFacesInfo);

    // Set initial field on cells.
    List<deltaData> cellDeltaData(mesh.nCells());

    forAll(geometricDelta, celli)
    {
        cellDeltaData[celli] = geometricDelta[celli];
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

    forAll(delta_, celli)
    {
        delta_[celli] = cellDeltaData[celli].delta();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LESModels::smoothDelta::smoothDelta
(
    const word& name,
    const momentumTransportModel& turbulence,
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
            dict.optionalSubDict(type() + "Coeffs")
        )
    ),
    maxDeltaRatio_
    (
        dict.optionalSubDict(type() + "Coeffs").lookup<scalar>("maxDeltaRatio")
    )
{
    calcDelta();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::LESModels::smoothDelta::read(const dictionary& dict)
{
    const dictionary& coeffsDict(dict.optionalSubDict(type() + "Coeffs"));

    geometricDelta_().read(coeffsDict);
    coeffsDict.lookup("maxDeltaRatio") >> maxDeltaRatio_;
    calcDelta();
}


void Foam::LESModels::smoothDelta::correct()
{
    geometricDelta_().correct();

    if (momentumTransportModel_.mesh().changing())
    {
        calcDelta();
    }
}


// ************************************************************************* //
