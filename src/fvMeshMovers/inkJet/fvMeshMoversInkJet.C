/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2022 OpenFOAM Foundation
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

#include "fvMeshMoversInkJet.H"
#include "volFields.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshMovers
{
    defineTypeNameAndDebug(inkJet, 0);
    addToRunTimeSelectionTable(fvMeshMover, inkJet, fvMesh);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshMovers::inkJet::inkJet(fvMesh& mesh)
:
    fvMeshMover(mesh),
    meshCoeffs_(dict()),
    amplitude_(meshCoeffs_.lookup<scalar>("amplitude")),
    frequency_(meshCoeffs_.lookup<scalar>("frequency")),
    refPlaneX_(meshCoeffs_.lookup<scalar>("refPlaneX")),
    stationaryPoints_
    (
        IOobject
        (
            "points",
            mesh.time().constant(),
            fvMesh::meshSubDir,
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    velocityMotionCorrection_(mesh, dict())
{
    Info<< "Performing a dynamic mesh calculation: " << endl
        << "amplitude: " << amplitude_
        << " frequency: " << frequency_
        << " refPlaneX: " << refPlaneX_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMeshMovers::inkJet::~inkJet()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fvMeshMovers::inkJet::update()
{
    const scalar scalingFunction =
        0.5*
        (
            cos(constant::mathematical::twoPi*frequency_*mesh().time().value())
          - 1.0
        );

    Info<< "Mesh scaling. Time = " << mesh().time().value() << " scaling: "
        << scalingFunction << endl;

    pointField newPoints = stationaryPoints_;

    newPoints.replace
    (
        vector::X,
        stationaryPoints_.component(vector::X)*
        (
            1.0
          + pos0
            (
              - (stationaryPoints_.component(vector::X))
              - refPlaneX_
            )*amplitude_*scalingFunction
        )
    );

    mesh().movePoints(newPoints);

    velocityMotionCorrection_.update();

    return true;
}


void Foam::fvMeshMovers::inkJet::topoChange(const polyTopoChangeMap&)
{
    NotImplemented;
}


void Foam::fvMeshMovers::inkJet::mapMesh(const polyMeshMap&)
{
    NotImplemented;
}


void Foam::fvMeshMovers::inkJet::distribute
(
    const polyDistributionMap&
)
{
    NotImplemented;
}


// ************************************************************************* //
