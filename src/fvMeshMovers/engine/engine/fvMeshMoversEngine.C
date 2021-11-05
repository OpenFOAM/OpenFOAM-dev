/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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

#include "fvMeshMoversEngine.H"
#include "engineTime.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshMovers
{
    defineTypeNameAndDebug(engine, 0);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshMovers::engine::engine(fvMesh& mesh)
:
    fvMeshMover(mesh),
    meshCoeffs_(dict()),
    rpm_
    (
        refCast<const userTimes::engine>(mesh.time().userTime()).rpm()
    ),
    conRodLength_("conRodLength", dimLength, meshCoeffs_),
    bore_("bore", dimLength, meshCoeffs_),
    stroke_("stroke", dimLength, meshCoeffs_),
    clearance_("clearance", dimLength, meshCoeffs_),
    pistonIndex_(-1),
    linerIndex_(-1),
    cylinderHeadIndex_(-1),
    deckHeight_("deckHeight", dimLength, great),
    pistonPosition_("pistonPosition", dimLength, -great)
{
    bool foundPiston = false;
    bool foundLiner = false;
    bool foundCylinderHead = false;

    forAll(mesh.boundary(), i)
    {
        if (mesh.boundary()[i].name() == "piston")
        {
            pistonIndex_ = i;
            foundPiston = true;
        }
        else if (mesh.boundary()[i].name() == "liner")
        {
            linerIndex_ = i;
            foundLiner = true;
        }
        else if (mesh.boundary()[i].name() == "cylinderHead")
        {
            cylinderHeadIndex_ = i;
            foundCylinderHead = true;
        }
    }

    reduce(foundPiston, orOp<bool>());
    reduce(foundLiner, orOp<bool>());
    reduce(foundCylinderHead, orOp<bool>());

    if (!foundPiston)
    {
        FatalErrorInFunction
            << "cannot find piston patch"
            << exit(FatalError);
    }

    if (!foundLiner)
    {
        FatalErrorInFunction
            << "cannot find liner patch"
            << exit(FatalError);
    }

    if (!foundCylinderHead)
    {
        FatalErrorInFunction
            << "cannot find cylinderHead patch"
            << exit(FatalError);
    }

    {
        if (pistonIndex_ != -1)
        {
            pistonPosition_.value() = -great;
            if (mesh.boundary()[pistonIndex_].patch().localPoints().size())
            {
                pistonPosition_.value() =
                    max(mesh.boundary()[pistonIndex_].patch().localPoints())
                   .z();
            }
        }
        reduce(pistonPosition_.value(), maxOp<scalar>());

        if (cylinderHeadIndex_ != -1)
        {
            deckHeight_.value() = great;
            if
            (
                mesh.boundary()[cylinderHeadIndex_].patch().localPoints().size()
            )
            {
                deckHeight_.value() = min
                (
                    mesh.boundary()[cylinderHeadIndex_].patch().localPoints()
                ).z();
            }
        }
        reduce(deckHeight_.value(), minOp<scalar>());

        Info<< "deckHeight: " << deckHeight_.value() << nl
            << "piston position: " << pistonPosition_.value() << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMeshMovers::engine::~engine()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::fvMeshMovers::engine::theta() const
{
    return mesh().time().userTimeValue();
}


Foam::scalar Foam::fvMeshMovers::engine::deltaTheta() const
{
    return mesh().time().timeToUserTime(mesh().time().deltaTValue());
}


Foam::scalar Foam::fvMeshMovers::engine::pistonPosition
(
    const scalar theta
) const
{
    return
    (
        conRodLength_.value()
      + stroke_.value()/2.0
      + clearance_.value()
    )
  - (
        stroke_.value()*::cos(degToRad(theta))/2.0
      + ::sqrt
        (
            sqr(conRodLength_.value())
          - sqr(stroke_.value()*::sin(degToRad(theta))/2.0)
        )
    );
}


Foam::dimensionedScalar Foam::fvMeshMovers::engine::pistonPosition() const
{
    return dimensionedScalar
    (
        "pistonPosition",
        dimLength,
        pistonPosition(theta())
    );
}


Foam::dimensionedScalar Foam::fvMeshMovers::engine::pistonDisplacement() const
{
    return dimensionedScalar
    (
        "pistonDisplacement",
        dimLength,
        pistonPosition(theta() - deltaTheta()) - pistonPosition().value()
    );
}


Foam::dimensionedScalar Foam::fvMeshMovers::engine::pistonSpeed() const
{
    return dimensionedScalar
    (
        "pistonSpeed",
        dimVelocity,
        pistonDisplacement().value()/(mesh().time().deltaTValue() + vSmall)
    );
}


// ************************************************************************* //
