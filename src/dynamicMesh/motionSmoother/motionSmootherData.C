/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#include "motionSmootherData.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::motionSmootherData::motionSmootherData
(
    const pointMesh& pMesh
)
:
    displacement_
    (
        IOobject
        (
            "displacement",
            pMesh.time().timeName(),
            pMesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        pMesh
    ),
    scale_
    (
        IOobject
        (
            "scale",
            pMesh.time().timeName(),
            pMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        pMesh,
        dimensionedScalar("scale", dimless, 1.0)
    ),
    oldPoints_(pMesh().points())
{}


Foam::motionSmootherData::motionSmootherData
(
    const pointVectorField& displacement
)
:
    displacement_
    (
        IOobject
        (
            "displacement",
            displacement.time().timeName(),
            displacement.mesh()(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        displacement
    ),
    scale_
    (
        IOobject
        (
            "scale",
            displacement.time().timeName(),
            displacement.mesh()(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        displacement.mesh(),
        dimensionedScalar("scale", dimless, 1.0)
    ),
    oldPoints_(displacement.mesh()().points())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::pointVectorField& Foam::motionSmootherData::displacement()
{
    return displacement_;
}


const Foam::pointVectorField& Foam::motionSmootherData::displacement() const
{
    return displacement_;
}


const Foam::pointScalarField& Foam::motionSmootherData::scale() const
{
    return scale_;
}


const Foam::pointField& Foam::motionSmootherData::oldPoints() const
{
    return oldPoints_;
}


// ************************************************************************* //
