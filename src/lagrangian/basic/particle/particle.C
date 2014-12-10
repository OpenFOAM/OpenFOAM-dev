/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "particle.H"
#include "transform.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::label Foam::particle::particleCount_ = 0;

const Foam::scalar Foam::particle::trackingCorrectionTol = 1e-5;

const Foam::scalar Foam::particle::lambdaDistanceToleranceCoeff = 1e3*SMALL;

const Foam::scalar Foam::particle::minStepFractionTol = 1e5*SMALL;

namespace Foam
{
    defineTypeNameAndDebug(particle, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::particle::particle
(
    const polyMesh& mesh,
    const vector& position,
    const label cellI,
    const label tetFaceI,
    const label tetPtI
)
:
    mesh_(mesh),
    position_(position),
    cellI_(cellI),
    faceI_(-1),
    stepFraction_(0.0),
    tetFaceI_(tetFaceI),
    tetPtI_(tetPtI),
    origProc_(Pstream::myProcNo()),
    origId_(getNewParticleID())
{}


Foam::particle::particle
(
    const polyMesh& mesh,
    const vector& position,
    const label cellI,
    bool doCellFacePt
)
:
    mesh_(mesh),
    position_(position),
    cellI_(cellI),
    faceI_(-1),
    stepFraction_(0.0),
    tetFaceI_(-1),
    tetPtI_(-1),
    origProc_(Pstream::myProcNo()),
    origId_(getNewParticleID())
{
    if (doCellFacePt)
    {
        initCellFacePt();
    }
}


Foam::particle::particle(const particle& p)
:
    mesh_(p.mesh_),
    position_(p.position_),
    cellI_(p.cellI_),
    faceI_(p.faceI_),
    stepFraction_(p.stepFraction_),
    tetFaceI_(p.tetFaceI_),
    tetPtI_(p.tetPtI_),
    origProc_(p.origProc_),
    origId_(p.origId_)
{}


Foam::particle::particle(const particle& p, const polyMesh& mesh)
:
    mesh_(mesh),
    position_(p.position_),
    cellI_(p.cellI_),
    faceI_(p.faceI_),
    stepFraction_(p.stepFraction_),
    tetFaceI_(p.tetFaceI_),
    tetPtI_(p.tetPtI_),
    origProc_(p.origProc_),
    origId_(p.origId_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::particle::transformProperties(const tensor&)
{}


void Foam::particle::transformProperties(const vector&)
{}


Foam::scalar Foam::particle::wallImpactDistance(const vector&) const
{
    Info<< "particle::wallImpactDistance" << endl;

    return 0.0;
}


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //

bool Foam::operator==(const particle& pA, const particle& pB)
{
    return (pA.origProc() == pB.origProc() && pA.origId() == pB.origId());
}


bool Foam::operator!=(const particle& pA, const particle& pB)
{
    return !(pA == pB);
}


// ************************************************************************* //
