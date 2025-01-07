/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "cloudFlux.H"
#include "cloud.H"
#include "CompactListList.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(cloudFlux, 0);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::cloudFlux::crossPatchFaces
(
    const LagrangianSubScalarSubField& fraction,
    const label sign
)
{
    const Foam::cloud& c = cloud();

    SubField<label> facei = fraction.mesh().sub(c.mesh().facei());

    const LagrangianSubScalarField dqdt(q(fraction)/time_.deltaT());

    const surfaceScalarField::Boundary& magSfb = mesh().magSf().boundaryField();

    forAll(fraction, subi)
    {
        const LagrangianState state =
            c.mesh().state(subi + fraction.mesh().start());

        if (sign < 0 && state == LagrangianState::toBeRemoved) continue;

        const label bFacei = facei[subi] - mesh().nInternalFaces();

        const labelUList patchis = mesh().polyBFacePatches()[bFacei];
        const labelUList patchFaceis = mesh().polyBFacePatchFaces()[bFacei];

        forAll(patchis, i)
        {
            phi_.boundaryFieldRef()[patchis[i]][patchFaceis[i]] +=
                sign
               *dqdt[subi]
               *magSfb[patchis[i]][patchFaceis[i]]
               /mesh().magFaceAreas()[facei[subi]];
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::cloudFlux::cloudFlux
(
    const word& name,
    const Time& runTime,
    const dictionary& dict,
    const word& qName,
    const dimensionSet& qDims
)
:
    cloudFvMeshFunctionObject(name, runTime, dict),
    phi_
    (
        IOobject
        (
            cloud().mesh().name()
          + ":phi"
          + (qName.size() == 1 ? qName : qName.capitalise()),
            runTime.name(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar(qDims/dimTime, scalar(0))
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::cloudFlux::~cloudFlux()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::functionObjects::cloudFlux::fields() const
{
    return wordList::null();
}


bool Foam::functionObjects::cloudFlux::executeAtStart() const
{
    return false;
}


void Foam::functionObjects::cloudFlux::preSolve()
{
    phi_ == Zero;
}


bool Foam::functionObjects::cloudFlux::execute()
{
    return true;
}


void Foam::functionObjects::cloudFlux::preCrossFaces
(
    const LagrangianSubScalarSubField& fraction
)
{
    if (fraction.mesh().group() != LagrangianGroup::inInternalMesh)
    {
        crossPatchFaces(fraction, +1);
    }
    else
    {
        const Foam::cloud& c = cloud();

        const labelList& owner = mesh().owner();

        SubField<label> celli = fraction.mesh().sub(c.mesh().celli());
        SubField<label> facei = fraction.mesh().sub(c.mesh().facei());

        const LagrangianSubScalarField dqdt(q(fraction)/time_.deltaT());

        forAll(fraction, subi)
        {
            phi_.internalFieldRef()[facei[subi]] +=
                (owner[facei[subi]] == celli[subi] ? +1 : -1)
               *dqdt[subi]
               /mesh().time().deltaTValue();
        }
    }
}


void Foam::functionObjects::cloudFlux::postCrossFaces
(
    const LagrangianSubScalarSubField& fraction
)
{
    if (fraction.mesh().group() != LagrangianGroup::inInternalMesh)
    {
        crossPatchFaces(fraction, -1);
    }
}


bool Foam::functionObjects::cloudFlux::write()
{
    return true;
}


bool Foam::functionObjects::cloudFlux::clear()
{
    return true;
}


// ************************************************************************* //
