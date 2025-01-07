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

#include "cloudBoundaryCollisionFlux.H"
#include "cloud.H"
#include "CompactListList.H"
#include "addToRunTimeSelectionTable.H"
#include "functionName.H"
#include "functionObject.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(cloudBoundaryCollisionFlux, 0);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::cloudBoundaryCollisionFlux::cloudBoundaryCollisionFlux
(
    const word& name,
    const Time& runTime,
    const dictionary& dict,
    const word& phiName,
    const dimensionSet& phiDims
)
:
    cloudFvMeshFunctionObject(name, runTime, dict),
    faceiPtr_(nullptr),
    qPtr_(nullptr),
    phiName_(phiName),
    phiDims_(phiDims),
    phib_
    (
        mesh().boundary(),
        surfaceScalarField::Internal::null(),
        calculatedFvsPatchField<scalar>::typeName
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::cloudBoundaryCollisionFlux::~cloudBoundaryCollisionFlux()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::functionObjects::cloudBoundaryCollisionFlux::fields() const
{
    return wordList::null();
}


bool Foam::functionObjects::cloudBoundaryCollisionFlux::executeAtStart() const
{
    return false;
}


bool Foam::functionObjects::cloudBoundaryCollisionFlux::execute()
{
    return true;
}


void Foam::functionObjects::cloudBoundaryCollisionFlux::preSolve()
{
    phib_ == Zero;
}


void Foam::functionObjects::cloudBoundaryCollisionFlux::preCrossFaces
(
    const LagrangianScalarInternalDynamicField& fraction
)
{
    faceiPtr_.set
    (
        new LagrangianLabelField
        (
            IOobject
            (
                name() + ":facei",
                time_.name(),
                cloud().mesh()
            ),
            cloud().mesh(),
            dimensioned<label>(dimless, -1)
        )
    );

    qPtr_.set
    (
        new LagrangianScalarField
        (
            IOobject
            (
                name() + ":q",
                time_.name(),
                cloud().mesh()
            ),
            cloud().mesh(),
            dimensionedScalar(phiDims_*dimTime, scalar(0))
        )
    );
}


void Foam::functionObjects::cloudBoundaryCollisionFlux::preCrossFaces
(
    const LagrangianSubScalarSubField& fraction
)
{
    const LagrangianSubMesh& subMesh = fraction.mesh();

    if
    (
        faceiPtr_.empty()
     || subMesh.group() == LagrangianGroup::none
     || subMesh.group() == LagrangianGroup::inInternalMesh
    ) return;

    const Foam::cloud& c = cloud();

    SubField<label> facei(subMesh.sub(faceiPtr_->primitiveFieldRef()));
    facei = subMesh.sub(c.mesh().facei());

    LagrangianSubScalarSubField q(subMesh.sub(qPtr_()));
    q += this->q(fraction, +1);
}


void Foam::functionObjects::cloudBoundaryCollisionFlux::postCrossFaces
(
    const LagrangianSubScalarSubField& fraction
)
{
    const LagrangianSubMesh& subMesh = fraction.mesh();

    if
    (
        faceiPtr_.empty()
     || subMesh.group() == LagrangianGroup::none
     || subMesh.group() == LagrangianGroup::inInternalMesh
    ) return;

    const Foam::cloud& c = cloud();

    SubField<label> facei0(subMesh.sub(faceiPtr_->primitiveFieldRef()));
    SubField<label> facei = subMesh.sub(c.mesh().facei());

    LagrangianSubScalarSubField q(subMesh.sub(qPtr_()));
    q += this->q(fraction, -1);

    const surfaceScalarField::Boundary& magSfb = mesh().magSf().boundaryField();

    forAll(subMesh, subi)
    {
        if (facei0[subi] != facei[subi]) continue;

        const LagrangianState state = c.mesh().state(subi + subMesh.start());

        if (state != LagrangianState::inCell) continue;

        const label bFacei = facei[subi] - mesh().nInternalFaces();

        const labelUList patchis = mesh().polyBFacePatches()[bFacei];
        const labelUList patchFaceis = mesh().polyBFacePatchFaces()[bFacei];

        forAll(patchis, i)
        {
            phib_[patchis[i]][patchFaceis[i]] +=
                q[subi]
               /time_.deltaTValue()
               *magSfb[patchis[i]][patchFaceis[i]]
               /mesh().magFaceAreas()[facei[subi]];
        }
    }
}


void Foam::functionObjects::cloudBoundaryCollisionFlux::postCrossFaces
(
    const LagrangianScalarInternalDynamicField& fraction
)
{
    faceiPtr_.clear();
    qPtr_.clear();
}


bool Foam::functionObjects::cloudBoundaryCollisionFlux::write()
{
    return
        surfaceScalarField
        (
            IOobject
            (
                cloud().mesh().name() + ":" + phiName_ + "Coll",
                time_.name(),
                mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh(),
            phiDims_,
            scalarField(mesh().nInternalFaces(), scalar(0)),
            phib_
        ).write();
}


bool Foam::functionObjects::cloudBoundaryCollisionFlux::clear()
{
    return true;
}


// ************************************************************************* //
