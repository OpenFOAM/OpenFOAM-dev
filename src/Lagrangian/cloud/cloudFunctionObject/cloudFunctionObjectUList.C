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

#include "cloudFunctionObjectUList.H"
#include "cloud.H"
#include "timeControlFunctionObject.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cloudFunctionObjectUList, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cloudFunctionObjectUList::cloudFunctionObjectUList
(
    const cloud& c,
    const bool inner
)
:
    UPtrList<functionObjects::cloudFunctionObject>
    (
        c.time().functionObjects().size()
    ),
    inner_(inner)
{
    const functionObjectList& functions = c.time().functionObjects();

    label cloudFunctioni = 0;

    forAll(functions, functioni)
    {
        const functionObject& fo =
            isA<functionObjects::timeControl>(functions[functioni])
          ? refCast<const functionObjects::timeControl>
            (
                functions[functioni]
            ).filter()
          : functions[functioni];

        if (!isA<functionObjects::cloudFunctionObject>(fo)) continue;

        const functionObjects::cloudFunctionObject& cfo =
            refCast<const functionObjects::cloudFunctionObject>(fo);

        if (&cfo.cloud() != &c) continue;

        UPtrList<functionObjects::cloudFunctionObject>::set
        (
            cloudFunctioni ++,
            const_cast<functionObjects::cloudFunctionObject*>(&cfo)
        );
    }

    UPtrList<functionObjects::cloudFunctionObject>::resize(cloudFunctioni);

    if (inner_) return;

    forAll(*this, cloudFunctioni)
    {
        this->operator[](cloudFunctioni).preSolve();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cloudFunctionObjectUList::~cloudFunctionObjectUList()
{
    if (inner_) return;

    forAll(*this, cloudFunctioni)
    {
        this->operator[](cloudFunctioni).postSolve();
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::cloudFunctionObjectUList::calculate
(
    const LagrangianSubScalarField& deltaT,
    const bool final
)
{
    forAll(*this, cloudFunctioni)
    {
        this->operator[](cloudFunctioni).calculate(deltaT, final);
    }
}


void Foam::cloudFunctionObjectUList::preCrossFaces
(
    const LagrangianScalarInternalDynamicField& fraction
)
{
    forAll(*this, cloudFunctioni)
    {
        this->operator[](cloudFunctioni).preCrossFaces(fraction);
    }

    if (size())
    {
        const Foam::cloud& cloud = first().cloud();

        preCrossFaces
        (
            cloud.mesh().sub(LagrangianGroup::inInternalMesh).sub(fraction)
        );

        forAll(cloud.mesh().boundary(), patchi)
        {
            preCrossFaces
            (
                cloud.mesh().boundary()[patchi].mesh().sub(fraction)
            );
        }
    }
}


void Foam::cloudFunctionObjectUList::preCrossFaces
(
    const LagrangianSubScalarSubField& fraction
)
{
    forAll(*this, cloudFunctioni)
    {
        this->operator[](cloudFunctioni).preCrossFaces(fraction);
    }
}


void Foam::cloudFunctionObjectUList::postCrossFaces
(
    const LagrangianSubScalarSubField& fraction
)
{
    forAll(*this, cloudFunctioni)
    {
        this->operator[](cloudFunctioni).postCrossFaces(fraction);
    }
}


void Foam::cloudFunctionObjectUList::postCrossFaces
(
    const LagrangianScalarInternalDynamicField& fraction
)
{
    if (size())
    {
        const Foam::cloud& cloud = first().cloud();

        postCrossFaces
        (
            cloud.mesh().sub(LagrangianGroup::inInternalMesh).sub(fraction)
        );

        forAll(cloud.mesh().boundary(), patchi)
        {
            postCrossFaces
            (
                cloud.mesh().boundary()[patchi].mesh().sub(fraction)
            );
        }
    }

    forAll(*this, cloudFunctioni)
    {
        this->operator[](cloudFunctioni).postCrossFaces(fraction);
    }
}


// ************************************************************************* //
