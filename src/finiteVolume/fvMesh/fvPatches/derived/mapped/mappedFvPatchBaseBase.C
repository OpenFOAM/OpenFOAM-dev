/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

#include "mappedFvPatchBaseBase.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mappedFvPatchBaseBase, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mappedFvPatchBaseBase::mappedFvPatchBaseBase(const fvPatch& patch)
:
    patch_(patch),
    mapper_(refCast<const mappedPatchBaseBase>(patch.patch()))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mappedFvPatchBaseBase::~mappedFvPatchBaseBase()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::mappedFvPatchBaseBase& Foam::mappedFvPatchBaseBase::getMap
(
    const fvPatch& patch
)
{
    if (!isA<mappedFvPatchBaseBase>(patch))
    {
        FatalErrorInFunction
            << "Patch " << patch.name() << " is not of type "
            << typeName << exit(FatalError);
    }

    return refCast<const mappedFvPatchBaseBase>(patch);
}


bool Foam::mappedFvPatchBaseBase::haveNbr() const
{
    const fvMesh& mesh = patch_.boundaryMesh().mesh();

    return mesh.time().foundObject<fvMesh>(mapper_.nbrRegionName());
}


const Foam::fvMesh& Foam::mappedFvPatchBaseBase::nbrMesh() const
{
    const fvMesh& mesh = patch_.boundaryMesh().mesh();

    return mesh.time().lookupObject<fvMesh>(mapper_.nbrRegionName());
}


const Foam::fvPatch& Foam::mappedFvPatchBaseBase::nbrFvPatch() const
{
    const fvMesh& nbrMesh = this->nbrMesh();

    const label patchi =
        nbrMesh.boundaryMesh().findIndex(mapper_.nbrPatchName());

    if (patchi == -1)
    {
        FatalErrorInFunction
            << "Cannot find patch " << mapper_.nbrPatchName()
            << " in region " << mapper_.nbrRegionName() << endl
            << "Valid patches are " << nbrMesh.boundaryMesh().names()
            << exit(FatalError);
    }

    return nbrMesh.boundary()[patchi];
}


// ************************************************************************* //
