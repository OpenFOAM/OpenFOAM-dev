/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2023 OpenFOAM Foundation
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

#include "nonConformalCoupledPolyPatch.H"
#include "nonConformalErrorPolyPatch.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(nonConformalCoupledPolyPatch, 0);
}


// * * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * //

void Foam::nonConformalCoupledPolyPatch::rename(const wordList& newNames)
{
    nonConformalPolyPatch::rename(newNames);

    errorPatchName_ = word::null;
    errorPatchIndex_ = -1;
}


void Foam::nonConformalCoupledPolyPatch::reorder
(
    const labelUList& newToOldIndex
)
{
    nonConformalPolyPatch::reorder(newToOldIndex);

    errorPatchName_ = word::null;
    errorPatchIndex_ = -1;
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::nonConformalCoupledPolyPatch::nonConformalCoupledPolyPatch
(
    const polyPatch& patch
)
:
    nonConformalPolyPatch(patch),
    patch_(refCast<const coupledPolyPatch>(patch)),
    errorPatchName_(word::null),
    errorPatchIndex_(-1)
{}


Foam::nonConformalCoupledPolyPatch::nonConformalCoupledPolyPatch
(
    const polyPatch& patch,
    const word& origPatchName
)
:
    nonConformalPolyPatch(patch, origPatchName),
    patch_(refCast<const coupledPolyPatch>(patch)),
    errorPatchName_(word::null),
    errorPatchIndex_(-1)
{}


Foam::nonConformalCoupledPolyPatch::nonConformalCoupledPolyPatch
(
    const polyPatch& patch,
    const dictionary& dict
)
:
    nonConformalPolyPatch(patch, dict),
    patch_(refCast<const coupledPolyPatch>(patch)),
    errorPatchName_(word::null),
    errorPatchIndex_(-1)
{}


Foam::nonConformalCoupledPolyPatch::nonConformalCoupledPolyPatch
(
    const polyPatch& patch,
    const nonConformalCoupledPolyPatch& nccPatch
)
:
    nonConformalPolyPatch(patch, nccPatch),
    patch_(refCast<const coupledPolyPatch>(patch)),
    errorPatchName_(word::null),
    errorPatchIndex_(-1)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nonConformalCoupledPolyPatch::~nonConformalCoupledPolyPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::nonConformalCoupledPolyPatch::owner() const
{
    return patch_.owner();
}


bool Foam::nonConformalCoupledPolyPatch::neighbour() const
{
    return patch_.neighbour();
}


const Foam::transformer& Foam::nonConformalCoupledPolyPatch::transform() const
{
    return patch_.transform();
}


const Foam::word& Foam::nonConformalCoupledPolyPatch::errorPatchName() const
{
    return patch_.boundaryMesh()[errorPatchIndex()].name();
}


Foam::label Foam::nonConformalCoupledPolyPatch::errorPatchIndex() const
{
    if (errorPatchIndex_ == -1)
    {
        forAll(patch_.boundaryMesh(), patchi)
        {
            const polyPatch& p = patch_.boundaryMesh()[patchi];

            if
            (
                isA<nonConformalErrorPolyPatch>(p)
             && refCast<const nonConformalErrorPolyPatch>(p).origPatchIndex()
             == origPatchIndex()
            )
            {
                errorPatchIndex_ = patchi;
                break;
            }
        }
    }

    if (errorPatchIndex_ == -1)
    {
        FatalErrorInFunction
            << "No error patch of type "
            << nonConformalErrorPolyPatch::typeName
            << " defined for patch " << origPatchName()
            << exit(FatalError);
    }

    return errorPatchIndex_;
}


const Foam::nonConformalErrorPolyPatch&
Foam::nonConformalCoupledPolyPatch::errorPatch() const
{
    return
        refCast<const nonConformalErrorPolyPatch>
        (
            patch_.boundaryMesh()[errorPatchIndex()]
        );
}


void Foam::nonConformalCoupledPolyPatch::write(Ostream& os) const
{
    nonConformalPolyPatch::write(os);
}


// ************************************************************************* //
