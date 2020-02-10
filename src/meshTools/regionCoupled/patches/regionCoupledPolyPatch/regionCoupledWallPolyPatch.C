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

#include "regionCoupledWallPolyPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "polyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(regionCoupledWallPolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, regionCoupledWallPolyPatch, word);
    addToRunTimeSelectionTable
    (
        polyPatch,
        regionCoupledWallPolyPatch,
        dictionary
    );
}


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::regionCoupledWallPolyPatch::regionCoupledWallPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    wallPolyPatch(name, size, start, index, bm, patchType),
    regionCoupledBase(static_cast<const polyPatch&>(*this))
{}


Foam::regionCoupledWallPolyPatch::regionCoupledWallPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    wallPolyPatch(name, dict, index, bm, patchType),
    regionCoupledBase(static_cast<const polyPatch&>(*this), dict)
{}


Foam::regionCoupledWallPolyPatch::regionCoupledWallPolyPatch
(
    const regionCoupledWallPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    wallPolyPatch(pp, bm),
    regionCoupledBase(*this, pp)
{}


Foam::regionCoupledWallPolyPatch::regionCoupledWallPolyPatch
(
    const regionCoupledWallPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    wallPolyPatch(pp, bm, index, newSize, newStart),
    regionCoupledBase(*this, pp)
{}


Foam::regionCoupledWallPolyPatch::regionCoupledWallPolyPatch
(
    const regionCoupledWallPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const label newStart
)
:
    wallPolyPatch(pp, bm, index, mapAddressing, newStart),
    regionCoupledBase(*this, pp)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionCoupledWallPolyPatch::~regionCoupledWallPolyPatch()
{
    regionCoupledBase::clearGeom();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionCoupledWallPolyPatch::initCalcGeometry(PstreamBuffers& pBufs)
{
    wallPolyPatch::initCalcGeometry(pBufs);
}


void Foam::regionCoupledWallPolyPatch::calcGeometry(PstreamBuffers& pBufs)
{
    wallPolyPatch::calcGeometry(pBufs);
    regionCoupledBase::clearGeom();
}


void Foam::regionCoupledWallPolyPatch::initMovePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    wallPolyPatch::initMovePoints(pBufs, p);
}


void Foam::regionCoupledWallPolyPatch::movePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    wallPolyPatch::movePoints(pBufs, p);
    regionCoupledBase::clearGeom();
}


void Foam::regionCoupledWallPolyPatch::initUpdateMesh(PstreamBuffers& pBufs)
{
    wallPolyPatch::initUpdateMesh(pBufs);
}


void Foam::regionCoupledWallPolyPatch::updateMesh(PstreamBuffers& pBufs)
{
    wallPolyPatch::updateMesh(pBufs);
    regionCoupledBase::clearGeom();
}


void Foam::regionCoupledWallPolyPatch::write(Ostream& os) const
{
    wallPolyPatch::write(os);
    regionCoupledBase::write(os);
}


// ************************************************************************* //
