/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "wedgePolyPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "SubField.H"
#include "transform.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(wedgePolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, wedgePolyPatch, word);
    addToRunTimeSelectionTable(polyPatch, wedgePolyPatch, dictionary);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::wedgePolyPatch::calcGeometry(PstreamBuffers&)
{
    if (axis_ != vector::rootMax)
    {
        return;
    }

    if (returnReduce(size(), sumOp<label>()))
    {
        const vectorField& nf(faceNormals());
        n_ = gAverage(nf);

        if (debug)
        {
            Info<< "Patch " << name() << " calculated average normal "
                << n_ << endl;
        }


        // Check the wedge is planar
        forAll(nf, facei)
        {
            if (magSqr(n_ - nf[facei]) > small)
            {
                // only issue warning instead of error so that the case can
                // still be read for post-processing
                WarningInFunction
                    << "Wedge patch '" << name() << "' is not planar." << nl
                    << "At local face at "
                    << primitivePatch::faceCentres()[facei]
                    << " the normal " << nf[facei]
                    << " differs from the average normal " << n_
                    << " by " << magSqr(n_ - nf[facei]) << nl
                    << "Either correct the patch or split it into planar parts"
                    << endl;
            }
        }

        centreNormal_ =
            vector
            (
                sign(n_.x())*(max(mag(n_.x()), 0.5) - 0.5),
                sign(n_.y())*(max(mag(n_.y()), 0.5) - 0.5),
                sign(n_.z())*(max(mag(n_.z()), 0.5) - 0.5)
            );
        centreNormal_ /= mag(centreNormal_);

        cosAngle_ = centreNormal_ & n_;

        const scalar cnCmptSum =
            centreNormal_.x() + centreNormal_.y() + centreNormal_.z();

        if (mag(cnCmptSum) < (1 - small))
        {
            FatalErrorInFunction
                << "wedge " << name()
                << " centre plane does not align with a coordinate plane by "
                << 1 - mag(cnCmptSum)
                << exit(FatalError);
        }

        axis_ = centreNormal_ ^ n_;
        scalar magAxis = mag(axis_);

        if (magAxis < small)
        {
            FatalErrorInFunction
                << "wedge " << name()
                << " plane aligns with a coordinate plane." << nl
                << "    The wedge plane should make a small angle (~2.5deg)"
                   " with the coordinate plane" << nl
                << "    and the pair of wedge planes should be symmetric"
                << " about the coordinate plane." << nl
                << "    Normal of wedge plane is " << n_
                << " , implied coordinate plane direction is " << centreNormal_
                << exit(FatalError);
        }

        axis_ /= magAxis;

        faceT_ = rotationTensor(centreNormal_, n_);
        cellT_ = faceT_ & faceT_;
    }
}


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::wedgePolyPatch::wedgePolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    polyPatch(name, size, start, index, bm, patchType),
    axis_(vector::rootMax),
    centreNormal_(vector::rootMax),
    n_(vector::rootMax),
    cosAngle_(0.0),
    faceT_(Zero),
    cellT_(Zero)
{}


Foam::wedgePolyPatch::wedgePolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    polyPatch(name, dict, index, bm, patchType),
    axis_(vector::rootMax),
    centreNormal_(vector::rootMax),
    n_(vector::rootMax),
    cosAngle_(0.0),
    faceT_(Zero),
    cellT_(Zero)
{}


Foam::wedgePolyPatch::wedgePolyPatch
(
    const wedgePolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    polyPatch(pp, bm),
    axis_(pp.axis_),
    centreNormal_(pp.centreNormal_),
    n_(pp.n_),
    cosAngle_(pp.cosAngle_),
    faceT_(pp.faceT_),
    cellT_(pp.cellT_)
{}


Foam::wedgePolyPatch::wedgePolyPatch
(
    const wedgePolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    polyPatch(pp, bm, index, newSize, newStart),
    axis_(pp.axis_),
    centreNormal_(pp.centreNormal_),
    n_(pp.n_),
    cosAngle_(pp.cosAngle_),
    faceT_(pp.faceT_),
    cellT_(pp.cellT_)
{}


Foam::wedgePolyPatch::wedgePolyPatch
(
    const wedgePolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const label newStart
)
:
    polyPatch(pp, bm, index, mapAddressing, newStart),
    axis_(pp.axis_),
    centreNormal_(pp.centreNormal_),
    n_(pp.n_),
    cosAngle_(pp.cosAngle_),
    faceT_(pp.faceT_),
    cellT_(pp.cellT_)
{}


// ************************************************************************* //
