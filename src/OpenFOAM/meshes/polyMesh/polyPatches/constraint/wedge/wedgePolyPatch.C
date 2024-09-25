/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
#include "polyMesh.H"
#include "polyBoundaryMesh.H"
#include "RemoteData.H"
#include "SubField.H"
#include "Tuple3.H"
#include "transform.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(wedgePolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, wedgePolyPatch, word);
    addToRunTimeSelectionTable(polyPatch, wedgePolyPatch, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::wedgePolyPatch::isQuadFace(const label facei) const
{
    const polyMesh& mesh = boundaryMesh().mesh();

    const face& f = mesh.faces()[facei];

    // Face must be of size 4
    if (f.size() != 4) return false;

    // Face must not have any repeated points
    forAll(f, fpi)
    {
        for (label fpj = fpi + 1; fpj < f.size(); ++ fpj)
        {
            if (f[fpi] == f[fpj]) return false;
        }
    }

    return true;
}


bool Foam::wedgePolyPatch::isWedgeFace(const label facei) const
{
    return
        facei >= boundaryMesh().mesh().nInternalFaces()
     && isA<wedgePolyPatch>
        (
            boundaryMesh()
            [
                boundaryMesh().patchIndices()
                [
                    facei - boundaryMesh().mesh().nInternalFaces()
                ]
            ]
        );
}


Foam::label Foam::wedgePolyPatch::oppositeWedgeFace(const label thisFacei) const
{
    const polyMesh& mesh = boundaryMesh().mesh();

    // Function to generate a generic error that the search for the opposite
    // wedge face has failed
    auto error = [&]()
    {
        FatalErrorInFunction
            << "Wedge face not found opposite face " << thisFacei;

        if (Pstream::parRun())
        {
            FatalErrorInFunction
                << " on processor " << Pstream::myProcNo();
        }

        FatalErrorInFunction
            << " at " << mesh.faceCentres()[thisFacei]
            << exit(FatalError);
    };

    const face& thisF = mesh.faces()[thisFacei];

    const cell& c = mesh.cells()[mesh.faceOwner()[thisFacei]];

    // Look for a "mid" face. This is any a quad face in the cell that is not
    // on a wedge patch and is edge-connected to this face.
    label midFacei = -1, midThisFaceEdgei = -1;
    for (label cfi = 0; cfi < c.size() && midFacei == -1; ++ cfi)
    {
        const label facei = c[cfi];
        const face& f = mesh.faces()[facei];

        if
        (
            thisFacei == facei
         || !isQuadFace(facei)
         || isWedgeFace(facei)
        ) continue;

        for
        (
            label thisFei = 0;
            thisFei < thisF.size() && midFacei == -1;
            ++ thisFei
        )
        {
            const edge thisE = thisF.faceEdge(thisFei);

            for (label fei = 0; fei < f.size() && midFacei == -1; ++ fei)
            {
                const edge e = f.faceEdge(fei);

                if (edge::compare(thisE, e) != 0)
                {
                    midFacei = facei;
                    midThisFaceEdgei = fei;
                }
            }
        }
    }

    if (midFacei == -1) error();

    // Get the edge on the opposite side of the "mid" face
    const label midOppFaceEdgei = (midThisFaceEdgei + 2) % 4;
    const edge midOppE = mesh.faces()[midFacei].faceEdge(midOppFaceEdgei);

    // Look for the opposite face. This is a face in the cell that is on a
    // wedge patch and is edge connected to the edge of the "mid" face opposite
    // the provided face.
    label oppFacei = -1;
    for (label cfi = 0; cfi < c.size() && oppFacei == -1; ++ cfi)
    {
        const label facei = c[cfi];
        const face& f = mesh.faces()[facei];

        if
        (
            facei == thisFacei
         || facei == midFacei
         || !isWedgeFace(facei)
        ) continue;

        for (label fei = 0; fei < f.size() && oppFacei == -1; ++ fei)
        {
            const edge e = f.faceEdge(fei);

            if (edge::compare(midOppE, e) != 0)
            {
                oppFacei = facei;
            }
        }
    }

    if (oppFacei == -1) error();

    return oppFacei;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::wedgePolyPatch::calcGeometry(PstreamBuffers&)
{
    if (axis_ != vector::rootMax)
    {
        return;
    }

    if (!returnReduce(size(), sumOp<label>()))
    {
        return;
    }

    const pointField& cf(primitivePatch::faceCentres());
    const vectorField& nf(faceNormals());
    n_ = gAverage(nf);

    if (debug)
    {
        Info<< "Patch " << name() << " calculated average normal "
            << n_ << endl;
    }

    // Get the largest variation from the average normal
    RemoteData<Tuple3<scalar, vector, vector>> maxDeltaN;
    if (size())
    {
        const scalarField deltaNSqr(magSqr(n_ - nf));
        const label i = findMax(deltaNSqr);
        maxDeltaN = {Pstream::myProcNo(), i, {deltaNSqr[i], cf[i], nf[i]}};
    }
    combineReduce
    (
        maxDeltaN,
        RemoteData<Tuple3<scalar, vector, vector>>::greatestFirstEqOp()
    );

    // Warn if the wedge is not planar
    if (maxDeltaN.data.first() > small)
    {
        Ostream& wos =
            WarningInFunction
                << "Wedge patch '" << name()
                << "' may not be sufficiently planar" << nl;

        wos << "At patch face #" << maxDeltaN.elementi;

        if (Pstream::parRun())
        {
            wos << " on processor #" << maxDeltaN.proci;
        }

        wos << " with centre " << maxDeltaN.data.second()
            << " the normal " << maxDeltaN.data.third()
            << " differs from the average normal " << n_
            << " by " << maxDeltaN.data.first() << nl << endl;
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
    cellT_(Zero),
    oppositePatchIndex_(-1)
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
    cellT_(Zero),
    oppositePatchIndex_(-1)
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
    cellT_(pp.cellT_),
    oppositePatchIndex_(-1)
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
    cellT_(pp.cellT_),
    oppositePatchIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::wedgePolyPatch::oppositePatchIndex() const
{
    if (oppositePatchIndex_ == -1)
    {
        if (size())
        {
            oppositePatchIndex_ =
                boundaryMesh().patchIndices()
                [
                    oppositeWedgeFace(start())
                  - boundaryMesh().mesh().nInternalFaces()
                ];
        }

        reduce(oppositePatchIndex_, maxOp<label>());
    }

    return oppositePatchIndex_;
}


const Foam::wedgePolyPatch& Foam::wedgePolyPatch::oppositePatch() const
{
    const polyPatch& pp = boundaryMesh()[oppositePatchIndex()];
    return refCast<const wedgePolyPatch>(pp);
}


// ************************************************************************* //
