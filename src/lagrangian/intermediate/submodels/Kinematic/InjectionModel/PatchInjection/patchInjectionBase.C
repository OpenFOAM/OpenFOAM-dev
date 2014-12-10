/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2014 OpenFOAM Foundation
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

#include "patchInjectionBase.H"
#include "polyMesh.H"
#include "SubField.H"
#include "cachedRandom.H"
#include "triPointRef.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchInjectionBase::patchInjectionBase
(
    const polyMesh& mesh,
    const word& patchName
)
:
    patchName_(patchName),
    patchId_(mesh.boundaryMesh().findPatchID(patchName_)),
    patchArea_(0.0),
    patchNormal_(),
    cellOwners_(),
    triFace_(),
    triToFace_(),
    triCumulativeMagSf_(),
    sumTriMagSf_(Pstream::nProcs() + 1, 0.0)
{
    if (patchId_ < 0)
    {
        FatalErrorIn
        (
            "Foam::patchInjectionBase::patchInjectionBase"
            "("
                "const polyMesh&, "
                "const word&"
            ")"
        )   << "Requested patch " << patchName_ << " not found" << nl
            << "Available patches are: " << mesh.boundaryMesh().names()
            << nl << exit(FatalError);
    }

    updateMesh(mesh);
}


Foam::patchInjectionBase::patchInjectionBase(const patchInjectionBase& pib)
:
    patchName_(pib.patchName_),
    patchId_(pib.patchId_),
    patchArea_(pib.patchArea_),
    patchNormal_(pib.patchNormal_),
    cellOwners_(pib.cellOwners_),
    triFace_(pib.triFace_),
    triToFace_(pib.triToFace_),
    triCumulativeMagSf_(pib.triCumulativeMagSf_),
    sumTriMagSf_(pib.sumTriMagSf_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::patchInjectionBase::~patchInjectionBase()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::patchInjectionBase::updateMesh(const polyMesh& mesh)
{
    // Set/cache the injector cells
    const polyPatch& patch = mesh.boundaryMesh()[patchId_];
    const pointField& points = patch.points();

    cellOwners_ = patch.faceCells();

    // triangulate the patch faces and create addressing
    DynamicList<label> triToFace(2*patch.size());
    DynamicList<scalar> triMagSf(2*patch.size());
    DynamicList<face> triFace(2*patch.size());
    DynamicList<face> tris(5);

    // set zero value at the start of the tri area list
    triMagSf.append(0.0);

    forAll(patch, faceI)
    {
        const face& f = patch[faceI];

        tris.clear();
        f.triangles(points, tris);

        forAll(tris, i)
        {
            triToFace.append(faceI);
            triFace.append(tris[i]);
            triMagSf.append(tris[i].mag(points));
        }
    }

    forAll(sumTriMagSf_, i)
    {
        sumTriMagSf_[i] = 0.0;
    }

    sumTriMagSf_[Pstream::myProcNo() + 1] = sum(triMagSf);

    Pstream::listCombineGather(sumTriMagSf_, maxEqOp<scalar>());
    Pstream::listCombineScatter(sumTriMagSf_);

    for (label i = 1; i < triMagSf.size(); i++)
    {
        triMagSf[i] += triMagSf[i-1];
    }

    // transfer to persistent storage
    triFace_.transfer(triFace);
    triToFace_.transfer(triToFace);
    triCumulativeMagSf_.transfer(triMagSf);

    // convert sumTriMagSf_ into cumulative sum of areas per proc
    for (label i = 1; i < sumTriMagSf_.size(); i++)
    {
        sumTriMagSf_[i] += sumTriMagSf_[i-1];
    }

    const scalarField magSf(mag(patch.faceAreas()));
    patchArea_ = sum(magSf);
    patchNormal_ = patch.faceAreas()/magSf;
    reduce(patchArea_, sumOp<scalar>());
}


void Foam::patchInjectionBase::setPositionAndCell
(
    const polyMesh& mesh,
    cachedRandom& rnd,
    vector& position,
    label& cellOwner,
    label& tetFaceI,
    label& tetPtI
)
{
    if (cellOwners_.size() > 0)
    {
        // determine which processor to inject from
        scalar areaFraction = rnd.position<scalar>(0, patchArea_);

        label procI = 0;
        forAllReverse(sumTriMagSf_, i)
        {
            if (areaFraction > sumTriMagSf_[i])
            {
                procI = i;
                break;
            }
        }

        if (Pstream::myProcNo() == procI)
        {
            // find corresponding decomposed face triangle
            label triI = 0;
            scalar offset = sumTriMagSf_[procI];
            forAllReverse(triCumulativeMagSf_, i)
            {
                if (areaFraction > triCumulativeMagSf_[i] + offset)
                {
                    triI = i;
                    break;
                }
            }

            // set cellOwner
            label faceI = triToFace_[triI];
            cellOwner = cellOwners_[faceI];

            // find random point in triangle
            const polyPatch& patch = mesh.boundaryMesh()[patchId_];
            const pointField& points = patch.points();
            const face& tf = triFace_[triI];
            const triPointRef tri(points[tf[0]], points[tf[1]], points[tf[2]]);
            const point pf(tri.randomPoint(rnd));

            // position perturbed away from face (into domain)
            const scalar a = rnd.position(scalar(0.1), scalar(0.5));
            const vector& pc = mesh.cellCentres()[cellOwner];
            const vector d = mag(pf - pc)*patchNormal_[faceI];

            position = pf - a*d;

            // The position is between the face and cell centre, which could
            // be in any tet of the decomposed cell, so arbitrarily choose the
            // first face of the cell as the tetFace and the first point after
            // the base point on the face as the tetPt.  The tracking will pick
            // the cell consistent with the motion in the first tracking step
            tetFaceI = mesh.cells()[cellOwner][0];
            tetPtI = 1;
        }
        else
        {
            cellOwner = -1;
            tetFaceI = -1;
            tetPtI = -1;

            // dummy position
            position = pTraits<vector>::max;
        }
    }
    else
    {
        cellOwner = -1;
        tetFaceI = -1;
        tetPtI = -1;

        // dummy position
        position = pTraits<vector>::max;
    }
}


// ************************************************************************* //
