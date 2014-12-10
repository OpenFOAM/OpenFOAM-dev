/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "twoDPointCorrector.H"
#include "polyMesh.H"
#include "wedgePolyPatch.H"
#include "emptyPolyPatch.H"
#include "SubField.H"
#include "meshTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::scalar Foam::twoDPointCorrector::edgeOrthogonalityTol = 1.0 - 1e-4;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::twoDPointCorrector::calcAddressing() const
{
    // Find geometry normal
    planeNormalPtr_ = new vector(0, 0, 0);
    vector& pn = *planeNormalPtr_;

    // Algorithm:
    // Attempt to find wedge patch and work out the normal from it.
    // If not found, find an empty patch with faces in it and use the
    // normal of the first face.  If neither is found, declare an
    // error.

    // Try and find a wedge patch
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    forAll(patches, patchI)
    {
        if (isA<wedgePolyPatch>(patches[patchI]))
        {
            isWedge_ = true;

            const wedgePolyPatch& wp =
                refCast<const wedgePolyPatch>(patches[patchI]);

            pn = wp.centreNormal();

            wedgeAxis_ = wp.axis();
            wedgeAngle_ = mag(acos(wp.cosAngle()));

            if (polyMesh::debug)
            {
                Pout<< "Found normal from wedge patch " << patchI;
            }

            break;
        }
    }

    // Try to find an empty patch with faces
    if (!isWedge_)
    {
        forAll(patches, patchI)
        {
            if (isA<emptyPolyPatch>(patches[patchI]) && patches[patchI].size())
            {
                pn = patches[patchI].faceAreas()[0];

                if (polyMesh::debug)
                {
                    Pout<< "Found normal from empty patch " << patchI;
                }

                break;
            }
        }
    }


    if (mag(pn) < VSMALL)
    {
        FatalErrorIn("twoDPointCorrector::calcAddressing()")
            << "Cannot determine normal vector from patches."
            << abort(FatalError);
    }
    else
    {
        pn /= mag(pn);
    }

    if (polyMesh::debug)
    {
        Pout<< " twoDPointCorrector normal: " << pn << endl;
    }

    // Select edges to be included in check.
    normalEdgeIndicesPtr_ = new labelList(mesh_.nEdges());
    labelList& neIndices = *normalEdgeIndicesPtr_;

    const edgeList& meshEdges = mesh_.edges();

    const pointField& meshPoints = mesh_.points();

    label nNormalEdges = 0;

    forAll(meshEdges, edgeI)
    {
        const edge& e = meshEdges[edgeI];

        vector edgeVector = e.vec(meshPoints)/(e.mag(meshPoints) + VSMALL);

        if (mag(edgeVector & pn) > edgeOrthogonalityTol)
        {
            // this edge is normal to the plane. Add it to the list
            neIndices[nNormalEdges++] = edgeI;
        }
    }

    neIndices.setSize(nNormalEdges);

    // Construction check: number of points in a read 2-D or wedge geometry
    // should be odd and the number of edges normal to the plane should be
    // exactly half the number of points
    if (!isWedge_)
    {
        if (meshPoints.size() % 2 != 0)
        {
            WarningIn("twoDPointCorrector::calcAddressing()")
                << "the number of vertices in the geometry "
                << "is odd - this should not be the case for a 2-D case. "
                << "Please check the geometry."
                << endl;
        }

        if (2*nNormalEdges != meshPoints.size())
        {
            WarningIn("twoDPointCorrector::calcAddressing()")
                << "The number of points in the mesh is "
                << "not equal to twice the number of edges normal to the plane "
                << "- this may be OK only for wedge geometries.\n"
                << "    Please check the geometry or adjust "
                << "the orthogonality tolerance.\n" << endl
                << "Number of normal edges: " << nNormalEdges
                << " number of points: " << meshPoints.size()
                << endl;
        }
    }
}


void Foam::twoDPointCorrector::clearAddressing() const
{
    deleteDemandDrivenData(planeNormalPtr_);
    deleteDemandDrivenData(normalEdgeIndicesPtr_);
}


void Foam::twoDPointCorrector::snapToWedge
(
    const vector& n,
    const point& A,
    point& p
) const
{
    scalar ADash = mag(A - wedgeAxis_*(wedgeAxis_ & A));
    vector pDash = ADash*tan(wedgeAngle_)*planeNormal();

    p = A + sign(n & p)*pDash;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoDPointCorrector::twoDPointCorrector(const polyMesh& mesh)
:
    MeshObject<polyMesh, Foam::UpdateableMeshObject, twoDPointCorrector>(mesh),
    required_(mesh_.nGeometricD() == 2),
    planeNormalPtr_(NULL),
    normalEdgeIndicesPtr_(NULL),
    isWedge_(false),
    wedgeAxis_(vector::zero),
    wedgeAngle_(0.0)
{}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::twoDPointCorrector::~twoDPointCorrector()
{
    clearAddressing();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::direction Foam::twoDPointCorrector::normalDir() const
{
    const vector& pn = planeNormal();

    if (mag(pn.x()) >= edgeOrthogonalityTol)
    {
        return vector::X;
    }
    else if (mag(pn.y()) >= edgeOrthogonalityTol)
    {
        return vector::Y;
    }
    else if (mag(pn.z()) >= edgeOrthogonalityTol)
    {
        return vector::Z;
    }
    else
    {
        FatalErrorIn("direction twoDPointCorrector::normalDir() const")
            << "Plane normal not aligned with the coordinate system" << nl
            << "    pn = " << pn
            << abort(FatalError);

        return vector::Z;
    }
}


const Foam::vector& Foam::twoDPointCorrector::planeNormal() const
{
    if (!planeNormalPtr_)
    {
        calcAddressing();
    }

    return *planeNormalPtr_;
}


const Foam::labelList& Foam::twoDPointCorrector::normalEdgeIndices() const
{
    if (!normalEdgeIndicesPtr_)
    {
        calcAddressing();
    }

    return *normalEdgeIndicesPtr_;
}


void Foam::twoDPointCorrector::correctPoints(pointField& p) const
{
    if (!required_) return;

    // Algorithm:
    // Loop through all edges. Calculate the average point position A for
    // the front and the back. Correct the position of point P (in two planes)
    // such that vectors AP and planeNormal are parallel

    // Get reference to edges
    const edgeList&  meshEdges = mesh_.edges();

    const labelList& neIndices = normalEdgeIndices();
    const vector& pn = planeNormal();

    forAll(neIndices, edgeI)
    {
        point& pStart = p[meshEdges[neIndices[edgeI]].start()];

        point& pEnd = p[meshEdges[neIndices[edgeI]].end()];

        // calculate average point position
        point A = 0.5*(pStart + pEnd);
        meshTools::constrainToMeshCentre(mesh_, A);

        if (isWedge_)
        {
            snapToWedge(pn, A, pStart);
            snapToWedge(pn, A, pEnd);
        }
        else
        {
            // correct point locations
            pStart = A + pn*(pn & (pStart - A));
            pEnd = A + pn*(pn & (pEnd - A));
        }
    }
}


void Foam::twoDPointCorrector::correctDisplacement
(
    const pointField& p,
    vectorField& disp
) const
{
    if (!required_) return;

    // Algorithm:
    // Loop through all edges. Calculate the average point position A for
    // the front and the back. Correct the position of point P (in two planes)
    // such that vectors AP and planeNormal are parallel

    // Get reference to edges
    const edgeList&  meshEdges = mesh_.edges();

    const labelList& neIndices = normalEdgeIndices();
    const vector& pn = planeNormal();

    forAll(neIndices, edgeI)
    {
        const edge& e = meshEdges[neIndices[edgeI]];

        label startPointI = e.start();
        point pStart = p[startPointI] + disp[startPointI];

        label endPointI = e.end();
        point pEnd = p[endPointI] + disp[endPointI];

        // calculate average point position
        point A = 0.5*(pStart + pEnd);
        meshTools::constrainToMeshCentre(mesh_, A);

        if (isWedge_)
        {
            snapToWedge(pn, A, pStart);
            snapToWedge(pn, A, pEnd);
            disp[startPointI] = pStart - p[startPointI];
            disp[endPointI] = pEnd - p[endPointI];
        }
        else
        {
            // correct point locations
            disp[startPointI] = (A + pn*(pn & (pStart - A))) - p[startPointI];
            disp[endPointI] = (A + pn*(pn & (pEnd - A))) - p[endPointI];
        }
    }
}


void Foam::twoDPointCorrector::updateMesh(const mapPolyMesh&)
{
    clearAddressing();
}


bool Foam::twoDPointCorrector::movePoints()
{
    return true;
}


// ************************************************************************* //
