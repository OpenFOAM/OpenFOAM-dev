/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "polyMesh.H"
#include "meshTools.H"
#include "hexMatcher.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null constructor
inline Foam::directionInfo::directionInfo()
:
    index_(-3),
    n_(vector::zero)
{}


// Construct from components
inline Foam::directionInfo::directionInfo
(
    const label index,
    const vector& n
)
:
    index_(index),
    n_(n)
{}


// Construct as copy
inline Foam::directionInfo::directionInfo(const directionInfo& w2)
:
    index_(w2.index()),
    n_(w2.n())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class TrackingData>
inline bool Foam::directionInfo::valid(TrackingData& td) const
{
    return index_ != -3;
}


// No geometric data so never any problem on cyclics
template<class TrackingData>
inline bool Foam::directionInfo::sameGeometry
(
    const polyMesh&,
    const directionInfo& w2,
    const scalar tol,
    TrackingData& td
)
 const
{
    return true;
}


// index_ is already offset in face
template<class TrackingData>
inline void Foam::directionInfo::leaveDomain
(
    const polyMesh&,
    const polyPatch& patch,
    const label patchFaceI,
    const point& faceCentre,
    TrackingData& td
)
{}


// index_ is offset in face on other side. So reverse it here.
// (Note: f[0] on other domain is connected to f[0] in this domain,
//        f[1]      ,,                         f[size-1]   ,,
// etc.)
template<class TrackingData>
inline void Foam::directionInfo::enterDomain
(
    const polyMesh&,
    const polyPatch& patch,
    const label patchFaceI,
    const point& faceCentre,
    TrackingData& td
)
{
    if (index_ >= 0)
    {
        const face& f = patch[patchFaceI];

        index_ = (f.size() - index_) % f.size();
    }
}


// No geometric data.
template<class TrackingData>
inline void Foam::directionInfo::transform
(
    const polyMesh&,
    const tensor& rotTensor,
    TrackingData& td
)
{}


// Update this cell with neighbouring face information
template<class TrackingData>
inline bool Foam::directionInfo::updateCell
(
    const polyMesh& mesh,
    const label thisCellI,
    const label neighbourFaceI,
    const directionInfo& neighbourInfo,
    const scalar,       // tol
    TrackingData& td
)
{
    if (index_ >= -2)
    {
        // Already determined.
        return false;
    }

    if (hexMatcher().isA(mesh, thisCellI))
    {
        const face& f = mesh.faces()[neighbourFaceI];

        if (neighbourInfo.index() == -2)
        {
            // Geometric information from neighbour
            index_ = -2;
        }
        else if (neighbourInfo.index() == -1)
        {
            // Cut tangential to face. Take any edge connected to face
            // but not used in face.

            // Get first edge on face.
            label edgeI = mesh.faceEdges()[neighbourFaceI][0];

            const edge& e = mesh.edges()[edgeI];

            // Find face connected to face through edgeI and on same cell.
            label faceI =
                meshTools::otherFace
                (
                    mesh,
                    thisCellI,
                    neighbourFaceI,
                    edgeI
                );

            // Find edge on faceI which is connected to e.start() but not edgeI.
            index_ =
                meshTools::otherEdge
                (
                    mesh,
                    mesh.faceEdges()[faceI],
                    edgeI,
                    e.start()
                );
        }
        else
        {
            // Index is a vertex on the face. Convert to mesh edge.

            // Get mesh edge between f[index_] and f[index_+1]
            label v0 = f[neighbourInfo.index()];
            label v1 = f[(neighbourInfo.index() + 1) % f.size()];

            index_ = findEdge(mesh, mesh.faceEdges()[neighbourFaceI], v0, v1);
        }
    }
    else
    {
        // Not a hex so mark this as geometric.
        index_ = -2;
    }


    n_ = neighbourInfo.n();

    return true;
}


// Update this face with neighbouring cell information
template<class TrackingData>
inline bool Foam::directionInfo::updateFace
(
    const polyMesh& mesh,
    const label thisFaceI,
    const label neighbourCellI,
    const directionInfo& neighbourInfo,
    const scalar,   // tol
    TrackingData& td
)
{
    // Handle special cases first

    if (index_ >= -2)
    {
        // Already determined
        return false;
    }

    // Handle normal cases where topological or geometrical info comes from
    // neighbouring cell

    if (neighbourInfo.index() >= 0)
    {
        // Neighbour has topological direction (and hence is hex). Find cut
        // edge on face.
        index_ =
            edgeToFaceIndex
            (
                mesh,
                neighbourCellI,
                thisFaceI,
                neighbourInfo.index()
            );
    }
    else
    {
        // Neighbour has geometric information. Use.
        index_ = -2;
    }


    n_ = neighbourInfo.n();

    return true;
}


// Merge this with information on same face
template<class TrackingData>
inline bool Foam::directionInfo::updateFace
(
    const polyMesh& mesh,
    const label,    // thisFaceI
    const directionInfo& neighbourInfo,
    const scalar,   // tol
    TrackingData& td
)
{
    if (index_ >= -2)
    {
        // Already visited.
        return false;
    }
    else
    {
        index_ = neighbourInfo.index();

        n_ = neighbourInfo.n();

        return true;
    }
}


template<class TrackingData>
inline bool Foam::directionInfo::equal
(
    const directionInfo& rhs,
    TrackingData& td
) const
{
    return operator==(rhs);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::directionInfo::operator==(const Foam::directionInfo& rhs)
 const
{
    return  index() == rhs.index() && n() == rhs.n();
}


inline bool Foam::directionInfo::operator!=(const Foam::directionInfo& rhs)
 const
{
    return !(*this == rhs);
}


// ************************************************************************* //
