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

#include "cellClassification.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Update this with w2 information
template<class TrackingData>
inline bool Foam::cellInfo::update
(
    const cellInfo& w2,
    const label thisFaceI,
    const label thisCellI,
    const label neighbourFaceI,
    const label neighbourCellI,
    TrackingData& td
)
{
    if
    (
        (w2.type() == cellClassification::NOTSET)
     || (w2.type() == cellClassification::CUT)
    )
    {
        FatalErrorIn("cellInfo::update(const cellInfo&)")
            << "Problem: trying to propagate NOTSET or CUT type:" << w2.type()
            << " into cell/face with type:" << type() << endl
            << "thisFaceI:" << thisFaceI
            << "  thisCellI:" << thisCellI
            << "  neighbourFaceI:" << neighbourFaceI
            << "  neighbourCellI:" << neighbourCellI
            << abort(FatalError);
        return false;
    }

    if (type() == cellClassification::NOTSET)
    {
        type_ = w2.type();

        return true;
    }

    if (type() == cellClassification::CUT)
    {
        // Reached boundary. Stop.
        return false;
    }

    if (type() == w2.type())
    {
        // Should never happen; already checked in meshWave
        return false;
    }

    // Two conflicting types
    FatalErrorIn("cellInfo::update(const cellInfo&)")
        << "Problem: trying to propagate conflicting types:" << w2.type()
        << " into cell/face with type:" << type() << endl
        << "thisFaceI:" << thisFaceI
        << "  thisCellI:" << thisCellI
        << "  neighbourFaceI:" << neighbourFaceI
        << "  neighbourCellI:" << neighbourCellI
        << abort(FatalError);

    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null constructor
inline Foam::cellInfo::cellInfo()
:
    type_(cellClassification::NOTSET)
{}


// Construct from components
inline Foam::cellInfo::cellInfo(const label type)
:
    type_(type)
{}


// Construct as copy
inline Foam::cellInfo::cellInfo(const cellInfo& w2)
:
    type_(w2.type())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class TrackingData>
inline bool Foam::cellInfo::valid(TrackingData& td) const
{
    return type_ != cellClassification::NOTSET;
}


// No geometric data so never any problem on cyclics
template<class TrackingData>
inline bool Foam::cellInfo::sameGeometry
(
    const polyMesh&,
    const cellInfo& w2,
    const scalar tol,
    TrackingData& td
)
 const
{
    return true;
}


// No geometric data.
template<class TrackingData>
inline void Foam::cellInfo::leaveDomain
(
    const polyMesh&,
    const polyPatch& patch,
    const label patchFaceI,
    const point& faceCentre,
    TrackingData& td
)
{}


// No geometric data.
template<class TrackingData>
inline void Foam::cellInfo::transform
(
    const polyMesh&,
    const tensor& rotTensor,
    TrackingData& td
)
{}


// No geometric data.
template<class TrackingData>
inline void Foam::cellInfo::enterDomain
(
    const polyMesh&,
    const polyPatch& patch,
    const label patchFaceI,
    const point& faceCentre,
    TrackingData& td
)
{}


// Update this with neighbour information
template<class TrackingData>
inline bool Foam::cellInfo::updateCell
(
    const polyMesh&,
    const label thisCellI,
    const label neighbourFaceI,
    const cellInfo& neighbourInfo,
    const scalar tol,
    TrackingData& td
)
{
    return update
    (
        neighbourInfo,
        -1,
        thisCellI,
        neighbourFaceI,
        -1,
        td
    );
}


// Update this with neighbour information
template<class TrackingData>
inline bool Foam::cellInfo::updateFace
(
    const polyMesh&,
    const label thisFaceI,
    const label neighbourCellI,
    const cellInfo& neighbourInfo,
    const scalar tol,
    TrackingData& td
)
{
    return update
    (
        neighbourInfo,
        thisFaceI,
        -1,
        -1,
        neighbourCellI,
        td
    );
}

// Update this with neighbour information
template<class TrackingData>
inline bool Foam::cellInfo::updateFace
(
    const polyMesh&,
    const label thisFaceI,
    const cellInfo& neighbourInfo,
    const scalar tol,
    TrackingData& td
)
{
    return update
    (
        neighbourInfo,
        thisFaceI,
        -1,
        -1,
        -1,
        td
    );
}


template<class TrackingData>
inline bool Foam::cellInfo::equal
(
    const cellInfo& rhs,
    TrackingData& td
) const
{
    return operator==(rhs);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::cellInfo::operator==(const Foam::cellInfo& rhs) const
{
    return type() == rhs.type();
}


inline bool Foam::cellInfo::operator!=(const Foam::cellInfo& rhs) const
{
    return !(*this == rhs);
}


// ************************************************************************* //
