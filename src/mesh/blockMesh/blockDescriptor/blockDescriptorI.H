/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::pointField& Foam::blockDescriptor::vertices() const
{
    return vertices_;
}


inline const Foam::blockFaceList& Foam::blockDescriptor::faces() const
{
    return faces_;
}


inline const Foam::cellShape& Foam::blockDescriptor::blockShape() const
{
    return blockShape_;
}


inline const Foam::Vector<Foam::label>& Foam::blockDescriptor::density() const
{
    return density_;
}


inline const Foam::word& Foam::blockDescriptor::zoneName() const
{
    return zoneName_;
}


inline Foam::label Foam::blockDescriptor::nPoints() const
{
    return
    (
        (density_.x() + 1)
      * (density_.y() + 1)
      * (density_.z() + 1)
    );
}


inline Foam::label Foam::blockDescriptor::nCells() const
{
    return
    (
        density_.x()
      * density_.y()
      * density_.z()
    );
}


inline const Foam::FixedList<Foam::label, 6>&
Foam::blockDescriptor::curvedFaces() const
{
    return curvedFaces_;
}


inline Foam::label Foam::blockDescriptor::nCurvedFaces() const
{
    return nCurvedFaces_;
}


inline const Foam::point& Foam::blockDescriptor::blockPoint(const label i) const
{
    return vertices_[blockShape_[i]];
}


inline Foam::label Foam::blockDescriptor::pointLabel
(
    const label i,
    const label j,
    const label k
) const
{
    return
    (
        i
      + j*(density_.x() + 1)
      + k*(density_.x() + 1)*(density_.y() + 1)
    );
}


inline Foam::label Foam::blockDescriptor::facePointLabel
(
    const label facei,
    const label i,
    const label j
) const
{
    if (facei == 0 || facei == 1)
    {
        return
        (
            i
          + j*(density_.y() + 1)
        );
    }
    else if (facei == 2 || facei == 3)
    {
        return
        (
            i
          + j*(density_.x() + 1)
        );
    }
    else
    {
        return
        (
            i
          + j*(density_.x() + 1)
        );
    }
}


inline bool Foam::blockDescriptor::vertex
(
    const label i, const label j, const label k
) const
{
    bool iEnd = (i == 0 || i == density_.x());
    bool jEnd = (j == 0 || j == density_.y());
    bool kEnd = (k == 0 || k == density_.z());

    return (iEnd && jEnd && kEnd);
}


inline bool Foam::blockDescriptor::edge
(
    const label i, const label j, const label k
) const
{
    bool iEnd = (i == 0 || i == density_.x());
    bool jEnd = (j == 0 || j == density_.y());
    bool kEnd = (k == 0 || k == density_.z());

    return (iEnd && jEnd) || (iEnd && kEnd) || (jEnd && kEnd);
}


inline bool Foam::blockDescriptor::flatFaceOrEdge
(
    const label i, const label j, const label k
) const
{
    if (i == 0 && curvedFaces_[0] == -1) return true;
    if (i == density_.x() && curvedFaces_[1] == -1) return true;
    if (j == 0 && curvedFaces_[2] == -1) return true;
    if (j == density_.y() && curvedFaces_[3] == -1) return true;
    if (k == 0 && curvedFaces_[4] == -1) return true;
    if (k == density_.z() && curvedFaces_[5] == -1) return true;

    bool iEnd = (i == 0 || i == density_.x());
    bool jEnd = (j == 0 || j == density_.y());
    bool kEnd = (k == 0 || k == density_.z());

    return (iEnd && jEnd) || (iEnd && kEnd) || (jEnd && kEnd);
}


// ************************************************************************* //
