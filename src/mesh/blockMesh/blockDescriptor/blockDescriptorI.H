/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

inline const Foam::cellShape& Foam::blockDescriptor::blockShape() const
{
    return blockShape_;
}


inline const Foam::List<Foam::List<Foam::point>>&
Foam::blockDescriptor::blockEdgePoints() const
{
    return edgePoints_;
}


inline const Foam::scalarListList&
Foam::blockDescriptor::blockEdgeWeights() const
{
    return edgeWeights_;
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


inline const Foam::point& Foam::blockDescriptor::blockPoint(const label i) const
{
    return vertices_[blockShape_[i]];
}


// ************************************************************************* //
