/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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

#include "meshToMesh.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

inline const Foam::polyMesh& Foam::meshToMesh::srcRegion() const
{
    return srcRegion_;
}


inline const Foam::polyMesh& Foam::meshToMesh::tgtRegion() const
{
    return tgtRegion_;
}


inline const Foam::labelListList&
Foam::meshToMesh::srcToTgtCellAddr() const
{
    return srcToTgtCellAddr_;
}


inline const Foam::labelListList&
Foam::meshToMesh::tgtToSrcCellAddr() const
{
    return tgtToSrcCellAddr_;
}


inline const Foam::scalarListList&
Foam::meshToMesh::srcToTgtCellWght() const
{
    return srcToTgtCellWght_;
}


inline const Foam::scalarListList&
Foam::meshToMesh::tgtToSrcCellWght() const
{
    return tgtToSrcCellWght_;
}


inline Foam::scalar Foam::meshToMesh::V() const
{
    return V_;
}


inline const Foam::PtrList<Foam::AMIInterpolation>&
Foam::meshToMesh::patchAMIs() const
{
    return patchAMIs_;
}


// ************************************************************************* //
