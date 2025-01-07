/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "LagrangianSubMesh.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

inline const Foam::LagrangianMesh& Foam::LagrangianSubMesh::mesh() const
{
    return mesh_;
}


inline Foam::LagrangianGroup Foam::LagrangianSubMesh::group() const
{
    return group_;
}


inline Foam::label Foam::LagrangianSubMesh::size() const
{
    return size_;
}


inline Foam::label Foam::LagrangianSubMesh::globalSize() const
{
    return returnReduce(size_, sumOp<label>());
}


inline Foam::label Foam::LagrangianSubMesh::size
(
    const LagrangianSubMesh& subMesh
)
{
    return subMesh.size_;
}


inline bool Foam::LagrangianSubMesh::empty() const
{
    return !size_;
}


inline Foam::label Foam::LagrangianSubMesh::start() const
{
    return start_;
}


inline Foam::label Foam::LagrangianSubMesh::end() const
{
    return start_ + size_;
}


inline Foam::label Foam::LagrangianSubMesh::index() const
{
    return index_;
}


// ************************************************************************* //
