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

#ifndef LagrangianMeshLocation_H
#define LagrangianMeshLocation_H

#include "LagrangianMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline Ostream& operator<<(Ostream& os, const LagrangianMesh::location l)
{
    return os << static_cast<label>(l);
}


inline Istream& operator>>(Istream& is, LagrangianMesh::location& l)
{
    const label i = readLabel(is);
    l = static_cast<LagrangianMesh::location>(i);
    return is;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace LagrangianMeshLocation
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

struct closestOp
{
    inline LagrangianMesh::location operator()
    (
        const LagrangianMesh::location a,
        const LagrangianMesh::location b
    ) const
    {
        return static_cast<label>(a) < static_cast<label>(b) ? a : b;
    }
};


struct furthestOp
{
    inline LagrangianMesh::location operator()
    (
        const LagrangianMesh::location a,
        const LagrangianMesh::location b
    ) const
    {
        return static_cast<label>(a) < static_cast<label>(b) ? a : b;
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LagrangianMeshLocation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
