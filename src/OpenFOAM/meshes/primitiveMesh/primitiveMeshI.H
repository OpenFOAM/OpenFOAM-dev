/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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

#include "primitiveMesh.H"

// * * * * * * * * * * Protected Static Member Functions * * * * * * * * * * //

template<class Class, class Type>
inline bool Foam::primitiveMesh::isMethod
(
    const Type& (Class::*method)() const,
    const Type& (Class::*otherMethod)() const
)
{
    return method == otherMethod;
}


template<class Class, class Type, class OtherType>
inline bool Foam::primitiveMesh::isMethod
(
    const Type& (Class::*method)() const,
    const OtherType& (Class::*otherMethod)() const
)
{
    return false;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::label Foam::primitiveMesh::nInternalPoints() const
{
    return nInternalPoints_;
}


inline Foam::label Foam::primitiveMesh::nPoints() const
{
    return nPoints_;
}


inline Foam::label Foam::primitiveMesh::nInternal0Edges() const
{
    // Force edge calculation
    (void)nEdges();
    return nInternal0Edges_;
}


inline Foam::label Foam::primitiveMesh::nInternal1Edges() const
{
    // Force edge calculation
    (void)nEdges();
    return nInternal1Edges_;
}


inline Foam::label Foam::primitiveMesh::nInternalEdges() const
{
    // Force edge calculation
    (void)nEdges();
    return nInternalEdges_;
}


inline Foam::label Foam::primitiveMesh::nEdges() const
{
    if (nEdges_ < 0)
    {
        nEdges_ = edges().size();
    }

    return nEdges_;
}


inline Foam::label Foam::primitiveMesh::nInternalFaces() const
{
    return nInternalFaces_;
}


inline Foam::label Foam::primitiveMesh::nFaces() const
{
    return nFaces_;
}


inline Foam::label Foam::primitiveMesh::nCells() const
{
    return nCells_;
}


inline bool Foam::primitiveMesh::isInternalFace(const label faceIndex) const
{
    return faceIndex < nInternalFaces();
}


template<class Type>
inline bool Foam::primitiveMesh::has
(
    const Type& (primitiveMesh::*method)() const
) const
{
    if (isMethod(method, &primitiveMesh::cellShapes)) return cellShapesPtr_;
    if (isMethod(method, &primitiveMesh::edges)) return edgesPtr_;
    if (isMethod(method, &primitiveMesh::cellCells)) return ccPtr_;
    if (isMethod(method, &primitiveMesh::edgeCells)) return ecPtr_;
    if (isMethod(method, &primitiveMesh::pointCells)) return pcPtr_;
    if (isMethod(method, &primitiveMesh::cells)) return cfPtr_;
    if (isMethod(method, &primitiveMesh::edgeFaces)) return efPtr_;
    if (isMethod(method, &primitiveMesh::pointFaces)) return pfPtr_;
    if (isMethod(method, &primitiveMesh::cellEdges)) return cePtr_;
    if (isMethod(method, &primitiveMesh::faceEdges)) return fePtr_;
    if (isMethod(method, &primitiveMesh::pointEdges)) return pePtr_;
    if (isMethod(method, &primitiveMesh::pointPoints)) return ppPtr_;
    if (isMethod(method, &primitiveMesh::cellPoints)) return cpPtr_;
    if (isMethod(method, &primitiveMesh::cellCentres)) return cellCentresPtr_;
    if (isMethod(method, &primitiveMesh::faceCentres)) return faceCentresPtr_;
    if (isMethod(method, &primitiveMesh::cellVolumes)) return cellVolumesPtr_;
    if (isMethod(method, &primitiveMesh::faceAreas)) return faceAreasPtr_;
    if (isMethod(method, &primitiveMesh::magFaceAreas)) return magFaceAreasPtr_;

    return false;
}


// ************************************************************************* //
