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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::label Foam::extendedEdgeMesh::convexStart() const
{
    return convexStart_;
}


inline Foam::label Foam::extendedEdgeMesh::concaveStart() const
{
    return concaveStart_;
}


inline Foam::label Foam::extendedEdgeMesh::mixedStart() const
{
    return mixedStart_;
}


inline Foam::label Foam::extendedEdgeMesh::nonFeatureStart() const
{
    return nonFeatureStart_;
}


inline Foam::label Foam::extendedEdgeMesh::externalStart() const
{
    return externalStart_;
}


inline Foam::label Foam::extendedEdgeMesh::internalStart() const
{
    return internalStart_;
}


inline Foam::label Foam::extendedEdgeMesh::flatStart() const
{
    return flatStart_;
}


inline Foam::label Foam::extendedEdgeMesh::openStart() const
{
    return openStart_;
}


inline Foam::label Foam::extendedEdgeMesh::multipleStart() const
{
    return multipleStart_;
}


inline bool Foam::extendedEdgeMesh::featurePoint(label ptI) const
{
    return ptI < nonFeatureStart_;
}


inline const Foam::vectorField& Foam::extendedEdgeMesh::normals() const
{
    return normals_;
}


inline const Foam::List<Foam::extendedEdgeMesh::sideVolumeType>&
Foam::extendedEdgeMesh::normalVolumeTypes() const
{
    return normalVolumeTypes_;
}


inline const Foam::vectorField& Foam::extendedEdgeMesh::edgeDirections()
const
{
    return edgeDirections_;
}


inline const Foam::labelListList&
Foam::extendedEdgeMesh::normalDirections() const
{
    return normalDirections_;
}


inline Foam::vector Foam::extendedEdgeMesh::edgeDirection
(
    label edgeI,
    label ptI
) const
{
    const edge& e = edges()[edgeI];

    if (ptI == e.start())
    {
        return edgeDirections()[edgeI];
    }
    else if (ptI == e.end())
    {
        return -edgeDirections()[edgeI];
    }
    else
    {
        FatalErrorIn("Foam::extendedEdgeMesh::edgeDirection")
            << "Requested ptI " << ptI << " is not a point on the requested "
            << "edgeI " << edgeI << ". edgeI start and end: "
            << e.start() << " " << e.end()
            << exit(FatalError);

        return vector::zero;
    }
}


inline const Foam::labelListList& Foam::extendedEdgeMesh::edgeNormals()
const
{
    return edgeNormals_;
}


inline Foam::vectorField Foam::extendedEdgeMesh::edgeNormals
(
    const labelList& edgeNormIs
) const
{
    vectorField norms(edgeNormIs.size());

    forAll(edgeNormIs, i)
    {
        norms[i] = normals_[edgeNormIs[i]];
    }

    return norms;
}


inline Foam::vectorField Foam::extendedEdgeMesh::edgeNormals(label edgeI)
const
{
    return edgeNormals(edgeNormals_[edgeI]);
}


inline const Foam::labelListList&
Foam::extendedEdgeMesh::featurePointNormals() const
{
    return featurePointNormals_;
}


inline Foam::vectorField Foam::extendedEdgeMesh::featurePointNormals
(
    label ptI
) const
{
    if (!featurePoint(ptI))
    {
        WarningIn("vectorField extendedEdgeMesh::featurePointNormals")
            << "Requesting the normals of a non-feature point. "
            << "Returned zero length vectorField."
            << endl;

        return vectorField(0);
    }

    labelList featPtNormIs(featurePointNormals_[ptI]);

    vectorField norms(featPtNormIs.size());

    forAll(featPtNormIs, i)
    {
        norms[i] = normals_[featPtNormIs[i]];
    }

    return norms;
}


inline const Foam::labelListList&
Foam::extendedEdgeMesh::featurePointEdges() const
{
    return featurePointEdges_;
}


inline const Foam::labelList& Foam::extendedEdgeMesh::regionEdges() const
{
    return regionEdges_;
}


inline Foam::extendedEdgeMesh::pointStatus
Foam::extendedEdgeMesh::getPointStatus(label ptI) const
{
    if (ptI < concaveStart_)
    {
        return CONVEX;
    }
    else if (ptI < mixedStart_)
    {
        return CONCAVE;
    }
    else if (ptI < nonFeatureStart_)
    {
        return MIXED;
    }
    else
    {
        return NONFEATURE;
    }
}


inline Foam::extendedEdgeMesh::edgeStatus
Foam::extendedEdgeMesh::getEdgeStatus(label edgeI) const
{
    if (edgeI < internalStart_)
    {
        return EXTERNAL;
    }
    else if (edgeI < flatStart_)
    {
        return INTERNAL;
    }
    else if (edgeI < openStart_)
    {
        return FLAT;
    }
    else if (edgeI < multipleStart_)
    {
        return OPEN;
    }
    else
    {
        return MULTIPLE;
    }
}


inline Foam::PackedList<2> Foam::extendedEdgeMesh::edgeBaffles
(
    label edgeI
) const
{
    const labelList& eNormals = edgeNormals_[edgeI];

    DynamicList<label> edgeBaffles(eNormals.size());

    forAll(eNormals, enI)
    {
        const label normI = eNormals[enI];

        if (normalVolumeTypes_[normI])
        {
            edgeBaffles.append(normI);
        }
    }

    return PackedList<2>(edgeBaffles);
}


// ************************************************************************* //
