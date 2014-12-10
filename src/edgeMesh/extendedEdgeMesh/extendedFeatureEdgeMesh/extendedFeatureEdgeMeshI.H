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

inline Foam::label Foam::extendedFeatureEdgeMesh::convexStart() const
{
    return convexStart_;
}


inline Foam::label Foam::extendedFeatureEdgeMesh::concaveStart() const
{
    return concaveStart_;
}


inline Foam::label Foam::extendedFeatureEdgeMesh::mixedStart() const
{
    return mixedStart_;
}


inline Foam::label Foam::extendedFeatureEdgeMesh::nonFeatureStart() const
{
    return nonFeatureStart_;
}


inline Foam::label Foam::extendedFeatureEdgeMesh::externalStart() const
{
    return externalStart_;
}


inline Foam::label Foam::extendedFeatureEdgeMesh::internalStart() const
{
    return internalStart_;
}


inline Foam::label Foam::extendedFeatureEdgeMesh::flatStart() const
{
    return flatStart_;
}


inline Foam::label Foam::extendedFeatureEdgeMesh::openStart() const
{
    return openStart_;
}


inline Foam::label Foam::extendedFeatureEdgeMesh::multipleStart() const
{
    return multipleStart_;
}


inline bool Foam::extendedFeatureEdgeMesh::featurePoint(label ptI) const
{
    return ptI < nonFeatureStart_;
}


inline const Foam::vectorField& Foam::extendedFeatureEdgeMesh::normals() const
{
    return normals_;
}

inline const Foam::List<Foam::extendedFeatureEdgeMesh::sideVolumeType>&
Foam::extendedFeatureEdgeMesh::normalVolumeTypes() const
{
    return normalVolumeTypes_;
}

inline const Foam::vectorField& Foam::extendedFeatureEdgeMesh::edgeDirections()
const
{
    return edgeDirections_;
}

inline const Foam::labelListList&
Foam::extendedFeatureEdgeMesh::normalDirections() const
{
    return normalDirections_;
}


inline Foam::vector Foam::extendedFeatureEdgeMesh::edgeDirection
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
        FatalErrorIn("Foam::extendedFeatureEdgeMesh::edgeDirection")
            << "Requested ptI " << ptI << " is not a point on the requested "
            << "edgeI " << edgeI << ". edgeI start and end: "
            << e.start() << " " << e.end()
            << exit(FatalError);

        return vector::zero;
    }
}


inline const Foam::labelListList& Foam::extendedFeatureEdgeMesh::edgeNormals()
const
{
    return edgeNormals_;
}


inline Foam::vectorField Foam::extendedFeatureEdgeMesh::edgeNormals
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


inline Foam::vectorField Foam::extendedFeatureEdgeMesh::edgeNormals(label edgeI)
const
{
    return edgeNormals(edgeNormals_[edgeI]);
}


inline const Foam::labelListList&
Foam::extendedFeatureEdgeMesh::featurePointNormals() const
{
    return featurePointNormals_;
}


inline Foam::vectorField Foam::extendedFeatureEdgeMesh::featurePointNormals
(
    label ptI
) const
{
    if (!featurePoint(ptI))
    {
        WarningIn("vectorField extendedFeatureEdgeMesh::featurePointNormals")
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
Foam::extendedFeatureEdgeMesh::featurePointEdges() const
{
    return featurePointEdges_;
}


inline const Foam::labelList& Foam::extendedFeatureEdgeMesh::regionEdges() const
{
    return regionEdges_;
}


inline Foam::extendedFeatureEdgeMesh::pointStatus
Foam::extendedFeatureEdgeMesh::getPointStatus(label ptI) const
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


inline Foam::extendedFeatureEdgeMesh::edgeStatus
Foam::extendedFeatureEdgeMesh::getEdgeStatus(label edgeI) const
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


inline Foam::PackedList<2> Foam::extendedFeatureEdgeMesh::edgeBaffles
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
