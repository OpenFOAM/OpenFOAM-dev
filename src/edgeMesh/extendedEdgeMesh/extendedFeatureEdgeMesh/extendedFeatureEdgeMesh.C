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

#include "extendedFeatureEdgeMesh.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(extendedFeatureEdgeMesh, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::extendedFeatureEdgeMesh::extendedFeatureEdgeMesh(const IOobject& io)
:
    regIOobject(io),
    extendedEdgeMesh()
{
    if
    (
        io.readOpt() == IOobject::MUST_READ
     || io.readOpt() == IOobject::MUST_READ_IF_MODIFIED
     || (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        if (readOpt() == IOobject::MUST_READ_IF_MODIFIED)
        {
            WarningIn
            (
                "extendedFeatureEdgeMesh::extendedFeatureEdgeMesh"
                "(const IOobject&)"
            )   << "Specified IOobject::MUST_READ_IF_MODIFIED but class"
                << " does not support automatic rereading."
                << endl;
        }

        readStream(typeName) >> *this;
        close();

        {
            // Calculate edgeDirections

            const edgeList& eds(edges());

            const pointField& pts(points());

            edgeDirections_.setSize(eds.size());

            forAll(eds, eI)
            {
                edgeDirections_[eI] = eds[eI].vec(pts);
            }

            edgeDirections_ /= (mag(edgeDirections_) + SMALL);
        }
    }

    if (debug)
    {
        Pout<< "extendedFeatureEdgeMesh::extendedFeatureEdgeMesh :"
            << " constructed from IOobject :"
            << " points:" << points().size()
            << " edges:" << edges().size()
            << endl;
    }
}


Foam::extendedFeatureEdgeMesh::extendedFeatureEdgeMesh
(
    const IOobject& io,
    const extendedEdgeMesh& em
)
:
    regIOobject(io),
    extendedEdgeMesh(em)
{}


Foam::extendedFeatureEdgeMesh::extendedFeatureEdgeMesh
(
    const surfaceFeatures& sFeat,
    const objectRegistry& obr,
    const fileName& sFeatFileName,
    const boolList& surfBaffleRegions
)
:
    regIOobject
    (
        IOobject
        (
            sFeatFileName,
            obr.time().constant(),
            "extendedFeatureEdgeMesh",
            obr,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    extendedEdgeMesh(sFeat, surfBaffleRegions)
{}


Foam::extendedFeatureEdgeMesh::extendedFeatureEdgeMesh
(
    const IOobject& io,
    const PrimitivePatch<face, List, pointField, point>& surf,
    const labelList& featureEdges,
    const labelList& regionFeatureEdges,
    const labelList& featurePoints
)
:
    regIOobject(io),
    extendedEdgeMesh(surf, featureEdges, regionFeatureEdges, featurePoints)
{}


Foam::extendedFeatureEdgeMesh::extendedFeatureEdgeMesh
(
    const IOobject& io,
    const pointField& pts,
    const edgeList& eds,
    label concaveStart,
    label mixedStart,
    label nonFeatureStart,
    label internalStart,
    label flatStart,
    label openStart,
    label multipleStart,
    const vectorField& normals,
    const List<sideVolumeType>& normalVolumeTypes,
    const vectorField& edgeDirections,
    const labelListList& normalDirections,
    const labelListList& edgeNormals,
    const labelListList& featurePointNormals,
    const labelListList& featurePointEdges,
    const labelList& regionEdges
)
:
    regIOobject(io),
    extendedEdgeMesh
    (
        pts,
        eds,
        concaveStart,
        mixedStart,
        nonFeatureStart,
        internalStart,
        flatStart,
        openStart,
        multipleStart,
        normals,
        normalVolumeTypes,
        edgeDirections,
        normalDirections,
        edgeNormals,
        featurePointNormals,
        featurePointEdges,
        regionEdges
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::extendedFeatureEdgeMesh::~extendedFeatureEdgeMesh()
{}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::extendedFeatureEdgeMesh::readData(Istream& is)
{
    is >> *this;
    return !is.bad();
}


bool Foam::extendedFeatureEdgeMesh::writeData(Ostream& os) const
{
    os << *this;

    return os.good();
}



//bool Foam::extendedFeatureEdgeMesh::writeData(Ostream& os) const
//{
//    os  << "// points" << nl
//        << points() << nl
//        << "// edges" << nl
//        << edges() << nl
//        << "// concaveStart mixedStart nonFeatureStart" << nl
//        << concaveStart_ << token::SPACE
//        << mixedStart_ << token::SPACE
//        << nonFeatureStart_ << nl
//        << "// internalStart flatStart openStart multipleStart" << nl
//        << internalStart_ << token::SPACE
//        << flatStart_ << token::SPACE
//        << openStart_ << token::SPACE
//        << multipleStart_ << nl
//        << "// normals" << nl
//        << normals_ << nl
//        << "// normal volume types" << nl
//        << normalVolumeTypes_ << nl
//        << "// normalDirections" << nl
//        << normalDirections_ << nl
//        << "// edgeNormals" << nl
//        << edgeNormals_ << nl
//        << "// featurePointNormals" << nl
//        << featurePointNormals_ << nl
//        << "// featurePointEdges" << nl
//        << featurePointEdges_ << nl
//        << "// regionEdges" << nl
//        << regionEdges_
//        << endl;
//
//    return os.good();
//}

//
//Foam::Istream& Foam::operator>>
//(
//    Istream& is,
//    Foam::extendedFeatureEdgeMesh::sideVolumeType& vt
//)
//{
//    label type;
//    is  >> type;
//
//    vt = static_cast<Foam::extendedFeatureEdgeMesh::sideVolumeType>(type);
//
//    // Check state of Istream
//    is.check("operator>>(Istream&, sideVolumeType&)");
//
//    return is;
//}
//
//
//Foam::Ostream& Foam::operator<<
//(
//    Ostream& os,
//    const Foam::extendedFeatureEdgeMesh::sideVolumeType& vt
//)
//{
//    os  << static_cast<label>(vt);
//
//    return os;
//}
//

// ************************************************************************* //
