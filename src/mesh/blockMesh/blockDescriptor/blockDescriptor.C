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

#include "blockDescriptor.H"
#include "blockMeshTools.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::blockDescriptor::check(const Istream& is)
{
    forAll(blockShape_, pi)
    {
        if (blockShape_[pi] < 0)
        {
            FatalIOErrorInFunction(is)
                << "Negative point label " << blockShape_[pi]
                << " in block " << *this
                << exit(FatalIOError);
        }
        else if (blockShape_[pi] >= vertices_.size())
        {
            FatalIOErrorInFunction(is)
                << "Point label " << blockShape_[pi]
                << " out of range 0.." << vertices_.size() - 1
                << " in block " << *this
                << exit(FatalIOError);
        }
    }

    const point blockCentre(blockShape_.centre(vertices_));
    const faceList faces(blockShape_.faces());

    // Check each face is outward-pointing with respect to the block centre
    label outwardFaceCount = 0;
    boolList correctFaces(faces.size(), true);

    forAll(faces, i)
    {
        point faceCentre(faces[i].centre(vertices_));
        vector faceArea(faces[i].area(vertices_));
        if (mag(faceArea) > small)
        {
            if (((faceCentre - blockCentre) & faceArea) > 0)
            {
                outwardFaceCount++;
            }
            else
            {
                correctFaces[i] = false;
            }
        }
        else
        {
            outwardFaceCount++;
        }
    }

    // If all faces are inward-pointing the block is inside-out
    if (outwardFaceCount == 0)
    {
        FatalIOErrorInFunction(is)
            << "Block " << *this << " is inside-out"
            << exit(FatalIOError);
    }
    else if (outwardFaceCount != faces.size())
    {
        FatalIOErrorInFunction(is)
            << "Block " << *this << " has inward-pointing faces"
            << nl << "    ";

        forAll(correctFaces, i)
        {
            if (!correctFaces[i])
            {
                FatalIOError<< faces[i] << token::SPACE;
            }
        }

        FatalIOError << exit(FatalIOError);
    }
}


void Foam::blockDescriptor::findCurvedFaces()
{
    const faceList blockFaces(blockShape().faces());

    forAll(blockFaces, blockFacei)
    {
        forAll(faces_, facei)
        {
            if
            (
                face::sameVertices
                (
                    faces_[facei].vertices(),
                    blockFaces[blockFacei]
                )
            )
            {
                curvedFaces_[blockFacei] = facei;
                nCurvedFaces_++;
                break;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blockDescriptor::blockDescriptor
(
    const cellShape& bshape,
    const pointField& vertices,
    const blockEdgeList& edges,
    const blockFaceList& faces,
    const Vector<label>& density,
    const UList<gradingDescriptors>& expand,
    const word& zoneName
)
:
    vertices_(vertices),
    edges_(edges),
    faces_(faces),
    blockShape_(bshape),
    density_(density),
    expand_(expand),
    zoneName_(zoneName),
    curvedFaces_(-1),
    nCurvedFaces_(0)
{
    if (expand_.size() != 12)
    {
        FatalErrorInFunction
            << "Unknown definition of expansion ratios"
            << exit(FatalError);
    }

    findCurvedFaces();
}


Foam::blockDescriptor::blockDescriptor
(
    const dictionary& dict,
    const label index,
    const pointField& vertices,
    const blockEdgeList& edges,
    const blockFaceList& faces,
    Istream& is
)
:
    vertices_(vertices),
    edges_(edges),
    faces_(faces),
    density_(),
    expand_(12, gradingDescriptors()),
    zoneName_(),
    curvedFaces_(-1),
    nCurvedFaces_(0)
{
    // Read cell model and list of vertices (potentially with variables)
    word model(is);
    blockShape_ = cellShape
    (
        model,
        blockMeshTools::read<label>
        (
            is,
            dict.subOrEmptyDict("namedVertices")
        )
    );

    // Examine next token
    token t(is);

    // Optional zone name
    if (t.isWord())
    {
        zoneName_ = t.wordToken();

        // Examine next token
        is >> t;
    }
    is.putBack(t);

    if (t.isPunctuation())
    {
        // New-style: read a list of 3 values
        if (t.pToken() == token::BEGIN_LIST)
        {
            is >> density_;
        }
        else
        {
            FatalIOErrorInFunction
            (
                is
            )   << "incorrect token while reading n, expected '(', found "
                << t.info()
                << exit(FatalIOError);
        }
    }
    else
    {
        // Old-style: read three labels
        is  >> density_.x()
            >> density_.y()
            >> density_.z();
    }

    is >> t;
    if (!t.isWord())
    {
        is.putBack(t);
    }

    List<gradingDescriptors> expRatios(is);

    if (expRatios.size() == 1)
    {
        // Identical in x/y/z-directions
        expand_ = expRatios[0];
    }
    else if (expRatios.size() == 3)
    {
        // x-direction
        expand_[0]  = expRatios[0];
        expand_[1]  = expRatios[0];
        expand_[2]  = expRatios[0];
        expand_[3]  = expRatios[0];

        // y-direction
        expand_[4]  = expRatios[1];
        expand_[5]  = expRatios[1];
        expand_[6]  = expRatios[1];
        expand_[7]  = expRatios[1];

        // z-direction
        expand_[8]  = expRatios[2];
        expand_[9]  = expRatios[2];
        expand_[10] = expRatios[2];
        expand_[11] = expRatios[2];
    }
    else if (expRatios.size() == 12)
    {
        expand_ = expRatios;
    }
    else
    {
        FatalIOErrorInFunction(is)
            << "Unknown definition of expansion ratios: " << expRatios
            << exit(FatalIOError);
    }

    check(is);

    findCurvedFaces();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::FixedList<Foam::pointField, 6>
Foam::blockDescriptor::facePoints(const pointField& points) const
{
    // Caches points for curvature correction
    FixedList<pointField, 6> facePoints;

    // Set local variables for mesh specification
    const label ni = density_.x();
    const label nj = density_.y();
    const label nk = density_.z();

    facePoints[0].setSize((nj + 1)*(nk + 1));
    facePoints[1].setSize((nj + 1)*(nk + 1));

    for (label j=0; j<=nj; j++)
    {
        for (label k=0; k<=nk; k++)
        {
            facePoints[0][facePointLabel(0, j, k)] =
                points[pointLabel(0, j, k)];
            facePoints[1][facePointLabel(1, j, k)] =
                points[pointLabel(ni, j, k)];
        }
    }

    facePoints[2].setSize((ni + 1)*(nk + 1));
    facePoints[3].setSize((ni + 1)*(nk + 1));

    for (label i=0; i<=ni; i++)
    {
        for (label k=0; k<=nk; k++)
        {
            facePoints[2][facePointLabel(2, i, k)] =
                points[pointLabel(i, 0, k)];
            facePoints[3][facePointLabel(3, i, k)] =
                points[pointLabel(i, nj, k)];
        }
    }

    facePoints[4].setSize((ni + 1)*(nj + 1));
    facePoints[5].setSize((ni + 1)*(nj + 1));

    for (label i=0; i<=ni; i++)
    {
        for (label j=0; j<=nj; j++)
        {
            facePoints[4][facePointLabel(4, i, j)] =
                points[pointLabel(i, j, 0)];
            facePoints[5][facePointLabel(5, i, j)] =
                points[pointLabel(i, j, nk)];
        }
    }

    return facePoints;
}


void Foam::blockDescriptor::correctFacePoints
(
    FixedList<pointField, 6>& facePoints
) const
{
    forAll(curvedFaces_, blockFacei)
    {
        if (curvedFaces_[blockFacei] != -1)
        {
            faces_[curvedFaces_[blockFacei]].project
            (
                *this,
                blockFacei,
                facePoints[blockFacei]
            );
        }
    }
}


void Foam::blockDescriptor::write
(
    Ostream& os,
    const label val,
    const dictionary& d
)
{
    const dictionary* varDictPtr = d.subDictPtr("namedBlocks");
    if (varDictPtr)
    {
        blockMeshTools::write(os, val, *varDictPtr);
    }
    else
    {
        os << val;
    }
}


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const blockDescriptor& bd)
{
    const cellShape& bshape = bd.blockShape();
    const labelList& blockLabels = bshape;

    os  << bshape.model().name() << " (";

    forAll(blockLabels, labelI)
    {
        if (labelI)
        {
            os  << ' ';
        }
        os  << blockLabels[labelI];
    }
    os  << ')';

    if (bd.zoneName().size())
    {
        os  << ' ' << bd.zoneName();
    }

    os  << ' '  << bd.density()
        << " simpleGrading (";


    const List<gradingDescriptors>& expand = bd.expand_;

    // Can we use a compact notation?
    if
    (
        // x-direction
        (
            expand[0] == expand[1]
         && expand[0] == expand[2]
         && expand[0] == expand[3]
        )
     && // y-direction
        (
            expand[4] == expand[5]
         && expand[4] == expand[6]
         && expand[4] == expand[7]
        )
     && // z-direction
        (
            expand[8] == expand[9]
         && expand[8] == expand[10]
         && expand[8] == expand[11]
        )
    )
    {
        os  << expand[0] << ' ' << expand[4] << ' ' << expand[8];
    }
    else
    {
        forAll(expand, edgei)
        {
            if (edgei)
            {
                os  << ' ';
            }
            os  << expand[edgei];
        }
    }


    os  << ")";

    return os;
}


// ************************************************************************* //
