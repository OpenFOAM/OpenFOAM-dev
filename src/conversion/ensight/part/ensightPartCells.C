/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "ensightPartCells.H"
#include "IOstream.H"
#include "IStringStream.H"
#include "dictionary.H"
#include "cellModeller.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ensightPartCells, 0);
    addToRunTimeSelectionTable(ensightPart, ensightPartCells, istream);
}

const Foam::polyMesh* const Foam::ensightPartCells::polyMeshNullPtr_ = nullptr;

const Foam::List<Foam::word> Foam::ensightPartCells::elemTypes_
(
    IStringStream
    (
        "(tetra4 pyramid5 penta6 hexa8 nfaced)"
    )()
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::ensightPartCells::classify
(
    const polyMesh& mesh,
    const labelUList& idList
)
{
    // References to cell shape models
    const cellModel& tet   = *(cellModeller::lookup("tet"));
    const cellModel& pyr   = *(cellModeller::lookup("pyr"));
    const cellModel& prism = *(cellModeller::lookup("prism"));
    const cellModel& hex   = *(cellModeller::lookup("hex"));

    const cellShapeList& cellShapes = mesh.cellShapes();

    offset_ = 0;
    size_ = mesh.nCells();

    bool limited = false;
    if (notNull(idList))
    {
        limited = true;
        size_ = idList.size();
    }

    // count the shapes
    label nTet   = 0;
    label nPyr   = 0;
    label nPrism = 0;
    label nHex   = 0;
    label nPoly  = 0;

    for (label listI = 0; listI < size_; ++listI)
    {
        label cellId = listI;
        if (limited)
        {
            cellId = idList[listI];
        }

        const cellShape& cellShape = cellShapes[cellId];
        const cellModel& cellModel = cellShape.model();

        if (cellModel == tet)
        {
            nTet++;
        }
        else if (cellModel == pyr)
        {
            nPyr++;
        }
        else if (cellModel == prism)
        {
            nPrism++;
        }
        else if (cellModel == hex)
        {
            nHex++;
        }
        else
        {
            nPoly++;
        }
    }


    // we can avoid double looping, but at the cost of allocation
    labelList tetCells(nTet);
    labelList pyramidCells(nPyr);
    labelList prismCells(nPrism);
    labelList hexCells(nHex);
    labelList polyCells(nPoly);

    nTet   = 0,
    nPyr   = 0;
    nPrism = 0;
    nHex   = 0;
    nPoly  = 0;

    // classify the shapes
    for (label listI = 0; listI < size_; ++listI)
    {
        label cellId = listI;
        if (limited)
        {
            cellId = idList[listI];
        }

        const cellShape& cellShape = cellShapes[cellId];
        const cellModel& cellModel = cellShape.model();

        if (cellModel == tet)
        {
            tetCells[nTet++] = cellId;
        }
        else if (cellModel == pyr)
        {
            pyramidCells[nPyr++] = cellId;
        }
        else if (cellModel == prism)
        {
            prismCells[nPrism++] = cellId;
        }
        else if (cellModel == hex)
        {
            hexCells[nHex++] = cellId;
        }
        else
        {
            polyCells[nPoly++] = cellId;
        }
    }


    // MUST match with elementTypes
    elemLists_.setSize(elementTypes().size());

    elemLists_[tetra4Elements].transfer(tetCells);
    elemLists_[pyramid5Elements].transfer(pyramidCells);
    elemLists_[penta6Elements].transfer(prismCells);
    elemLists_[hexa8Elements].transfer(hexCells);
    elemLists_[nfacedElements].transfer(polyCells);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ensightPartCells::ensightPartCells
(
    label partNumber,
    const string& partDescription
)
:
    ensightPart(partNumber, partDescription),
    mesh_(*polyMeshNullPtr_)
{}


Foam::ensightPartCells::ensightPartCells
(
    label partNumber,
    const polyMesh& mesh
)
:
    ensightPart(partNumber, "cells", mesh.points()),
    mesh_(mesh)
{
    classify(mesh);
}


Foam::ensightPartCells::ensightPartCells
(
    label partNumber,
    const polyMesh& mesh,
    const labelUList& idList
)
:
    ensightPart(partNumber, "cells", mesh.points()),
    mesh_(mesh)
{
    classify(mesh, idList);
}


Foam::ensightPartCells::ensightPartCells
(
    label partNumber,
    const polyMesh& mesh,
    const cellZone& cZone
)
:
    ensightPart(partNumber, cZone.name(), mesh.points()),
    mesh_(mesh)
{
    classify(mesh, cZone);
}


Foam::ensightPartCells::ensightPartCells(const ensightPartCells& part)
:
    ensightPart(part),
    mesh_(part.mesh_)
{}


Foam::ensightPartCells::ensightPartCells(Istream& is)
:
    ensightPart(),
    mesh_(*polyMeshNullPtr_)
{
    reconstruct(is);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ensightPartCells::~ensightPartCells()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::ensightPart::localPoints Foam::ensightPartCells::calcLocalPoints() const
{
    localPoints ptList(points_);
    labelList& usedPoints = ptList.list;
    label nPoints = 0;

    forAll(elemLists_, typeI)
    {
        const labelUList& idList = elemLists_[typeI];

        // add all points from cells
        forAll(idList, i)
        {
            const label id = idList[i] + offset_;
            const labelUList& cFaces = mesh_.cells()[id];

            forAll(cFaces, cFacei)
            {
                const face& f = mesh_.faces()[cFaces[cFacei]];

                forAll(f, fp)
                {
                    if (usedPoints[f[fp]] == -1)
                    {
                        usedPoints[f[fp]] = nPoints++;
                    }
                }
            }
        }
    }

    // this is not absolutely necessary, but renumber anyhow
    nPoints = 0;
    forAll(usedPoints, ptI)
    {
        if (usedPoints[ptI] > -1)
        {
            usedPoints[ptI] = nPoints++;
        }
    }

    ptList.nPoints = nPoints;
    return ptList;
}


void Foam::ensightPartCells::writeConnectivity
(
    ensightGeoFile& os,
    const word& key,
    const labelUList& idList,
    const labelUList& pointMap
) const
{
    os.writeKeyword(key);
    os.write(idList.size());
    os.newline();

    // write polyhedral
    if (key == "nfaced")
    {
        const faceList& meshFaces = mesh_.faces();
        const labelUList& owner = mesh_.faceOwner();

        // write the number of faces per element
        forAll(idList, i)
        {
            const label id = idList[i] + offset_;
            const labelUList& cFace = mesh_.cells()[id];

            os.write(cFace.size());
            os.newline();
        }

        // write the number of points per element face
        forAll(idList, i)
        {
            const label id = idList[i] + offset_;
            const labelUList& cFace = mesh_.cells()[id];

            forAll(cFace, facei)
            {
                const face& cf = meshFaces[cFace[facei]];

                os.write(cf.size());
                os.newline();
            }
        }

        // write the points describing each element face
        forAll(idList, i)
        {
            const label id = idList[i] + offset_;
            const labelUList& cFace = mesh_.cells()[id];

            forAll(cFace, cFacei)
            {
                const label faceId = cFace[cFacei];
                const face& cf = meshFaces[faceId];

                // convert global -> local index
                // (note: Ensight indices start with 1)

                // ensight >= 9 needs consistently oriented nfaced cells
                if (id == owner[faceId])
                {
                    forAll(cf, ptI)
                    {
                        os.write(pointMap[cf[ptI]] + 1);
                    }
                }
                else
                {
                    // as per face::reverseFace(), but without copying

                    os.write(pointMap[cf[0]] + 1);
                    for (label ptI = cf.size()-1; ptI > 0; --ptI)
                    {
                        os.write(pointMap[cf[ptI]] + 1);
                    }
                }

                os.newline();
            }
        }
    }
    else
    {
        // write primitive
        const cellShapeList& cellShapes = mesh_.cellShapes();

        forAll(idList, i)
        {
            const label id = idList[i] + offset_;
            const cellShape& cellPoints = cellShapes[id];

            // convert global -> local index
            // (note: Ensight indices start with 1)
            forAll(cellPoints, ptI)
            {
                os.write(pointMap[cellPoints[ptI]] + 1);
            }
            os.newline();
        }
    }
}


void Foam::ensightPartCells::writeGeometry(ensightGeoFile& os) const
{
    ensightPart::writeGeometry(os, points_);
}


// ************************************************************************* //
