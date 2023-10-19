/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "ensightParts.H"
#include "processorPolyPatch.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ensightParts::ensightParts(const polyMesh& mesh)
:
    partsList_()
{
    recalculate(mesh);
}


Foam::ensightParts::ensightParts(const IOobject& ioObj)
:
    partsList_()
{
    IOPtrList<ensightPart> ioList(ioObj);
    partsList_.transfer(ioList);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ensightParts::~ensightParts()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ensightParts::recalculate(const polyMesh& mesh)
{
    partsList_.clear();

    // extra space for unzoned cells
    label nPart =
    (
        mesh.cellZones().size()
      + mesh.boundaryMesh().size()
      + 1
    );

    partsList_.setSize(nPart);
    nPart = 0;

    label nZoneCells = 0;

    // do cell zones
    forAll(mesh.cellZones(), zoneI)
    {
        const cellZone& cZone = mesh.cellZones()[zoneI];
        nZoneCells += cZone.size();

        if (cZone.size())
        {
            partsList_.set
            (
                nPart,
                new ensightPartCells(nPart, mesh, cZone)
            );

            nPart++;
        }
    }

    // collect unzoned cells

    // special case: no zones at all - do entire mesh
    if (nZoneCells == 0)
    {
        partsList_.set
        (
            nPart,
            new ensightPartCells(nPart, mesh)
        );

        nPart++;
    }
    else if (mesh.nCells() > nZoneCells)
    {
        // determine which cells are not in a cellZone
        labelList unzoned(mesh.nCells(), -1);

        forAll(mesh.cellZones(), zoneI)
        {
            const labelUList& idList = mesh.cellZones()[zoneI];

            forAll(idList, i)
            {
                unzoned[idList[i]] = idList[i];
            }
        }

        label nUnzoned = 0;
        forAll(unzoned, i)
        {
            if (unzoned[i] < 0)
            {
                unzoned[nUnzoned] = i;
                nUnzoned++;
            }
        }
        unzoned.setSize(nUnzoned);

        if (unzoned.size())
        {
            partsList_.set
            (
                nPart,
                new ensightPartCells(nPart, mesh, unzoned)
            );

            nPart++;
        }
    }


    // do boundaries, skipping empty and processor patches
    forAll(mesh.boundaryMesh(), patchi)
    {
        const polyPatch& patch = mesh.boundaryMesh()[patchi];
        if (patch.size() && !isA<processorPolyPatch>(patch))
        {
            partsList_.set
            (
                nPart,
                new ensightPartFaces(nPart, mesh, patch)
            );

            nPart++;
        }
    }

    // truncate to correct size
    partsList_.setSize(nPart);
}


void Foam::ensightParts::renumber
(
    const labelUList& origCellId,
    const labelUList& origFaceId
)
{
    forAll(partsList_, partI)
    {
        if (partsList_[partI].isCellData())
        {
            partsList_[partI].renumber(origCellId);
        }
        else
        {
            partsList_[partI].renumber(origFaceId);
        }
    }
}


void Foam::ensightParts::writeGeometry(ensightGeoFile& os) const
{
    // with some feedback
    Info<< "write geometry part:" << nl << flush;

    forAll(partsList_, partI)
    {
        Info<< " " << partI << flush;
        partsList_[partI].writeGeometry(os);
    }
}


bool Foam::ensightParts::writeSummary(Ostream& os) const
{
    forAll(partsList_, partI)
    {
        partsList_[partI].writeSummary(os);
    }

    return true;
}


void Foam::ensightParts::writeData(Ostream& os) const
{
    // Begin write list
    os  << nl << partsList_.size()
        << nl << token::BEGIN_LIST;

    // Write list contents
    forAll(partsList_, i)
    {
        os  << nl << partsList_[i];
    }

    // End write list
    os  << nl << token::END_LIST << nl;

    // Check state of IOstream
    os.check("ensightParts::writeData(Ostream&)");
}


void Foam::ensightParts::writeScalarField
(
    ensightFile& os,
    const List<scalar>& field,
    const bool useFaceData,
    const bool perNode
) const
{
    forAll(partsList_, partI)
    {
        if
        (
            useFaceData
          ? partsList_[partI].isFaceData()
          : partsList_[partI].isCellData()
        )
        {
            partsList_[partI].writeScalarField(os, field, perNode);
        }
    }
}


void Foam::ensightParts::writeVectorField
(
    ensightFile& os,
    const List<scalar>& field0,
    const List<scalar>& field1,
    const List<scalar>& field2,
    const bool useFaceData,
    const bool perNode
) const
{
    forAll(partsList_, partI)
    {
        if
        (
            useFaceData
          ? partsList_[partI].isFaceData()
          : partsList_[partI].isCellData()
        )
        {
            partsList_[partI].writeVectorField
            (
                os,
                field0, field1, field2,
                perNode
            );
        }
    }
}


// * * * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * //

Foam::ensightGeoFile& Foam::operator<<
(
    ensightGeoFile& os,
    const ensightParts& parts
)
{
    parts.writeGeometry(os);
    return os;
}


// ************************************************************************* //
