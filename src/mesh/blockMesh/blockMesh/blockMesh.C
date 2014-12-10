/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "blockMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

bool Foam::blockMesh::blockMesh::verboseOutput(false);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blockMesh::blockMesh(const IOdictionary& dict, const word& regionName)
:
    blockPointField_(dict.lookup("vertices")),
    scaleFactor_(1.0),
    topologyPtr_(createTopology(dict, regionName))
{
    calcMergeInfo();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::blockMesh::~blockMesh()
{
    delete topologyPtr_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::blockMesh::verbose(const bool on)
{
    verboseOutput = on;
}


const Foam::pointField& Foam::blockMesh::blockPointField() const
{
    return blockPointField_;
}


const Foam::polyMesh& Foam::blockMesh::topology() const
{
    if (!topologyPtr_)
    {
        FatalErrorIn("blockMesh::topology() const")
            << "topologyPtr_ not allocated"
            << exit(FatalError);
    }

    return *topologyPtr_;
}


Foam::PtrList<Foam::dictionary> Foam::blockMesh::patchDicts() const
{
    const polyPatchList& patchTopologies = topology().boundaryMesh();

    PtrList<dictionary> patchDicts(patchTopologies.size());

    forAll(patchTopologies, patchI)
    {
        OStringStream os;
        patchTopologies[patchI].write(os);
        IStringStream is(os.str());
        patchDicts.set(patchI, new dictionary(is));
    }
    return patchDicts;
}


Foam::scalar Foam::blockMesh::scaleFactor() const
{
    return scaleFactor_;
}


const Foam::pointField& Foam::blockMesh::points() const
{
    if (points_.empty())
    {
        createPoints();
    }

    return points_;
}


const Foam::cellShapeList& Foam::blockMesh::cells() const
{
    if (cells_.empty())
    {
        createCells();
    }

    return cells_;
}


const Foam::faceListList& Foam::blockMesh::patches() const
{
    if (patches_.empty())
    {
        createPatches();
    }

    return patches_;
}


Foam::wordList Foam::blockMesh::patchNames() const
{
    return topology().boundaryMesh().names();
}


//Foam::wordList Foam::blockMesh::patchTypes() const
//{
//    return topology().boundaryMesh().types();
//}
//
//
//Foam::wordList Foam::blockMesh::patchPhysicalTypes() const
//{
//    return topology().boundaryMesh().physicalTypes();
//}


Foam::label Foam::blockMesh::numZonedBlocks() const
{
    label num = 0;

    forAll(*this, blockI)
    {
        if (operator[](blockI).zoneName().size())
        {
            num++;
        }
    }

    return num;
}


void Foam::blockMesh::writeTopology(Ostream& os) const
{
    const pointField& pts = topology().points();

    forAll(pts, pI)
    {
        const point& pt = pts[pI];

        os << "v " << pt.x() << ' ' << pt.y() << ' ' << pt.z() << endl;
    }

    const edgeList& edges = topology().edges();

    forAll(edges, eI)
    {
        const edge& e = edges[eI];

        os << "l " << e.start() + 1 << ' ' << e.end() + 1 << endl;
    }
}

// ************************************************************************* //
