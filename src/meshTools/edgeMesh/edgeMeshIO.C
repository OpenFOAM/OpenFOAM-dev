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

#include "edgeMesh.H"
#include "boundBox.H"
#include "edgeMeshFormat.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::edgeMesh::edgeMesh
(
    const fileName& name,
    const word& ext
)
:
    points_(0),
    edges_(0),
    pointEdgesPtr_(nullptr)
{
    read(name, ext);
}


Foam::edgeMesh::edgeMesh(const fileName& name)
:
    points_(0),
    edges_(0),
    pointEdgesPtr_(nullptr)
{
    read(name);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::edgeMesh::read(const fileName& name)
{
    word ext = name.ext();
    if (ext == "gz")
    {
        fileName unzipName = name.lessExt();
        return read(unzipName, unzipName.ext());
    }
    else
    {
        return read(name, ext);
    }
}


// Read from file in given format
bool Foam::edgeMesh::read
(
    const fileName& name,
    const word& ext
)
{
    // read via selector mechanism
    transfer(New(name, ext)());
    return true;
}


void Foam::edgeMesh::write
(
    const fileName& name,
    const edgeMesh& mesh
)
{
    if (debug)
    {
        InfoInFunction << "Writing to " << name << endl;
    }

    const word ext = name.ext();

    writefileExtensionMemberFunctionTable::iterator mfIter =
        writefileExtensionMemberFunctionTablePtr_->find(ext);

    if (mfIter == writefileExtensionMemberFunctionTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown file extension " << ext << nl << nl
            << "Valid types are :" << endl
            << writefileExtensionMemberFunctionTablePtr_->sortedToc()
            << exit(FatalError);
    }
    else
    {
        mfIter()(name, mesh);
    }
}


void Foam::edgeMesh::writeStats(Ostream& os) const
{
    os  << indent << "points      : " << points().size() << nl;
    os  << indent << "edges       : " << edges().size() << nl;
    os  << indent << "boundingBox : " << boundBox(this->points()) << endl;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const edgeMesh& em)
{
    fileFormats::edgeMeshFormat::write(os, em.points_, em.edges_);

    // Check state of Ostream
    os.check("Ostream& operator<<(Ostream&, const edgeMesh&)");

    return os;
}


Foam::Istream& Foam::operator>>(Istream& is, edgeMesh& em)
{
    fileFormats::edgeMeshFormat::read(is, em.points_, em.edges_);

    em.pointEdgesPtr_.clear();

    // Check state of Istream
    is.check("Istream& operator>>(Istream&, edgeMesh&)");

    return is;
}


// ************************************************************************* //
