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

#include "featureEdgeMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(featureEdgeMesh, 0);

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::featureEdgeMesh::featureEdgeMesh(const IOobject& io)
:
    regIOobject(io),
    edgeMesh(pointField(0), edgeList(0))
{
    if
    (
        io.readOpt() == IOobject::MUST_READ
     || io.readOpt() == IOobject::MUST_READ_IF_MODIFIED
     || (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        readStream(typeName) >> *this;
        close();
    }

    if (debug)
    {
        Pout<< "featureEdgeMesh::featureEdgeMesh :"
            << " constructed from IOobject :"
            << " points:" << points().size()
            << " edges:" << edges().size()
            << endl;
    }
}


//- Construct from components
Foam::featureEdgeMesh::featureEdgeMesh
(
    const IOobject& io,
    const pointField& points,
    const edgeList& edges
)
:
    regIOobject(io),
    edgeMesh(points, edges)
{}


// Construct as copy
Foam::featureEdgeMesh::featureEdgeMesh
(
    const IOobject& io,
    const featureEdgeMesh& em
)
:
    regIOobject(io),
    edgeMesh(em)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::featureEdgeMesh::readData(Istream& is)
{
    is >> *this;
    return !is.bad();
}


bool Foam::featureEdgeMesh::writeData(Ostream& os) const
{
    os << *this;

    return os.good();
}


// ************************************************************************* //
