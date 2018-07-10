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

#include "block.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(block, 0);
    defineRunTimeSelectionTable(block, Istream);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::block::block
(
    const dictionary& dict,
    const label index,
    const pointField& vertices,
    const blockEdgeList& edges,
    const blockFaceList& faces,
    Istream& is
)
:
    blockDescriptor(dict, index, vertices, edges, faces, is)
{
    createPoints();
    createBoundary();
}


Foam::block::block(const blockDescriptor& blockDesc)
:
    blockDescriptor(blockDesc)
{
    createPoints();
    createBoundary();
}


Foam::autoPtr<Foam::block> Foam::block::New
(
    const dictionary& dict,
    const label index,
    const pointField& points,
    const blockEdgeList& edges,
    const blockFaceList& faces,
    Istream& is
)
{
    if (debug)
    {
        InfoInFunction << "Constructing block" << endl;
    }

    const word blockOrCellShapeType(is);

    IstreamConstructorTable::iterator cstrIter =
        IstreamConstructorTablePtr_->find(blockOrCellShapeType);

    if (cstrIter == IstreamConstructorTablePtr_->end())
    {
        is.putBack(blockOrCellShapeType);
        return autoPtr<block>(new block(dict, index, points, edges, faces, is));
    }
    else
    {
        return autoPtr<block>
        (
            cstrIter()
            (
                dict,
                index,
                points,
                edges,
                faces,
                is
            )
        );
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const block& b)
{
    os << b.points() << nl
       << b.cells() << nl
       << b.boundaryPatches() << endl;

    return os;
}


// ************************************************************************* //
