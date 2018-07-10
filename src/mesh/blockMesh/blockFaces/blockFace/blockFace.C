/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
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

#include "blockFace.H"
#include "blockMeshTools.H"
#include "blockVertex.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(blockFace, 0);
    defineRunTimeSelectionTable(blockFace, Istream);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blockFace::blockFace(const face& vertices)
:
    vertices_(vertices)
{}


Foam::blockFace::blockFace
(
    const dictionary& dict,
    const label index,
    Istream& is
)
:
    vertices_
    (
        blockMeshTools::read<label>
        (
            is,
            dict.subOrEmptyDict("namedVertices")
        )
    )
{}


Foam::autoPtr<Foam::blockFace> Foam::blockFace::clone() const
{
    NotImplemented;
    return autoPtr<blockFace>(nullptr);
}


Foam::autoPtr<Foam::blockFace> Foam::blockFace::New
(
    const dictionary& dict,
    const label index,
    const searchableSurfaces& geometry,
    Istream& is
)
{
    if (debug)
    {
        InfoInFunction << "Constructing blockFace" << endl;
    }

    const word faceType(is);

    IstreamConstructorTable::iterator cstrIter =
        IstreamConstructorTablePtr_->find(faceType);

    if (cstrIter == IstreamConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown blockFace type "
            << faceType << nl << nl
            << "Valid blockFace types are" << endl
            << IstreamConstructorTablePtr_->sortedToc()
            << abort(FatalError);
    }

    return autoPtr<blockFace>(cstrIter()(dict, index, geometry, is));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::blockFace::write(Ostream& os, const dictionary& d) const
{
    // Write size and start delimiter
    os << vertices_.size() << token::BEGIN_LIST;

    // Write contents
    forAll(vertices_, i)
    {
        if (i > 0) os << token::SPACE;
        blockVertex::write(os, vertices_[i], d);
    }

    // Write end delimiter
    os << token::END_LIST;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const blockFace& p)
{
    os << p.vertices_ << endl;

    return os;
}


// ************************************************************************* //
