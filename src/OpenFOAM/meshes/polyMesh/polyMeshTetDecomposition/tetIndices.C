/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "tetIndices.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::tetIndices::tetIndices()
:
    celli_(-1),
    facei_(-1),
    faceBasePtI_(-1),
    facePtAI_(-1),
    facePtBI_(-1),
    tetPti_(-1)
{}


Foam::tetIndices::tetIndices
(
    label celli,
    label facei,
    label faceBasePtI,
    label facePtAI,
    label facePtBI,
    label tetPtI
)
:
    celli_(celli),
    facei_(facei),
    faceBasePtI_(faceBasePtI),
    facePtAI_(facePtAI),
    facePtBI_(facePtBI),
    tetPti_(tetPtI)
{}


Foam::tetIndices::tetIndices
(
    label celli,
    label facei,
    label tetPtI,
    const polyMesh& mesh
)
:
    celli_(celli),
    facei_(facei),
    faceBasePtI_(-1),
    facePtAI_(-1),
    facePtBI_(-1),
    tetPti_(tetPtI)
{
    const faceList& pFaces = mesh.faces();
    const labelList& pOwner = mesh.faceOwner();

    const Foam::face& f = pFaces[facei_];

    bool own = (pOwner[facei_] == celli_);

    faceBasePtI_ = mesh.tetBasePtIs()[facei_];

    label facePtI = (tetPti_ + faceBasePtI_) % f.size();
    label otherFacePtI = f.fcIndex(facePtI);

    if (own)
    {
        facePtAI_ = facePtI;
        facePtBI_ = otherFacePtI;
    }
    else
    {
        facePtAI_ = otherFacePtI;
        facePtBI_ = facePtI;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::tetIndices::~tetIndices()
{}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, tetIndices& tI)
{
    is  >> tI.cell()
        >> tI.face()
        >> tI.faceBasePt()
        >> tI.facePtA()
        >> tI.facePtB()
        >> tI.tetPt();

    // Check state of Istream
    is.check
    (
        "Foam::Istream& Foam::operator>>(Foam::Istream&, Foam::tetIndices&)"
    );

    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const tetIndices& tI)
{
    os  << tI.cell() << token::SPACE
        << tI.face() << token::SPACE
        << tI.faceBasePt() << token::SPACE
        << tI.facePtA() << token::SPACE
        << tI.facePtB() << token::SPACE
        << tI.tetPt() << token::SPACE
        << endl;

    // Check state of Ostream
    os.check
    (
        "Foam::Ostream& Foam::operator<<(Foam::Ostream&, "
        "const Foam::tetIndices&)"
    );

    return os;
}


// ************************************************************************* //
