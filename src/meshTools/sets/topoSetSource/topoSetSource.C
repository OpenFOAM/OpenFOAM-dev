/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "topoSetSource.H"
#include "polyMesh.H"
#include "topoSet.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(topoSetSource, 0);
    defineRunTimeSelectionTable(topoSetSource, word);

    template<>
    const char* Foam::NamedEnum
    <
        Foam::topoSetSource::setAction,
        8
    >::names[] =
    {
        "clear",
        "new",
        "invert",
        "add",
        "delete",
        "subset",
        "list",
        "remove"
    };
}


const Foam::NamedEnum<Foam::topoSetSource::setAction, 8>
    Foam::topoSetSource::actionNames_;


const Foam::string Foam::topoSetSource::illegalSource_
(
    "Illegal topoSetSource name"
);



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::topoSetSource> Foam::topoSetSource::New
(
    const word& topoSetSourceType,
    const polyMesh& mesh,
    const dictionary& dict
)
{
    wordConstructorTable::iterator cstrIter =
        wordConstructorTablePtr_->find(topoSetSourceType);

    if (cstrIter == wordConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown topoSetSource type " << topoSetSourceType
            << endl << endl
            << "Valid topoSetSource types : " << endl
            << wordConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<topoSetSource>(cstrIter()(mesh, dict));
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::topoSetSource::addOrDelete
(
    topoSet& set,
    const label celli,
    const bool add
) const
{
    if (add)
    {
        set.insert(celli);
    }
    else
    {
        set.erase(celli);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::topoSetSource::topoSetSource(const polyMesh& mesh)
:
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::topoSetSource::~topoSetSource()
{}


// ************************************************************************* //
