/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2024 OpenFOAM Foundation
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

#include "list_fvMeshTopoChanger.H"
#include "polyTopoChangeMap.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshTopoChangers
{
    defineTypeNameAndDebug(list, 0);
    addToRunTimeSelectionTable(fvMeshTopoChanger, list, fvMesh);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshTopoChangers::list::list(fvMesh& mesh, const dictionary& dict)
:
    fvMeshTopoChanger(mesh)
{
    const dictionary& solversDict = dict.subDict("topoChangers");

    forAllConstIter(dictionary, solversDict, iter)
    {
        if (iter().isDict())
        {
            const word& name = iter().keyword();
            const dictionary& dict = iter().dict();

            list_.insert
            (
                name,
                fvMeshTopoChanger::New(mesh, dict).ptr()
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMeshTopoChangers::list::~list()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fvMeshTopoChangers::list::update()
{
    bool updated = false;

    forAllIter(PtrDictionary<fvMeshTopoChanger>, list_, iter)
    {
        updated = iter().update() || updated;
        mesh().topoChanged_ = updated;
    }

    return updated;
}


void Foam::fvMeshTopoChangers::list::topoChange(const polyTopoChangeMap& map)
{
    forAllIter(PtrDictionary<fvMeshTopoChanger>, list_, iter)
    {
        iter().topoChange(map);
    }
}


void Foam::fvMeshTopoChangers::list::mapMesh(const polyMeshMap& map)
{
    forAllIter(PtrDictionary<fvMeshTopoChanger>, list_, iter)
    {
        iter().mapMesh(map);
    }
}


void Foam::fvMeshTopoChangers::list::distribute
(
    const polyDistributionMap& map
)
{
    forAllIter(PtrDictionary<fvMeshTopoChanger>, list_, iter)
    {
        iter().distribute(map);
    }
}


// ************************************************************************* //
