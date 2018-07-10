/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2018 OpenFOAM Foundation
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

#include "decompositionModel.H"
#include "polyMesh.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(decompositionModel, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::decompositionModel::decompositionModel
(
    const polyMesh& mesh,
    const fileName& decompDictFile
)
:
    MeshObject
    <
        polyMesh,
        Foam::UpdateableMeshObject,
        decompositionModel
    >(mesh),
    IOdictionary
    (
        selectIO
        (
            IOobject
            (
                "decomposeParDict",
                mesh.time().system(),
                mesh.local(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false   // io.registerObject()
            ),
            decompDictFile
        )
    )
{}


Foam::decompositionModel::decompositionModel
(
    const polyMesh& mesh,
    const dictionary& dict,
    const fileName& decompDictFile
)
:
    MeshObject
    <
        polyMesh,
        Foam::UpdateableMeshObject,
        decompositionModel
    >(mesh),
    IOdictionary
    (
        selectIO
        (
            IOobject
            (
                "decomposeParDict",
                mesh.time().system(),
                mesh.local(),
                mesh,
                (dict.size() ? IOobject::NO_READ : IOobject::MUST_READ),
                IOobject::NO_WRITE,
                false   // io.registerObject()
            ),
            decompDictFile
        ),
        dict
    )
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

const Foam::decompositionModel& Foam::decompositionModel::New
(
    const polyMesh& mesh,
    const fileName& decompDictFile
)
{
    return
        MeshObject
        <
            polyMesh,
            Foam::UpdateableMeshObject,
            decompositionModel
        >::New(mesh, decompDictFile);
}


const Foam::decompositionModel& Foam::decompositionModel::New
(
    const polyMesh& mesh,
    const dictionary& dict,
    const fileName& decompDictFile
)
{
    return
        MeshObject
        <
            polyMesh,
            Foam::UpdateableMeshObject,
            decompositionModel
        >::New(mesh, dict, decompDictFile);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::IOobject Foam::decompositionModel::selectIO
(
    const IOobject& io,
    const fileName& f
)
{
    return
    (
        f.size()
      ? IOobject        // construct from filePath instead
        (
            fileName(f).toAbsolute(),
            io.db(),
            io.readOpt(),
            io.writeOpt(),
            io.registerObject()
        )
     :  io
    );
}


// ************************************************************************* //
