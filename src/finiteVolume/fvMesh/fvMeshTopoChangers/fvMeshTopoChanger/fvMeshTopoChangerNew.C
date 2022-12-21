/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2022 OpenFOAM Foundation
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

#include "fvMeshTopoChangersNone.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::fvMeshTopoChanger> Foam::fvMeshTopoChanger::New
(
    fvMesh& mesh,
    const dictionary& dict
)
{
    const word fvMeshTopoChangerTypeName
    (
        dict.lookup("type")
    );

    Info<< "Selecting fvMeshTopoChanger "
        << fvMeshTopoChangerTypeName << endl;

    libs.open
    (
        dict,
        "libs",
        fvMeshConstructorTablePtr_
    );

    if (!fvMeshConstructorTablePtr_)
    {
        FatalErrorInFunction
            << "fvMeshTopoChangers table is empty"
            << exit(FatalError);
    }

    fvMeshConstructorTable::iterator cstrIter =
        fvMeshConstructorTablePtr_->find(fvMeshTopoChangerTypeName);

    if (cstrIter == fvMeshConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown fvMeshTopoChanger type "
            << fvMeshTopoChangerTypeName << nl << nl
            << "Valid fvMeshTopoChangers are :" << endl
            << fvMeshConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<fvMeshTopoChanger>
    (
        cstrIter()(mesh, dict)
    );
}


Foam::autoPtr<Foam::fvMeshTopoChanger> Foam::fvMeshTopoChanger::New
(
    fvMesh& mesh
)
{
    typeIOobject<IOdictionary> dictHeader
    (
        IOobject
        (
            "dynamicMeshDict",
            mesh.time().constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        )
    );

    if (dictHeader.headerOk())
    {
        IOdictionary dict(dictHeader);

        if (dict.found("topoChanger"))
        {
            return New(mesh, dict.subDict("topoChanger"));
        }
    }

    return autoPtr<fvMeshTopoChanger>(new fvMeshTopoChangers::none(mesh));
}


// ************************************************************************* //
