/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "dynamicFvMesh.H"
#include "Time.H"
#include "dlLibraryTable.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::dynamicFvMesh> Foam::dynamicFvMesh::New(const IOobject& io)
{
    // Note: - do not register the dictionary since dynamicFvMeshes themselves
    // do this.
    // - defaultRegion (region0) gets loaded from constant, other ones
    //   get loaded from constant/<regionname>. Normally we'd use
    //   polyMesh::dbDir() but we haven't got a polyMesh yet ...
    IOdictionary dict
    (
        IOobject
        (
            "dynamicMeshDict",
            io.time().constant(),
            (io.name() == polyMesh::defaultRegion ? "" : io.name()),
            io.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false
        )
    );

    const word dynamicFvMeshTypeName(dict.lookup("dynamicFvMesh"));

    Info<< "Selecting dynamicFvMesh " << dynamicFvMeshTypeName << endl;

    const_cast<Time&>(io.time()).libs().open
    (
        dict,
        "dynamicFvMeshLibs",
        IOobjectConstructorTablePtr_
    );

    if (!IOobjectConstructorTablePtr_)
    {
        FatalErrorIn
        (
            "dynamicFvMesh::New(const IOobject&)"
        )   << "dynamicFvMesh table is empty"
            << exit(FatalError);
    }

    IOobjectConstructorTable::iterator cstrIter =
        IOobjectConstructorTablePtr_->find(dynamicFvMeshTypeName);

    if (cstrIter == IOobjectConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "dynamicFvMesh::New(const IOobject&)"
        )   << "Unknown dynamicFvMesh type "
            << dynamicFvMeshTypeName << nl << nl
            << "Valid dynamicFvMesh types are :" << endl
            << IOobjectConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<dynamicFvMesh>(cstrIter()(io));
}


// ************************************************************************* //
