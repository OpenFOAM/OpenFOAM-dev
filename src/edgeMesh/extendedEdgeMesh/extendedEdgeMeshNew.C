/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#include "extendedEdgeMesh.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineRunTimeSelectionTable(extendedEdgeMesh, fileExtension);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::extendedEdgeMesh> Foam::extendedEdgeMesh::New
(
    const fileName& name,
    const word& ext
)
{
    fileExtensionConstructorTable::iterator cstrIter =
        fileExtensionConstructorTablePtr_->find(ext);

    if (cstrIter == fileExtensionConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "extendedEdgeMesh::New(const fileName&, const word&) : "
            "constructing extendedEdgeMesh"
        )   << "Unknown file extension " << ext
            << " for file " << name << nl << nl
            << "Valid extensions are :" << nl
            << fileExtensionConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<extendedEdgeMesh>(cstrIter()(name));
}


Foam::autoPtr<Foam::extendedEdgeMesh> Foam::extendedEdgeMesh::New
(
    const fileName& name
)
{
    word ext = name.ext();
    if (ext == "gz")
    {
        ext = name.lessExt().ext();
    }
    return New(name, ext);
}


// ************************************************************************* //
