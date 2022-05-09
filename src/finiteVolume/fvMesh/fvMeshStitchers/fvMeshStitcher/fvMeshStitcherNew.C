/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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

#include "fvMeshStitcher.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::fvMeshStitcher> Foam::fvMeshStitcher::New
(
    fvMesh& mesh,
    const bool changing
)
{
    typeIOobject<IOdictionary> dictHeader
    (
        IOobject
        (
            "dynamicMeshDict",
            mesh.time().constant(),
            mesh.dbDir(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        )
    );

    if (changing && dictHeader.headerOk())
    {
        IOdictionary dict(dictHeader);

        forAllConstIter(dictionary, dict, iter)
        {
            if (!iter().isDict()) continue;

            const dictionary& changerDict = iter().dict();

            if (!changerDict.found("libs")) continue;

            libs.open
            (
                changerDict,
                "libs",
                fvMeshConstructorTablePtr_
            );
        }
    }

    forAllConstIter
    (
        fvMeshConstructorTable,
        *fvMeshConstructorTablePtr_,
        cstrIter
    )
    {
        autoPtr<fvMeshStitcher> stitcherPtr(cstrIter()(mesh));

        if (stitcherPtr->changing() == (changing && dictHeader.headerOk()))
        {
            return stitcherPtr;
        }
    }

    FatalErrorInFunction
        << typeName << " for " << (changing ? "" : "non-")
        << "changing mesh not found " << nl << nl
        << "Valid " << typeName << "s are : " << endl
        << fvMeshConstructorTablePtr_->toc()
        << exit(FatalError);

    return autoPtr<fvMeshStitcher>(nullptr);
}


// ************************************************************************* //
