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
    bool changing
)
{
    // Determine if the mesh is actually changing and load the
    // fvMeshStitchers library if so
    if (changing)
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

        changing = changing && dictHeader.headerOk();

        if (changing)
        {
            libs.open("lib" + typeName + "s.so");
        }
    }

    // Select a non-changing or changing stitcher as appropriate
    forAllConstIter
    (
        fvMeshConstructorTable,
        *fvMeshConstructorTablePtr_,
        cstrIter
    )
    {
        autoPtr<fvMeshStitcher> stitcherPtr(cstrIter()(mesh));

        if (stitcherPtr->changing() == changing)
        {
            return stitcherPtr;
        }
    }

    // Error if an appropriate stitcher was not found
    FatalErrorInFunction
        << typeName << " for " << (changing ? "" : "non-")
        << "changing mesh not found " << nl << nl
        << "Valid " << typeName << "s are : " << endl
        << fvMeshConstructorTablePtr_->toc()
        << exit(FatalError);

    return autoPtr<fvMeshStitcher>(nullptr);
}


// ************************************************************************* //
