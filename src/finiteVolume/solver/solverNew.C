/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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

#include "solver.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

void Foam::solver::load(const word& solverName)
{
    libs.open("lib" + solverName + ".so");
}


Foam::autoPtr<Foam::solver> Foam::solver::New
(
    const word& solverName,
    fvMesh& mesh
)
{
    Info<< "Selecting solver " << solverName << endl;

    load(solverName);

    if (!fvMeshConstructorTablePtr_)
    {
        FatalErrorInFunction
            << "solvers table is empty"
            << exit(FatalError);
    }

    fvMeshConstructorTable::iterator cstrIter =
        fvMeshConstructorTablePtr_->find(solverName);

    if (cstrIter == fvMeshConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown solver type "
            << solverName << nl << nl
                << "Valid solvers are :" << endl
                << fvMeshConstructorTablePtr_->sortedToc()
                << exit(FatalError);
    }

    autoPtr<solver> solverPtr(cstrIter()(mesh));

    // Ensure fvModels and fvConstraints are constructed
    // before time is incremented
    solverPtr->fvModels();
    solverPtr->fvConstraints();

    return solverPtr;
}


// ************************************************************************* //
