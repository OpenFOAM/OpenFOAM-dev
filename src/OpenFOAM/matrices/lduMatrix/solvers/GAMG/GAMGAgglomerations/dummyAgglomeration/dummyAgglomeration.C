/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2019 OpenFOAM Foundation
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

#include "dummyAgglomeration.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dummyAgglomeration, 0);

    addToRunTimeSelectionTable
    (
        GAMGAgglomeration,
        dummyAgglomeration,
        lduMesh
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dummyAgglomeration::dummyAgglomeration
(
    const lduMesh& mesh,
    const dictionary& controlDict
)
:
    GAMGAgglomeration(mesh, controlDict),
    nLevels_(controlDict.lookup<label>("nLevels"))
{
    const label nCoarseCells = mesh.lduAddr().size();

    for
    (
        label nCreatedLevels = 0;
        nCreatedLevels < nLevels_;
        nCreatedLevels++
    )
    {
        nCells_[nCreatedLevels] = nCoarseCells;
        restrictAddressing_.set
        (
            nCreatedLevels,
            new labelField(identity(nCoarseCells))
        );

        agglomerateLduAddressing(nCreatedLevels);
    }

    // Shrink the storage of the levels to those created
    compactLevels(nLevels_);
}


// ************************************************************************* //
