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

#include "noSoot.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{
namespace sootModels
{
    defineTypeNameAndDebug(noSoot, 0);

    addToRunTimeSelectionTable
    (
        sootModel,
        noSoot,
        dictionary
    );
}
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiationModels::sootModels::noSoot::noSoot
(
    const dictionary& dict,
    const fvMesh& mesh,
    const word& modelType
)
:
    sootModel(dict, mesh, modelType)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiationModels::sootModels::noSoot::~noSoot()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::radiationModels::sootModels::noSoot::correct()
{}


const Foam::volScalarField&
Foam::radiationModels::sootModels::noSoot::soot() const
{
    NotImplemented;
    return tmp<volScalarField>(nullptr);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
