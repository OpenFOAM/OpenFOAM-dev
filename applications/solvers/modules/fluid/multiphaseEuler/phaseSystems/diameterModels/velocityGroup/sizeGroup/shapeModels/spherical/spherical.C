/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2022 OpenFOAM Foundation
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

#include "spherical.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace shapeModels
{
    defineTypeNameAndDebug(spherical, 0);
    addToRunTimeSelectionTable
    (
        shapeModel,
        spherical,
        dictionary
    );
}
}
}

using Foam::constant::mathematical::pi;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::shapeModels::spherical::
spherical
(
    const dictionary& dict,
    const sizeGroup& group
)
:
    shapeModel(dict, group)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diameterModels::shapeModels::spherical::~spherical()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::tmp<Foam::volScalarField>
Foam::diameterModels::shapeModels::spherical::a() const
{
    return tmp<volScalarField>
    (
        volScalarField::New
        (
            "a",
            sizeGroup_.mesh(),
            6/sizeGroup_.dSph()*sizeGroup_.x()
        )
    );
}


const Foam::tmp<Foam::volScalarField>
Foam::diameterModels::shapeModels::spherical::d() const
{
    return tmp<volScalarField>
    (
        volScalarField::New("d", sizeGroup_.mesh(), sizeGroup_.dSph())
    );
}


// ************************************************************************* //
