/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2022 OpenFOAM Foundation
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

#include "noVirtualMass.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace virtualMassModels
{
    defineTypeNameAndDebug(noVirtualMass, 0);
    addToRunTimeSelectionTable(virtualMassModel, noVirtualMass, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::virtualMassModels::noVirtualMass::noVirtualMass
(
    const dictionary& dict,
    const phaseInterface& interface,
    const bool registerObject
)
:
    virtualMassModel(dict, interface, registerObject),
    interface_(interface)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::virtualMassModels::noVirtualMass::~noVirtualMass()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::virtualMassModels::noVirtualMass::K() const
{
    return volScalarField::New
    (
        "K",
        interface_.mesh(),
        dimensionedScalar(dimK, 0)
    );
}


Foam::tmp<Foam::surfaceScalarField>
Foam::virtualMassModels::noVirtualMass::Kf() const
{
    return surfaceScalarField::New
    (
        "Kf",
        interface_.mesh(),
        dimensionedScalar(dimK, 0)
    );
}


// ************************************************************************* //
