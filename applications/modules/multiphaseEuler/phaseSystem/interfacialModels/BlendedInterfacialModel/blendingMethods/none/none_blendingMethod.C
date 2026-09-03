/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2026 OpenFOAM Foundation
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

#include "none_blendingMethod.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace blendingMethods
{
    defineTypeNameAndDebug(none, 0);
    addToRunTimeSelectionTable(blendingMethod, none, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * /

Foam::tmp<Foam::volScalarField> Foam::blendingMethods::none::fContinuous
(
    const UPtrList<const volScalarField>& alphas,
    const label index
) const
{
    return constant(alphas, 0);
}


Foam::tmp<Foam::volScalarField> Foam::blendingMethods::none::fDisplaced
(
    const UPtrList<const volScalarField>& alphas,
    const label displacingPhasei
) const
{
    return constant(alphas, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blendingMethods::none::none
(
    const dictionary& dict,
    const phaseInterface& interface,
    const bool allowDisplaced
)
:
    blendingMethod(interface)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::blendingMethods::none::~none()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::blendingMethods::none::canBeContinuous(const label index) const
{
    return false;
}


bool Foam::blendingMethods::none::canSegregate() const
{
    return false;
}


bool Foam::blendingMethods::none::isDisplacedBy
(
    const label displacingPhasei
) const
{
    return false;
}


bool Foam::blendingMethods::none::functionOfAlphas() const
{
    return true;
}


// ************************************************************************* //
