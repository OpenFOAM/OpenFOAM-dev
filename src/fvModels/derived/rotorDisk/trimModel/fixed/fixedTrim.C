/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2024 OpenFOAM Foundation
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

#include "fixedTrim.H"
#include "rotorDisk.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fixedTrim, 0);

    addToRunTimeSelectionTable(trimModel, fixedTrim, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fixedTrim::fixedTrim
(
    const fv::rotorDisk& rotor,
    const dictionary& dict
)
:
    trimModel(rotor, dict, typeName),
    thetag_(rotor.set().nCells(), 0.0)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

Foam::fixedTrim::~fixedTrim()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fixedTrim::read(const dictionary& dict)
{
    trimModel::read(dict);

    scalar theta0 = coeffs_.lookup<scalar>("theta0", unitDegrees);
    scalar theta1c = coeffs_.lookup<scalar>("theta1c", unitDegrees);
    scalar theta1s = coeffs_.lookup<scalar>("theta1s", unitDegrees);

    const List<point>& x = rotor_.x();
    forAll(thetag_, i)
    {
        scalar psi = x[i].y();
        thetag_[i] = theta0 + theta1c*cos(psi) + theta1s*sin(psi);
    }
}


Foam::tmp<Foam::scalarField> Foam::fixedTrim::thetag() const
{
    return tmp<scalarField>(thetag_);
}


void Foam::fixedTrim::correct
(
    const vectorField& U,
    vectorField& force
) const
{}


void Foam::fixedTrim::correct
(
    const volScalarField rho,
    const vectorField& U,
    vectorField& force
) const
{}


// ************************************************************************* //
