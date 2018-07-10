/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "dampedCoulomb.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace pairPotentials
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(dampedCoulomb, 0);

addToRunTimeSelectionTable
(
    pairPotential,
    dampedCoulomb,
    dictionary
);

scalar dampedCoulomb::oneOverFourPiEps0 =
    1.0/(4.0*constant::mathematical::pi*8.854187817e-12);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

dampedCoulomb::dampedCoulomb
(
    const word& name,
    const dictionary& pairPotentialProperties
)
:
    pairPotential(name, pairPotentialProperties),
    dampedCoulombCoeffs_
    (
        pairPotentialProperties.subDict(typeName + "Coeffs")
    ),
    alpha_(readScalar(dampedCoulombCoeffs_.lookup("alpha")))
{
    setLookupTables();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

scalar dampedCoulomb::unscaledEnergy(const scalar r) const
{
    return oneOverFourPiEps0*erfc(alpha_*r)/r;
}


bool dampedCoulomb::read(const dictionary& pairPotentialProperties)
{
    pairPotential::read(pairPotentialProperties);

    dampedCoulombCoeffs_ =
        pairPotentialProperties.subDict(typeName + "Coeffs");

    dampedCoulombCoeffs_.lookup("alpha") >> alpha_;

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace pairPotentials
} // End namespace Foam

// ************************************************************************* //
