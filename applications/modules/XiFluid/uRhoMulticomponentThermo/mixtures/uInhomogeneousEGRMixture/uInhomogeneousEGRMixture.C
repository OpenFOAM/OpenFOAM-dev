/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025-2026 OpenFOAM Foundation
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

#include "uInhomogeneousEGRMixture.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(uInhomogeneousEGRMixture, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uInhomogeneousEGRMixture::uInhomogeneousEGRMixture
(
    const dictionary& dict
)
:
    species_({"fu", "egr"}),
    stoicRatio_(dict.lookup<scalar>("stoichiometricAirFuelMassRatio")),
    active_(2, true)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::uInhomogeneousEGRMixture::Phi
(
    const scalarFieldListSlice& Yu
) const
{
    const scalar ft = Yu[FU] + Yu[EGR]/(stoicRatio_ + 1);
    return stoicRatio_*ft/max(1 - ft, small);
}


void Foam::uInhomogeneousEGRMixture::read(const dictionary& dict)
{
    stoicRatio_ = dict.lookup<scalar>("stoichiometricAirFuelMassRatio");
}


// ************************************************************************* //
