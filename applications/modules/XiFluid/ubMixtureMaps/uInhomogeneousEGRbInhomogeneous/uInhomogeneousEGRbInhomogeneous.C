/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "uInhomogeneousEGRbInhomogeneous.H"
#include "uInhomogeneousEGRMixture.H"
#include "bInhomogeneousMixture.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace ubMixtureMaps
{
    defineTypeNameAndDebug(uInhomogeneousEGRbInhomogeneous, 0);
    addToRunTimeSelectionTable
    (
        ubMixtureMap,
        uInhomogeneousEGRbInhomogeneous,
        thermo
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ubMixtureMaps::uInhomogeneousEGRbInhomogeneous::
uInhomogeneousEGRbInhomogeneous
(
    const uRhoMulticomponentThermo& uThermo,
    const bRhoMulticomponentThermo& bThermo
)
:
    ubMixtureMap(uThermo, bThermo)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ubMixtureMaps::uInhomogeneousEGRbInhomogeneous::
~uInhomogeneousEGRbInhomogeneous()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::PtrList<Foam::volScalarField::Internal>
Foam::ubMixtureMaps::uInhomogeneousEGRbInhomogeneous::prompt
(
    const PtrList<volScalarField>& Yu
) const
{
    PtrList<volScalarField::Internal> Yp(1);
    Yp.set(bInhomogeneousMixture::FT, Yu[uInhomogeneousEGRMixture::FU]());

    return Yp;
}


void Foam::ubMixtureMaps::uInhomogeneousEGRbInhomogeneous::reset
(
    const volScalarField& b,
    UPtrList<volScalarField>& Yu,
    const volScalarField& c,
    const UPtrList<const volScalarField>& Yb
) const
{
    const uInhomogeneousEGRMixture& uIEGR =
        uMixtureCast<uInhomogeneousEGRMixture>();

    volScalarField& fuu = Yu[uInhomogeneousEGRMixture::FU];
    volScalarField& egr = Yu[uInhomogeneousEGRMixture::EGR];
    const volScalarField& ftb = Yb[bInhomogeneousMixture::FT];

    const volScalarField fub
    (
        max
        (
            ftb - (scalar(1) - ftb)/uIEGR.stoicRatio(),
            scalar(0)
        )
    );
    const volScalarField oxb
    (
        1 - ftb - (ftb - fub)*uIEGR.stoicRatio()
    );

    fuu = b*fuu + c*fub;
    egr = b*egr + c*(1 - fub - oxb);
}


// ************************************************************************* //
