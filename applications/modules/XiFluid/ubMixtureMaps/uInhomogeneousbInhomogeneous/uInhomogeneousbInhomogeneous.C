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

#include "uInhomogeneousbInhomogeneous.H"
#include "uInhomogeneousMixture.H"
#include "bInhomogeneousMixture.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace ubMixtureMaps
{
    defineTypeNameAndDebug(uInhomogeneousbInhomogeneous, 0);
    addToRunTimeSelectionTable
    (
        ubMixtureMap,
        uInhomogeneousbInhomogeneous,
        thermo
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ubMixtureMaps::uInhomogeneousbInhomogeneous::
uInhomogeneousbInhomogeneous
(
    const uRhoMulticomponentThermo& uThermo,
    const bRhoMulticomponentThermo& bThermo
)
:
    ubMixtureMap(uThermo, bThermo)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ubMixtureMaps::uInhomogeneousbInhomogeneous::
~uInhomogeneousbInhomogeneous()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::PtrList<Foam::volScalarField::Internal>
Foam::ubMixtureMaps::uInhomogeneousbInhomogeneous::prompt
(
    const PtrList<volScalarField>& Yu
) const
{
    PtrList<volScalarField::Internal> Yp(1);
    Yp.set(bInhomogeneousMixture::FT, Yu[uInhomogeneousMixture::FU]());

    return Yp;
}


// ************************************************************************* //
