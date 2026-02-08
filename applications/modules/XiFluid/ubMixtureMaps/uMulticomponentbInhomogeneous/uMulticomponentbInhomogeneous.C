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

#include "uMulticomponentbInhomogeneous.H"
#include "uMulticomponentMixture.H"
#include "bInhomogeneousMixture.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace ubMixtureMaps
{
    defineTypeNameAndDebug(uMulticomponentbInhomogeneous, 0);
    addToRunTimeSelectionTable
    (
        ubMixtureMap,
        uMulticomponentbInhomogeneous,
        thermo
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ubMixtureMaps::uMulticomponentbInhomogeneous::
uMulticomponentbInhomogeneous
(
    const uRhoMulticomponentThermo& uThermo,
    const bRhoMulticomponentThermo& bThermo
)
:
    ubMixtureMap(uThermo, bThermo),
    fu_(uMixtureCast<uMulticomponentMixture>().fu())
{
    const dictionary& uDict = uThermo.properties();
    const speciesTable& uSpecies = uThermo.species();

    if (uDict.found("oxidantSpecies"))
    {
        const dictionary& oxidantSpecies(uDict.subDict("oxidantSpecies"));
        ox_.setSize(oxidantSpecies.size());

        label i = 0;
        forAllConstIter(dictionary, oxidantSpecies, iter)
        {
            const word specieName(iter().keyword());
            const scalar massFraction(iter().stream()[0].scalarToken());

            ox_[i++] = {uSpecies[specieName], massFraction};
        }
    }

    if (uDict.found("productSpecies"))
    {
        const dictionary& productSpecies(uDict.subDict("productSpecies"));
        pr_.setSize(productSpecies.size());

        label i = 0;
        forAllConstIter(dictionary, productSpecies, iter)
        {
            const word specieName(iter().keyword());
            const scalar massFraction(iter().stream()[0].scalarToken());

            pr_[i++] = {uSpecies[specieName], massFraction};
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ubMixtureMaps::uMulticomponentbInhomogeneous::
~uMulticomponentbInhomogeneous()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::PtrList<Foam::volScalarField::Internal>
Foam::ubMixtureMaps::uMulticomponentbInhomogeneous::prompt
(
    const PtrList<volScalarField>& Yu
) const
{
    PtrList<volScalarField::Internal> Yp(1);
    Yp.set(bInhomogeneousMixture::FT, Yu[fu_]());

    return Yp;
}


void Foam::ubMixtureMaps::uMulticomponentbInhomogeneous::reset
(
    const volScalarField& b,
    UPtrList<volScalarField>& Yu,
    const volScalarField& c,
    const UPtrList<const volScalarField>& Yb
) const
{
    if (!ox_.size() || !pr_.size())
    {
        FatalErrorInFunction
            << "oxidantSpecies or productSpecies not specified"
            << exit(FatalError);
    }

    const uMulticomponentMixture& um = uMixtureCast<uMulticomponentMixture>();

    volScalarField& fuu = Yu[fu_];

    const volScalarField& ftb = Yb[bInhomogeneousMixture::FT];

    const volScalarField fub
    (
        max
        (
            ftb - (scalar(1) - ftb)/um.stoicRatio(),
            scalar(0)
        )
    );

    fuu = b*fuu + c*fub;

    const volScalarField oxb(1 - ftb - (ftb - fub)*um.stoicRatio());

    forAll(ox_, i)
    {
        Yu[ox_[i].first()] = b*Yu[ox_[i].first()] + c*ox_[i].second()*oxb;
    }

    const volScalarField prb(1 - fub - oxb);

    forAll(pr_, i)
    {
        Yu[pr_[i].first()] = b*Yu[pr_[i].first()] + c*pr_[i].second()*prb;
    }
}


// ************************************************************************* //
