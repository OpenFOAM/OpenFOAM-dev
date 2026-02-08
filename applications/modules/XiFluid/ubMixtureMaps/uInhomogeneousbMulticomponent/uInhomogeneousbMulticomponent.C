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

#include "uInhomogeneousbMulticomponent.H"
#include "uInhomogeneousMixture.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace ubMixtureMaps
{
    defineTypeNameAndDebug(uInhomogeneousbMulticomponent, 0);
    addToRunTimeSelectionTable
    (
        ubMixtureMap,
        uInhomogeneousbMulticomponent,
        thermo
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ubMixtureMaps::uInhomogeneousbMulticomponent::
uInhomogeneousbMulticomponent
(
    const uRhoMulticomponentThermo& uThermo,
    const bRhoMulticomponentThermo& bThermo
)
:
    ubMixtureMap(uThermo, bThermo)
{
    const dictionary& bDict = bThermo.properties();
    const speciesTable& bSpecies = bThermo.species();

    fu_ = bSpecies[bDict.lookup<word>("fuelSpecie")];

    const dictionary& oxidantSpecies(bDict.subDict("oxidantSpecies"));
    ox_.setSize(oxidantSpecies.size());

    label i = 0;
    forAllConstIter(dictionary, oxidantSpecies, iter)
    {
        const word specieName(iter().keyword());
        const scalar massFraction(iter().stream()[0].scalarToken());

        ox_[i++] = {bSpecies[specieName], massFraction};
    }

    const dictionary& productSpecies(bDict.subDict("productSpecies"));
    pr_.setSize(productSpecies.size());

    i = 0;
    forAllConstIter(dictionary, productSpecies, iter)
    {
        const word specieName(iter().keyword());
        const scalar massFraction(iter().stream()[0].scalarToken());

        pr_[i++] = {bSpecies[specieName], massFraction};
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ubMixtureMaps::uInhomogeneousbMulticomponent::
~uInhomogeneousbMulticomponent()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::PtrList<Foam::volScalarField::Internal>
Foam::ubMixtureMaps::uInhomogeneousbMulticomponent::prompt
(
    const PtrList<volScalarField>& Yu
) const
{
    const PtrList<volScalarField>& Yb = bThermo_.Y();
    PtrList<volScalarField::Internal> Yp(Yb.size());

    const uInhomogeneousMixture& um = uMixtureCast<uInhomogeneousMixture>();

    const volScalarField::Internal& ftb = Yu[uInhomogeneousMixture::FU];
    const volScalarField::Internal fub
    (
        max
        (
            ftb - (scalar(1) - ftb)/um.stoicRatio(),
            scalar(0)
        )
    );

    Yp.set(fu_, fub);

    const volScalarField::Internal oxb(1 - ftb - (ftb - fub)*um.stoicRatio());

    forAll(ox_, i)
    {
        Yp.set(ox_[i].first(), ox_[i].second()*oxb);
    }

    const volScalarField::Internal prb(1 - fub - oxb);

    forAll(pr_, i)
    {
        Yp.set(pr_[i].first(), pr_[i].second()*prb);
    }

    forAll(Yp, i)
    {
        if (!Yp.set(i))
        {
            Yp.set
            (
                i,
                volScalarField::Internal::New
                (
                    Yb[i].name(),
                    ftb.mesh(),
                    dimensionedScalar(dimless, 0)
                )
            );
        }
    }

    return Yp;
}


// ************************************************************************* //
