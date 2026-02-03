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

#include "uMulticomponentMixture.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::uMulticomponentMixture<ThermoType>::uMulticomponentMixture
(
    const dictionary& dict
)
:
    coefficientMulticomponentMixture<ThermoType>(dict),
    fu_(this->species()[dict.lookup<word>("fuelSpecie")]),
    stoicRatio_(dict.lookup<scalar>("stoichiometricAirFuelMassRatio"))
{
    if (dict.found("oxidantSpecies"))
    {
        const dictionary& oxidantSpecies(dict.subDict("oxidantSpecies"));
        ox_.setSize(oxidantSpecies.size());

        label i = 0;
        forAllConstIter(dictionary, oxidantSpecies, iter)
        {
            const word specieName(iter().keyword());
            const scalar massFraction(iter().stream()[0].scalarToken());

            ox_[i++] = {this->species()[specieName], massFraction};
        }
    }

    if (dict.found("productSpecies"))
    {
        const dictionary& productSpecies(dict.subDict("productSpecies"));
        pr_.setSize(productSpecies.size());

        label i = 0;
        forAllConstIter(dictionary, productSpecies, iter)
        {
            const word specieName(iter().keyword());
            const scalar massFraction(iter().stream()[0].scalarToken());

            pr_[i++] = {this->species()[specieName], massFraction};
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
Foam::scalar Foam::uMulticomponentMixture<ThermoType>::Phi
(
    const scalarFieldListSlice& Yu
) const
{
    return stoicRatio_*Yu[fu_]/max(scalar(1) - Yu[fu_], small);
}


template<class ThermoType>
Foam::PtrList<Foam::volScalarField::Internal>
Foam::uMulticomponentMixture<ThermoType>::prompt
(
    const PtrList<volScalarField>& Yu
) const
{
    PtrList<volScalarField::Internal> Yp(1);
    Yp.set(bInhomogeneousMixture<ThermoType>::FT, Yu[fu_]());

    return Yp;
}


template<class ThermoType>
void Foam::uMulticomponentMixture<ThermoType>::reset
(
    const volScalarField& b,
    PtrList<volScalarField>& Yu,
    const volScalarField& c,
    const PtrList<volScalarField>& Yb
) const
{
    if (!ox_.size() || !pr_.size())
    {
        FatalErrorInFunction
            << "oxidantSpecies or productSpecies not specified"
            << exit(FatalError);
    }

    volScalarField& fuu = Yu[fu_];

    const volScalarField& ftb = Yb[bInhomogeneousMixture<ThermoType>::FT];

    for (label t=0; t<=fuu.nOldTimes(); t++)
    {
        const volScalarField fub
        (
            max
            (
                ftb.oldTime(t) - (scalar(1) - ftb.oldTime(t))/stoicRatio_,
                scalar(0)
            )
        );

        fuu.oldTimeRef(t) = b.oldTime(t)*fuu.oldTime(t) + c.oldTime(t)*fub;

        const volScalarField oxb
        (
            1 - ftb.oldTime(t) - (ftb.oldTime(t) - fub)*stoicRatio_
        );

        forAll(ox_, i)
        {
            Yu[ox_[i].first()].oldTimeRef(t) =
                b.oldTime(t)*Yu[ox_[i].first()].oldTime(t)
              + c.oldTime(t)*ox_[i].second()*oxb;
        }

        const volScalarField prb(1 - fub - oxb);

        forAll(pr_, i)
        {
            Yu[pr_[i].first()].oldTimeRef(t) =
                b.oldTime(t)*Yu[pr_[i].first()].oldTime(t)
              + c.oldTime(t)*pr_[i].second()*prb;
        }
    }
}


// ************************************************************************* //
