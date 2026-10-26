/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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

#include "NonEquilibriumReversible_reaction.H"
#include "delimitDictionary.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Rate>
Foam::reactions::NonEquilibriumReversible<Rate>::NonEquilibriumReversible
(
    const reactionCoeffs& r,
    const specieThermosRef& specieThermos,
    const Rate& kf,
    const Rate& kr
)
:
    reaction(r),
    kf_(kf),
    kr_(kr)
{}


template<class Rate>
Foam::reactions::NonEquilibriumReversible<Rate>::NonEquilibriumReversible
(
    const speciesTable& species,
    const specieThermosRef& specieThermos,
    const dictionary& dict
)
:
    reaction(species, dict),
    kf_(species, this->kDimensions(this->lhs()), dict.subDict("forward")),
    kr_(species, this->kDimensions(this->rhs()), dict.subDict("reverse"))
{}


template<class Rate>
Foam::reactions::NonEquilibriumReversible<Rate>::NonEquilibriumReversible
(
    const speciesTable& species,
    const specieThermosRef& specieThermos,
    const objectRegistry& ob,
    const dictionary& dict
)
:
    reaction(species, dict),
    kf_(species, ob, this->kDimensions(this->lhs()), dict.subDict("forward")),
    kr_(species, ob, this->kDimensions(this->rhs()), dict.subDict("reverse"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Rate>
Foam::reaction::jacobianDensity
Foam::reactions::NonEquilibriumReversible<Rate>::density() const
{
    return
        Rate::hasDdc
      ? reaction::jacobianDensity::all
      : reaction::jacobianDensity::reactantsAndProducts;
}


template<class Rate>
void Foam::reactions::NonEquilibriumReversible<Rate>::preEvaluate() const
{
    kf_.preEvaluate();
    kr_.preEvaluate();
}


template<class Rate>
void Foam::reactions::NonEquilibriumReversible<Rate>::postEvaluate() const
{
    kf_.postEvaluate();
    kr_.postEvaluate();
}


template<class Rate>
Foam::scalar Foam::reactions::NonEquilibriumReversible<Rate>::omega
(
    const state& s,
    scalar& omegaf,
    scalar& omegar
) const
{
    omegaf = kf_.k(s)*this->cPowExponentProduct(this->lhs(), s);
    omegar = kr_.k(s)*this->cPowExponentProduct(this->rhs(), s);
    return omegaf - omegar;
}


template<class Rate>
Foam::scalar Foam::reactions::NonEquilibriumReversible<Rate>::omega
(
    const state& s,
    scalar& dOmegaDT,
    scalarField& dOmegaDc
) const
{
    // Rate constants
    const scalar kf = kf_.k(s);
    const scalar kr = kr_.k(s);
    const scalar dkfdT = kf_.ddT(s);
    const scalar dkrdT = kr_.ddT(s);
    scalarField& dkfdc = s.specieWork_.elmt<1>(); // !!! elmt<0> is dOmegaDc
    scalarField& dkrdc = s.specieWork_.elmt<2>();
    if (Rate::hasDdc)
    {
        kf_.ddc(s, dkfdc);
        kr_.ddc(s, dkrdc);
    }

    // Concentration products
    scalarField& dCfdc = this->lWork_, & dCrdc = this->rWork_;
    const scalar Cf = this->cPowExponentProduct(this->lhs(), s, dCfdc);
    const scalar Cr = this->cPowExponentProduct(this->rhs(), s, dCrdc);

    // Rate
    const scalar omega = kf*Cf - kr*Cr;

    // Derivative w.r.t temperature
    dOmegaDT = dkfdT*Cf - dkrdT*Cr;

    // Derivative w.r.t concentration
    if (Rate::hasDdc)
    {
        forAll(s.c, i)
        {
            dOmegaDc[i] = dkfdc[i]*Cf - dkrdc[i]*Cr;
        }
    }
    else
    {
        forAll(specieIndices(), i)
        {
            dOmegaDc[specieIndices()[i]] = 0;
        }
    }
    forAll(this->lhs(), i)
    {
        dOmegaDc[this->lhs()[i].index] += kf*dCfdc[i];
    }
    forAll(this->rhs(), i)
    {
        dOmegaDc[this->rhs()[i].index] -= kr*dCrdc[i];
    }

    return omega;
}


template<class Rate>
void Foam::reactions::NonEquilibriumReversible<Rate>::write(Ostream& os) const
{
    reaction::write(os);
    {
        const delimitDictionary delimitForward(os, "forward");
        kf_.write(os);
    }
    {
        const delimitDictionary delimitReverse(os, "reverse");
        kr_.write(os);
    }
}


// ************************************************************************* //
