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

#include "Irreversible_reaction.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Rate>
Foam::reactions::Irreversible<Rate>::Irreversible
(
    const reactionCoeffs& r,
    const specieThermosRef& specieThermos,
    const Rate& k
)
:
    reaction(r),
    k_(k)
{}


template<class Rate>
Foam::reactions::Irreversible<Rate>::Irreversible
(
    const speciesTable& species,
    const specieThermosRef& specieThermos,
    const dictionary& dict
)
:
    reaction(species, dict),
    k_(species, this->kDimensions(this->lhs()), dict)
{}


template<class Rate>
Foam::reactions::Irreversible<Rate>::Irreversible
(
    const speciesTable& species,
    const specieThermosRef& specieThermos,
    const objectRegistry& ob,
    const dictionary& dict
)
:
    reaction(species, dict),
    k_(species, ob, this->kDimensions(this->lhs()), dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Rate>
Foam::reaction::jacobianDensity
Foam::reactions::Irreversible<Rate>::density() const
{
    return
        Rate::hasDdc
      ? reaction::jacobianDensity::all
      : reaction::jacobianDensity::reactants;
}


template<class Rate>
void Foam::reactions::Irreversible<Rate>::preEvaluate() const
{
    k_.preEvaluate();
}


template<class Rate>
void Foam::reactions::Irreversible<Rate>::postEvaluate() const
{
    k_.postEvaluate();
}


template<class Rate>
Foam::scalar Foam::reactions::Irreversible<Rate>::omega
(
    const state& s,
    scalar& omegaf,
    scalar& omegar
) const
{
    omegaf = k_.k(s)*this->cPowExponentProduct(this->lhs(), s);
    omegar = 0;
    return omegaf;
}


template<class Rate>
Foam::scalar Foam::reactions::Irreversible<Rate>::omega
(
    const state& s,
    scalar& dOmegaDT,
    scalarField& dOmegaDc
) const
{
    // Rate constant
    const scalar kf = k_.k(s);
    const scalar dkfdT = k_.ddT(s);
    scalarField& dkfdc = s.specieWork_.elmt<1>(); // !!! elmt<0> is dOmegaDc
    if (Rate::hasDdc) k_.ddc(s, dkfdc);

    // Concentration product
    scalarField& dCfdc = this->lWork_;
    const scalar Cf = this->cPowExponentProduct(this->lhs(), s, dCfdc);

    // Rate
    const scalar omega = kf*Cf;

    // Derivative w.r.t temperature
    dOmegaDT = dkfdT*Cf;

    // Derivative w.r.t concentration
    if (Rate::hasDdc)
    {
        forAll(s.c, i)
        {
            dOmegaDc[i] = dkfdc[i]*Cf;
        }
    }
    else
    {
        forAll(this->lhs(), i)
        {
            dOmegaDc[this->lhs()[i].index] = 0;
        }
    }
    forAll(this->lhs(), i)
    {
        dOmegaDc[this->lhs()[i].index] += kf*dCfdc[i];
    }

    return omega;
}


template<class Rate>
void Foam::reactions::Irreversible<Rate>::write(Ostream& os) const
{
    reaction::write(os);
    k_.write(os);
}


// ************************************************************************* //
