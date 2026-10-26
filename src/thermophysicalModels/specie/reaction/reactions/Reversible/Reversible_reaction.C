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

#include "Reversible_reaction.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Rate, class ThermoType>
typename ThermoType::thermoType
Foam::reactions::Reversible<Rate, ThermoType>::thermo
(
    const reactionCoeffs& r,
    const SpecieThermosRef<ThermoType>& specieThermos
)
{
    typename ThermoType::thermoType rhsThermo
    (
        r.rhs()[0].stoichCoeff
       *specieThermos[r.rhs()[0].index].W()
       *specieThermos[r.rhs()[0].index]
    );

    for (label i=1; i<r.rhs().size(); ++i)
    {
        rhsThermo +=
            r.rhs()[i].stoichCoeff
           *specieThermos[r.rhs()[i].index].W()
           *specieThermos[r.rhs()[i].index];
    }

    typename ThermoType::thermoType lhsThermo
    (
        r.lhs()[0].stoichCoeff
       *specieThermos[r.lhs()[0].index].W()
       *specieThermos[r.lhs()[0].index]
    );

    for (label i=1; i<r.lhs().size(); ++i)
    {
        lhsThermo +=
            r.lhs()[i].stoichCoeff
           *specieThermos[r.lhs()[i].index].W()
           *specieThermos[r.lhs()[i].index];
    }

    // Check for mass imbalance in the reaction. A value of 1 corresponds to an
    // error of 1 H atom in the reaction; i.e. 1 kg/kmol.
    if (mag(lhsThermo.Y() - rhsThermo.Y()) > 0.1)
    {
        FatalErrorInFunction
            << "Mass imbalance for reaction " << r.name() << ": "
            << mag(lhsThermo.Y() - rhsThermo.Y()) << " kg/kmol"
            << exit(FatalError);
    }

    return lhsThermo == rhsThermo;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Rate, class ThermoType>
Foam::reactions::Reversible<Rate, ThermoType>::Reversible
(
    const reactionCoeffs& r,
    const specieThermosRef& specieThermos,
    const Rate& k
)
:
    reaction(r),
    k_(k),
    thermo_(thermo(*this, specieThermos.cast<ThermoType>()))
{}


template<class Rate, class ThermoType>
Foam::reactions::Reversible<Rate, ThermoType>::Reversible
(
    const speciesTable& species,
    const specieThermosRef& specieThermos,
    const dictionary& dict
)
:
    reaction(species, dict),
    k_(species, this->kDimensions(this->lhs()), dict),
    thermo_(thermo(*this, specieThermos.cast<ThermoType>()))
{}


template<class Rate, class ThermoType>
Foam::reactions::Reversible<Rate, ThermoType>::Reversible
(
    const speciesTable& species,
    const specieThermosRef& specieThermos,
    const objectRegistry& ob,
    const dictionary& dict
)
:
    reaction(species, dict),
    k_(species, ob, this->kDimensions(this->lhs()), dict),
    thermo_(thermo(*this, specieThermos.cast<ThermoType>()))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Rate, class ThermoType>
Foam::word Foam::reactions::Reversible<Rate, ThermoType>::selectionType() const
{
    // Don't include thermo in the name
    return word(typeName_(), '<', Rate::type(), '>');
}


template<class Rate, class ThermoType>
Foam::reaction::jacobianDensity
Foam::reactions::Reversible<Rate, ThermoType>::density() const
{
    return
        Rate::hasDdc
      ? reaction::jacobianDensity::all
      : reaction::jacobianDensity::reactantsAndProducts;
}


template<class Rate, class ThermoType>
void Foam::reactions::Reversible<Rate, ThermoType>::preEvaluate() const
{
    k_.preEvaluate();
}


template<class Rate, class ThermoType>
void Foam::reactions::Reversible<Rate, ThermoType>::postEvaluate() const
{
    k_.postEvaluate();
}


template<class Rate, class ThermoType>
Foam::scalar Foam::reactions::Reversible<Rate, ThermoType>::omega
(
    const state& s,
    scalar& omegaf,
    scalar& omegar
) const
{
    const scalar kf = k_.k(s);
    const scalar Kc = thermo_.Kc(s.p, s.T);
    omegaf = kf*this->cPowExponentProduct(this->lhs(), s);
    omegar = kf/Kc*this->cPowExponentProduct(this->rhs(), s);
    return omegaf - omegar;
}


template<class Rate, class ThermoType>
Foam::scalar Foam::reactions::Reversible<Rate, ThermoType>::omega
(
    const state& s,
    scalar& dOmegaDT,
    scalarField& dOmegaDc
) const
{
    // Rate constants
    const scalar kf = k_.k(s);
    const scalar Kc = thermo_.Kc(s.p, s.T);
    const scalar kr = kf/Kc;
    const scalar dkfdT = k_.ddT(s);
    const scalar dkrdT = dkfdT/Kc - kr*thermo_.dKcdTbyKc(s.p, s.T);
    scalarField& dkfdc = s.specieWork_.elmt<1>(); // !!! elmt<0> is dOmegaDc
    if (Rate::hasDdc) k_.ddc(s, dkfdc);

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
            dOmegaDc[i] = dkfdc[i]*Cf - dkfdc[i]/Kc*Cr;
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


template<class Rate, class ThermoType>
void Foam::reactions::Reversible<Rate, ThermoType>::write(Ostream& os) const
{
    reaction::write(os);
    k_.write(os);
}


// ************************************************************************* //
