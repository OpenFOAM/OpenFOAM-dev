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

#include "pairPotential.H"
#include "energyScalingFunction.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pairPotential, 0);
    defineRunTimeSelectionTable(pairPotential, dictionary);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::pairPotential::scaleEnergy(scalar& e, const scalar r) const
{
    if (!esfPtr_)
    {
        esfPtr_ = energyScalingFunction::New
        (
            name_, pairPotentialProperties_, *this
        ).ptr();
    }

    esfPtr_->scaleEnergy(e, r);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pairPotential::pairPotential
(
    const word& name,
    const dictionary& pairPotentialProperties
)
:
    name_(name),
    pairPotentialProperties_(pairPotentialProperties),
    rCut_(readScalar(pairPotentialProperties_.lookup("rCut"))),
    rCutSqr_(rCut_*rCut_),
    rMin_(readScalar(pairPotentialProperties_.lookup("rMin"))),
    dr_(readScalar(pairPotentialProperties_.lookup("dr"))),
    forceLookup_(0),
    energyLookup_(0),
    esfPtr_(nullptr),
    writeTables_(Switch(pairPotentialProperties_.lookup("writeTables")))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::pairPotential::setLookupTables()
{
    label N = label((rCut_ - rMin_)/dr_) + 1;

    forceLookup_.setSize(N);

    energyLookup_.setSize(N);

    forAll(forceLookup_, k)
    {
        energyLookup_[k] = scaledEnergy(k*dr_ + rMin_);

        forceLookup_[k] = -energyDerivative((k*dr_ + rMin_), true);
    }
}


Foam::scalar Foam::pairPotential::force(const scalar r) const
{
    scalar k_rIJ = (r - rMin_)/dr_;

    label k = label(k_rIJ);

    if (k < 0)
    {
        FatalErrorInFunction
            << "r less than rMin in pair potential " << name_ << nl
            << abort(FatalError);
    }

    scalar f =
        (k_rIJ - k)*forceLookup_[k+1]
      + (k + 1 - k_rIJ)*forceLookup_[k];

    return f;
}


Foam::List<Foam::Pair<Foam::scalar>>
Foam::pairPotential::forceTable() const
{
    List<Pair<scalar>> forceTab(forceLookup_.size());

    forAll(forceLookup_,k)
    {
        forceTab[k].first() = rMin_ + k*dr_;

        forceTab[k].second() = forceLookup_[k];
    }

    return forceTab;
}


Foam::scalar Foam::pairPotential::energy(const scalar r) const
{
    scalar k_rIJ = (r - rMin_)/dr_;

    label k = label(k_rIJ);

    if (k < 0)
    {
        FatalErrorInFunction
            << "r less than rMin in pair potential " << name_ << nl
            << abort(FatalError);
    }

    scalar e =
        (k_rIJ - k)*energyLookup_[k+1]
      + (k + 1 - k_rIJ)*energyLookup_[k];

    return e;
}


Foam::List<Foam::Pair<Foam::scalar>>
    Foam::pairPotential::energyTable() const
{
    List<Pair<scalar>> energyTab(energyLookup_.size());

    forAll(energyLookup_,k)
    {
        energyTab[k].first() = rMin_ + k*dr_;

        energyTab[k].second() = energyLookup_[k];
    }

    return energyTab;
}


Foam::scalar Foam::pairPotential::scaledEnergy(const scalar r) const
{
    scalar e = unscaledEnergy(r);

    scaleEnergy(e, r);

    return e;
}


Foam::scalar Foam::pairPotential::energyDerivative
(
    const scalar r,
    const bool scaledEnergyDerivative
) const
{
    // Local quadratic fit to energy: E = a0 + a1*r + a2*r^2
    // Differentiate to give f = -dE/dr = -a1 - 2*a2*r

    scalar ra = r - dr_;
    scalar rf = r;
    scalar rb = r + dr_;

    scalar Ea, Ef, Eb;

    if (scaledEnergyDerivative)
    {
        Ea = scaledEnergy(ra);
        Ef = scaledEnergy(rf);
        Eb = scaledEnergy(rb);
    }
    else
    {
        Ea = unscaledEnergy(ra);
        Ef = unscaledEnergy(rf);
        Eb = unscaledEnergy(rb);
    }

    scalar denominator = (ra - rf)*(ra - rb)*(rf - rb);

    scalar a1 =
    (
        rb*rb*(Ea - Ef) + ra*ra*(Ef - Eb) + rf*rf*(Eb - Ea)
    ) / denominator;

    scalar a2 =
    (
        rb*(Ef - Ea) + rf*(Ea - Eb) + ra*(Eb - Ef)
    ) / denominator;

    return a1 + 2.0*a2*r;
}


bool Foam::pairPotential::read(const dictionary& pairPotentialProperties)
{
    pairPotentialProperties_ = pairPotentialProperties;

    return true;
}


// ************************************************************************* //
