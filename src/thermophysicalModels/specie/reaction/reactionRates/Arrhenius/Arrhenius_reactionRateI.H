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

#include "Arrhenius_reactionRate.H"
#include "physicoChemicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::reactionRates::Arrhenius::Arrhenius
(
    const scalar A,
    const scalar beta,
    const scalar Ta
)
:
    beta_(beta),
    A_(A),
    Ta_(Ta)
{}


inline Foam::reactionRates::Arrhenius::Arrhenius
(
    const speciesTable&,
    const dimensionSet& dims,
    const dictionary& dict
)
:
    beta_(dict.lookup<scalar>("beta", dimless)),
    A_(dict.lookup<scalar>("A", dims/pow(dimensions::temperature, beta_))),
    Ta_
    (
        dict.found("Ta") || !dict.found("Ea")
      ? dict.lookup<scalar>("Ta", dimensions::temperature)
      : dict.lookup<scalar>("Ea", dimensions::energy/dimensions::moles)
       /constant::physicoChemical::RR.value()
    )
{}


inline Foam::reactionRates::Arrhenius::Arrhenius
(
    const speciesTable& species,
    const objectRegistry& ob,
    const dimensionSet& dims,
    const dictionary& dict
)
:
    Arrhenius(species, dims, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline void Foam::reactionRates::Arrhenius::preEvaluate() const
{}


inline void Foam::reactionRates::Arrhenius::postEvaluate() const
{}


inline Foam::scalar Foam::reactionRates::Arrhenius::k
(
    const reaction::state& s
) const
{
    return A_*exp(beta_*s.lnT - Ta_*s.rT);
}


inline Foam::scalar Foam::reactionRates::Arrhenius::ddT
(
    const reaction::state& s
) const
{
    return A_*exp(beta_*s.lnT - Ta_*s.rT)*(beta_ + Ta_*s.rT)*s.rT;
}


inline void Foam::reactionRates::Arrhenius::ddc
(
    const reaction::state& s,
    scalarField& ddc
) const
{
    NotImplemented;
}


inline void Foam::reactionRates::Arrhenius::write(Ostream& os) const
{
    writeEntry(os, "A", A_);
    writeEntry(os, "beta", beta_);
    writeEntry(os, "Ta", Ta_);
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const reactionRates::Arrhenius& rr
)
{
    rr.write(os);
    return os;
}


// ************************************************************************* //
