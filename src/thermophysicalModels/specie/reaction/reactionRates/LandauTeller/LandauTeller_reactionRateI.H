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

#include "LandauTeller_reactionRate.H"
#include "physicoChemicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::reactionRates::LandauTeller::LandauTeller
(
    const scalar A,
    const scalar beta,
    const scalar Ta,
    const scalar B,
    const scalar C
)
:
    beta_(beta),
    A_(A),
    Ta_(Ta),
    B_(B),
    C_(C)
{}


inline Foam::reactionRates::LandauTeller::LandauTeller
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
    ),
    B_(dict.lookup<scalar>("B", cbrt(dimensions::temperature))),
    C_(dict.lookup<scalar>("C", sqr(cbrt(dimensions::temperature))))
{}


inline Foam::reactionRates::LandauTeller::LandauTeller
(
    const speciesTable& species,
    const objectRegistry& ob,
    const dimensionSet& dims,
    const dictionary& dict
)
:
    LandauTeller(species, dims, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline void Foam::reactionRates::LandauTeller::preEvaluate() const
{}


inline void Foam::reactionRates::LandauTeller::postEvaluate() const
{}


inline Foam::scalar Foam::reactionRates::LandauTeller::k
(
    const reaction::state& s
) const
{
    const scalar cbrtT = cbrt(s.T);
    return A_*exp(beta_*s.lnT - Ta_*s.rT + B_/cbrtT + C_/sqr(cbrtT));
}


inline Foam::scalar Foam::reactionRates::LandauTeller::ddT
(
    const reaction::state& s
) const
{
    const scalar cbrtT = cbrt(s.T);
    const scalar k = A_*exp(beta_*s.lnT - Ta_*s.rT + B_/cbrtT + C_/sqr(cbrtT));
    const scalar deriv = beta_ + Ta_*s.rT - B_/(3*cbrtT) - 2*C_/(3*sqr(cbrtT));
    return k*deriv*s.rT;
}


inline void Foam::reactionRates::LandauTeller::ddc
(
    const reaction::state& s,
    scalarField& ddc
) const
{
    NotImplemented;
}


inline void Foam::reactionRates::LandauTeller::write(Ostream& os) const
{
    writeEntry(os, "A", A_);
    writeEntry(os, "beta", beta_);
    writeEntry(os, "Ta", Ta_);
    writeEntry(os, "B", B_);
    writeEntry(os, "C", C_);
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const reactionRates::LandauTeller& rr
)
{
    rr.write(os);
    return os;
}


// ************************************************************************* //
