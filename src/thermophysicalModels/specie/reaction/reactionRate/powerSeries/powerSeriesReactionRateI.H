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

#include "powerSeriesReactionRate.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::powerSeriesReactionRate::powerSeriesReactionRate
(
    const scalar A,
    const scalar beta,
    const FixedList<scalar, nCoeff_> coeffs
)
:
    beta_(beta),
    A_(A),
    coeffs_(coeffs)
{}


inline Foam::powerSeriesReactionRate::powerSeriesReactionRate
(
    const speciesTable&,
    const dimensionSet& dims,
    const dictionary& dict
)
:
    beta_(dict.lookup<scalar>("beta", dimless)),
    A_(dict.lookup<scalar>("A", dims/pow(dimensions::temperature, beta_))),
    coeffs_(dict.lookup("coeffs"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline void Foam::powerSeriesReactionRate::preEvaluate() const
{}


inline void Foam::powerSeriesReactionRate::postEvaluate() const
{}


inline Foam::scalar Foam::powerSeriesReactionRate::k
(
    const reaction::state& s
) const
{
    scalar expArg = beta_*s.lnT;

    scalar rTn = 1;
    forAll(coeffs_, n)
    {
        rTn *= s.rT;
        expArg += coeffs_[n]*rTn;
    }

    return A_*exp(expArg);
}


inline Foam::scalar Foam::powerSeriesReactionRate::ddT
(
    const reaction::state& s
) const
{
    scalar expArg = beta_*s.lnT;
    scalar deriv = beta_;

    scalar rTn = 1;
    forAll(coeffs_, n)
    {
        rTn *= s.rT;
        expArg += coeffs_[n]*rTn;
        deriv -= (n + 1)*coeffs_[n]*rTn;
    }

    return A_*exp(expArg)*deriv*s.rT;
}


inline bool Foam::powerSeriesReactionRate::hasDdc() const
{
    return false;
}


inline void Foam::powerSeriesReactionRate::ddc
(
    const reaction::state& s,
    scalarField& ddc
) const
{
    NotImplemented;
}


inline void Foam::powerSeriesReactionRate::write(Ostream& os) const
{
    writeEntry(os, "A", A_);
    writeEntry(os, "beta", beta_);
    writeEntry(os, "coeffs", coeffs_);
}


inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const powerSeriesReactionRate& psrr
)
{
    psrr.write(os);
    return os;
}


// ************************************************************************* //
