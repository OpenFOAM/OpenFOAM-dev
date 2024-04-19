/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "JanevReactionRate.H"
#include "physicoChemicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::JanevReactionRate::JanevReactionRate
(
    const scalar A,
    const scalar beta,
    const scalar Ta,
    const FixedList<scalar, nb_> b
)
:
    beta_(beta),
    A_(A),
    Ta_(Ta),
    b_(b)
{}


inline Foam::JanevReactionRate::JanevReactionRate
(
    const speciesTable&,
    const dimensionSet& dims,
    const dictionary& dict
)
:
    beta_(dict.lookup<scalar>("beta", dimless)),
    A_(dict.lookup<scalar>("A", dims/pow(dimTemperature, beta_))),
    Ta_
    (
        dict.found("Ta") || !dict.found("Ea")
      ? dict.lookup<scalar>("Ta", dimTemperature)
      : dict.lookup<scalar>("Ea", dimEnergy/dimMoles)
       /constant::physicoChemical::RR.value()
    ),
    b_(dict.lookup("b"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline void Foam::JanevReactionRate::preEvaluate() const
{}


inline void Foam::JanevReactionRate::postEvaluate() const
{}


inline Foam::scalar Foam::JanevReactionRate::operator()
(
    const scalar p,
    const scalar T,
    const scalarField&,
    const label
) const
{
    scalar lta = A_;

    if (mag(beta_) > vSmall)
    {
        lta *= pow(T, beta_);
    }

    scalar expArg = 0;

    if (mag(Ta_) > vSmall)
    {
        expArg -= Ta_/T;
    }

    scalar lnT = log(T);

    for (int n=0; n<nb_; n++)
    {
        expArg += b_[n]*pow(lnT, n);
    }

    lta *= exp(expArg);

    return lta;
}


inline Foam::scalar Foam::JanevReactionRate::ddT
(
    const scalar p,
    const scalar T,
    const scalarField&,
    const label
) const
{
    scalar lta = A_;

    if (mag(beta_) > vSmall)
    {
        lta *= pow(T, beta_);
    }

    scalar expArg = 0;

    if (mag(Ta_) > vSmall)
    {
        expArg -= Ta_/T;
    }

    scalar lnT = log(T);

    for (int n=0; n<nb_; n++)
    {
        expArg += b_[n]*pow(lnT, n);
    }

    scalar deriv = b_[1];

    for (int n=2; n<nb_; n++)
    {
        deriv += n*b_[n]*pow(lnT, n-1);
    }

    lta *= exp(expArg);

    return lta*(beta_+Ta_/T+deriv)/T;
}


inline bool Foam::JanevReactionRate::hasDdc() const
{
    return false;
}


inline void Foam::JanevReactionRate::ddc
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li,
    scalarField& ddc
) const
{
    ddc = 0;
}


inline void Foam::JanevReactionRate::write(Ostream& os) const
{
    writeKeyword(os, "A") << A_ << nl;
    writeKeyword(os, "beta") << beta_ << nl;
    writeKeyword(os, "Ta") << Ta_ << nl;
    writeKeyword(os, "b") << b_ << nl;
}


inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const JanevReactionRate& jrr
)
{
    jrr.write(os);
    return os;
}


// ************************************************************************* //
