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

#include "LandauTellerReactionRate.H"
#include "physicoChemicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::LandauTellerReactionRate::LandauTellerReactionRate
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


inline Foam::LandauTellerReactionRate::LandauTellerReactionRate
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
    B_(dict.lookup<scalar>("B", cbrt(dimTemperature))),
    C_(dict.lookup<scalar>("C", sqr(cbrt(dimTemperature))))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline void Foam::LandauTellerReactionRate::preEvaluate() const
{}


inline void Foam::LandauTellerReactionRate::postEvaluate() const
{}


inline Foam::scalar Foam::LandauTellerReactionRate::operator()
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

    if (mag(B_) > vSmall)
    {
        expArg += B_/cbrt(T);
    }

    if (mag(C_) > vSmall)
    {
        expArg += C_/sqr(cbrt(T));
    }

    if (mag(expArg) > vSmall)
    {
        lta *= exp(expArg);
    }

    return lta;
}


inline Foam::scalar Foam::LandauTellerReactionRate::ddT
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
    scalar deriv = 0;

    if (mag(Ta_) > vSmall)
    {
        scalar TaT = Ta_/T;
        expArg -= TaT;
        deriv += TaT;
    }

    if (mag(B_) > vSmall)
    {
        scalar BT = B_/cbrt(T);
        expArg += BT;
        deriv -= BT/3;
    }

    if (mag(C_) > vSmall)
    {
        scalar CT = C_/sqr(cbrt(T));
        expArg += CT;
        deriv -= 2*CT/3;
    }

    if (mag(expArg) > vSmall)
    {
        lta *= exp(expArg);
    }

    return lta*(beta_ + deriv)/T;
}


inline bool Foam::LandauTellerReactionRate::hasDdc() const
{
    return false;
}


inline void Foam::LandauTellerReactionRate::ddc
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


inline void Foam::LandauTellerReactionRate::write(Ostream& os) const
{
    writeEntry(os, "A", A_);
    writeEntry(os, "beta", beta_);
    writeEntry(os, "Ta", Ta_);
    writeEntry(os, "B", B_);
    writeEntry(os, "C", C_);
}


inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const LandauTellerReactionRate& arr
)
{
    arr.write(os);
    return os;
}


// ************************************************************************* //
