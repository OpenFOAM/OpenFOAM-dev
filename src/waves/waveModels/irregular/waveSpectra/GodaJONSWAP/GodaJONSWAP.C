/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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

#include "GodaJONSWAP.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace waveSpectra
{
    defineTypeNameAndDebug(GodaJONSWAP, 0);
    addToRunTimeSelectionTable(waveSpectrum, GodaJONSWAP, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::waveSpectra::GodaJONSWAP::GodaJONSWAP
(
    const GodaJONSWAP& spectrum
)
:
    waveSpectrum(spectrum),
    Hs_(spectrum.Hs_),
    Ts_(spectrum.Ts_),
    gamma_(spectrum.gamma_)
{}


Foam::waveSpectra::GodaJONSWAP::GodaJONSWAP
(
    const dictionary& dict,
    const scalar g
)
:
    waveSpectrum(dict, g),
    Hs_(dict.lookup<scalar>("Hs")),
    Ts_(dict.lookup<scalar>("Ts")),
    gamma_(dict.lookupOrDefault<scalar>("gamma", 3.3))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::waveSpectra::GodaJONSWAP::~GodaJONSWAP()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::waveSpectra::GodaJONSWAP::S
(
    const scalarField& f
) const
{
    const scalar betaJ =
        0.0624
       /(0.230 + 0.0336*gamma_ - 0.185/(1.9 + gamma_))
       *(1.094 - 0.01915*log(gamma_));
    const scalar Tp = Ts_/(1 - 0.132*pow(gamma_ + 0.2, -0.559));
    const scalarField sigma(0.07 + 0.02*pos(f - 1/Tp));
    const scalarField r(exp(- sqr(Tp*f - 1)/(2*sqr(sigma))));
    return betaJ*sqr(Hs_)/f/pow4(Tp*f)*exp(- 1.25/pow4(Tp*f))*pow(gamma_, r);
}


Foam::scalar Foam::waveSpectra::GodaJONSWAP::fFraction
(
    const scalar fraction
) const
{
    const scalar Tp = Ts_/(1 - 0.132*pow(gamma_ + 0.2, -0.559));
    const scalar f1 = 1/Tp/pow025(rootSmall/1.25);
    return waveSpectrum::fFraction(fraction, f1);
}


void Foam::waveSpectra::GodaJONSWAP::write(Ostream& os) const
{
    waveSpectrum::write(os);

    writeEntry(os, "Hs", Hs_);
    writeEntry(os, "Ts", Ts_);
    writeEntryIfDifferent(os, "gamma", scalar(3.3), gamma_);
}


// ************************************************************************* //
