/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

#include "PiersonMoskowitz.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace waveSpectra
{
    defineTypeNameAndDebug(PiersonMoskowitz, 0);
    addToRunTimeSelectionTable(waveSpectrum, PiersonMoskowitz, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::waveSpectra::PiersonMoskowitz::PiersonMoskowitz
(
    const PiersonMoskowitz& spectrum
)
:
    waveSpectrum(spectrum),
    U19_5_(spectrum.U19_5_),
    alpha_(spectrum.alpha_),
    beta_(spectrum.beta_)
{}


Foam::waveSpectra::PiersonMoskowitz::PiersonMoskowitz
(
    const dictionary& dict,
    const scalar g
)
:
    waveSpectrum(dict, g),
    U19_5_(dict.lookup<scalar>("U19_5")),
    alpha_(dict.lookupOrDefault<scalar>("alpha", 8.1e-3)),
    beta_(dict.lookupOrDefault<scalar>("beta", 0.74))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::waveSpectra::PiersonMoskowitz::~PiersonMoskowitz()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::waveSpectra::PiersonMoskowitz::S
(
    const scalarField& f
) const
{
    const scalarField w(twoPi*f);
    const scalar w0 = g()/U19_5_;
    return twoPi*alpha_*sqr(g())/pow5(w)*exp(- beta_*pow4(w0/w));
}


Foam::tmp<Foam::scalarField> Foam::waveSpectra::PiersonMoskowitz::integralS
(
    const scalarField& f
) const
{
    const scalarField w(twoPi*f);
    const scalar w0 = g()/U19_5_;
    const scalar integralSInf = twoPi*alpha_*sqr(g())/(4*beta_*pow4(w0))/twoPi;
    return integralSInf*exp(- beta_*pow4(w0/w));
}


Foam::scalar Foam::waveSpectra::PiersonMoskowitz::fFraction
(
    const scalar fraction
) const
{
    const scalar w0 = g()/U19_5_;
    return w0/pow025(- log(fraction)/beta_)/twoPi;
}


void Foam::waveSpectra::PiersonMoskowitz::write(Ostream& os) const
{
    waveSpectrum::write(os);

    writeEntry(os, "U19_5", U19_5_);
    writeEntryIfDifferent<scalar>(os, "alpha", 8.1e-3, alpha_);
    writeEntryIfDifferent<scalar>(os, "beta", 0.74, beta_);
}


// ************************************************************************* //
