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

#include "JONSWAP.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace waveSpectra
{
    defineTypeNameAndDebug(JONSWAP, 0);
    addToRunTimeSelectionTable(waveSpectrum, JONSWAP, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::waveSpectra::JONSWAP::JONSWAP
(
    const JONSWAP& spectrum
)
:
    waveSpectrum(spectrum),
    U10_(spectrum.U10_),
    F_(spectrum.F_),
    gamma_(spectrum.gamma_)
{}


Foam::waveSpectra::JONSWAP::JONSWAP
(
    const dictionary& dict,
    const scalar g
)
:
    waveSpectrum(dict, g),
    U10_(dict.lookup<scalar>("U10")),
    F_(dict.lookup<scalar>("F")),
    gamma_(dict.lookupOrDefault<scalar>("gamma", 3.3))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::waveSpectra::JONSWAP::~JONSWAP()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::waveSpectra::JONSWAP::S
(
    const scalarField& f
) const
{
    const scalar alpha = 0.076*pow(sqr(U10_)/(F_*g()), 0.22);
    const scalarField w(twoPi*f);
    const scalar wp = 22*pow(sqr(g())/(U10_*F_), 1.0/3.0);
    const scalarField sigma(0.07 + 0.02*pos(w - wp));
    const scalarField r(exp(- sqr(w - wp)/(2*sqr(sigma)*sqr(wp))));
    return twoPi*alpha*sqr(g())/pow5(w)*exp(- 1.2*pow4(wp/w))*pow(gamma_, r);
}


Foam::scalar Foam::waveSpectra::JONSWAP::fFraction(const scalar fraction) const
{
    const scalar wp = 22*pow(sqr(g())/(U10_*F_), 1.0/3.0);
    const scalar f1 = wp/pow025(rootSmall/1.2)/twoPi;
    return waveSpectrum::fFraction(fraction, f1);
}


void Foam::waveSpectra::JONSWAP::write(Ostream& os) const
{
    waveSpectrum::write(os);

    writeEntry(os, "U10", U10_);
    writeEntry(os, "F", F_);
    writeEntryIfDifferent(os, "gamma", scalar(3.3), gamma_);
}


// ************************************************************************* //
