/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2021 OpenFOAM Foundation
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

#include "solitary.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace waveModels
{
    defineTypeNameAndDebug(solitary, 0);
    addToRunTimeSelectionTable(waveModel, solitary, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::waveModels::solitary::k() const
{
    return sqrt(0.75*amplitude()/pow3(depth()));
}


Foam::scalar Foam::waveModels::solitary::alpha() const
{
    return amplitude()/depth();
}


Foam::scalar Foam::waveModels::solitary::celerity() const
{
    return sqrt(g()*depth()/(1 - alpha()));
}


Foam::tmp<Foam::scalarField> Foam::waveModels::solitary::parameter
(
    const scalar t,
    const scalarField& x
) const
{
    return k()*(x - offset() - celerity()*t);
}


Foam::tmp<Foam::scalarField> Foam::waveModels::solitary::Pi
(
    const scalar t,
    const scalarField& x
) const
{
    const scalar clip = log(great);

    return 1/sqr(cosh(max(- clip, min(clip, parameter(t, x)))));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::waveModels::solitary::solitary(const solitary& wave)
:
    waveModel(wave),
    depth_(wave.depth_),
    amplitude_(wave.amplitude_),
    offset_(wave.offset_)
{}


Foam::waveModels::solitary::solitary
(
    const dictionary& dict,
    const scalar g
)
:
    waveModel(dict, g),
    depth_(dict.lookup<scalar>("depth")),
    amplitude_(dict.lookup<scalar>("amplitude")),
    offset_(dict.lookup<scalar>("offset"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::waveModels::solitary::~solitary()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::waveModels::solitary::elevation
(
    const scalar t,
    const scalarField& x
) const
{
    return amplitude()*Pi(t, x);
}


Foam::tmp<Foam::vector2DField> Foam::waveModels::solitary::velocity
(
    const scalar t,
    const vector2DField& xz
) const
{
    const scalar A = alpha();
    const scalarField Z(max(scalar(0), 1 + xz.component(1)/depth()));
    const scalarField P(Pi(t, xz.component(0)));

    return
        celerity()
       *zip
        (
            A/4
           *(
                (4 + 2*A - 6*A*sqr(Z))*P
              + (- 7*A + 9*A*sqr(Z))*sqr(P)
            ),
            A*Z*depth()*k()*tanh(parameter(t, xz.component(0)))
           *(
                (2 + A - A*sqr(Z))*P
              + (- 7*A + 3*A*sqr(Z))*sqr(P)
            )
        );
}


void Foam::waveModels::solitary::write(Ostream& os) const
{
    waveModel::write(os);

    writeEntry(os, "offset", offset_);
    writeEntry(os, "amplitude", amplitude_);
    writeEntry(os, "depth", depth_);
}


// ************************************************************************* //
