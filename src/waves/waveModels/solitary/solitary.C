/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
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
    addToRunTimeSelectionTable(waveModel, solitary, objectRegistry);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::waveModels::solitary::k(const scalar t) const
{
    return sqrt(0.75*amplitude(t)/pow3(depth()));
}


Foam::scalar Foam::waveModels::solitary::alpha(const scalar t) const
{
    return amplitude(t)/depth();
}


Foam::scalar Foam::waveModels::solitary::celerity(const scalar t) const
{
    return sqrt(g()*depth()/(1 - alpha(t)));
}


Foam::tmp<Foam::scalarField> Foam::waveModels::solitary::parameter
(
    const scalar t,
    const scalar u,
    const scalarField& x
) const
{
    return k(t)*(x - offset_ - (u + celerity(t))*t);
}


Foam::tmp<Foam::scalarField> Foam::waveModels::solitary::Pi
(
    const scalar t,
    const scalar u,
    const scalarField& x
) const
{
    const scalar clip = log(great);

    return 1/sqr(cosh(max(- clip, min(clip, parameter(t, u, x)))));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::waveModels::solitary::solitary(const solitary& wave)
:
    waveModel(wave),
    offset_(wave.offset_),
    depth_(wave.depth_)
{}


Foam::waveModels::solitary::solitary
(
    const objectRegistry& db,
    const dictionary& dict
)
:
    waveModel(db, dict),
    offset_(readScalar(dict.lookup("offset"))),
    depth_(readScalar(dict.lookup("depth")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::waveModels::solitary::~solitary()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::waveModels::solitary::elevation
(
    const scalar t,
    const scalar u,
    const scalarField& x
) const
{
    return amplitude(t)*Pi(t, u, x);
}


Foam::tmp<Foam::vector2DField> Foam::waveModels::solitary::velocity
(
    const scalar t,
    const scalar u,
    const vector2DField& xz
) const
{
    const scalar A = alpha(t);
    const scalarField Z(max(scalar(0), 1 - mag(xz.component(1)/depth())));
    const scalarField P(Pi(t, u, xz.component(0)));

    return
        celerity(t)
       *zip
        (
            A/4
           *(
                (4 + 2*A - 6*A*sqr(Z))*P
              + (- 7*A + 9*A*sqr(Z))*sqr(P)
            ),
          - A*Z*depth()*k(t)*tanh(parameter(t, u, xz.component(0)))
           *(
                (2 + A - A*sqr(Z))*P
              + (- 7*A + 3*A*sqr(Z))*sqr(P)
            )
        );
}


Foam::tmp<Foam::scalarField> Foam::waveModels::solitary::pressure
(
    const scalar t,
    const scalar u,
    const vector2DField& xz
) const
{
    NotImplemented;
    return tmp<scalarField>(nullptr);
}


void Foam::waveModels::solitary::write(Ostream& os) const
{
    waveModel::write(os);

    os.writeKeyword("offset") << offset_ << token::END_STATEMENT << nl;
    os.writeKeyword("depth") << depth_ << token::END_STATEMENT << nl;
}


// ************************************************************************* //
