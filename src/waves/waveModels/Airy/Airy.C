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

#include "Airy.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace waveModels
{
    defineTypeNameAndDebug(Airy, 0);
    addToRunTimeSelectionTable(waveModel, Airy, objectRegistry);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::waveModels::Airy::k() const
{
    return 2*Foam::constant::mathematical::pi/length_;
}


Foam::scalar Foam::waveModels::Airy::celerity() const
{
    return sqrt(g()/k()*tanh(k()*depth()));
}


Foam::tmp<Foam::scalarField> Foam::waveModels::Airy::angle
(
    const scalar t,
    const scalar u,
    const scalarField& x
) const
{
    return phase_ + k()*(x - (u + celerity())*t);
}


bool Foam::waveModels::Airy::deep() const
{
    return k()*depth() > log(great);
}


Foam::tmp<Foam::vector2DField> Foam::waveModels::Airy::vi
(
    const label i,
    const scalar t,
    const scalar u,
    const vector2DField& xz
) const
{
    const scalarField x(xz.component(0));
    const scalarField z(xz.component(1));

    const scalarField phi(angle(t, u, x));
    const scalarField kz(- k()*z);

    if (deep())
    {
        return i*exp(- mag(kz))*zip(cos(i*phi), sin(i*phi));
    }
    else
    {
        const scalar kd = k()*depth();
        const scalarField kdz(max(scalar(0), kd - mag(kz)));
        return i*zip(cosh(i*kdz)*cos(i*phi), sinh(i*kdz)*sin(i*phi))/sinh(kd);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::waveModels::Airy::Airy(const Airy& wave)
:
    waveModel(wave),
    length_(wave.length_),
    phase_(wave.phase_),
    depth_(wave.depth_)
{}


Foam::waveModels::Airy::Airy
(
    const objectRegistry& db,
    const dictionary& dict
)
:
    waveModel(db, dict),
    length_(readScalar(dict.lookup("length"))),
    phase_(readScalar(dict.lookup("phase"))),
    depth_(dict.lookupOrDefault<scalar>("depth", log(2*great)/k()))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::waveModels::Airy::~Airy()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::waveModels::Airy::elevation
(
    const scalar t,
    const scalar u,
    const scalarField& x
) const
{
    return amplitude(t)*cos(angle(t, u, x));
}


Foam::tmp<Foam::vector2DField> Foam::waveModels::Airy::velocity
(
    const scalar t,
    const scalar u,
    const vector2DField& xz
) const
{
    const scalar ka = k()*amplitude(t);

    return celerity()*ka*vi(1, t, u, xz);
}


Foam::tmp<Foam::scalarField> Foam::waveModels::Airy::pressure
(
    const scalar t,
    const scalar u,
    const vector2DField& xz
) const
{
    // It is a fluke of the formulation that the time derivative of the velocity
    // potential equals the x-derivative multiplied by the celerity. This allows
    // for this shortcut in evaluating the unsteady pressure.
    return celerity()*velocity(t, u, xz)->component(0);
}


void Foam::waveModels::Airy::write(Ostream& os) const
{
    waveModel::write(os);

    os.writeKeyword("length") << length_ << token::END_STATEMENT << nl;
    os.writeKeyword("phase") << phase_ << token::END_STATEMENT << nl;
    if (!deep())
    {
        os.writeKeyword("depth") << depth_ << token::END_STATEMENT << nl;
    }
}


// ************************************************************************* //
