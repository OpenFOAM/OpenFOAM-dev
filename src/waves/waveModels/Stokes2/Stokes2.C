/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenFOAM Foundation
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

#include "Stokes2.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace waveModels
{
    defineTypeNameAndDebug(Stokes2, 0);
    addToRunTimeSelectionTable(waveModel, Stokes2, objectRegistry);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::waveModels::Stokes2::Stokes2
(
    const objectRegistry& db,
    const dictionary& dict
)
:
    Airy(db, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::waveModels::Stokes2::~Stokes2()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::waveModels::Stokes2::elevation
(
    const scalar t,
    const scalar u,
    const scalarField& x
) const
{
    const scalar kh = k()*depth();
    const scalarField correction(k()*sqr(amplitude(t))*cos(2*angle(t, u, x))/2);

    if (shallow())
    {
        const scalar factor((3/sqr(tanh(kh)) - 1)/tanh(kh)/2);
        return Airy::elevation(t, u, x) + factor*correction;
    }
    else
    {
        return Airy::elevation(t, u, x) + correction;
    }
}


Foam::tmp<Foam::vector2DField> Foam::waveModels::Stokes2::velocity
(
    const scalar t,
    const scalar u,
    const vector2DField& xz
) const
{
    if (shallow())
    {
        const scalarField x(xz.component(0));
        const scalarField z(xz.component(1));

        const scalarField phi(angle(t, u, x));
        const scalarField kz(- k()*mag(z - elevation(t, u, x)));
        const scalar kh = k()*depth();
        const scalarField khz((kh + kz)*pos0(kh + kz));
        const scalar kwaa = k()*omega(u)*sqr(amplitude(t));

        return
            Airy::velocity(t, u, xz)
          + 6*kwaa*zip
            (
                cosh(2*khz)*cos(2*phi),
                sinh(2*khz)*sin(2*phi)
            )/pow4(sinh(kh));
    }
    else
    {
        return Airy::velocity(t, u, xz);
    }
}


// ************************************************************************* //
