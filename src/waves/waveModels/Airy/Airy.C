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

#include "Airy.H"
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::waveModels::Airy::Airy
(
    const objectRegistry& db,
    const dictionary& dict
)
:
    waveModel(db, dict)
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
    const scalarField x(xz.component(0));
    const scalarField z(xz.component(1));

    const scalarField phi(angle(t, u, x));
    const scalarField kz(- k()*mag(z - elevation(t, u, x)));
    const scalar wa = omega(u)*amplitude(t);

    if (shallow())
    {
        const scalar kh = k()*depth();
        const scalarField khz((kh + kz)*pos0(kh + kz));
        return wa*zip(cosh(khz)*cos(phi), sinh(khz)*sin(phi))/sinh(kh);
    }
    else
    {
        return wa*exp(kz)*zip(cos(phi), sin(phi));
    }
}


// ************************************************************************* //
