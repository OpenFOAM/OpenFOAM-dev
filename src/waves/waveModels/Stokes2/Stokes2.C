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
    const scalar kd = k()*depth(), ka = k()*amplitude(t);

    const scalar T = deep() ? 1 : tanh(kd);

    const scalar B22 = (3/sqr(T) - 1)/T/4;

    if (debug)
    {
        Info<< "B22 = " << B22 << endl;
    }

    return
        Airy::elevation(t, u, x)
      + (1/k())*sqr(ka)*B22*cos(2*angle(t, u, x));
}


Foam::tmp<Foam::vector2DField> Foam::waveModels::Stokes2::velocity
(
    const scalar t,
    const scalar u,
    const vector2DField& xz
) const
{
    const scalar kd = k()*depth(), ka = k()*amplitude(t);

    const scalar A22ByA11 = deep() ? 0 : 0.375/pow3(sinh(kd));

    if (debug)
    {
        const scalar A11 = 1/sinh(kd);
        Info<< "A22 = " << A22ByA11*A11 << endl;
    }

    return
        Airy::velocity(t, u, xz)
      + celerity()*sqr(ka)*A22ByA11*vi(2, t, u, xz);
}


// ************************************************************************* //
