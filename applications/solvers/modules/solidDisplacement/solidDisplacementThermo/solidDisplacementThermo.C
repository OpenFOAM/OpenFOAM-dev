/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2023 OpenFOAM Foundation
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

#include "solidDisplacementThermo.H"

/* * * * * * * * * * * * * * * Private Static Data * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(solidDisplacementThermo, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidDisplacementThermo::solidDisplacementThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    constSolidThermo(mesh, phaseName),
    planeStress_(lookup("planeStress")),
    thermalStress_(lookup("thermalStress")),
    E_(readProperty<scalar>("E", dimPressure)),
    nu_(readProperty<scalar>("nu", dimless)),
    alphav_(readProperty<scalar>("alphav", dimless/dimTemperature))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidDisplacementThermo::~solidDisplacementThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::volScalarField& Foam::solidDisplacementThermo::E() const
{
    return E_;
}


const Foam::scalarField& Foam::solidDisplacementThermo::E
(
    const label patchi
) const
{
    return E_.boundaryField()[patchi];
}


const Foam::volScalarField& Foam::solidDisplacementThermo::nu() const
{
    return nu_;
}


const Foam::scalarField& Foam::solidDisplacementThermo::nu
(
    const label patchi
) const
{
    return nu_.boundaryField()[patchi];
}


const Foam::volScalarField& Foam::solidDisplacementThermo::alphav() const
{
    return alphav_;
}


const Foam::scalarField& Foam::solidDisplacementThermo::alphav
(
    const label patchi
) const
{
    return alphav_.boundaryField()[patchi];
}


// ************************************************************************* //
