/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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

#include "incompressibleTurbulenceModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(incompressibleTurbulenceModel, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::incompressibleTurbulenceModel::incompressibleTurbulenceModel
(
    const geometricOneField&,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const word& propertiesName
)
:
    turbulenceModel
    (
        U,
        alphaRhoPhi,
        phi,
        propertiesName
    )
{}


Foam::tmp<Foam::volScalarField>
Foam::incompressibleTurbulenceModel::mu() const
{
    return nu();
}


Foam::tmp<Foam::scalarField>
Foam::incompressibleTurbulenceModel::mu(const label patchi) const
{
    return nu(patchi);
}


Foam::tmp<Foam::volScalarField>
Foam::incompressibleTurbulenceModel::mut() const
{
    return nut();
}


Foam::tmp<Foam::scalarField>
Foam::incompressibleTurbulenceModel::mut(const label patchi) const
{
    return nut(patchi);
}


Foam::tmp<Foam::volScalarField>
Foam::incompressibleTurbulenceModel::muEff() const
{
    return nuEff();
}


Foam::tmp<Foam::scalarField>
Foam::incompressibleTurbulenceModel::muEff(const label patchi) const
{
    return nuEff(patchi);
}


// ************************************************************************* //
