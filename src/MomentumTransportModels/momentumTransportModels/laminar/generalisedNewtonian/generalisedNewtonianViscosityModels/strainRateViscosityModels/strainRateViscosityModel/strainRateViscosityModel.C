/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2024 OpenFOAM Foundation
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

#include "strainRateViscosityModel.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarModels
{
namespace generalisedNewtonianViscosityModels
{
    defineTypeNameAndDebug(strainRateViscosityModel, 0);
}
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::laminarModels::generalisedNewtonianViscosityModels::
strainRateViscosityModel::strainRate() const
{
    return sqrt(2.0)*mag(symm(fvc::grad(U_)));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::laminarModels::generalisedNewtonianViscosityModels::
strainRateViscosityModel::strainRateViscosityModel
(
    const dictionary& viscosityProperties,
    const viscosity& viscosity,
    const volVectorField& U
)
:
    generalisedNewtonianViscosityModel(viscosityProperties, viscosity, U),
    nu_
    (
        IOobject
        (
            IOobject::groupName
            (
                typedName("nu"),
                U.group()
            ),
            U.time().name(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar(dimKinematicViscosity, 0)
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::laminarModels::generalisedNewtonianViscosityModels::
strainRateViscosityModel::correct()
{
    nu_ = nu(viscosity_.nu(), strainRate());
}


// ************************************************************************* //
