/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "LaplacePressureLagrangianScalarFieldSource.H"
#include "coupledToThermalFluid.H"
#include "spherical.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LaplacePressureLagrangianScalarFieldSource::
LaplacePressureLagrangianScalarFieldSource
(
    const regIOobject& iIo,
    const dictionary& dict
)
:
    LagrangianScalarFieldSource(iIo, dict),
    cloudLagrangianFieldSource(*this),
    sigma_("sigma", dimForce/dimLength, dict)
{}


Foam::LaplacePressureLagrangianScalarFieldSource::
LaplacePressureLagrangianScalarFieldSource
(
    const LaplacePressureLagrangianScalarFieldSource& field,
    const regIOobject& iIo
)
:
    LagrangianScalarFieldSource(field, iIo),
    cloudLagrangianFieldSource(*this),
    sigma_(field.sigma_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::LaplacePressureLagrangianScalarFieldSource::
~LaplacePressureLagrangianScalarFieldSource()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::LagrangianSubScalarField>
Foam::LaplacePressureLagrangianScalarFieldSource::value
(
    const LagrangianInjection& injection,
    const LagrangianSubMesh& subMesh
) const
{
    return
        cloud<clouds::coupledToThermalFluid>().pc(subMesh)
      + 2*sigma_/cloud<clouds::spherical>().d(injection, subMesh);
}


Foam::tmp<Foam::LagrangianSubScalarField>
Foam::LaplacePressureLagrangianScalarFieldSource::value
(
    const LagrangianSubMesh& subMesh
) const
{
    return
        cloud<clouds::coupledToThermalFluid>().pc(subMesh)
      + 2*sigma_/cloud<clouds::spherical>().d(subMesh);
}


void Foam::LaplacePressureLagrangianScalarFieldSource::write(Ostream& os) const
{
    LagrangianScalarFieldSource::write(os);

    writeEntry(os, sigma_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeLagrangianTypeFieldSource
    (
        LagrangianScalarFieldSource,
        LaplacePressureLagrangianScalarFieldSource
    );
}

// ************************************************************************* //
