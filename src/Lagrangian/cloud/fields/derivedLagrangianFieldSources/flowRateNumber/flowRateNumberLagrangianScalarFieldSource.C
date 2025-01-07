/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "flowRateNumberLagrangianScalarFieldSource.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::flowRateNumberLagrangianScalarFieldSource::
flowRateNumberLagrangianScalarFieldSource
(
    const regIOobject& iIo,
    const dictionary& dict
)
:
    uniformSizeNumberLagrangianScalarFieldSource(iIo, dict),
    volumetricFlowRate_
    (
        dict.found("volumetricFlowRate")
      ? Function1<scalar>::New
        (
            "volumetricFlowRate",
            db().time().userUnits(),
            dimVolumetricFlux,
            dict
        ).ptr()
      : nullptr
    ),
    massFlowRate_
    (
        dict.found("massFlowRate")
      ? Function1<scalar>::New
        (
            "massFlowRate",
            db().time().userUnits(),
            dimMassFlux,
            dict
        ).ptr()
      : nullptr
    )
{
    if (volumetricFlowRate_.valid() == massFlowRate_.valid())
    {
        FatalIOErrorInFunction(dict)
            << "keywords volumetricFlowRate and massFlowRate both "
            << (volumetricFlowRate_.valid() ? "" : "un") << "defined in "
            << "dictionary " << dict.name()
            << exit(FatalIOError);
    }
}


Foam::flowRateNumberLagrangianScalarFieldSource::
flowRateNumberLagrangianScalarFieldSource
(
    const flowRateNumberLagrangianScalarFieldSource& field,
    const regIOobject& iIo
)
:
    uniformSizeNumberLagrangianScalarFieldSource(field, iIo),
    volumetricFlowRate_(field.volumetricFlowRate_, false),
    massFlowRate_(field.massFlowRate_, false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::flowRateNumberLagrangianScalarFieldSource::
~flowRateNumberLagrangianScalarFieldSource()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::LagrangianSubScalarField>
Foam::flowRateNumberLagrangianScalarFieldSource::value
(
    const LagrangianInjection& injection,
    const LagrangianSubMesh& subMesh
) const
{
    // Get the range of the time-step
    const scalar t1 = db().time().value();
    const scalar t0 = t1 - db().time().deltaTValue();

    // Calculate the necessary sizes
    tmp<LagrangianSubScalarField> size, v, m;
    calcSizes
    (
        injection, subMesh,
        size,
        volumetricFlowRate_.valid(), v,
        massFlowRate_.valid(), m
    );

    // Return the numbers that equalise the sizes on all parcels and recover
    // the specified flow rate
    if (volumetricFlowRate_.valid())
    {
        const dimensionedScalar V
        (
            dimVolume,
            volumetricFlowRate_->integral(t0, t1)
        );

        return V/size()/sum(v()/size());
    }
    else
    {
        const dimensionedScalar M
        (
            dimMass,
            massFlowRate_->integral(t0, t1)
        );

        return M/size()/sum(m()/size());
    }
}


void Foam::flowRateNumberLagrangianScalarFieldSource::write(Ostream& os) const
{
    uniformSizeNumberLagrangianScalarFieldSource::write(os);

    if (volumetricFlowRate_.valid())
    {
        writeEntry(os, volumetricFlowRate_());
    }
    else
    {
        writeEntry(os, massFlowRate_());
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeLagrangianTypeFieldSource
    (
        LagrangianScalarFieldSource,
        flowRateNumberLagrangianScalarFieldSource
    );
}

// ************************************************************************* //
