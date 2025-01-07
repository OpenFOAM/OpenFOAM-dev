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

#include "totalPressureVelocityLagrangianVectorFieldSource.H"
#include "coupledToIncompressibleFluid.H"
#include "massive.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::totalPressureVelocityLagrangianVectorFieldSource::
totalPressureVelocityLagrangianVectorFieldSource
(
    const regIOobject& iIo,
    const dictionary& dict
)
:
    LagrangianVectorFieldSource(iIo, dict),
    Function1LagrangianFieldSource<vector>(*this),
    CloudLagrangianFieldSource<vector>(*this),
    dict_(new dictionary(dict)),
    direction_
    (
        Function1<vector>::New
        (
            "direction",
            iIo.time().userUnits(),
            dimless,
            dict
        )
    ),
    p0_(nullptr),
    rhocName_
    (
        dict.lookupOrDefault<word>
        (
            clouds::coupled::carrierName("rho"),
            "rho"
        )
    ),
    pcName_
    (
        dict.lookupOrDefault<word>
        (
            clouds::coupled::carrierName("p"),
            "p"
        )
    )
{}


Foam::totalPressureVelocityLagrangianVectorFieldSource::
totalPressureVelocityLagrangianVectorFieldSource
(
    const totalPressureVelocityLagrangianVectorFieldSource& field,
    const regIOobject& iIo
)
:
    LagrangianVectorFieldSource(field, iIo),
    Function1LagrangianFieldSource<vector>(*this),
    CloudLagrangianFieldSource<vector>(*this),
    dict_(field.dict_, false),
    direction_(field.direction_, false),
    p0_(field.p0_, false),
    rhocName_(field.rhocName_),
    pcName_(field.pcName_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::totalPressureVelocityLagrangianVectorFieldSource::
~totalPressureVelocityLagrangianVectorFieldSource()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::LagrangianSubVectorField>
Foam::totalPressureVelocityLagrangianVectorFieldSource::value
(
    const LagrangianInjection& injection,
    const LagrangianSubMesh& subMesh
) const
{
    // Get the carrier pressure
    const volScalarField& pcVf =
        subMesh.mesh().mesh().lookupObject<volScalarField>(pcName_);
    const CarrierField<scalar>& pc =
        cloud<clouds::coupled>(injection, subMesh).carrierField(pcVf);

    // Construct the total pressure function now we know the dimensions
    if (dict_.valid())
    {
        p0_ =
            Function1<scalar>::New
            (
                "p0",
                subMesh.mesh().time().userUnits(),
                pcVf.dimensions(),
                dict_()
            ).ptr();

        dict_.clear();
    }

    // Evaluate the direction the pressure drop
    const LagrangianSubVectorField direction
    (
        normalised(value(injection, subMesh, dimless, direction_()))
    );
    const LagrangianSubScalarField deltaP
    (
        max
        (
            value(injection, subMesh, pcVf.dimensions(), p0_())
          - pc(subMesh),
            dimensionedScalar(pcVf.dimensions(), scalar(0))
        )
    );

    // Return the direction multiplied by the velocity magnitude associated
    // with the dynamic head
    if (pcVf.dimensions() == dimKinematicPressure)
    {
        const clouds::coupledToIncompressibleFluid& ctifCloud =
            cloud<clouds::coupledToIncompressibleFluid>(injection, subMesh);

        return direction*sqrt(2*deltaP/ctifCloud.rhoByRhoc);
    }
    else if (pcVf.dimensions() == dimPressure)
    {
        assertCloud
        <
            clouds::coupledToIncompressibleFluid,
            clouds::massive
        >(injection, subMesh);

        if (isCloud<clouds::coupledToIncompressibleFluid>())
        {
            const clouds::coupledToIncompressibleFluid& ctifCloud =
                cloud<clouds::coupledToIncompressibleFluid>(injection, subMesh);

            // Get the carrier density
            const volScalarField& rhocVf =
                subMesh.mesh().mesh().lookupObject<volScalarField>(rhocName_);
            const CarrierField<scalar>& rhoc =
                cloud<clouds::coupled>(injection, subMesh).carrierField(rhocVf);

            return direction*sqrt(2*deltaP/(ctifCloud.rhoByRhoc*rhoc(subMesh)));
        }
        else // if (isCloud<clouds::massive>())
        {
            const clouds::massive& mCloud =
                cloud<clouds::massive>(injection, subMesh);

            return direction*sqrt(2*deltaP/mCloud.rho(subMesh));
        }
    }
    else
    {
        FatalErrorInFunction
            << "Dimensions of field " << pcVf.name()
            << " not recognised as pressure"
            << exit(FatalError);

        return tmp<LagrangianSubVectorField>(nullptr);
    }
}


void Foam::totalPressureVelocityLagrangianVectorFieldSource::write
(
    Ostream& os
) const
{
    if (dict_.valid())
    {
        dict_->write(os, false);
    }
    else
    {
        LagrangianVectorFieldSource::write(os);

        writeEntry(os, db().time().userUnits(), dimless, direction_());
        writeEntry(os, db().time().userUnits(), dimPressure, p0_());
        writeEntryIfDifferent<word>
        (
            os,
            clouds::coupled::carrierName("rho"),
            "rho",
            rhocName_
        );
        writeEntryIfDifferent<word>
        (
            os,
            clouds::coupled::carrierName("p"),
            "p",
            pcName_
        );
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeLagrangianTypeFieldSource
    (
        LagrangianVectorFieldSource,
        totalPressureVelocityLagrangianVectorFieldSource
    );
}

// ************************************************************************* //
