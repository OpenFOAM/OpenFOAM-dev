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

#include "totalPressureVelocityMagnitudeLagrangianScalarFieldSource.H"
#include "cloudLagrangianFieldSource.H"
#include "coupledToIncompressibleFluid.H"
#include "Function1LagrangianFieldSource.H"
#include "LagrangianFieldSource.H"
#include "massive.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::totalPressureVelocityMagnitudeLagrangianScalarFieldSource::
totalPressureVelocityMagnitudeLagrangianScalarFieldSource
(
    const LagrangianFieldSourceBase& field,
    const cloudLagrangianFieldSource& cloudField,
    const dictionary& dict
)
:
    field_(field),
    cloudField_(cloudField),
    p0Entry_(dict.lookupEntry("p0", false, true).clone()),
    p0_(nullptr),
    rhocName_
    (
        dict.lookupOrDefault<word>
        (
            clouds::carried::carrierName("rho"),
            "rho"
        )
    ),
    pcName_
    (
        dict.lookupOrDefault<word>
        (
            clouds::carried::carrierName("p"),
            "p"
        )
    )
{}


Foam::totalPressureVelocityMagnitudeLagrangianScalarFieldSource::
totalPressureVelocityMagnitudeLagrangianScalarFieldSource
(
    const totalPressureVelocityMagnitudeLagrangianScalarFieldSource& tpvmlsfs,
    const LagrangianFieldSourceBase& field,
    const cloudLagrangianFieldSource& cloudField
)
:
    field_(field),
    cloudField_(cloudField),
    p0Entry_(tpvmlsfs.p0Entry_, false),
    p0_(tpvmlsfs.p0_, false),
    rhocName_(tpvmlsfs.rhocName_),
    pcName_(tpvmlsfs.pcName_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::totalPressureVelocityMagnitudeLagrangianScalarFieldSource::
~totalPressureVelocityMagnitudeLagrangianScalarFieldSource()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::LagrangianSubScalarField>
Foam::totalPressureVelocityMagnitudeLagrangianScalarFieldSource::Umag
(
    const LagrangianInjection& injection,
    const LagrangianSubMesh& subMesh
) const
{
    // Get the carrier pressure
    const volScalarField& pcVf =
        subMesh.mesh().mesh().lookupObject<volScalarField>(pcName_);
    const CarrierField<scalar>& pc =
        cloudField_
       .cloud<clouds::carried>(injection, subMesh)
       .carrierField(pcVf);

    // Construct the total pressure function now we know the dimensions
    if (p0Entry_.valid())
    {
        p0_ =
            Function1<scalar>::New
            (
                field_.db().time().userUnits(),
                pcVf.dimensions(),
                p0Entry_()
            ).ptr();

        p0Entry_.clear();
    }

    // Evaluate the pressure drop
    const LagrangianSubScalarField deltaP
    (
        max
        (
            Function1LagrangianFieldSource::value
            (
                injection,
                subMesh,
                pcVf.dimensions(),
                p0_()
            ) - pc(subMesh),
            dimensionedScalar(pcVf.dimensions(), scalar(0))
        )
    );

    // Return the direction multiplied by the velocity magnitude associated
    // with the dynamic head
    if (pcVf.dimensions() == dimKinematicPressure)
    {
        const clouds::coupledToIncompressibleFluid& ctifCloud =
            cloudField_
           .cloud<clouds::coupledToIncompressibleFluid>(injection, subMesh);

        return sqrt(2*deltaP/ctifCloud.rhoByRhoc);
    }
    else if (pcVf.dimensions() == dimPressure)
    {
        cloudField_.assertCloud
        <
            clouds::coupledToIncompressibleFluid,
            clouds::massive
        >(injection, subMesh);

        if (cloudField_.isCloud<clouds::coupledToIncompressibleFluid>())
        {
            const clouds::coupledToIncompressibleFluid& ctifCloud =
                cloudField_
               .cloud<clouds::coupledToIncompressibleFluid>(injection, subMesh);

            // Get the carrier density
            const volScalarField& rhocVf =
                subMesh.mesh().mesh().lookupObject<volScalarField>(rhocName_);
            const CarrierField<scalar>& rhoc =
                cloudField_
               .cloud<clouds::carried>(injection, subMesh)
               .carrierField(rhocVf);

            return sqrt(2*deltaP/(ctifCloud.rhoByRhoc*rhoc(subMesh)));
        }
        else // if (isCloud<clouds::massive>())
        {
            const clouds::massive& mCloud =
                cloudField_.cloud<clouds::massive>(injection, subMesh);

            return sqrt(2*deltaP/mCloud.rho(subMesh));
        }
    }
    else
    {
        FatalErrorInFunction
            << "Dimensions of field " << pcVf.name()
            << " not recognised as pressure"
            << exit(FatalError);

        return tmp<LagrangianSubScalarField>(nullptr);
    }
}


void Foam::totalPressureVelocityMagnitudeLagrangianScalarFieldSource::write
(
    Ostream& os
) const
{
    if (p0Entry_.valid())
    {
        p0Entry_->write(os);
    }
    else
    {
        writeEntry(os, field_.db().time().userUnits(), unitAny, p0_());
    }

    writeEntryIfDifferent<word>
    (
        os,
        clouds::carried::carrierName("rho"),
        "rho",
        rhocName_
    );

    writeEntryIfDifferent<word>
    (
        os,
        clouds::carried::carrierName("p"),
        "p",
        pcName_
    );
}


// ************************************************************************* //
