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

#include "coupledToConstantDensityFluid.H"
#include "dimensionSet.H"
#include "dimensionedScalar.H"
#include "physicalProperties.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace clouds
{
    defineTypeNameAndDebug(coupledToConstantDensityFluid, 0);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::physicalProperties&
Foam::clouds::coupledToConstantDensityFluid::physicalProperties() const
{
    if (!physicalPropertiesPtr_.valid())
    {
        physicalPropertiesPtr_.set
        (
            new Foam::physicalProperties(mesh_, word::null)
        );
    }

    return physicalPropertiesPtr_();
}


Foam::tmp<Foam::LagrangianSubScalarField>
Foam::clouds::coupledToConstantDensityFluid::one
(
    const LagrangianSubMesh& subMesh,
    const word& phaseName
)
{
    return
        LagrangianSubScalarField::New
        (
            subMesh.sub(IOobject::groupName(nameToCarrierName("1"), phaseName)),
            subMesh,
            dimensionedScalar(dimless, scalar(1))
        );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::clouds::coupledToConstantDensityFluid::coupledToConstantDensityFluid
(
    const cloud& c,
    const dictionary& dict
)
:
    coupled(c, dict),
    mesh_(c.mesh()),
    physicalPropertiesPtr_(nullptr),
    rho_
    (
        phaseName() == word::null
      ? NullObjectRef<dimensionedScalar>()
      : mesh_.mesh().lookupObject<uniformDimensionedScalarField>
        (
            IOobject::groupName("rho", phaseName())
        )
    ),
    rhoc_
    (
        carrierPhaseName() == word::null
      ? NullObjectRef<dimensionedScalar>()
      : mesh_.mesh().lookupObject<uniformDimensionedScalarField>
        (
            IOobject::groupName("rho", carrierPhaseName())
        )
    ),
    onec
    (
        c.derivedField<scalar>
        (
            [&](const LagrangianModelRef&, const LagrangianSubMesh& subMesh)
            {
                return one(subMesh, carrierPhaseName());
            }
        )
    ),
    onecPhase
    (
        c.derivedField<scalar>
        (
            [&](const LagrangianModelRef&, const LagrangianSubMesh& subMesh)
            {
                return one(subMesh, phaseName());
            }
        )
    ),
    rhoByRhoc
    (
        isNull(rho_) || isNull(rhoc_)
      ? dimensionedScalar("rhoByRhoc", dimless, physicalProperties())
      : rho_/rhoc_
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::clouds::coupledToConstantDensityFluid::~coupledToConstantDensityFluid()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::dimensionedScalar&
Foam::clouds::coupledToConstantDensityFluid::rho() const
{
    if (isNull(rho_))
    {
        FatalErrorInFunction
            << "Constant cloud density (rho) requested for cloud "
            << mesh_.name() << ", but only the density ratio (rhoByRhoc) is "
            << "defined. Constant cloud density is only available when the "
            << "Eulerian system is multiphase" << exit(FatalError);
    }

    return rho_;
}


const Foam::dimensionedScalar&
Foam::clouds::coupledToConstantDensityFluid::rhoc() const
{
    if (isNull(rhoc_))
    {
        FatalErrorInFunction
            << "Constant carrier density (rhoc) requested for cloud "
            << mesh_.name() << ", but only the density ratio (rhoByRhoc) is "
            << "defined. Constant carrier density is only available when the "
            << "Eulerian system is multiphase" << exit(FatalError);
    }

    return rhoc_;
}


const Foam::dimensionedScalar&
Foam::clouds::coupledToConstantDensityFluid::rhocPhase() const
{
    if (isNull(rho_))
    {
        FatalErrorInFunction
            << "Constant corresponding Eulerian phase density (rhocPhase) "
            << "requested for cloud " << mesh_.name() << ", but only the "
            << "density ratio (rhoByRhoc) is defined. Constant corresponding "
            << "Eulerian phase density is only available when the Eulerian "
            << "system is multiphase" << exit(FatalError);
    }

    return rho_;
}


// ************************************************************************* //
