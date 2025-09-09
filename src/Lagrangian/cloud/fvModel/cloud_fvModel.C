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

#include "cloud_fvModel.H"
#include "coupledToConstantDensityFluid.H"
#include "fvmSup.H"
#include "pimpleNoLoopControl.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(cloud, 0);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::fv::cloud::fail(const fvMatrix<Type>& eqn) const
{
    FatalErrorInFunction
        << "Could not add a source for conservation of "
        << "1"
        << " in the Lagrangian cloud " << cloud_.name()
        << " to the finite-volume equation for " << eqn.psi().name()
        << exit(FatalError);
}


template<class Type, class ... AlphaRhoFieldTypes>
void Foam::fv::cloud::fail
(
    const fvMatrix<Type>& eqn,
    const AlphaRhoFieldTypes& ... alphaRhoFields
) const
{
    FatalErrorInFunction
        << "Could not add a source for conservation of "
        << LagrangianModel::fieldsName(alphaRhoFields ...)
        << " in the Lagrangian cloud " << cloud_.name()
        << " to the finite-volume equation for " << eqn.psi().name()
        << exit(FatalError);
}


Foam::tmp<Foam::volScalarField::Internal> Foam::fv::cloud::S
(
    const word& phaseName,
    const dimensionSet& dims
) const
{
    typedef HashTable<const CarrierEqn<scalar>*> carrierEqnTable;

    const bool isPhase = phaseName != word::null;
    const bool isMultiphase = coupledCloud_.carrierPhaseName() != word::null;

    const clouds::coupledToConstantDensityFluid& ctcdfCloud =
        refCastNull<const clouds::coupledToConstantDensityFluid>(cloud_);

    // Volume source
    if (!isPhase && !isMultiphase && dims == dimless && notNull(ctcdfCloud))
    {
        return
            ctcdfCloud.rhoByRhoc
           *coupledCloud_.carrierEqn<scalar>("1").residual(dimless/dimTime);
    }
    if (!isPhase && isMultiphase && dims == dimless && notNull(ctcdfCloud))
    {
        tmp<volScalarField::Internal> tS =
            volScalarField::Internal::New
            (
                "S",
                mesh(),
                dimensionedScalar(dimless/dimTime, scalar(0))
            );
        volScalarField::Internal& S = tS.ref();

        const carrierEqnTable carrierEqns =
            coupledCloud_.carrierEqns<scalar>("1");

        forAllConstIter(carrierEqnTable, carrierEqns, iter)
        {
            const uniformDimensionedScalarField& rhoPhase =
                mesh().lookupObject<uniformDimensionedScalarField>
                (
                    IOobject::groupName("rho", iter.key())
                );

            S +=
                ctcdfCloud.rho()/rhoPhase
               *iter()->residual(dimless/dimTime);
        }

        return tS;
    }

    // Phase volume source
    if (isPhase && dims == dimless && notNull(ctcdfCloud))
    {
        const uniformDimensionedScalarField& rhoPhase =
            mesh().lookupObject<uniformDimensionedScalarField>
            (
                IOobject::groupName("rho", phaseName)
            );

        return
            ctcdfCloud.rho()/rhoPhase
           *coupledCloud_
           .carrierEqn<scalar>(IOobject::groupName("1", phaseName))
           .residual(dimless/dimTime);
    }

    // Mass source
    if (!isPhase && dims == dimDensity)
    {
        tmp<volScalarField::Internal> tS =
            volScalarField::Internal::New
            (
                "S",
                mesh(),
                dimensionedScalar(dimDensity/dimTime, scalar(0))
            );
        volScalarField::Internal& S = tS.ref();

        const carrierEqnTable carrierEqns =
            coupledCloud_.carrierEqns<scalar>("rho");

        forAllConstIter(carrierEqnTable, carrierEqns, iter)
        {
            S += iter()->residual(dimDensity/dimTime);
        }

        return tS;
    }

    // Phase mass source
    if (isPhase && dims == dimDensity)
    {
        return
            coupledCloud_
           .carrierEqn<scalar>(IOobject::groupName("rho", phaseName))
           .residual(dimDensity/dimTime);
    }

    FatalError
        << "Could not create material source "
        << (isPhase ? "for phase " + phaseName : "").c_str()
        << " with dimensions " << dims << exit(FatalError);

    return tmp<volScalarField::Internal>(nullptr);
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fv::cloud::Sfield
(
    const VolField<Type>& field,
    const dimensionSet& dims
) const
{
    const volScalarField::Internal S(this->S(field.group(), dims));

    tmp<typename VolField<Type>::Internal> sourceCoeff =
        field.sources()[name()].sourceCoeff(*this, S);
    tmp<typename volScalarField::Internal> internalCoeff =
        field.sources()[name()].internalCoeff(*this, S);

    return S*sourceCoeff + fvm::Sp(S*internalCoeff, field);
}


template<class Type>
void Foam::fv::cloud::addSupType
(
    const VolField<Type>& field,
    fvMatrix<Type>& eqn
) const
{
    const word phaseName = field.group();

    const bool isPhase = phaseName != word::null;
    const bool isMultiphase = coupledCloud_.carrierPhaseName() != word::null;

    const bool hasCarrierEqn = coupledCloud_.hasCarrierEqn(field);

    const clouds::coupledToConstantDensityFluid& ctcdfCloud =
        refCastNull<const clouds::coupledToConstantDensityFluid>(cloud_);

    DebugInFunction
        << "field=" << field.name()
        << ", eqnField=" << eqn.psi().name()
        << ", hasCarrierEqn(" << field.name() << ")=" << hasCarrierEqn << endl;

    // Volume-weighted cloud source into a volume-weighted single-phase
    // property equation
    if (!isPhase && !isMultiphase && hasCarrierEqn && notNull(ctcdfCloud))
    {
        eqn += ctcdfCloud.rhoByRhoc*coupledCloud_.carrierEqn(field);
    }
    // Volume-weighted cloud source into a volume-weighted multiphase mixture
    // property equation
    else if (!isPhase && isMultiphase && hasCarrierEqn && notNull(ctcdfCloud))
    {
        // This is not possible. The parts of the source for the mixture
        // property that relate to different phases have already been combined.
        // To do this correctly would require separating them again and then
        // scaling them by their individual cloud-to-phase density ratios.
        fail(eqn, field);
    }
    // Generic material source into a volume-weighted property equation
    else if (!isPhase && !hasCarrierEqn)
    {
        eqn += Sfield(field, dimless);
    }
    // Not recognised
    else
    {
        fail(eqn, field);
    }
}


void Foam::fv::cloud::addSupType
(
    const volScalarField& alphaRhoOrField,
    fvMatrix<scalar>& eqn
) const
{
    const word phaseName = alphaRhoOrField.group();

    const bool isPhase = phaseName != word::null;
    const bool isAlpha =
        alphaRhoOrField.member() == "alpha"
     && alphaRhoOrField.dimensions() == dimless;
    const bool isRho =
        alphaRhoOrField.member() == "rho"
     && alphaRhoOrField.dimensions() == dimDensity;

    const word oneName = IOobject::groupName("1", phaseName);
    const word rhoName = IOobject::groupName("rho", phaseName);

    const clouds::coupledToConstantDensityFluid& ctcdfCloud =
        refCastNull<const clouds::coupledToConstantDensityFluid>(cloud_);

    DebugInFunction
        << "alphaRhoOrField=" << alphaRhoOrField.name()
        << ", eqnField=" << eqn.psi().name()
        << ", hasCarrierEqn(" << oneName << ")="
        << coupledCloud_.hasCarrierEqn<scalar>(oneName) << endl;

    // Mass continuity equation
    if (!isPhase && isRho)
    {
        eqn += coupledCloud_.carrierEqn<scalar>(rhoName);
    }
    // Phase volume continuity equation
    else if (isPhase && isAlpha && notNull(ctcdfCloud))
    {
        eqn +=
            ctcdfCloud.rho()
           /mesh().lookupObject<uniformDimensionedScalarField>(rhoName)
           *coupledCloud_.carrierEqn<scalar>(oneName);
    }
    // Try property equations
    else
    {
        addSupType<scalar>(alphaRhoOrField, eqn);
    }
}


template<class Type>
void Foam::fv::cloud::addSupType
(
    const volScalarField& alphaOrRho,
    const VolField<Type>& field,
    fvMatrix<Type>& eqn
) const
{
    const word phaseName = alphaOrRho.group();

    const bool isPhase = phaseName != word::null;
    const bool isAlpha =
        alphaOrRho.member() == "alpha"
     && alphaOrRho.dimensions() == dimless;
    const bool isRho =
        alphaOrRho.member() == "rho"
     && alphaOrRho.dimensions() == dimDensity;

    const bool hasCarrierEqn = coupledCloud_.hasCarrierEqn(field);

    const word oneName = IOobject::groupName("1", phaseName);
    const word rhoName = IOobject::groupName("rho", phaseName);

    const clouds::coupledToConstantDensityFluid& ctcdfCloud =
        refCastNull<const clouds::coupledToConstantDensityFluid>(cloud_);

    DebugInFunction
        << "alphaOrRho=" << alphaOrRho.name()
        << ", field=" << field.name()
        << ", eqnField=" << eqn.psi().name()
        << ", hasCarrierEqn(" << field.name() << ")=" << hasCarrierEqn << endl;

    // Mass-weighted cloud source into a mass-weighted single-phase property
    // equation
    if (!isPhase && isRho && hasCarrierEqn && isNull(ctcdfCloud))
    {
        eqn += coupledCloud_.carrierEqn(field);
    }
    // Volume-weighted cloud source into a mass-weighted multiphase mixture
    // property equation
    else if (!isPhase && isRho && hasCarrierEqn && notNull(ctcdfCloud))
    {
        eqn += ctcdfCloud.rho()*coupledCloud_.carrierEqn(field);
    }
    // Generic material source into a mass-weighted property equation
    else if (!isPhase && isRho && !hasCarrierEqn)
    {
        eqn += Sfield(field, dimDensity);
    }
    // Volume-weighted cloud source into a volume-weighted phase property
    // equation
    else if (isPhase && isAlpha && hasCarrierEqn && notNull(ctcdfCloud))
    {
        eqn +=
            ctcdfCloud.rho()
           /mesh().lookupObject<uniformDimensionedScalarField>(rhoName)
           *coupledCloud_.carrierEqn(field);
    }
    // Generic material source into a volume-weighted phase property equation
    else if (isPhase && isAlpha && !hasCarrierEqn && notNull(ctcdfCloud))
    {
        eqn += Sfield(field, dimless);
    }
    // Not recognised
    else
    {
        fail(eqn, alphaOrRho, field);
    }
}


void Foam::fv::cloud::addSupType
(
    const volScalarField& alphaOrRho,
    const volScalarField& rhoOrField,
    fvMatrix<scalar>& eqn
) const
{
    const word phaseName = alphaOrRho.group();

    const bool isPhase = phaseName != word::null;
    const bool isAlpha =
        alphaOrRho.member() == "alpha"
     && alphaOrRho.dimensions() == dimless;
    const bool isRhoField =
        rhoOrField.member() == "rho"
     && rhoOrField.dimensions() == dimDensity;

    DebugInFunction
        << "alphaOrRho=" << alphaOrRho.name()
        << ", rhoOrField=" << rhoOrField.name()
        << ", eqnField=" << eqn.psi().name()
        << ", hasCarrierEqn(" << rhoOrField.name() << ")="
        << coupledCloud_.hasCarrierEqn(rhoOrField) << endl;

    // Phase mass continuity equation
    if (isPhase && isAlpha && isRhoField)
    {
        eqn += coupledCloud_.carrierEqn(rhoOrField);
    }
    // Try property equations
    else
    {
        addSupType<scalar>(alphaOrRho, rhoOrField, eqn);
    }
}


template<class Type>
void Foam::fv::cloud::addSupType
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const VolField<Type>& field,
    fvMatrix<Type>& eqn
) const
{
    const bool hasCarrierEqn = coupledCloud_.hasCarrierEqn(field);

    DebugInFunction
        << "alpha=" << alpha.name()
        << ", rho=" << rho.name()
        << ", field=" << field.name()
        << ", eqnField=" << eqn.psi().name()
        << ", hasCarrierEqn(" << field.name() << ")=" << hasCarrierEqn << endl;

    // Mass-weighted cloud source into a mass-weighted phase property equation
    if (hasCarrierEqn)
    {
        eqn += coupledCloud_.carrierEqn(field);
    }
    // Generic material source into a mass-weighted phase property equation
    else
    {
        eqn += Sfield(field, dimDensity);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::cloud::~cloud()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::fv::cloud::addsSupToField(const word& fieldName) const
{
    const LagrangianModels& models = cloud_.LagrangianModels();

    const word phaseName = IOobject::group(fieldName);

    #define hasCarrierEqnType(Type, nullArg) \
        || coupledCloud_.hasCarrierEqn<Type>(fieldName)

    return
        models.addsSupToField
        (
            word::null,
            IOobject::groupName
            (
                clouds::carried::nameToCarrierName("1"),
                phaseName
            )
        )
     || models.addsSupToField
        (
            word::null,
            IOobject::groupName
            (
                clouds::carried::nameToCarrierName("rho"),
                phaseName
            )
        )
        FOR_ALL_FIELD_TYPES(hasCarrierEqnType);
}


void Foam::fv::cloud::addSup(fvMatrix<scalar>& eqn) const
{
    const word oneName = "1";

    const clouds::coupledToConstantDensityFluid& ctcdfCloud =
        refCastNull<const clouds::coupledToConstantDensityFluid>(cloud_);

    DebugInFunction
        << ", eqnField=" << eqn.psi().name()
        << ", hasCarrierEqn(" << oneName << ")="
        << coupledCloud_.hasCarrierEqn<scalar>(oneName) << endl;

    // Volume continuity equation
    if (notNull(ctcdfCloud))
    {
        eqn += ctcdfCloud.rhoByRhoc*coupledCloud_.carrierEqn<scalar>(oneName);
    }
    // Not recognised
    else
    {
        fail(eqn);
    }
}


FOR_ALL_FIELD_TYPES(IMPLEMENT_FV_MODEL_ADD_FIELD_SUP, fv::cloud)


FOR_ALL_FIELD_TYPES(IMPLEMENT_FV_MODEL_ADD_RHO_FIELD_SUP, fv::cloud)


FOR_ALL_FIELD_TYPES(IMPLEMENT_FV_MODEL_ADD_ALPHA_RHO_FIELD_SUP, fv::cloud)


void Foam::fv::cloud::correct()
{
    if (mesh().foundObject<pimpleNoLoopControl>(solutionControl::typeName))
    {
        const pimpleNoLoopControl& pimple =
            mesh().lookupObject<pimpleNoLoopControl>(solutionControl::typeName);

        const bool outerCorrectors =
            cloud_.mesh().solution().lookup<bool>("outerCorrectors");

        if (pimple.firstIter() || outerCorrectors)
        {
            cloudPtr_->solve(pimple.firstIter(), pimple.finalIter());
        }
    }
    else
    {
        cloudPtr_->solve(true, true);
    }
}


void Foam::fv::cloud::preUpdateMesh()
{
    cloudPtr_->storePosition();
}


bool Foam::fv::cloud::movePoints()
{
    return true;
}


void Foam::fv::cloud::topoChange(const polyTopoChangeMap& map)
{
    cloudPtr_->topoChange(map);
}


void Foam::fv::cloud::mapMesh(const polyMeshMap& map)
{
    cloudPtr_->mapMesh(map);
}


void Foam::fv::cloud::distribute(const polyDistributionMap& map)
{
    cloudPtr_->distribute(map);
}


// ************************************************************************* //
