/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2019 OpenFOAM Foundation
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

#include "StationaryPhaseModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::StationaryPhaseModel<BasePhaseModel>::StationaryPhaseModel
(
    const phaseSystem& fluid,
    const word& phaseName,
    const label index
)
:
    BasePhaseModel(fluid, phaseName, index)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::StationaryPhaseModel<BasePhaseModel>::~StationaryPhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasePhaseModel>
bool Foam::StationaryPhaseModel<BasePhaseModel>::stationary() const
{
    return true;
}


template<class BasePhaseModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::StationaryPhaseModel<BasePhaseModel>::UEqn()
{
    FatalErrorInFunction
        << "Cannot construct a momentum equation for a stationary phase"
        << exit(FatalError);

    return tmp<fvVectorMatrix>();
}


template<class BasePhaseModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::StationaryPhaseModel<BasePhaseModel>::UfEqn()
{
    FatalErrorInFunction
        << "Cannot construct a momentum equation for a stationary phase"
        << exit(FatalError);

    return tmp<fvVectorMatrix>();
}


template<class BasePhaseModel>
Foam::tmp<Foam::volVectorField>
Foam::StationaryPhaseModel<BasePhaseModel>::U() const
{
    return zeroVolField<vector>(*this, "U", dimVelocity);
}


template<class BasePhaseModel>
Foam::volVectorField&
Foam::StationaryPhaseModel<BasePhaseModel>::URef()
{
    FatalErrorInFunction
        << "Cannot access the velocity of a stationary phase"
        << exit(FatalError);

    return const_cast<volVectorField&>(volVectorField::null());
}


template<class BasePhaseModel>
Foam::tmp<Foam::surfaceScalarField>
Foam::StationaryPhaseModel<BasePhaseModel>::phi() const
{
    return zeroSurfaceField<scalar>(*this, "phi", dimVolume/dimTime);
}


template<class BasePhaseModel>
Foam::surfaceScalarField&
Foam::StationaryPhaseModel<BasePhaseModel>::phiRef()
{
    FatalErrorInFunction
        << "Cannot access the flux of a stationary phase"
        << exit(FatalError);

    return const_cast<surfaceScalarField&>(surfaceScalarField::null());
}


template<class BasePhaseModel>
Foam::tmp<Foam::surfaceScalarField>
Foam::StationaryPhaseModel<BasePhaseModel>::alphaPhi() const
{
    return zeroSurfaceField<scalar>(*this, "alphaPhi", dimVolume/dimTime);
}


template<class BasePhaseModel>
Foam::surfaceScalarField&
Foam::StationaryPhaseModel<BasePhaseModel>::alphaPhiRef()
{
    FatalErrorInFunction
        << "Cannot access the volumetric flux of a stationary phase"
        << exit(FatalError);

    return const_cast<surfaceScalarField&>(surfaceScalarField::null());
}


template<class BasePhaseModel>
Foam::tmp<Foam::surfaceScalarField>
Foam::StationaryPhaseModel<BasePhaseModel>::alphaRhoPhi() const
{
    return zeroSurfaceField<scalar>(*this, "alphaRhoPhi", dimMass/dimTime);
}


template<class BasePhaseModel>
Foam::surfaceScalarField&
Foam::StationaryPhaseModel<BasePhaseModel>::alphaRhoPhiRef()
{
    FatalErrorInFunction
        << "Cannot access the mass flux of a stationary phase"
        << exit(FatalError);

    return const_cast<surfaceScalarField&>(surfaceScalarField::null());
}


template<class BasePhaseModel>
Foam::tmp<Foam::volVectorField>
Foam::StationaryPhaseModel<BasePhaseModel>::DUDt() const
{
    return zeroVolField<vector>(*this, "DUDt", dimVelocity/dimTime);
}


template<class BasePhaseModel>
Foam::tmp<Foam::surfaceScalarField>
Foam::StationaryPhaseModel<BasePhaseModel>::DUDtf() const
{
    return zeroSurfaceField<scalar>(*this, "DUDtf", dimVolume/sqr(dimTime));
}


template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::StationaryPhaseModel<BasePhaseModel>::continuityError() const
{
    return zeroVolField<scalar>(*this, "contErr", dimDensity/dimTime);
}


template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::StationaryPhaseModel<BasePhaseModel>::continuityErrorFlow() const
{
    return zeroVolField<scalar>(*this, "contErrFlow", dimDensity/dimTime);
}


template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::StationaryPhaseModel<BasePhaseModel>::continuityErrorSources() const
{
    return zeroVolField<scalar>(*this, "contErrSources", dimDensity/dimTime);
}


template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::StationaryPhaseModel<BasePhaseModel>::K() const
{
    return zeroVolField<scalar>(*this, "K", sqr(dimVelocity));
}


template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::StationaryPhaseModel<BasePhaseModel>::divU() const
{
    return tmp<volScalarField>();
}


template<class BasePhaseModel>
void Foam::StationaryPhaseModel<BasePhaseModel>::divU
(
    tmp<volScalarField> divU
)
{
    FatalErrorInFunction
        << "Cannot set the dilatation rate of a stationary phase"
        << exit(FatalError);
}


template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::StationaryPhaseModel<BasePhaseModel>::mut() const
{
    return zeroVolField<scalar>(*this, "mut", dimDynamicViscosity);
}


template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::StationaryPhaseModel<BasePhaseModel>::muEff() const
{
    return this->thermo().mu();
}


template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::StationaryPhaseModel<BasePhaseModel>::nut() const
{
    return zeroVolField<scalar>(*this, "nut", dimViscosity);
}


template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::StationaryPhaseModel<BasePhaseModel>::nuEff() const
{
    return this->thermo().nu();
}


template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::StationaryPhaseModel<BasePhaseModel>::kappaEff() const
{
    return this->thermo().kappa();
}


template<class BasePhaseModel>
Foam::tmp<Foam::scalarField>
Foam::StationaryPhaseModel<BasePhaseModel>::kappaEff(const label patchi) const
{
    return this->thermo().kappa(patchi);
}


template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::StationaryPhaseModel<BasePhaseModel>::alphaEff() const
{
    return this->thermo().alpha();
}


template<class BasePhaseModel>
Foam::tmp<Foam::scalarField>
Foam::StationaryPhaseModel<BasePhaseModel>::alphaEff(const label patchi) const
{
    return this->thermo().alpha(patchi);
}


template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::StationaryPhaseModel<BasePhaseModel>::k() const
{
    return zeroVolField<scalar>(*this, "k", sqr(dimVelocity));
}


template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::StationaryPhaseModel<BasePhaseModel>::pPrime() const
{
    return zeroVolField<scalar>(*this, "pPrime", dimPressure);
}


// ************************************************************************* //
