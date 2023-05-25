/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2023 OpenFOAM Foundation
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
#include "fvcLaplacian.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::StationaryPhaseModel<BasePhaseModel>::StationaryPhaseModel
(
    const phaseSystem& fluid,
    const word& phaseName,
    const bool referencePhase,
    const label index
)
:
    BasePhaseModel(fluid, phaseName, referencePhase, index)
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
    return volVectorField::New
    (
        IOobject::groupName("U", this->name()),
        this->mesh(),
        dimensionedVector(dimVelocity, Zero)
    );
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
const Foam::volVectorField&
Foam::StationaryPhaseModel<BasePhaseModel>::URef() const
{
    FatalErrorInFunction
        << "Cannot access the velocity of a stationary phase"
        << exit(FatalError);

    return volVectorField::null();
}


template<class BasePhaseModel>
Foam::tmp<Foam::surfaceScalarField>
Foam::StationaryPhaseModel<BasePhaseModel>::phi() const
{
    return surfaceScalarField::New
    (
        IOobject::groupName("phi", this->name()),
        this->mesh(),
        dimensionedScalar(dimVolume/dimTime, 0)
    );
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
const Foam::surfaceScalarField&
Foam::StationaryPhaseModel<BasePhaseModel>::phiRef() const
{
    FatalErrorInFunction
        << "Cannot access the flux of a stationary phase"
        << exit(FatalError);

    return surfaceScalarField::null();
}


template<class BasePhaseModel>
const Foam::autoPtr<Foam::surfaceVectorField>&
Foam::StationaryPhaseModel<BasePhaseModel>::Uf() const
{
    FatalErrorInFunction
        << "Cannot access the face velocity of a stationary phase"
        << abort(FatalError);

    static autoPtr<Foam::surfaceVectorField> Uf_;
    return Uf_;
}


template<class BasePhaseModel>
Foam::surfaceVectorField&
Foam::StationaryPhaseModel<BasePhaseModel>::UfRef()
{
    FatalErrorInFunction
        << "Cannot access the face velocity of a stationary phase"
        << exit(FatalError);

    return const_cast<surfaceVectorField&>(surfaceVectorField::null());
}


template<class BasePhaseModel>
const Foam::surfaceVectorField&
Foam::StationaryPhaseModel<BasePhaseModel>::UfRef() const
{
    FatalErrorInFunction
        << "Cannot access the face velocity of a stationary phase"
        << exit(FatalError);

    return surfaceVectorField::null();
}


template<class BasePhaseModel>
Foam::tmp<Foam::surfaceScalarField>
Foam::StationaryPhaseModel<BasePhaseModel>::alphaPhi() const
{
    return surfaceScalarField::New
    (
        IOobject::groupName("alphaPhi", this->name()),
        this->mesh(),
        dimensionedScalar(dimVolume/dimTime, 0)
    );
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
const Foam::surfaceScalarField&
Foam::StationaryPhaseModel<BasePhaseModel>::alphaPhiRef() const
{
    FatalErrorInFunction
        << "Cannot access the volumetric flux of a stationary phase"
        << exit(FatalError);

    return surfaceScalarField::null();
}


template<class BasePhaseModel>
Foam::tmp<Foam::surfaceScalarField>
Foam::StationaryPhaseModel<BasePhaseModel>::alphaRhoPhi() const
{
    return surfaceScalarField::New
    (
        IOobject::groupName("alphaRhoPhi", this->name()),
        this->mesh(),
        dimensionedScalar(dimMass/dimTime, 0)
    );
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
const Foam::surfaceScalarField&
Foam::StationaryPhaseModel<BasePhaseModel>::alphaRhoPhiRef() const
{
    FatalErrorInFunction
        << "Cannot access the mass flux of a stationary phase"
        << exit(FatalError);

    return surfaceScalarField::null();
}


template<class BasePhaseModel>
Foam::tmp<Foam::volVectorField>
Foam::StationaryPhaseModel<BasePhaseModel>::DUDt() const
{
    return volVectorField::New
    (
        IOobject::groupName("DUDt", this->name()),
        this->mesh(),
        dimensionedVector(dimVelocity/dimTime, Zero)
    );
}


template<class BasePhaseModel>
Foam::tmp<Foam::surfaceScalarField>
Foam::StationaryPhaseModel<BasePhaseModel>::DUDtf() const
{
    return surfaceScalarField::New
    (
        IOobject::groupName("DUDtf", this->name()),
        this->mesh(),
        dimensionedScalar(dimVolume/sqr(dimTime), 0)
    );
}


template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::StationaryPhaseModel<BasePhaseModel>::continuityError() const
{
    return volScalarField::New
    (
        IOobject::groupName("continuityError", this->name()),
        this->mesh(),
        dimensionedScalar(dimDensity/dimTime, 0)
    );
}


template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::StationaryPhaseModel<BasePhaseModel>::K() const
{
    return volScalarField::New
    (
        IOobject::groupName("K", this->name()),
        this->mesh(),
        dimensionedScalar(sqr(dimVelocity), 0)
    );
}


template<class BasePhaseModel>
const Foam::autoPtr<Foam::volScalarField>&
Foam::StationaryPhaseModel<BasePhaseModel>::divU() const
{
    static autoPtr<volScalarField> divU_;
    return divU_;
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
Foam::tmp<Foam::scalarField>
Foam::StationaryPhaseModel<BasePhaseModel>::kappaEff(const label patchi) const
{
    return this->thermo().kappa().boundaryField()[patchi];
}


template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::StationaryPhaseModel<BasePhaseModel>::k() const
{
    return volScalarField::New
    (
        IOobject::groupName("k", this->name()),
        this->mesh(),
        dimensionedScalar(sqr(dimVelocity), 0)
    );
}


template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::StationaryPhaseModel<BasePhaseModel>::pPrime() const
{
    return volScalarField::New
    (
        IOobject::groupName("pPrime", this->name()),
        this->mesh(),
        dimensionedScalar(dimPressure, 0)
    );
}


template<class BasePhaseModel>
Foam::tmp<Foam::fvScalarMatrix>
Foam::StationaryPhaseModel<BasePhaseModel>::divq(volScalarField& he) const
{
    const volScalarField& alpha = *this;

    const surfaceScalarField alphaKappa
    (
        alpha.name() + '*' + this->thermo().kappa().name(),
        fvc::interpolate(alpha)*fvc::interpolate(this->thermo().kappa())
    );

    // Return heat flux source as an implicit energy correction
    // to the temperature gradient flux
    return
       -fvc::laplacian(alphaKappa, this->thermo().T())
       -fvm::laplacianCorrection
        (
            alphaKappa/fvc::interpolate(this->thermo().Cpv()),
            he
        );
}


// ************************************************************************* //
