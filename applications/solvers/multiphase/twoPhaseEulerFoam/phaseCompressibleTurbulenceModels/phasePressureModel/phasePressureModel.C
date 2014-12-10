/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2014 OpenFOAM Foundation
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

#include "phasePressureModel.H"
#include "twoPhaseSystem.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RASModels::phasePressureModel::phasePressureModel
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& phase,
    const word& propertiesName,
    const word& type
)
:
    eddyViscosity<RASModel<PhaseCompressibleTurbulenceModel<phaseModel> > >
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        phase,
        propertiesName
    ),

    phase_(phase),

    alphaMax_(readScalar(this->coeffDict_.lookup("alphaMax"))),
    preAlphaExp_(readScalar(this->coeffDict_.lookup("preAlphaExp"))),
    expMax_(readScalar(this->coeffDict_.lookup("expMax"))),
    g0_
    (
        "g0",
        dimensionSet(1, -1, -2, 0, 0),
        this->coeffDict_.lookup("g0")
    )
{
    this->nut_ == dimensionedScalar("zero", this->nut_.dimensions(), 0.0);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::RASModels::phasePressureModel::~phasePressureModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::RASModels::phasePressureModel::read()
{
    if
    (
        eddyViscosity
        <
            RASModel<PhaseCompressibleTurbulenceModel<phaseModel> >
        >::read()
    )
    {
        this->coeffDict().lookup("alphaMax") >> alphaMax_;
        this->coeffDict().lookup("preAlphaExp") >> preAlphaExp_;
        this->coeffDict().lookup("expMax") >> expMax_;
        g0_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


Foam::tmp<Foam::volScalarField>
Foam::RASModels::phasePressureModel::k() const
{
    notImplemented("phasePressureModel::k()");
    return nut_;
}


Foam::tmp<Foam::volScalarField>
Foam::RASModels::phasePressureModel::epsilon() const
{
    notImplemented("phasePressureModel::epsilon()");
    return nut_;
}


Foam::tmp<Foam::volSymmTensorField>
Foam::RASModels::phasePressureModel::R() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                IOobject::groupName("R", this->U_.group()),
                this->runTime_.timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh_,
            dimensioned<symmTensor>
            (
                "R",
                dimensionSet(0, 2, -2, 0, 0),
                symmTensor::zero
            )
        )
    );
}


Foam::tmp<Foam::volScalarField>
Foam::RASModels::phasePressureModel::pPrime() const
{
    return
        g0_
       *min
        (
            exp(preAlphaExp_*(this->alpha_ - alphaMax_)),
            expMax_
        );
}


Foam::tmp<Foam::surfaceScalarField>
Foam::RASModels::phasePressureModel::pPrimef() const
{
    return
        g0_
       *min
        (
            exp(preAlphaExp_*(fvc::interpolate(this->alpha_) - alphaMax_)),
            expMax_
        );
}


Foam::tmp<Foam::volSymmTensorField>
Foam::RASModels::phasePressureModel::devRhoReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                IOobject::groupName("devRhoReff", this->U_.group()),
                this->runTime_.timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh_,
            dimensioned<symmTensor>
            (
                "R",
                this->rho_.dimensions()*dimensionSet(0, 2, -2, 0, 0),
                symmTensor::zero
            )
        )
    );
}


Foam::tmp<Foam::fvVectorMatrix>
Foam::RASModels::phasePressureModel::divDevRhoReff
(
    volVectorField& U
) const
{
    return tmp<fvVectorMatrix>
    (
        new fvVectorMatrix
        (
            U,
            this->rho_.dimensions()*dimensionSet(0, 4, -2, 0, 0)
        )
    );
}


void Foam::RASModels::phasePressureModel::correct()
{}


// ************************************************************************* //
