/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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

#include "phaseModel.H"
#include "phaseSystem.H"
#include "diameterModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(phaseModel, 0);
    defineRunTimeSelectionTable(phaseModel, phaseSystem);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseModel::phaseModel
(
    const phaseSystem& fluid,
    const word& phaseName,
    const label index
)
:
    volScalarField
    (
        IOobject
        (
            IOobject::groupName("alpha", phaseName),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fluid.mesh(),
        dimensionedScalar("alpha", dimless, 0)
    ),

    fluid_(fluid),
    name_(phaseName),
    index_(index),
    residualAlpha_
    (
        "residualAlpha",
        dimless,
        fluid.subDict(phaseName).lookup("residualAlpha")
    ),
    alphaMax_(fluid.subDict(phaseName).lookupOrDefault("alphaMax", 1.0))
{
    diameterModel_ = diameterModel::New(fluid.subDict(phaseName), *this);
}


Foam::autoPtr<Foam::phaseModel> Foam::phaseModel::clone() const
{
    NotImplemented;
    return autoPtr<phaseModel>(nullptr);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseModel::~phaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::word& Foam::phaseModel::name() const
{
    return name_;
}


const Foam::word& Foam::phaseModel::keyword() const
{
    return name_;
}


Foam::label Foam::phaseModel::index() const
{
    return index_;
}


const Foam::phaseSystem& Foam::phaseModel::fluid() const
{
    return fluid_;
}


const Foam::dimensionedScalar& Foam::phaseModel::residualAlpha() const
{
    return residualAlpha_;
}


Foam::scalar Foam::phaseModel::alphaMax() const
{
    return alphaMax_;
}


Foam::tmp<Foam::volScalarField> Foam::phaseModel::d() const
{
    return diameterModel_().d();
}


const Foam::autoPtr<Foam::diameterModel>& Foam::phaseModel::dPtr() const
{
    return diameterModel_;
}


void Foam::phaseModel::correct()
{
    diameterModel_->correct();
}


void Foam::phaseModel::correctKinematics()
{}


void Foam::phaseModel::correctThermo()
{}


void Foam::phaseModel::correctTurbulence()
{}


void Foam::phaseModel::correctEnergyTransport()
{}


bool Foam::phaseModel::read()
{
    return diameterModel_->read(fluid_.subDict(name_));
}


bool Foam::phaseModel::compressible() const
{
    return false;
}


void Foam::phaseModel::correctInflowOutflow(surfaceScalarField& alphaPhi) const
{
    surfaceScalarField::Boundary& alphaPhiBf = alphaPhi.boundaryFieldRef();
    const volScalarField::Boundary& alphaBf = boundaryField();
    const surfaceScalarField::Boundary& phiBf = phi()().boundaryField();

    forAll(alphaPhiBf, patchi)
    {
        fvsPatchScalarField& alphaPhip = alphaPhiBf[patchi];

        if (!alphaPhip.coupled())
        {
            alphaPhip = phiBf[patchi]*alphaBf[patchi];
        }
    }
}


const Foam::tmp<Foam::volScalarField>& Foam::phaseModel::divU() const
{
    NotImplemented;
    static tmp<Foam::volScalarField> divU_(nullptr);
    return divU_;
}


void Foam::phaseModel::divU(const tmp<volScalarField>& divU)
{
    WarningInFunction
        << "Attempt to set the dilatation rate of an incompressible phase"
        << endl;
}


const Foam::volScalarField& Foam::phaseModel::K() const
{
    NotImplemented;
    return volScalarField::null();
}


// ************************************************************************* //
