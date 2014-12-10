/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "radiationModel.H"
#include "absorptionEmissionModel.H"
#include "scatterModel.H"
#include "sootModel.H"
#include "fvmSup.H"
#include "fluidThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(radiationModel, 0);
        defineRunTimeSelectionTable(radiationModel, T);
        defineRunTimeSelectionTable(radiationModel, dictionary);
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::IOobject Foam::radiation::radiationModel::createIOobject
(
    const fvMesh& mesh
) const
{
    IOobject io
    (
        "radiationProperties",
        mesh.time().constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (io.headerOk())
    {
        io.readOpt() = IOobject::MUST_READ_IF_MODIFIED;
        return io;
    }
    else
    {
        io.readOpt() = IOobject::NO_READ;
        return io;
    }
}


void Foam::radiation::radiationModel::initialise()
{
    if (radiation_)
    {
        solverFreq_ = max(1, lookupOrDefault<label>("solverFreq", 1));

        absorptionEmission_.reset
        (
            absorptionEmissionModel::New(*this, mesh_).ptr()
        );

        scatter_.reset(scatterModel::New(*this, mesh_).ptr());

        soot_.reset(sootModel::New(*this, mesh_).ptr());
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::radiationModel::radiationModel(const volScalarField& T)
:
    IOdictionary
    (
        IOobject
        (
            "radiationProperties",
            T.time().constant(),
            T.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(T.mesh()),
    time_(T.time()),
    T_(T),
    radiation_(false),
    coeffs_(dictionary::null),
    solverFreq_(0),
    firstIter_(true),
    absorptionEmission_(NULL),
    scatter_(NULL),
    soot_(NULL)
{}


Foam::radiation::radiationModel::radiationModel
(
    const word& type,
    const volScalarField& T
)
:
    IOdictionary(createIOobject(T.mesh())),
    mesh_(T.mesh()),
    time_(T.time()),
    T_(T),
    radiation_(lookupOrDefault("radiation", true)),
    coeffs_(subOrEmptyDict(type + "Coeffs")),
    solverFreq_(1),
    firstIter_(true),
    absorptionEmission_(NULL),
    scatter_(NULL),
    soot_(NULL)
{
    if (readOpt() == IOobject::NO_READ)
    {
        radiation_ = false;
    }

    initialise();
}


Foam::radiation::radiationModel::radiationModel
(
    const word& type,
    const dictionary& dict,
    const volScalarField& T
)
:
    IOdictionary
    (
        IOobject
        (
            "radiationProperties",
            T.time().constant(),
            T.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dict
    ),
    mesh_(T.mesh()),
    time_(T.time()),
    T_(T),
    radiation_(lookupOrDefault("radiation", true)),
    coeffs_(subOrEmptyDict(type + "Coeffs")),
    solverFreq_(1),
    firstIter_(true),
    absorptionEmission_(NULL),
    scatter_(NULL),
    soot_(NULL)
{
    initialise();
}


// * * * * * * * * * * * * * * * * Destructor    * * * * * * * * * * * * * * //

Foam::radiation::radiationModel::~radiationModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::radiation::radiationModel::read()
{
    if (regIOobject::read())
    {
        lookup("radiation") >> radiation_;
        coeffs_ = subOrEmptyDict(type() + "Coeffs");

        solverFreq_ = lookupOrDefault<label>("solverFreq", 1);
        solverFreq_ = max(1, solverFreq_);

        return true;
    }
    else
    {
        return false;
    }
}


void Foam::radiation::radiationModel::correct()
{
    if (!radiation_)
    {
        return;
    }

    if (firstIter_ || (time_.timeIndex() % solverFreq_ == 0))
    {
        calculate();
        firstIter_ = false;
    }

    if (!soot_.empty())
    {
        soot_->correct();
    }
}


Foam::tmp<Foam::fvScalarMatrix> Foam::radiation::radiationModel::Sh
(
    fluidThermo& thermo
) const
{
    volScalarField& he = thermo.he();
    const volScalarField Cpv(thermo.Cpv());
    const volScalarField T3(pow3(T_));

    return
    (
        Ru()
      - fvm::Sp(4.0*Rp()*T3/Cpv, he)
      - Rp()*T3*(T_ - 4.0*he/Cpv)
    );
}


Foam::tmp<Foam::fvScalarMatrix> Foam::radiation::radiationModel::ST
(
    const dimensionedScalar& rhoCp,
    volScalarField& T
) const
{
    return
    (
        Ru()/rhoCp
      - fvm::Sp(Rp()*pow3(T)/rhoCp, T)
    );
}


const Foam::radiation::absorptionEmissionModel&
Foam::radiation::radiationModel::absorptionEmission() const
{
    if (!absorptionEmission_.valid())
    {
        FatalErrorIn
        (
            "const Foam::radiation::absorptionEmissionModel&"
            "Foam::radiation::radiationModel::absorptionEmission() const"
        )
            << "Requested radiation absorptionEmission model, but model is "
            << "not activate" << abort(FatalError);
    }

    return absorptionEmission_();
}


const Foam::radiation::sootModel&
Foam::radiation::radiationModel::soot() const
{
    if (!soot_.valid())
    {
        FatalErrorIn
        (
            "const Foam::radiation::sootModel&"
            "Foam::radiation::radiationModel::soot() const"
        )
            << "Requested radiation sootModel model, but model is "
            << "not activate" << abort(FatalError);
    }

    return soot_();
}


// ************************************************************************* //
