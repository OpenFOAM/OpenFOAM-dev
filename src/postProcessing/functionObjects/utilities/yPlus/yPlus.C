/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2016 OpenFOAM Foundation
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

#include "yPlus.H"
#include "volFields.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(yPlus, 0);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::yPlus::writeFileHeader(const label i)
{
    writeHeader(file(), "y+ ()");

    writeCommented(file(), "Time");
    writeTabbed(file(), "patch");
    writeTabbed(file(), "min");
    writeTabbed(file(), "max");
    writeTabbed(file(), "average");
    file() << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::yPlus::yPlus
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    functionObjectFiles(obr, name, typeName),
    name_(name),
    obr_(obr),
    log_(true),
    phiName_("phi")
{
    if (!isA<fvMesh>(obr))
    {
        FatalErrorInFunction
            << "objectRegistry is not an fvMesh" << exit(FatalError);
    }

    const fvMesh& mesh = refCast<const fvMesh>(obr_);

    volScalarField* yPlusPtr
    (
        new volScalarField
        (
            IOobject
            (
                type(),
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("0", dimless, 0.0)
        )
    );

    mesh.objectRegistry::store(yPlusPtr);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::yPlus::~yPlus()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::functionObjects::yPlus::read(const dictionary& dict)
{
    log_ = dict.lookupOrDefault<Switch>("log", true);
    phiName_ = dict.lookupOrDefault<word>("phiName", "phi");
}


void Foam::functionObjects::yPlus::execute()
{
    typedef compressible::turbulenceModel cmpModel;
    typedef incompressible::turbulenceModel icoModel;

    functionObjectFiles::write();

    const fvMesh& mesh = refCast<const fvMesh>(obr_);

    volScalarField& yPlus =
        const_cast<volScalarField&>
        (
            mesh.lookupObject<volScalarField>(type())
        );

    if (log_) Info<< type() << " " << name_ << " output:" << nl;

    tmp<volSymmTensorField> Reff;
    if (mesh.foundObject<cmpModel>(turbulenceModel::propertiesName))
    {
        const cmpModel& model =
            mesh.lookupObject<cmpModel>(turbulenceModel::propertiesName);

        calcYPlus(model, mesh, yPlus);
    }
    else if (mesh.foundObject<icoModel>(turbulenceModel::propertiesName))
    {
        const icoModel& model =
            mesh.lookupObject<icoModel>(turbulenceModel::propertiesName);

        calcYPlus(model, mesh, yPlus);
    }
    else
    {
        FatalErrorInFunction
            << "Unable to find turbulence model in the "
            << "database" << exit(FatalError);
    }
}


void Foam::functionObjects::yPlus::end()
{
    execute();
}


void Foam::functionObjects::yPlus::timeSet()
{}


void Foam::functionObjects::yPlus::write()
{
    functionObjectFiles::write();

    const volScalarField& yPlus =
        obr_.lookupObject<volScalarField>(type());

    if (log_) Info<< "    writing field " << yPlus.name() << nl << endl;

    yPlus.write();
}


// ************************************************************************* //
