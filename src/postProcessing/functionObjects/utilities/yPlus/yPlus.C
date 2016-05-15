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
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(yPlus, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        yPlus,
        dictionary
    );
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
    const Time& runTime,
    const dictionary& dict
)
:
    writeFiles(name, runTime, dict, name),
    phiName_("phi")
{
    if (!isA<fvMesh>(obr_))
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

    resetName(typeName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::yPlus::~yPlus()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::yPlus::read(const dictionary& dict)
{
    writeFiles::read(dict);
    phiName_ = dict.lookupOrDefault<word>("phiName", "phi");

    return true;
}


bool Foam::functionObjects::yPlus::execute(const bool postProcess)
{
    typedef compressible::turbulenceModel cmpModel;
    typedef incompressible::turbulenceModel icoModel;

    writeFiles::write();

    const fvMesh& mesh = refCast<const fvMesh>(obr_);

    volScalarField& yPlus =
        const_cast<volScalarField&>
        (
            mesh.lookupObject<volScalarField>(type())
        );

    if (log_) Info<< type() << " " << name() << " output:" << nl;

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

    return true;
}


bool Foam::functionObjects::yPlus::write(const bool postProcess)
{
    writeFiles::write();

    const volScalarField& yPlus =
        obr_.lookupObject<volScalarField>(type());

    if (log_) Info<< "    writing field " << yPlus.name() << nl << endl;

    yPlus.write();

    return true;
}


// ************************************************************************* //
