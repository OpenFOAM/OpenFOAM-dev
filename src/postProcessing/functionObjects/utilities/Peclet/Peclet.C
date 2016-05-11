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

#include "Peclet.H"
#include "volFields.H"
#include "dictionary.H"
#include "surfaceFields.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(Peclet, 0);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::Peclet::Peclet
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    name_(name),
    obr_(obr),
    phiName_("phi"),
    rhoName_("rho")
{
    if (!isA<fvMesh>(obr))
    {
        FatalErrorInFunction
            << "objectRegistry is not an fvMesh" << exit(FatalError);
    }

    read(dict);

    const fvMesh& mesh = refCast<const fvMesh>(obr_);

    surfaceScalarField* PecletPtr
    (
        new surfaceScalarField
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

    mesh.objectRegistry::store(PecletPtr);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::Peclet::~Peclet()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::functionObjects::Peclet::read(const dictionary& dict)
{
    phiName_ = dict.lookupOrDefault<word>("phiName", "phi");
    rhoName_ = dict.lookupOrDefault<word>("rhoName", "rho");
}


void Foam::functionObjects::Peclet::execute()
{
    typedef compressible::turbulenceModel cmpTurbModel;
    typedef incompressible::turbulenceModel icoTurbModel;

    const fvMesh& mesh = refCast<const fvMesh>(obr_);

    tmp<volScalarField> nuEff;
    if (mesh.foundObject<cmpTurbModel>(turbulenceModel::propertiesName))
    {
        const cmpTurbModel& model =
            mesh.lookupObject<cmpTurbModel>
            (
                turbulenceModel::propertiesName
            );

        const volScalarField& rho =
            mesh.lookupObject<volScalarField>(rhoName_);

        nuEff = model.muEff()/rho;
    }
    else if
    (
        mesh.foundObject<icoTurbModel>(turbulenceModel::propertiesName)
    )
    {
        const icoTurbModel& model =
            mesh.lookupObject<icoTurbModel>
            (
                turbulenceModel::propertiesName
            );

        nuEff = model.nuEff();
    }
    else if (mesh.foundObject<dictionary>("transportProperties"))
    {
        const dictionary& model =
            mesh.lookupObject<dictionary>("transportProperties");

        nuEff =
            tmp<volScalarField>
            (
                new volScalarField
                (
                    IOobject
                    (
                        "nuEff",
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh,
                    dimensionedScalar(model.lookup("nu"))
                )
            );
    }
    else
    {
        FatalErrorInFunction
            << "Unable to determine the viscosity"
            << exit(FatalError);
    }

    const surfaceScalarField& phi =
        mesh.lookupObject<surfaceScalarField>(phiName_);

    surfaceScalarField& Peclet =
        const_cast<surfaceScalarField&>
        (
            mesh.lookupObject<surfaceScalarField>(type())
        );

    Peclet =
        mag(phi)
       /(
            mesh.magSf()
           *mesh.surfaceInterpolation::deltaCoeffs()
           *fvc::interpolate(nuEff)
        );
}


void Foam::functionObjects::Peclet::end()
{
    execute();
}

void Foam::functionObjects::Peclet::timeSet()
{}


void Foam::functionObjects::Peclet::write()
{
    const surfaceScalarField& Peclet =
        obr_.lookupObject<surfaceScalarField>(type());

    Info<< type() << " " << name_ << " output:" << nl
        << "    writing field " << Peclet.name() << nl
        << endl;

    Peclet.write();
}


// ************************************************************************* //
