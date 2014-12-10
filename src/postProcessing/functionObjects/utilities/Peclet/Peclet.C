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

#include "Peclet.H"
#include "volFields.H"
#include "dictionary.H"
#include "surfaceFields.H"
#include "incompressible/turbulenceModel/turbulenceModel.H"
#include "compressible/turbulenceModel/turbulenceModel.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(Peclet, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Peclet::Peclet
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    name_(name),
    obr_(obr),
    active_(true),
    phiName_("phi"),
    rhoName_("rho")
{
    // Check if the available mesh is an fvMesh, otherwise deactivate
    if (!isA<fvMesh>(obr_))
    {
        active_ = false;
        WarningIn
        (
            "Peclet::Peclet"
            "("
                "const word&, "
                "const objectRegistry&, "
                "const dictionary&, "
                "const bool"
            ")"
        )   << "No fvMesh available, deactivating " << name_ << nl
            << endl;
    }

    read(dict);

    if (active_)
    {
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
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Peclet::~Peclet()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::Peclet::read(const dictionary& dict)
{
    if (active_)
    {
        phiName_ = dict.lookupOrDefault<word>("phiName", "phi");
        rhoName_ = dict.lookupOrDefault<word>("rhoName", "rho");
    }
}


void Foam::Peclet::execute()
{
    typedef compressible::turbulenceModel cmpTurbModel;
    typedef incompressible::turbulenceModel icoTurbModel;

    if (active_)
    {
        const fvMesh& mesh = refCast<const fvMesh>(obr_);

        tmp<volScalarField> nuEff;
        if (mesh.foundObject<cmpTurbModel>("turbulenceModel"))
        {
            const cmpTurbModel& model =
                mesh.lookupObject<cmpTurbModel>("turbulenceModel");

            const volScalarField& rho =
                mesh.lookupObject<volScalarField>(rhoName_);

            nuEff = model.muEff()/rho;
        }
        else if (mesh.foundObject<icoTurbModel>("turbulenceModel"))
        {
            const icoTurbModel& model =
                mesh.lookupObject<icoTurbModel>("turbulenceModel");

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
            FatalErrorIn("void Foam::Peclet::write()")
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
}


void Foam::Peclet::end()
{
    if (active_)
    {
        execute();
    }
}

void Foam::Peclet::timeSet()
{
    // Do nothing
}


void Foam::Peclet::write()
{
    if (active_)
    {
        const surfaceScalarField& Peclet =
            obr_.lookupObject<surfaceScalarField>(type());

        Info<< type() << " " << name_ << " output:" << nl
            << "    writing field " << Peclet.name() << nl
            << endl;

        Peclet.write();
    }
}


// ************************************************************************* //
