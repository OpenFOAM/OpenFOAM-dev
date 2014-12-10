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

#include "DESModelRegions.H"
#include "volFields.H"
#include "compressible/turbulenceModel/turbulenceModel.H"
#include "compressible/LES/DESModel/DESModel.H"
#include "incompressible/turbulenceModel/turbulenceModel.H"
#include "incompressible/LES/DESModel/DESModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(DESModelRegions, 0);
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::DESModelRegions::writeFileHeader(const label i)
{
    writeHeader(file(), "DES model region coverage (% volume)");

    writeCommented(file(), "Time");
    writeTabbed(file(), "LES");
    writeTabbed(file(), "RAS");
    file() << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::DESModelRegions::DESModelRegions
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    functionObjectFile(obr, name, typeName),
    name_(name),
    obr_(obr),
    active_(true),
    log_(false)
{
    // Check if the available mesh is an fvMesh, otherwise deactivate
    if (!isA<fvMesh>(obr_))
    {
        active_ = false;
        WarningIn
        (
            "DESModelRegions::DESModelRegions"
            "("
                "const word&, "
                "const objectRegistry&, "
                "const dictionary&, "
                "const bool"
            ")"
        )   << "No fvMesh available, deactivating " << name_ << nl
            << endl;
    }

    if (active_)
    {
        const fvMesh& mesh = refCast<const fvMesh>(obr_);

        volScalarField* DESModelRegionsPtr
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

        mesh.objectRegistry::store(DESModelRegionsPtr);
    }

    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::DESModelRegions::~DESModelRegions()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::DESModelRegions::read(const dictionary& dict)
{
    if (active_)
    {
        log_ = dict.lookupOrDefault<Switch>("log", false);
    }
}


void Foam::DESModelRegions::execute()
{
    typedef incompressible::turbulenceModel icoModel;
    typedef incompressible::DESModel icoDESModel;

    typedef compressible::turbulenceModel cmpModel;
    typedef compressible::DESModel cmpDESModel;

    if (active_)
    {
        functionObjectFile::write();

        const fvMesh& mesh = refCast<const fvMesh>(obr_);

        Info(log_)<< type() << " " << name_ <<  " output:" << nl;

        volScalarField& DESModelRegions =
            const_cast<volScalarField&>
            (
                mesh.lookupObject<volScalarField>(type())
            );


        label DESpresent = false;
        if (mesh.foundObject<icoModel>("turbulenceModel"))
        {
            const icoModel& model =
                mesh.lookupObject<icoModel>("turbulenceModel");

            if (isA<icoDESModel>(model))
            {
                const icoDESModel& des =
                    dynamic_cast<const icoDESModel&>(model);
                DESModelRegions == des.LESRegion();
                DESpresent = true;
            }
        }
        else if (mesh.foundObject<cmpModel>("turbulenceModel"))
        {
            const cmpModel& model =
                mesh.lookupObject<cmpModel>("turbulenceModel");

            if (isA<cmpDESModel>(model))
            {
                const cmpDESModel& des =
                    dynamic_cast<const cmpDESModel&>(model);
                DESModelRegions == des.LESRegion();
                DESpresent = true;
            }
        }

        if (DESpresent)
        {
            scalar prc =
                gSum(DESModelRegions.internalField()*mesh.V())
               /gSum(mesh.V())*100.0;

            if (Pstream::master() && log_)
            {
                file() << obr_.time().value()
                    << token::TAB << prc
                    << token::TAB << 100.0 - prc
                    << endl;
            }

            Info(log_)
                << "    LES = " << prc << " % (volume)" << nl
                << "    RAS = " << 100.0 - prc << " % (volume)" << endl;
        }
        else
        {
            Info(log_)<< "    No DES turbulence model found in database" << nl
                << endl;
        }
    }
}


void Foam::DESModelRegions::end()
{
    if (active_)
    {
        execute();
    }
}


void Foam::DESModelRegions::timeSet()
{
    // Do nothing
}


void Foam::DESModelRegions::write()
{
    if (log_)
    {
        const volScalarField& DESModelRegions =
            obr_.lookupObject<volScalarField>(type());


        Info<< type() << " " << name_ <<  " output:" << nl
            << "    writing field " << DESModelRegions.name() << nl
            << endl;

        DESModelRegions.write();
    }
}


// ************************************************************************* //
