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

#include "wallShearStress.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"
#include "wallPolyPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(wallShearStress, 0);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::functionObjects::wallShearStress::writeFileHeader(const label i)
{
    // Add headers to output data
    writeHeader(file(), "Wall shear stress");
    writeCommented(file(), "Time");
    writeTabbed(file(), "patch");
    writeTabbed(file(), "min");
    writeTabbed(file(), "max");
    file() << endl;
}


void Foam::functionObjects::wallShearStress::calcShearStress
(
    const fvMesh& mesh,
    const volSymmTensorField& Reff,
    volVectorField& shearStress
)
{
    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        label patchi = iter.key();
        const polyPatch& pp = mesh.boundaryMesh()[patchi];

        vectorField& ssp = shearStress.boundaryFieldRef()[patchi];
        const vectorField& Sfp = mesh.Sf().boundaryField()[patchi];
        const scalarField& magSfp = mesh.magSf().boundaryField()[patchi];
        const symmTensorField& Reffp = Reff.boundaryField()[patchi];

        ssp = (-Sfp/magSfp) & Reffp;

        vector minSsp = gMin(ssp);
        vector maxSsp = gMax(ssp);

        if (Pstream::master())
        {
            file() << mesh.time().value()
                << token::TAB << pp.name()
                << token::TAB << minSsp
                << token::TAB << maxSsp
                << endl;
        }

        if (log_) Info<< "    min/max(" << pp.name() << ") = "
            << minSsp << ", " << maxSsp << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::wallShearStress::wallShearStress
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
    patchSet_()
{
    if (!isA<fvMesh>(obr))
    {
        FatalErrorInFunction
            << "objectRegistry is not an fvMesh" << exit(FatalError);
    }

    const fvMesh& mesh = refCast<const fvMesh>(obr_);

    volVectorField* wallShearStressPtr
    (
        new volVectorField
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
            dimensionedVector
            (
                "0",
                sqr(dimLength)/sqr(dimTime),
                Zero
            )
        )
    );

    mesh.objectRegistry::store(wallShearStressPtr);

    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::wallShearStress::~wallShearStress()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::functionObjects::wallShearStress::read(const dictionary& dict)
{
    log_ = dict.lookupOrDefault<Switch>("log", true);

    const fvMesh& mesh = refCast<const fvMesh>(obr_);
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    patchSet_ =
        mesh.boundaryMesh().patchSet
        (
            wordReList(dict.lookupOrDefault("patches", wordReList()))
        );

    Info<< type() << " " << name_ << ":" << nl;

    if (patchSet_.empty())
    {
        forAll(pbm, patchi)
        {
            if (isA<wallPolyPatch>(pbm[patchi]))
            {
                patchSet_.insert(patchi);
            }
        }

        Info<< "    processing all wall patches" << nl << endl;
    }
    else
    {
        Info<< "    processing wall patches: " << nl;
        labelHashSet filteredPatchSet;
        forAllConstIter(labelHashSet, patchSet_, iter)
        {
            label patchi = iter.key();
            if (isA<wallPolyPatch>(pbm[patchi]))
            {
                filteredPatchSet.insert(patchi);
                Info<< "        " << pbm[patchi].name() << endl;
            }
            else
            {
                WarningInFunction
                    << "Requested wall shear stress on non-wall boundary "
                    << "type patch: " << pbm[patchi].name() << endl;
            }
        }

        Info<< endl;

        patchSet_ = filteredPatchSet;
    }
}


void Foam::functionObjects::wallShearStress::execute()
{
    typedef compressible::turbulenceModel cmpModel;
    typedef incompressible::turbulenceModel icoModel;

    functionObjectFiles::write();

    const fvMesh& mesh = refCast<const fvMesh>(obr_);

    volVectorField& wallShearStress =
        const_cast<volVectorField&>
        (
            mesh.lookupObject<volVectorField>(type())
        );

    if (log_) Info<< type() << " " << name_ << " output:" << nl;


    tmp<volSymmTensorField> Reff;
    if (mesh.foundObject<cmpModel>(turbulenceModel::propertiesName))
    {
        const cmpModel& model =
            mesh.lookupObject<cmpModel>(turbulenceModel::propertiesName);

        Reff = model.devRhoReff();
    }
    else if (mesh.foundObject<icoModel>(turbulenceModel::propertiesName))
    {
        const icoModel& model =
            mesh.lookupObject<icoModel>(turbulenceModel::propertiesName);

        Reff = model.devReff();
    }
    else
    {
        FatalErrorInFunction
            << "Unable to find turbulence model in the "
            << "database" << exit(FatalError);
    }

    calcShearStress(mesh, Reff(), wallShearStress);
}


void Foam::functionObjects::wallShearStress::end()
{
    execute();
}


void Foam::functionObjects::wallShearStress::timeSet()
{}


void Foam::functionObjects::wallShearStress::write()
{
    functionObjectFiles::write();

    const volVectorField& wallShearStress =
        obr_.lookupObject<volVectorField>(type());

    if (log_) Info<< type() << " " << name_ << " output:" << nl
        << "    writing field " << wallShearStress.name() << nl
        << endl;

    wallShearStress.write();
}


// ************************************************************************* //
