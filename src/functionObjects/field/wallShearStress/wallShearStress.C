/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2024 OpenFOAM Foundation
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
#include "incompressibleMomentumTransportModel.H"
#include "compressibleMomentumTransportModel.H"
#include "wallPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(wallShearStress, 0);
    addToRunTimeSelectionTable(functionObject, wallShearStress, dictionary);
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


Foam::tmp<Foam::volVectorField>
Foam::functionObjects::wallShearStress::calcShearStress
(
    const surfaceVectorField& tau
)
{
    tmp<volVectorField> twallShearStress
    (
        volVectorField::New
        (
            type(),
            mesh_,
            dimensionedVector(tau.dimensions(), Zero)
        )
    );

    volVectorField::Boundary& wallShearStressBf =
        twallShearStress.ref().boundaryFieldRef();

    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        label patchi = iter.key();

        const vectorField& Sfp = mesh_.Sf().boundaryField()[patchi];
        const scalarField& magSfp = mesh_.magSf().boundaryField()[patchi];

        wallShearStressBf[patchi] = -tau.boundaryField()[patchi];
    }

    return twallShearStress;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::wallShearStress::wallShearStress
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    logFiles(obr_, name),
    writeLocalObjects(obr_, log),
    phaseName_(word::null),
    patchSet_()
{
    read(dict);
    resetLocalObjectName(IOobject::groupName(type(), phaseName_));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::wallShearStress::~wallShearStress()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::wallShearStress::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeLocalObjects::read(dict);

    phaseName_ = dict.lookupOrDefault<word>("phase", word::null);

    patchSet_ = patchSet(dict, true);

    Info<< type() << " " << name() << ":" << nl;

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

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
            const label patchi = iter.key();
            filteredPatchSet.insert(patchi);

            if (isA<wallPolyPatch>(pbm[patchi]))
            {
                Info<< "        " << pbm[patchi].name() << endl;
            }
            else
            {
                Info<< "        " << pbm[patchi].name()
                    << "    type " << pbm[patchi].type() << endl;
            }
        }

        Info<< endl;

        patchSet_ = filteredPatchSet;
    }

    resetName(typeName);
    resetLocalObjectName(typeName);

    return true;
}


bool Foam::functionObjects::wallShearStress::execute()
{
    const word fieldName(IOobject::groupName(type(), phaseName_));

    typedef compressible::momentumTransportModel cmpModel;
    typedef incompressible::momentumTransportModel icoModel;

    const word momentumTransportModelName
    (
        IOobject::groupName(momentumTransportModel::typeName, phaseName_)
    );

    if (mesh_.foundObject<cmpModel>(momentumTransportModelName))
    {
        const cmpModel& model =
            mesh_.lookupObject<cmpModel>(momentumTransportModelName);

        return store(fieldName, calcShearStress(model.devTau()));
    }
    else if (mesh_.foundObject<icoModel>(momentumTransportModelName))
    {
        const icoModel& model =
            mesh_.lookupObject<icoModel>(momentumTransportModelName);

        return store(fieldName, calcShearStress(model.devSigma()));
    }
    else
    {
        FatalErrorInFunction
            << "Unable to find turbulence model in the "
            << "database" << exit(FatalError);

        return false;
    }
}


bool Foam::functionObjects::wallShearStress::write()
{
    Log << type() << " " << name() << " write:" << nl;

    writeLocalObjects::write();

    logFiles::write();

    const volVectorField& wallShearStress = obr_.lookupObject<volVectorField>
    (
        IOobject::groupName(type(), phaseName_)
    );

    const fvPatchList& patches = mesh_.boundary();

    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        label patchi = iter.key();
        const fvPatch& pp = patches[patchi];

        const vectorField& ssp = wallShearStress.boundaryField()[patchi];

        vector minSsp = gMin(ssp);
        vector maxSsp = gMax(ssp);

        if (Pstream::master())
        {
            file() << time_.userTimeValue()
                << tab << pp.name()
                << tab << minSsp
                << tab << maxSsp
                << endl;
        }

        Log << "    min/max(" << pp.name() << ") = "
            << minSsp << ", " << maxSsp << endl;
    }

    Log << endl;

    return true;
}


// ************************************************************************* //
