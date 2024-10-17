/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2024 OpenFOAM Foundation
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

#include "wallHeatFlux.H"
#include "thermophysicalTransportModel.H"
#include "surfaceInterpolate.H"
#include "fvcGrad.H"
#include "wallPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(wallHeatFlux, 0);
    addToRunTimeSelectionTable(functionObject, wallHeatFlux, dictionary);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::functionObjects::wallHeatFlux::writeFileHeader(const label i)
{
    // Add headers to output data
    writeHeader(file(), "Wall heat-flux");
    writeCommented(file(), "Time");
    writeTabbed(file(), "patch");
    writeTabbed(file(), "min [W/m^2]");
    writeTabbed(file(), "max [W/m^2]");
    writeTabbed(file(), "Q [W]");
    writeTabbed(file(), "q [W/m^2]");
    file() << endl;
}


Foam::tmp<Foam::volScalarField>
Foam::functionObjects::wallHeatFlux::calcWallHeatFlux
(
    const surfaceScalarField& q
)
{
    tmp<volScalarField> twallHeatFlux
    (
        volScalarField::New
        (
            type(),
            mesh_,
            dimensionedScalar(dimMass/pow3(dimTime), 0)
        )
    );

    volScalarField::Boundary& wallHeatFluxBf =
        twallHeatFlux.ref().boundaryFieldRef();

    const surfaceScalarField::Boundary& qBf = q.boundaryField();

    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        const label patchi = iter.key();

        wallHeatFluxBf[patchi] = -qBf[patchi];
    }

    if (foundObject<volScalarField>("qr"))
    {
        const volScalarField& qr = lookupObject<volScalarField>("qr");

        const volScalarField::Boundary& radHeatFluxBf = qr.boundaryField();

        forAllConstIter(labelHashSet, patchSet_, iter)
        {
            const label patchi = iter.key();

            wallHeatFluxBf[patchi] -= radHeatFluxBf[patchi];
        }
    }

    return twallHeatFlux;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::wallHeatFlux::wallHeatFlux
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

Foam::functionObjects::wallHeatFlux::~wallHeatFlux()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::wallHeatFlux::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeLocalObjects::read(dict);

    phaseName_ = dict.lookupOrDefault<word>("phase", word::null);

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    patchSet_ = patchSet(dict, true);

    Info<< type() << " " << name() << ":" << nl;

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
                    << "Requested wall heat-flux on non-wall boundary "
                    << "type patch: " << pbm[patchi].name() << endl;
            }
        }

        Info<< endl;

        patchSet_ = filteredPatchSet;
    }

    resetName(typeName);
    resetLocalObjectName(typeName);

    return true;
}


bool Foam::functionObjects::wallHeatFlux::execute()
{
    const word fieldName(IOobject::groupName(type(), phaseName_));

    const word thermophysicalTransportModelName
    (
        IOobject::groupName(thermophysicalTransportModel::typeName, phaseName_)
    );

    if
    (
        foundObject<thermophysicalTransportModel>
        (
            thermophysicalTransportModelName
        )
    )
    {
        const thermophysicalTransportModel& ttm =
            lookupObject<thermophysicalTransportModel>
            (
                thermophysicalTransportModelName
            );

        store(fieldName, calcWallHeatFlux(ttm.q()));

        return true;
    }
    else
    {
        FatalErrorInFunction
            << "Unable to find thermophysicalTransportModel "
            << thermophysicalTransportModelName
            << " in the database"
            << exit(FatalError);
    }

    return true;
}


bool Foam::functionObjects::wallHeatFlux::write()
{
    Log << type() << " " << name() << " write:" << nl;

    writeLocalObjects::write();

    logFiles::write();

    const volScalarField& wallHeatFlux = obr_.lookupObject<volScalarField>
    (
        IOobject::groupName(type(), phaseName_)
    );

    const fvPatchList& patches = mesh_.boundary();

    const surfaceScalarField::Boundary& magSf = mesh_.magSf().boundaryField();

    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        label patchi = iter.key();
        const fvPatch& pp = patches[patchi];

        const scalarField& qp = wallHeatFlux.boundaryField()[patchi];

        const scalar minqp = gMin(qp);
        const scalar maxqp = gMax(qp);
        const scalar Q = gSum(magSf[patchi]*qp);
        const scalar q = Q/gSum(magSf[patchi]);

        if (Pstream::master())
        {
            file()
                << time_.userTimeValue()
                << tab << pp.name()
                << tab << minqp
                << tab << maxqp
                << tab << Q
                << tab << q
                << endl;
        }

        Log << "    min, max, Q [W], q [W/m^2] for patch " << pp.name() << " = "
            << minqp << ", " << maxqp << ", " << Q << ", " << q << endl;
    }

    Log << endl;

    return true;
}


// ************************************************************************* //
