/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2023 OpenFOAM Foundation
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

#include "wallHeatTransferCoeff.H"
#include "incompressibleMomentumTransportModel.H"
#include "compressibleMomentumTransportModel.H"
#include "fvsPatchField.H"
#include "basicThermo.H"
#include "wallPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(wallHeatTransferCoeff, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        wallHeatTransferCoeff,
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::wallHeatTransferCoeff::writeFileHeader
(
    const label i
)
{
    // Add headers to output data
    writeHeader(file(), "Wall heat transfer coefficient");
    writeCommented(file(), "Time");
    writeTabbed(file(), "patch");
    writeTabbed(file(), "min");
    writeTabbed(file(), "max");
    writeTabbed(file(), "average");
    file() << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::wallHeatTransferCoeff::wallHeatTransferCoeff
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    logFiles(obr_, name),
    writeLocalObjects(obr_, log),
    coeffModel_(wallHeatTransferCoeffModel::New(dict.name(), mesh_, dict)),
    rho_("rho", dimDensity, Zero),
    Cp_("Cp", dimArea/sqr(dimTime)/dimTemperature, Zero),
    runTime_(runTime),
    patchSet_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::wallHeatTransferCoeff::~wallHeatTransferCoeff()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::wallHeatTransferCoeff::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeLocalObjects::read(dict);

    if (!foundObject<basicThermo>(physicalProperties::typeName))
    {
        rho_.read(dict);
        Cp_.read(dict);
    }

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    patchSet_ =
        pbm.patchSet
        (
            wordReList(dict.lookupOrDefault("patches", wordReList()))
        );

    Info<< type() << ":" << nl;

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
                    << "Requested wall heat-transferCoeff on non-wall boundary"
                    << " type patch: " << pbm[patchi].name() << nl << endl;
            }
        }

        Info<< endl;

        patchSet_ = filteredPatchSet;
    }

    coeffModel_->read(dict);

    resetName(typeName);
    resetLocalObjectName(typeName);

    return true;
}


bool Foam::functionObjects::wallHeatTransferCoeff::execute()
{
    const momentumTransportModel& mmtm =
        lookupObject<momentumTransportModel>
        (
            momentumTransportModel::typeName
        );

    tmp<volScalarField> thtc;
    thtc = coeffModel_->htcByRhoCp(mmtm, patchSet_);

    if (!foundObject<basicThermo>(physicalProperties::typeName))
    {
        thtc.ref() *= rho_*Cp_;
    }
    else
    {
        const basicThermo& thermo =
            lookupObject<basicThermo>(physicalProperties::typeName);

        thtc.ref() *= thermo.rho()*thermo.Cp();
    }

    store("wallHeatTransferCoeff", thtc);

    return true;
}


bool Foam::functionObjects::wallHeatTransferCoeff::write()
{
    Log << name() << " write:" << nl;

    writeLocalObjects::write();
    logFiles::write();

    const volScalarField& htc = obr_.lookupObject<volScalarField>(type());

    const fvPatchList& patches = mesh_.boundary();

    const surfaceScalarField::Boundary& magSf =
        mesh_.magSf().boundaryField();

    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        label patchi = iter.key();
        const fvPatch& pp = patches[patchi];

        const scalarField& hfp = htc.boundaryField()[patchi];

        const scalar minHtcp = gMin(hfp);
        const scalar maxHtcp = gMax(hfp);
        const scalar averageHtcp =
            gSum(magSf[patchi]*hfp)/gSum(magSf[patchi]);

        if (Pstream::master())
        {
            file()
                << time_.userTimeValue()
                << tab << pp.name()
                << tab << minHtcp
                << tab << maxHtcp
                << tab << averageHtcp
                << endl;
        }

        Log << "    min/max/average(" << pp.name() << ") = "
            << minHtcp << ", " << maxHtcp << ", " << averageHtcp << endl;
    }

    Log << endl;

    return true;
}


// ************************************************************************* //
