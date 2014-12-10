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

Application
    execFlowFunctionObjects

Description
    Execute the set of functionObjects specified in the selected dictionary
    (which defaults to system/controlDict) for the selected set of times.
    Alternative dictionaries should be placed in the system/ directory.

    The flow (p-U) and optionally turbulence fields are available for the
    function objects to operate on allowing forces and other related properties
    to be calculated in addition to cutting planes etc.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"

#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "uniformDimensionedFields.H"
#include "ReadFields.H"
#include "fvIOoptionList.H"

#include "incompressible/singlePhaseTransportModel/singlePhaseTransportModel.H"

#include "incompressible/RAS/RASModel/RASModel.H"
#include "incompressible/LES/LESModel/LESModel.H"

#include "fluidThermo.H"
#include "compressible/RAS/RASModel/RASModel.H"
#include "compressible/LES/LESModel/LESModel.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Read all fields of type. Returns names of fields read. Guarantees all
// processors to read fields in same order.
template<class GeoField>
wordList ReadUniformFields
(
    const IOobjectList& objects,
    PtrList<GeoField>& fields,
    const bool syncPar
)
{
    // Search list of objects for wanted type
    IOobjectList fieldObjects(objects.lookupClass(GeoField::typeName));

    wordList masterNames(fieldObjects.names());

    if (syncPar && Pstream::parRun())
    {
        // Check that I have the same fields as the master
        const wordList localNames(masterNames);
        Pstream::scatter(masterNames);

        HashSet<word> localNamesSet(localNames);

        forAll(masterNames, i)
        {
            const word& masterFld = masterNames[i];

            HashSet<word>::iterator iter = localNamesSet.find(masterFld);

            if (iter == localNamesSet.end())
            {
                FatalErrorIn
                (
                    "ReadFields<class GeoField>"
                    "(const IOobjectList&, PtrList<GeoField>&"
                    ", const bool)"
                )   << "Fields not synchronised across processors." << endl
                    << "Master has fields " << masterNames
                    << "  processor " << Pstream::myProcNo()
                    << " has fields " << localNames << exit(FatalError);
            }
            else
            {
                localNamesSet.erase(iter);
            }
        }

        forAllConstIter(HashSet<word>, localNamesSet, iter)
        {
            FatalErrorIn
            (
                "ReadFields<class GeoField>"
                "(const IOobjectList&, PtrList<GeoField>&"
                ", const bool)"
            )   << "Fields not synchronised across processors." << endl
                << "Master has fields " << masterNames
                << "  processor " << Pstream::myProcNo()
                << " has fields " << localNames << exit(FatalError);
        }
    }


    fields.setSize(masterNames.size());

    // Make sure to read in masterNames order.

    forAll(masterNames, i)
    {
        Info<< "Reading " << GeoField::typeName << ' ' << masterNames[i]
            << endl;

        const IOobject& io = *fieldObjects[masterNames[i]];

        fields.set
        (
            i,
            new GeoField
            (
                IOobject
                (
                    io.name(),
                    io.instance(),
                    io.local(),
                    io.db(),
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE,
                    io.registerObject()
                )
            )
        );
    }
    return masterNames;
}


void calc
(
    const argList& args,
    const Time& runTime,
    const fvMesh& mesh,
    functionObjectList& fol
)
{
    if (args.optionFound("noFlow"))
    {
        Info<< "    Operating in no-flow mode; no models will be loaded."
            << " All vol, surface and point fields will be loaded." << endl;

        // Read objects in time directory
        IOobjectList objects(mesh, runTime.timeName());

        // Read vol fields.

        PtrList<volScalarField> vsFlds;
        ReadFields(mesh, objects, vsFlds);

        PtrList<volVectorField> vvFlds;
        ReadFields(mesh, objects, vvFlds);

        PtrList<volSphericalTensorField> vstFlds;
        ReadFields(mesh, objects, vstFlds);

        PtrList<volSymmTensorField> vsymtFlds;
        ReadFields(mesh, objects, vsymtFlds);

        PtrList<volTensorField> vtFlds;
        ReadFields(mesh, objects, vtFlds);

        // Read vol-internal fields.

        PtrList<volScalarField::DimensionedInternalField> vsiFlds;
        ReadFields(mesh, objects, vsiFlds);

        PtrList<volVectorField::DimensionedInternalField> vviFlds;
        ReadFields(mesh, objects, vviFlds);

        PtrList<volSphericalTensorField::DimensionedInternalField> vstiFlds;
        ReadFields(mesh, objects, vstiFlds);

        PtrList<volSymmTensorField::DimensionedInternalField> vsymtiFlds;
        ReadFields(mesh, objects, vsymtiFlds);

        PtrList<volTensorField::DimensionedInternalField> vtiFlds;
        ReadFields(mesh, objects, vtiFlds);

        // Read surface fields.

        PtrList<surfaceScalarField> ssFlds;
        ReadFields(mesh, objects, ssFlds);

        PtrList<surfaceVectorField> svFlds;
        ReadFields(mesh, objects, svFlds);

        PtrList<surfaceSphericalTensorField> sstFlds;
        ReadFields(mesh, objects, sstFlds);

        PtrList<surfaceSymmTensorField> ssymtFlds;
        ReadFields(mesh, objects, ssymtFlds);

        PtrList<surfaceTensorField> stFlds;
        ReadFields(mesh, objects, stFlds);

        // Read point fields.
        const pointMesh& pMesh = pointMesh::New(mesh);

        PtrList<pointScalarField> psFlds;
        ReadFields(pMesh, objects, psFlds);

        PtrList<pointVectorField> pvFlds;
        ReadFields(pMesh, objects, pvFlds);

        PtrList<pointSphericalTensorField> pstFlds;
        ReadFields(pMesh, objects, pstFlds);

        PtrList<pointSymmTensorField> psymtFlds;
        ReadFields(pMesh, objects, psymtFlds);

        PtrList<pointTensorField> ptFlds;
        ReadFields(pMesh, objects, ptFlds);

        // Read uniform dimensioned fields
        IOobjectList constantObjects(mesh, runTime.constant());

        PtrList<uniformDimensionedScalarField> usFlds;
        ReadUniformFields(constantObjects, usFlds, true);

        PtrList<uniformDimensionedVectorField> uvFlds;
        ReadUniformFields(constantObjects, uvFlds, true);

        PtrList<uniformDimensionedSphericalTensorField> ustFlds;
        ReadUniformFields(constantObjects, ustFlds, true);

        PtrList<uniformDimensionedSymmTensorField> usymmtFlds;
        ReadUniformFields(constantObjects, usymmtFlds, true);

        PtrList<uniformDimensionedTensorField> utFlds;
        ReadUniformFields(constantObjects, utFlds, true);

        fol.execute(true);
    }
    else
    {
        Info<< "    Reading phi" << endl;
        surfaceScalarField phi
        (
            IOobject
            (
                "phi",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ
            ),
            mesh
        );

        Info<< "    Reading U" << endl;
        volVectorField U
        (
            IOobject
            (
                "U",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ
            ),
            mesh
        );

        Info<< "    Reading p" << endl;
        volScalarField p
        (
            IOobject
            (
                "p",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ
            ),
            mesh
        );

        #include "createFvOptions.H"

        if (phi.dimensions() == dimVolume/dimTime)
        {
            IOobject RASPropertiesHeader
            (
                "RASProperties",
                runTime.constant(),
                mesh,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            );

            IOobject LESPropertiesHeader
            (
                "LESProperties",
                runTime.constant(),
                mesh,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            );

            if (RASPropertiesHeader.headerOk())
            {
                IOdictionary RASProperties(RASPropertiesHeader);

                singlePhaseTransportModel laminarTransport(U, phi);

                autoPtr<incompressible::RASModel> RASModel
                (
                    incompressible::RASModel::New
                    (
                        U,
                        phi,
                        laminarTransport
                    )
                );

                fol.execute(true);
            }
            else if (LESPropertiesHeader.headerOk())
            {
                IOdictionary LESProperties(LESPropertiesHeader);

                singlePhaseTransportModel laminarTransport(U, phi);

                autoPtr<incompressible::LESModel> sgsModel
                (
                    incompressible::LESModel::New(U, phi, laminarTransport)
                );

                fol.execute(true);
            }
            else
            {
                IOdictionary transportProperties
                (
                    IOobject
                    (
                        "transportProperties",
                        runTime.constant(),
                        mesh,
                        IOobject::MUST_READ_IF_MODIFIED,
                        IOobject::NO_WRITE
                    )
                );

                fol.execute(true);
            }
        }
        else if (phi.dimensions() == dimMass/dimTime)
        {
            autoPtr<fluidThermo> thermo(fluidThermo::New(mesh));

            volScalarField rho
            (
                IOobject
                (
                    "rho",
                    runTime.timeName(),
                    mesh
                ),
                thermo->rho()
            );

            IOobject RASPropertiesHeader
            (
                "RASProperties",
                runTime.constant(),
                mesh,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            );

            IOobject LESPropertiesHeader
            (
                "LESProperties",
                runTime.constant(),
                mesh,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            );

            if (RASPropertiesHeader.headerOk())
            {
                IOdictionary RASProperties(RASPropertiesHeader);

                autoPtr<compressible::RASModel> RASModel
                (
                    compressible::RASModel::New
                    (
                        rho,
                        U,
                        phi,
                        thermo()
                    )
                );

                fol.execute(true);
            }
            else if (LESPropertiesHeader.headerOk())
            {
                IOdictionary LESProperties(LESPropertiesHeader);

                autoPtr<compressible::LESModel> sgsModel
                (
                    compressible::LESModel::New(rho, U, phi, thermo())
                );

                fol.execute(true);
            }
            else
            {
                IOdictionary transportProperties
                (
                    IOobject
                    (
                        "transportProperties",
                        runTime.constant(),
                        mesh,
                        IOobject::MUST_READ_IF_MODIFIED,
                        IOobject::NO_WRITE
                    )
                );

                fol.execute(true);
            }
        }
        else
        {
            FatalErrorIn(args.executable())
                << "Incorrect dimensions of phi: " << phi.dimensions()
                << nl << exit(FatalError);
        }
    }
}


int main(int argc, char *argv[])
{
    Foam::timeSelector::addOptions();
    #include "addRegionOption.H"
    Foam::argList::addBoolOption
    (
        "noFlow",
        "suppress creating flow models"
    );
    #include "addDictOption.H"

    #include "setRootCase.H"
    #include "createTime.H"
    Foam::instantList timeDirs = Foam::timeSelector::select0(runTime, args);
    #include "createNamedMesh.H"

    // Construct functionObjectList

    autoPtr<functionObjectList> folPtr;
    // Externally stored dictionary for if fol constructed not from runTime
    dictionary folDict;

    if (args.optionFound("dict"))
    {
        folDict = IOdictionary
        (
            IOobject
            (
                args["dict"],
                runTime,
                IOobject::MUST_READ_IF_MODIFIED
            )
        );
        folPtr.reset(new functionObjectList(runTime, folDict));
    }
    else
    {
        folPtr.reset(new functionObjectList(runTime));
    }
    folPtr->start();


    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        Info<< "Time = " << runTime.timeName() << endl;

        mesh.readUpdate();

        FatalIOError.throwExceptions();

        try
        {
            calc(args, runTime, mesh, folPtr());
        }
        catch (IOerror& err)
        {
            Warning<< err << endl;
        }

        Info<< endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
