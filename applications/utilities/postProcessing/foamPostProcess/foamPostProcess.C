/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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
    foamPostProcess

Description
    Execute the set of functionObjects specified in the selected dictionary
    (which defaults to system/controlDict) or on the command-line for the
    selected set of times on the selected set of fields.

    The functionObjects are either executed directly or for the solver
    optionally specified as a command-line argument.

Usage
    \b foamPostProcess [OPTION]
      - \par -dict <file>
        Read control dictionary from specified location

      - \par -solver <name>
        Solver name

      - \par -libs '(\"lib1.so\" ... \"libN.so\")'
        Specify the additional libraries loaded

      -\par -region <name>
        Specify the region

      - \par -func <name>
        Specify the name of the functionObject to execute, e.g. Q

      - \par -funcs <list>
        Specify the names of the functionObjects to execute, e.g. '(Q div(U))'

      - \par -field <name>
        Specify the name of the field to be processed, e.g. U

      - \par -fields <list>
        Specify a list of fields to be processed,
        e.g. '(U T p)' - regular expressions not currently supported

      - \par -time <ranges>
        comma-separated time ranges - eg, ':10,20,40:70,1000:'

      - \par -latestTime
        Select the latest time

      - \par -list
        List the available configured functionObjects

    Example usage:
      - Print the list of available configured functionObjects:
        \verbatim
            foamPostProcess -list
        \endverbatim

      - Execute the functionObjects specified in the controlDict of the
        fluid region for all the available times:
        \verbatim
            foamPostProcess -region fluid
        \endverbatim

      - Execute the functionObjects specified in the controlDict
        for the 'fluid' solver in the 'cooling' region for the latest time only:
        \verbatim
            foamPostProcess -solver fluid -region cooling -latestTime
        \endverbatim

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "solver.H"
#include "ReadFields.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "uniformDimensionedFields.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define ReadFields(GeoFieldType)                                               \
    readFields<GeoFieldType>(mesh, objects, requiredFields, storedObjects);

#define ReadPointFields(GeoFieldType)                                          \
    readFields<GeoFieldType>(pMesh, objects, requiredFields, storedObjects);

#define ReadUniformFields(FieldType)                                           \
    readUniformFields<FieldType>                                               \
    (constantObjects, requiredFields, storedObjects);

void executeFunctionObjects
(
    const argList& args,
    const Time& runTime,
    fvMesh& mesh,
    const HashSet<word>& requiredFields0,
    functionObjectList& functions,
    bool lastTime
)
{
    Info<< nl << "Reading fields:" << endl;

    // Maintain a stack of the stored objects to clear after executing
    // the functionObjects
    LIFOStack<regIOobject*> storedObjects;

    // Read objects in time directory
    IOobjectList objects(mesh, runTime.name());

    HashSet<word> requiredFields(requiredFields0);
    forAll(functions, i)
    {
        requiredFields.insert(functions[i].fields());
    }

    // Read volFields
    ReadFields(volScalarField);
    ReadFields(volVectorField);
    ReadFields(volSphericalTensorField);
    ReadFields(volSymmTensorField);
    ReadFields(volTensorField);

    // Read internal fields
    ReadFields(volScalarField::Internal);
    ReadFields(volVectorField::Internal);
    ReadFields(volSphericalTensorField::Internal);
    ReadFields(volSymmTensorField::Internal);
    ReadFields(volTensorField::Internal);

    // Read surface fields
    ReadFields(surfaceScalarField);
    ReadFields(surfaceVectorField);
    ReadFields(surfaceSphericalTensorField);
    ReadFields(surfaceSymmTensorField);
    ReadFields(surfaceTensorField);

    // Read point fields.
    const pointMesh& pMesh = pointMesh::New(mesh);

    ReadPointFields(pointScalarField)
    ReadPointFields(pointVectorField);
    ReadPointFields(pointSphericalTensorField);
    ReadPointFields(pointSymmTensorField);
    ReadPointFields(pointTensorField);

    // Read uniform dimensioned fields
    IOobjectList constantObjects(mesh, runTime.constant());

    ReadUniformFields(uniformDimensionedScalarField);
    ReadUniformFields(uniformDimensionedVectorField);
    ReadUniformFields(uniformDimensionedSphericalTensorField);
    ReadUniformFields(uniformDimensionedSymmTensorField);
    ReadUniformFields(uniformDimensionedTensorField);

    Info<< nl << "Executing functionObjects" << endl;

    // Execute the functionObjects in post-processing mode
    functions.execute();

    // Execute the functionObject 'end()' function for the last time
    if (lastTime)
    {
        functions.end();
    }

    while (!storedObjects.empty())
    {
        storedObjects.pop()->checkOut();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addOption
    (
        "solver",
        "name",
        "Solver name"
    );

    timeSelector::addOptions();
    #include "addRegionOption.H"
    #include "addFunctionObjectOptions.H"

    // Set functionObject post-processing mode
    functionObject::postProcess = true;

    #include "setRootCase.H"

    if (args.optionFound("list"))
    {
        Info<< nl
            << "Available configured functionObjects:"
            << listAllConfigFiles
               (
                   functionEntries::includeFuncEntry::functionObjectDictPath
               )
            << endl;
        return 0;
    }

    #include "createTime.H"

    const instantList timeDirs = timeSelector::select0(runTime, args);

    word regionName = fvMesh::defaultRegion;

    if (args.optionReadIfPresent("region", regionName))
    {
        Info
            << "Create mesh " << regionName << " for time = "
            << runTime.name() << nl << endl;
    }
    else
    {
        Info
            << "Create mesh for time = "
            << runTime.name() << nl << endl;
    }

    fvMesh mesh
    (
        IOobject
        (
            regionName,
            runTime.name(),
            runTime,
            IOobject::MUST_READ
        )
    );

    // Either the solver name is specified...
    word solverName;

    // ...or the fields are specified on the command-line
    // or later inferred from the function arguments
    HashSet<word> requiredFields;

    if (args.optionReadIfPresent("solver", solverName))
    {
        libs.open("lib" + solverName + ".so");
    }
    else
    {
        // Initialise the set of selected fields from the command-line options
        if (args.optionFound("fields"))
        {
            args.optionLookup("fields")() >> requiredFields;
        }
        if (args.optionFound("field"))
        {
            requiredFields.insert(args.optionLookup("field")());
        }
    }

    // Externally stored dictionary for functionObjectList
    // if not constructed from runTime
    dictionary functionsControlDict("controlDict");

    // Construct functionObjectList
    autoPtr<functionObjectList> functionsPtr
    (
        functionObjectList::New
        (
            args,
            runTime,
            functionsControlDict
        )
    );

    forAll(timeDirs, timei)
    {
        runTime.setTime(timeDirs[timei], timei);

        Info<< "Time = " << runTime.userTimeName() << endl;

        if (mesh.readUpdate() != polyMesh::UNCHANGED)
        {
            // Update functionObjectList if mesh changes
            functionsPtr = functionObjectList::New
            (
                args,
                runTime,
                functionsControlDict
            );
        }

        FatalIOError.throwExceptions();

        try
        {
            if (solverName != word::null)
            {
                // Optionally instantiate the selected solver
                autoPtr<solver> solverPtr;

                solverPtr = solver::New(solverName, mesh);

                functionsPtr->execute();

                // Clear the objects owned by the mesh
                mesh.objectRegistry::clear();
            }
            else
            {
                executeFunctionObjects
                (
                    args,
                    runTime,
                    mesh,
                    requiredFields,
                    functionsPtr(),
                    timei == timeDirs.size()-1
                );
            }
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
