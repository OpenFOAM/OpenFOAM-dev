/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2022 OpenFOAM Foundation
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

#include "writeVTK.H"
#include "Time.H"
#include "vtkMesh.H"
#include "internalWriter.H"
#include "OSspecific.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(writeVTK, 0);
    addToRunTimeSelectionTable(functionObject, writeVTK, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::writeVTK::writeVTK
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    objectNames_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::writeVTK::~writeVTK()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::writeVTK::read(const dictionary& dict)
{
    dict.lookup("objects") >> objectNames_;

    return true;
}


Foam::wordList Foam::functionObjects::writeVTK::fields() const
{
    return objectNames_;
}


bool Foam::functionObjects::writeVTK::execute()
{
    return true;
}


bool Foam::functionObjects::writeVTK::write()
{
    Info<< type() << " " << name() << " output:" << nl;

    Info<< "Time: " << time_.name() << endl;

    word timeDesc = time_.name();

    // VTK/ directory in the case
    fileName fvPath(time_.path()/"VTK");

    mkDir(fvPath);

    string vtkName = time_.caseName();

    if (Pstream::parRun())
    {
        // Strip off leading casename, leaving just processor_DDD ending.
        string::size_type i = vtkName.rfind("processor");

        if (i != string::npos)
        {
            vtkName = vtkName.substr(i);
        }
    }

    // Create file and write header
    fileName vtkFileName
    (
        fvPath/vtkName
      + "_"
      + timeDesc
      + ".vtk"
    );

    Info<< "    Internal  : " << vtkFileName << endl;

    vtkMesh vMesh(const_cast<fvMesh&>(mesh_));

    // Write mesh
    internalWriter writer(vMesh, false, vtkFileName);

    // Declare UPtrLists to the volFields that are to be written
    #define DeclareTypeFields(Type, nullArg) \
        UPtrList<const VolField<Type>> Type##Fields;
    FOR_ALL_FIELD_TYPES(DeclareTypeFields);
    #undef DeclareTypeFields

    // Look up the volFields and store the pointers
    label nFields = 0;
    forAll(objectNames_, i)
    {
        bool objectFound = false;

        #define SetTypeFields(Type, nullarg)                            \
        {                                                               \
            if (obr_.foundObject<VolField<Type>>(objectNames_[i]))      \
            {                                                           \
                const VolField<Type>& field =                           \
                    obr_.lookupObject<VolField<Type>>(objectNames_[i]); \
                                                                        \
                Type##Fields.resize(Type##Fields.size() + 1);           \
                Type##Fields.set(Type##Fields.size() - 1, &field);      \
                                                                        \
                Info<< "    Writing " << VolField<Type>::typeName       \
                    << " field " << field.name() << endl;               \
                                                                        \
                nFields ++;                                             \
                objectFound = true;                                     \
            }                                                           \
        }
        FOR_ALL_FIELD_TYPES(SetTypeFields);
        #undef SetTypeFields

        if (!objectFound)
        {
            cannotFindObject(objectNames_[i]);
        }
    }

    // Write header for cellID and volFields
    vtkWriteOps::writeCellDataHeader
    (
        writer.os(),
        vMesh.nFieldCells(),
        1 + nFields
    );

    // Write cellID field
    writer.writeCellIDs();

    // Write volFields
    #define WriteTypeFields(Type, nullArg) \
        writer.write(Type##Fields);
    FOR_ALL_FIELD_TYPES(WriteTypeFields);
    #undef WriteTypeFields

    return true;
}


// ************************************************************************* //
