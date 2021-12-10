/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "vtkSetWriter.H"
#include "coordSet.H"
#include "faceList.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "vtkWritePolyData.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(vtkSetWriter, 0);
    addToRunTimeSelectionTable(setWriter, vtkSetWriter, word);
    addToRunTimeSelectionTable(setWriter, vtkSetWriter, dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::vtkSetWriter::~vtkSetWriter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::vtkSetWriter::write
(
    const fileName& outputDir,
    const fileName& setName,
    const coordSet& set,
    const wordList& valueSetNames
    #define TypeValueSetsConstArg(Type, nullArg) \
        , const UPtrList<const Field<Type>>& Type##ValueSets
    FOR_ALL_FIELD_TYPES(TypeValueSetsConstArg)
    #undef TypeValueSetsConstArg
) const
{
    if (!set.hasPointAxis())
    {
        FatalErrorInFunction
            << "Cannot write " << setName << " in " << typeName
            << " format as it does not have a point axis"
            << exit(FatalError);
    }

    if (!isDir(outputDir))
    {
        mkDir(outputDir);
    }

    // Write
    vtkWritePolyData::write
    (
        outputDir/setName + ".vtk",
        setName,
        writeFormat_ == IOstream::BINARY,
        set.pointCoords()(),
        set.vertices(),
        set.lines(),
        faceList(),
        valueSetNames,
        boolList(valueSetNames.size(), true),
        UPtrList<const Field<label>>(valueSetNames.size())
        #define TypeValueSetsParameter(Type, nullArg) , Type##ValueSets
        FOR_ALL_FIELD_TYPES(TypeValueSetsParameter)
        #undef TypeValueSetsParameter
    );
}


// ************************************************************************* //
