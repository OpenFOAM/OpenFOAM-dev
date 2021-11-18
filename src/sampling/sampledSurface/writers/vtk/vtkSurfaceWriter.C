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

#include "vtkSurfaceWriter.H"
#include "vtkWritePolyData.H"
#include "OSspecific.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(vtkSurfaceWriter, 0);
    addToRunTimeSelectionTable(surfaceWriter, vtkSurfaceWriter, word);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vtkSurfaceWriter::vtkSurfaceWriter
(
    const IOstream::streamFormat writeFormat
)
:
    surfaceWriter(writeFormat)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::vtkSurfaceWriter::~vtkSurfaceWriter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::vtkSurfaceWriter::write
(
    const fileName& outputDir,
    const fileName& surfaceName,
    const pointField& points,
    const faceList& faces,
    const wordList& fieldNames,
    const bool writePointValues
    #define FieldTypeValuesConstArg(Type, nullArg) \
        , const UPtrList<const Field<Type>>& field##Type##Values
    FOR_ALL_FIELD_TYPES(FieldTypeValuesConstArg)
    #undef FieldTypeValuesConstArg
) const
{
    if (!isDir(outputDir))
    {
        mkDir(outputDir);
    }

    vtkWritePolyData::write
    (
        outputDir/surfaceName + ".vtk",
        "sampleSurface",
        writeFormat_ == IOstream::BINARY,
        points,
        labelList(),
        edgeList(),
        faces,
        fieldNames,
        boolList(fieldNames.size(), writePointValues),
        UPtrList<const Field<label>>(fieldNames.size())
        #define FieldTypeValuesParameter(Type, nullArg) , field##Type##Values
        FOR_ALL_FIELD_TYPES(FieldTypeValuesParameter)
        #undef FieldTypeValuesParameter
    );
}


// ************************************************************************* //
