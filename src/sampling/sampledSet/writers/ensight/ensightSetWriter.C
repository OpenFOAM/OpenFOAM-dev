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

#include "ensightSetWriter.H"
#include "coordSet.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "IOmanip.H"
#include "ensightPart.H"
#include "ensightPTraits.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ensightSetWriter, 0);
    addToRunTimeSelectionTable(setWriter, ensightSetWriter, word);
    addToRunTimeSelectionTable(setWriter, ensightSetWriter, dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ensightSetWriter::~ensightSetWriter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ensightSetWriter::write
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

    // const scalar timeValue = Foam::name(this->mesh().time().timeValue());
    const scalar timeValue = 0.0;

    OFstream osCase(outputDir/setName + ".case");
    ensightGeoFile osGeom
    (
        outputDir/setName + ".000.mesh",
        writeFormat_
    );

    osCase
        << "FORMAT" << nl
        << "type: ensight gold" << nl
        << nl;

    osCase
        << "GEOMETRY" << nl
        << "model:        1     " << osGeom.name().name() << nl
        << nl;

    osCase
        << "VARIABLE" << nl;
    forAll(valueSetNames, fieldi)
    {
        #define WriteTypeCase(Type, nullArg)                               \
            if (Type##ValueSets.set(fieldi))                               \
            {                                                              \
                osCase                                                     \
                    << ensightPTraits<Type>::typeName << " per node:"      \
                    << setw(10) << 1 << "       " << valueSetNames[fieldi] \
                    << "       " << setName.c_str() << ".***."             \
                    << valueSetNames[fieldi] << nl;                        \
            }
        FOR_ALL_FIELD_TYPES(WriteTypeCase);
        #undef WriteTypeCase
    }
    osCase
        << nl;

    osCase
        << "TIME" << nl
        << "time set:                      1" << nl
        << "number of steps:               1" << nl
        << "filename start number:         0" << nl
        << "filename increment:            1" << nl
        << "time values:" << nl
        << timeValue << nl
        << nl;

    // Write mesh file
    {
        osGeom
            << "part" << nl
            << setw(10) << 1 << nl
            << "internalMesh" << nl
            << "coordinates" << nl
            << setw(10) << set.size() << nl;

        const tmp<pointField> tPoints = set.pointCoords();
        const pointField& points = tPoints();

        for (direction cmpt = 0; cmpt < vector::nComponents; cmpt++)
        {
            forAll(set, pointi)
            {
                if (mag(points[pointi][cmpt]) >= scalar(floatScalarVSmall))
                {
                    osGeom  << setw(12) << points[pointi][cmpt] << nl;
                }
                else
                {
                    osGeom  << setw(12) << scalar(0) << nl;
                }
            }
        }

        labelList vertices(set.vertices());
        osGeom
            << "point" << nl
            << setw(10) << vertices.size() << nl;
        forAll(vertices, vertexi)
        {
            osGeom
                << setw(10) << vertices[vertexi] + 1 << nl;
        }

        labelPairList edges(set.edges());
        osGeom
            << "bar2" << nl
            << setw(10) << edges.size() << nl;
        forAll(edges, edgei)
        {
            osGeom
                << setw(10) << edges[edgei].first() + 1
                << setw(10) << edges[edgei].second() + 1 << nl;
        }
    }

    // Write field files
    forAll(valueSetNames, fieldi)
    {
        ensightFile osField
        (
            outputDir/setName
          + ".000."
          + valueSetNames[fieldi],
            writeFormat_
        );

        #define WriteTypeValues(Type, nullArg)                        \
            if (Type##ValueSets.set(fieldi))                          \
            {                                                         \
                osField                                               \
                    << pTraits<Type>::typeName << nl                  \
                    << "part" << nl                                   \
                    << setw(10) << 1 << nl                            \
                    << "coordinates" << nl;                           \
                for                                                   \
                (                                                     \
                    direction cmpt = 0;                               \
                    cmpt < pTraits<Type>::nComponents;                \
                    cmpt++                                            \
                )                                                     \
                {                                                     \
                    const scalarField fld                             \
                    (                                                 \
                        Type##ValueSets[fieldi].component(cmpt)       \
                    );                                                \
                    forAll(fld, i)                                    \
                    {                                                 \
                        if (mag(fld[i]) >= scalar(floatScalarVSmall)) \
                        {                                             \
                            osField << setw(12) << fld[i] << nl;      \
                        }                                             \
                        else                                          \
                        {                                             \
                            osField << setw(12) << scalar(0) << nl;   \
                        }                                             \
                    }                                                 \
                }                                                     \
            }
        FOR_ALL_FIELD_TYPES(WriteTypeValues);
        #undef WriteTypeValues
    }
}


// ************************************************************************* //
