/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "sampledSets.H"
#include "dictionary.H"
#include "Time.H"
#include "volFields.H"
#include "ListListOps.H"
#include "SortableList.H"
#include "volPointInterpolation.H"
#include "polyTopoChangeMap.H"
#include "polyMeshMap.H"
#include "polyDistributionMap.H"
#include "writeFile.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(sampledSets, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        sampledSets,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::sampledSets::combineSampledSets()
{
    masterSets_.setSize(size());
    masterSetOrders_.setSize(size());

    forAll(*this, seti)
    {
        const sampledSet& s = operator[](seti);

        Tuple2<coordSet, labelList> g = s.gather();

        masterSets_.set(seti, new coordSet(g.first()));
        masterSetOrders_[seti] = g.second();

        if (Pstream::master() && masterSets_[seti].size() == 0)
        {
            WarningInFunction
                << "Sample set " << s.name()
                << " has zero points." << endl;
        }
    }
}


void Foam::functionObjects::sampledSets::correct()
{
    bool setsFound = dict_.found("sets");
    if (setsFound)
    {
        searchEngine_.correct();

        PtrList<sampledSet> newList
        (
            dict_.lookup("sets"),
            sampledSet::iNew(mesh_, searchEngine_)
        );
        transfer(newList);
        combineSampledSets();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::sampledSets::sampledSets
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, t, dict),
    PtrList<sampledSet>(),
    outputPath_
    (
        mesh_.time().globalPath()
       /writeFile::outputPrefix
       /(mesh_.name() != polyMesh::defaultRegion ? mesh_.name() : word())
       /name
    ),
    searchEngine_(mesh_),
    interpolationScheme_(word::null),
    formatter_(nullptr)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::sampledSets::~sampledSets()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::sampledSets::read(const dictionary& dict)
{
    dict_ = dict;

    bool setsFound = dict_.found("sets");
    if (setsFound)
    {
        dict_.lookup("fields") >> fields_;

        dict.lookup("interpolationScheme") >> interpolationScheme_;

        const word writeType(dict.lookup("setFormat"));

        // Define the set formatter
        formatter_ = setWriter::New(writeType, dict);

        PtrList<sampledSet> newList
        (
            dict_.lookup("sets"),
            sampledSet::iNew(mesh_, searchEngine_)
        );
        transfer(newList);
        combineSampledSets();

        if (this->size())
        {
            Info<< "Reading set description:" << nl;
            forAll(*this, seti)
            {
                Info<< "    " << operator[](seti).name() << nl;
            }
            Info<< endl;
        }
    }

    return true;
}


Foam::wordList Foam::functionObjects::sampledSets::fields() const
{
    return fields_;
}


bool Foam::functionObjects::sampledSets::execute()
{
    return true;
}


bool Foam::functionObjects::sampledSets::write()
{
    if (size())
    {
        if (Pstream::master())
        {
            if (debug)
            {
                Pout<< "Creating directory "
                    << outputPath_/mesh_.time().name() << nl << endl;
            }

            mkDir(outputPath_/mesh_.time().name());
        }

        // Create a list of names of fields that are actually available
        wordList fieldNames;
        forAll(fields_, fieldi)
        {
            #define FoundFieldType(Type, nullArg)             \
              || foundObject<VolField<Type>>(fields_[fieldi])
            if (false FOR_ALL_FIELD_TYPES(FoundFieldType))
            {
                fieldNames.append(fields_[fieldi]);
            }
            else
            {
                cannotFindObject(fields_[fieldi]);
            }
            #undef FoundFieldType
        }

        // Create table of cached interpolations, to prevent unnecessary work
        // when interpolating fields over multiple surfaces
        #define DeclareInterpolations(Type, nullArg) \
            HashPtrTable<interpolation<Type>> interpolation##Type##s;
        FOR_ALL_FIELD_TYPES(DeclareInterpolations);
        #undef DeclareInterpolations

        // Sample and write the sets
        forAll(*this, seti)
        {
            #define GenerateFieldTypeValues(Type, nullArg) \
                PtrList<Field<Type>> field##Type##Values = \
                    sampleType<Type>(seti, fieldNames, interpolation##Type##s);
            FOR_ALL_FIELD_TYPES(GenerateFieldTypeValues);
            #undef GenerateFieldTypeValues

            if (Pstream::parRun())
            {
                if (Pstream::master() && masterSets_[seti].size())
                {
                    formatter_->write
                    (
                        outputPath_/mesh_.time().name(),
                        operator[](seti).name(),
                        masterSets_[seti],
                        fieldNames
                        #define FieldTypeValuesParameter(Type, nullArg) \
                            , field##Type##Values
                        FOR_ALL_FIELD_TYPES(FieldTypeValuesParameter)
                        #undef FieldTypeValuesParameter
                    );

                }
            }
            else
            {
                if (operator[](seti).size())
                {
                    formatter_->write
                    (
                        outputPath_/mesh_.time().name(),
                        operator[](seti).name(),
                        operator[](seti),
                        fieldNames
                        #define FieldTypeValuesParameter(Type, nullArg) \
                            , field##Type##Values
                        FOR_ALL_FIELD_TYPES(FieldTypeValuesParameter)
                        #undef FieldTypeValuesParameter
                    );
                }
            }
        }
    }

    return true;
}


void Foam::functionObjects::sampledSets::movePoints(const polyMesh& mesh)
{
    if (&mesh == &mesh_)
    {
        correct();
    }
}



void Foam::functionObjects::sampledSets::topoChange
(
    const polyTopoChangeMap& map
)
{
    if (&map.mesh() == &mesh_)
    {
        correct();
    }
}


void Foam::functionObjects::sampledSets::mapMesh(const polyMeshMap& map)
{
    if (&map.mesh() == &mesh_)
    {
        correct();
    }
}


void Foam::functionObjects::sampledSets::distribute
(
    const polyDistributionMap& map
)
{
    if (&map.mesh() == &mesh_)
    {
        correct();
    }
}


// ************************************************************************* //
