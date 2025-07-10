/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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

void Foam::functionObjects::sampledSets::updateMasterSets()
{
    masterSets_.setSize(size());
    masterSetOrders_.setSize(size());

    forAll(*this, seti)
    {
        const sampledSet& s = operator[](seti);

        Tuple2<coordSet, labelList> g = s.coords().gather();

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
    if (!dict.found("sets")) return true;

    dict.lookup("fields") >> fields_;

    dict.lookup("interpolationScheme") >> interpolationScheme_;

    formatter_ = setWriter::New(dict.lookup("setFormat"), dict);

    if (dict.isDict("sets"))
    {
        const dictionary& setsDict = dict.subDict("sets");

        setSize(setsDict.size());

        label i = 0;

        forAllConstIter(dictionary, setsDict, iter)
        {
            set
            (
                i++,
                sampledSet::New
                (
                    iter().keyword(),
                    mesh_,
                    iter().dict()
                )
            );
        }
    }
    else
    {
        PtrList<sampledSet> newList
        (
            dict.lookup("sets"),
            sampledSet::iNew(mesh_)
        );
        transfer(newList);
    }

    if (this->size())
    {
        Info<< "Reading set description:" << nl;
        forAll(*this, seti)
        {
            Info<< "    " << operator[](seti).name() << nl;
        }
        Info<< endl;
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
        updateMasterSets();

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
        forAll(*this, seti)
        {
            operator[](seti).movePoints();
        }

        masterSets_.clear();
        masterSetOrders_.clear();
    }
}


void Foam::functionObjects::sampledSets::topoChange
(
    const polyTopoChangeMap& map
)
{
    if (&map.mesh() == &mesh_)
    {
        forAll(*this, seti)
        {
            operator[](seti).topoChange(map);
        }

        masterSets_.clear();
        masterSetOrders_.clear();
    }
}


void Foam::functionObjects::sampledSets::mapMesh(const polyMeshMap& map)
{
    if (&map.mesh() == &mesh_)
    {
        forAll(*this, seti)
        {
            operator[](seti).mapMesh(map);
        }

        masterSets_.clear();
        masterSetOrders_.clear();
    }
}


void Foam::functionObjects::sampledSets::distribute
(
    const polyDistributionMap& map
)
{
    if (&map.mesh() == &mesh_)
    {
        forAll(*this, seti)
        {
            operator[](seti).distribute(map);
        }

        masterSets_.clear();
        masterSetOrders_.clear();
    }
}


// ************************************************************************* //
