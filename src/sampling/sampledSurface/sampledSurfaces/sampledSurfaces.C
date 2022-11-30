/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "sampledSurfaces.H"
#include "PatchTools.H"
#include "polyTopoChangeMap.H"
#include "polyMeshMap.H"
#include "polyDistributionMap.H"
#include "OSspecific.H"
#include "writeFile.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(sampledSurfaces, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        sampledSurfaces,
        dictionary
    );
}
}

Foam::scalar Foam::functionObjects::sampledSurfaces::mergeTol_ = 1e-10;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::sampledSurfaces::needsUpdate() const
{
    forAll(*this, si)
    {
        if (operator[](si).needsUpdate())
        {
            return true;
        }
    }

    return false;
}


bool Foam::functionObjects::sampledSurfaces::expire()
{
    bool justExpired = false;

    forAll(*this, si)
    {
        if (operator[](si).expire())
        {
            justExpired = true;
        }

        // Clear merge information
        if (Pstream::parRun())
        {
            mergeList_[si].clear();
        }
    }

    // true if any surfaces just expired
    return justExpired;
}


bool Foam::functionObjects::sampledSurfaces::update()
{
    bool updated = false;

    if (!needsUpdate())
    {
        return updated;
    }

    // Serial: quick and easy, no merging required
    if (!Pstream::parRun())
    {
        forAll(*this, si)
        {
            if (operator[](si).update())
            {
                updated = true;
            }
        }

        return updated;
    }

    // Dimension as fraction of mesh bounding box
    scalar mergeDim = mergeTol_ * mesh_.bounds().mag();

    if (Pstream::master() && debug)
    {
        Pout<< nl << "Merging all points within "
            << mergeDim << " metre" << endl;
    }

    forAll(*this, si)
    {
        sampledSurface& s = operator[](si);

        if (s.update())
        {
            updated = true;
        }
        else
        {
            continue;
        }

        PatchTools::gatherAndMerge
        (
            mergeDim,
            primitivePatch
            (
                SubList<face>(s.faces(), s.faces().size()),
                s.points()
            ),
            mergeList_[si].points,
            mergeList_[si].faces,
            mergeList_[si].pointsMap
        );
    }

    return updated;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::sampledSurfaces::sampledSurfaces
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, t, dict),
    PtrList<sampledSurface>(),
    outputPath_
    (
        mesh_.time().globalPath()
       /writeFile::outputPrefix
       /(mesh_.name() != polyMesh::defaultRegion ? mesh_.name() : word())
       /name
    ),
    fields_(),
    interpolationScheme_(word::null),
    writeEmpty_(false),
    mergeList_(),
    formatter_(nullptr)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::sampledSurfaces::~sampledSurfaces()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::sampledSurfaces::read(const dictionary& dict)
{
    bool surfacesFound = dict.found("surfaces");

    if (surfacesFound)
    {
        dict.lookup("fields") >> fields_;

        dict.lookup("interpolationScheme") >> interpolationScheme_;

        dict.readIfPresent("writeEmpty", writeEmpty_);

        const word writeType(dict.lookup("surfaceFormat"));

        // Define the surface formatter
        formatter_ = surfaceWriter::New(writeType, dict);

        PtrList<sampledSurface> newList
        (
            dict.lookup("surfaces"),
            sampledSurface::iNew(mesh_)
        );
        transfer(newList);

        if (Pstream::parRun())
        {
            mergeList_.setSize(size());
        }

        // Ensure all surfaces and merge information are expired
        expire();

        if (this->size())
        {
            Info<< "Reading surface description:" << nl;
            forAll(*this, si)
            {
                Info<< "    " << operator[](si).name() << nl;
            }
            Info<< endl;
        }
    }

    if (Pstream::master() && debug)
    {
        Pout<< "sample fields:" << fields_ << nl
            << "sample surfaces:" << nl << "(" << nl;

        forAll(*this, si)
        {
            Pout<< "  " << operator[](si) << endl;
        }
        Pout<< ")" << endl;
    }

    return true;
}


Foam::wordList Foam::functionObjects::sampledSurfaces::fields() const
{
    wordList fields(fields_);

    forAll(*this, si)
    {
        fields.append(operator[](si).fields());
    }

    return fields;
}


bool Foam::functionObjects::sampledSurfaces::execute()
{
    return true;
}


bool Foam::functionObjects::sampledSurfaces::write()
{
    if (size())
    {
        // Finalise surfaces, merge points etc.
        update();

        // Create the output directory
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
              || foundObject<VolField<Type>>(fields_[fieldi]) \
              || foundObject<SurfaceField<Type>>(fields_[fieldi])
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

        // Sample and write the surfaces
        forAll(*this, surfi)
        {
            const sampledSurface& s = operator[](surfi);

            #define GenerateFieldTypeValues(Type, nullArg) \
                PtrList<Field<Type>> field##Type##Values = \
                    sampleType<Type>(surfi, fieldNames, interpolation##Type##s);
            FOR_ALL_FIELD_TYPES(GenerateFieldTypeValues);
            #undef GenerateFieldTypeValues

            if (Pstream::parRun())
            {
                if
                (
                    Pstream::master()
                 && (mergeList_[surfi].faces.size() || writeEmpty_)
                )
                {
                    formatter_->write
                    (
                        outputPath_/mesh_.time().name(),
                        s.name(),
                        mergeList_[surfi].points,
                        mergeList_[surfi].faces,
                        fieldNames,
                        s.interpolate()
                        #define FieldTypeValuesParameter(Type, nullArg) \
                            , field##Type##Values
                        FOR_ALL_FIELD_TYPES(FieldTypeValuesParameter)
                        #undef FieldTypeValuesParameter
                    );
                }
            }
            else
            {
                if (s.faces().size() || writeEmpty_)
                {
                    formatter_->write
                    (
                        outputPath_/mesh_.time().name(),
                        s.name(),
                        s.points(),
                        s.faces(),
                        fieldNames,
                        s.interpolate()
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


void Foam::functionObjects::sampledSurfaces::movePoints(const polyMesh& mesh)
{
    if (&mesh == &mesh_)
    {
        expire();
    }
}


void Foam::functionObjects::sampledSurfaces::topoChange
(
    const polyTopoChangeMap& map
)
{
    if (&map.mesh() == &mesh_)
    {
        expire();
    }
}


void Foam::functionObjects::sampledSurfaces::mapMesh
(
    const polyMeshMap& map
)
{
    if (&map.mesh() == &mesh_)
    {
        expire();
    }
}


void Foam::functionObjects::sampledSurfaces::distribute
(
    const polyDistributionMap& map
)
{
    if (&map.mesh() == &mesh_)
    {
        expire();
    }
}


// ************************************************************************* //
