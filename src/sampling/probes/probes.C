/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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

#include "probes.H"
#include "volFields.H"
#include "OSspecific.H"
#include "writeFile.H"
#include "meshSearch.H"
#include "polyTopoChangeMap.H"
#include "polyMeshMap.H"
#include "polyDistributionMap.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(probes, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        probes,
        dictionary
    );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::probes::findElements(const fvMesh& mesh)
{
    if (debug)
    {
        Info<< "probes: resetting sample locations" << endl;
    }

    const meshSearch& searchEngine = meshSearch::New(mesh_);

    cellList_.clear();
    cellList_.setSize(locations_.size());

    faceList_.clear();
    faceList_.setSize(locations_.size());

    forAll(locations_, probei)
    {
        const vector& location = locations_[probei];

        const label celli = searchEngine.findCell(location);

        cellList_[probei] = celli;

        if (celli != -1)
        {
            const labelList& cellFaces = mesh.cells()[celli];
            const vector& cellCentre = mesh.cellCentres()[celli];
            scalar minDistance = great;
            label minFaceID = -1;
            forAll(cellFaces, i)
            {
                label facei = cellFaces[i];
                vector dist = mesh.faceCentres()[facei] - cellCentre;
                if (mag(dist) < minDistance)
                {
                    minDistance = mag(dist);
                    minFaceID = facei;
                }
            }
            faceList_[probei] = minFaceID;
        }
        else
        {
            faceList_[probei] = -1;
        }

        if (debug && (cellList_[probei] != -1 || faceList_[probei] != -1))
        {
            Pout<< "probes : found point " << location
                << " in cell " << cellList_[probei]
                << " and face " << faceList_[probei] << endl;
        }
    }


    // Check if all probes have been found.
    forAll(cellList_, probei)
    {
        const vector& location = locations_[probei];
        label celli = cellList_[probei];
        label facei = faceList_[probei];

        // Check at least one processor with cell.
        reduce(celli, maxOp<label>());
        reduce(facei, maxOp<label>());

        if (celli == -1)
        {
            if (Pstream::master())
            {
                WarningInFunction
                    << "Did not find location " << location
                    << " in any cell. Skipping location." << endl;
            }
        }
        else if (facei == -1)
        {
            if (Pstream::master())
            {
                WarningInFunction
                    << "Did not find location " << location
                    << " in any face. Skipping location." << endl;
            }
        }
        else
        {
            // Make sure location not on two domains.
            if (cellList_[probei] != -1 && cellList_[probei] != celli)
            {
                WarningInFunction
                    << "Location " << location
                    << " seems to be on multiple domains:"
                    << " cell " << cellList_[probei]
                    << " on my domain " << Pstream::myProcNo()
                    << " and cell " << celli << " on some other domain." << endl
                    << "This might happen if the probe location is on"
                    << " a processor patch. Change the location slightly"
                    << " to prevent this." << endl;
            }

            if (faceList_[probei] != -1 && faceList_[probei] != facei)
            {
                WarningInFunction
                    << "Location " << location
                    << " seems to be on multiple domains:"
                    << " cell " << faceList_[probei]
                    << " on my domain " << Pstream::myProcNo()
                    << " and face " << facei << " on some other domain." << endl
                    << "This might happen if the probe location is on"
                    << " a processor patch. Change the location slightly"
                    << " to prevent this." << endl;
            }
        }
    }
}


Foam::label Foam::probes::prepare()
{
    const label nFields = classifyFields();

    // adjust file streams
    if (Pstream::master())
    {
        wordHashSet currentFields;

        currentFields.insert(scalarFields_);
        currentFields.insert(vectorFields_);
        currentFields.insert(sphericalTensorFields_);
        currentFields.insert(symmTensorFields_);
        currentFields.insert(tensorFields_);

        currentFields.insert(surfaceScalarFields_);
        currentFields.insert(surfaceVectorFields_);
        currentFields.insert(surfaceSphericalTensorFields_);
        currentFields.insert(surfaceSymmTensorFields_);
        currentFields.insert(surfaceTensorFields_);

        if (debug)
        {
            Info<< "Probing fields: " << currentFields << nl
                << "Probing locations: " << locations_ << nl
                << endl;
        }

        const fileName probeDir =
            mesh_.time().globalPath()
           /functionObjects::writeFile::outputPrefix
           /(mesh_.name() != polyMesh::defaultRegion ? mesh_.name() : word())
           /name()
           /mesh_.time().name();

        // ignore known fields, close streams for fields that no longer exist
        forAllIter(HashPtrTable<OFstream>, probeFilePtrs_, iter)
        {
            if (!currentFields.erase(iter.key()))
            {
                if (debug)
                {
                    Info<< "close probe stream: " << iter()->name() << endl;
                }

                delete probeFilePtrs_.remove(iter);
            }
        }

        // currentFields now just has the new fields - open streams for them
        forAllConstIter(wordHashSet, currentFields, iter)
        {
            const word& fieldName = iter.key();

            // Create directory if does not exist.
            mkDir(probeDir);

            OFstream* fPtr = new OFstream(probeDir/fieldName);
            OFstream& os = *fPtr;

            if (debug)
            {
                Info<< "open probe stream: " << os.name() << endl;
            }

            probeFilePtrs_.insert(fieldName, fPtr);

            const unsigned int w = IOstream::defaultPrecision() + 7;
            os << setf(ios_base::left);

            forAll(locations_, probei)
            {
                os<< "# Probe " << probei << ' ' << locations_[probei] << endl;
            }

            os  << setw(w) << "# Time";

            forAll(locations_, probei)
            {
                os<< ' ' << setw(w) << probei;
            }
            os<< endl;
        }
    }

    return nFields;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::probes::probes
(
    const word& name,
    const Time& t,
    const dictionary& dict,
    const bool initialise
)
:
    functionObject(name, t),
    mesh_
    (
        refCast<const fvMesh>
        (
            t.lookupObject<objectRegistry>
            (
                dict.lookupOrDefault("region", polyMesh::defaultRegion)
            )
        )
    ),
    fields_(),
    fixedLocations_(true),
    interpolationScheme_("cell")
{
    read(dict, initialise);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::probes::~probes()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::probes::read(const dictionary& dict, const bool initialise)
{
    dict.lookup("probeLocations") >> locations_;

    dict.lookup("fields") >> fields_;

    dict.readIfPresent("fixedLocations", fixedLocations_);
    if (dict.readIfPresent("interpolationScheme", interpolationScheme_))
    {
        if (!fixedLocations_ && interpolationScheme_ != "cell")
        {
            WarningInFunction
                << "Only cell interpolation can be applied when "
                << "not using fixedLocations.  InterpolationScheme "
                << "entry will be ignored";
        }
    }

    if (initialise)
    {
        // Initialise cells to sample from supplied locations
        findElements(mesh_);

        prepare();
    }

    return true;
}


bool Foam::probes::read(const dictionary& dict)
{
    return read(dict, true);
}


Foam::wordList Foam::probes::fields() const
{
    return fields_;
}


bool Foam::probes::execute()
{
    return true;
}


bool Foam::probes::write()
{
    if (locations_.size() && prepare())
    {
        sampleAndWrite(scalarFields_);
        sampleAndWrite(vectorFields_);
        sampleAndWrite(sphericalTensorFields_);
        sampleAndWrite(symmTensorFields_);
        sampleAndWrite(tensorFields_);

        sampleAndWriteSurfaceFields(surfaceScalarFields_);
        sampleAndWriteSurfaceFields(surfaceVectorFields_);
        sampleAndWriteSurfaceFields(surfaceSphericalTensorFields_);
        sampleAndWriteSurfaceFields(surfaceSymmTensorFields_);
        sampleAndWriteSurfaceFields(surfaceTensorFields_);
    }

    return true;
}


void Foam::probes::movePoints(const polyMesh& mesh)
{
    DebugInfo<< "probes: movePoints" << endl;

    if (&mesh != &mesh_) return;

    if (fixedLocations_)
    {
        findElements(mesh_);
    }
}


void Foam::probes::topoChange(const polyTopoChangeMap& map)
{
    DebugInfo<< "probes: topoChange" << endl;

    if (&map.mesh() != &mesh_) return;

    if (fixedLocations_)
    {
        findElements(mesh_);
    }
    else
    {
        if (debug)
        {
            Info<< "probes: remapping sample locations" << endl;
        }

        // 1. Update cells
        if (!map.reverseCellMap().empty())
        {
            DynamicList<label> elems(cellList_.size());

            const labelList& reverseMap = map.reverseCellMap();
            forAll(cellList_, i)
            {
                const label celli = cellList_[i];
                const label newCelli = reverseMap[celli];

                if (newCelli == -1)
                {
                    // cell removed
                }
                else if (newCelli < -1)
                {
                    // cell merged
                    elems.append(-newCelli - 2);
                }
                else
                {
                    // valid new cell
                    elems.append(newCelli);
                }
            }

            cellList_.transfer(elems);
        }

        // 2. Update faces
        if (!map.reverseFaceMap().empty())
        {
            DynamicList<label> elems(faceList_.size());

            const labelList& reverseMap = map.reverseFaceMap();
            forAll(faceList_, i)
            {
                const label facei = faceList_[i];
                const label newFacei = reverseMap[facei];

                if (newFacei == -1)
                {
                    // face removed
                }
                else if (newFacei < -1)
                {
                    // face merged
                    elems.append(-newFacei - 2);
                }
                else
                {
                    // valid new face
                    elems.append(newFacei);
                }
            }

            faceList_.transfer(elems);
        }
    }
}


void Foam::probes::mapMesh(const polyMeshMap& map)
{
    DebugInfo<< "probes: mapMesh" << endl;

    if (&map.mesh() != &mesh_) return;

    findElements(mesh_);
}


void Foam::probes::distribute(const polyDistributionMap& map)
{
    DebugInfo<< "probes: distribute" << endl;

    if (&map.mesh() != &mesh_) return;

    // Distribute the list of cells and faces. There might be a cheaper way of
    // doing this than distributing a full cell/face list. But this way is easy
    // and readable and uses the polyDistributionMap at a high level without
    // needing to think about the actual communication or addressing. And
    // run-distribution shouldn't be happening too often. So it's fine.
    auto distribute = []
    (
        const label nOldElements,
        const distributionMap& elementMap,
        labelList& probeElements
    )
    {
        labelList elementProbes(nOldElements, -1);
        forAll(probeElements, probei)
        {
            if (probeElements[probei] != -1)
            {
                elementProbes[probeElements[probei]] = probei;
            }
        }

        elementMap.distribute(elementProbes);

        probeElements = -1;
        forAll(elementProbes, elementi)
        {
            if (elementProbes[elementi] != -1)
            {
                probeElements[elementProbes[elementi]] = elementi;
            }
        }
    };

    distribute(map.nOldCells(), map.cellMap(), cellList_);
    distribute(map.nOldFaces(), map.faceMap(), faceList_);
}


// ************************************************************************* //
