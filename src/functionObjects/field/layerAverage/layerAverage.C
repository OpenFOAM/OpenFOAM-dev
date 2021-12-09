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

#include "layerAverage.H"
#include "regionSplit.H"
#include "syncTools.H"
#include "volFields.H"
#include "writeFile.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace functionObjects
    {
        defineTypeNameAndDebug(layerAverage, 0);
        addToRunTimeSelectionTable(functionObject, layerAverage, dictionary);
    }

    template<>
    const char*
        Foam::NamedEnum<Foam::vector::components, 3>::names[] =
        {"x", "y", "z"};
}

const Foam::NamedEnum<Foam::vector::components, 3>
    Foam::functionObjects::layerAverage::axisNames_;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::layerAverage::walkOppositeFaces
(
    const labelList& startFaces,
    const boolList& startFaceIntoOwners,
    boolList& blockedFace,
    boolList& cellIsLayer
)
{
    // Initialise the front
    DynamicList<label> frontFaces(startFaces);
    DynamicList<bool> frontFaceIntoOwners(startFaceIntoOwners);
    DynamicList<label> newFrontFaces(frontFaces.size());
    DynamicList<bool> newFrontFaceIntoOwners(frontFaceIntoOwners.size());

    // Set the start and end faces as blocked
    UIndirectList<bool>(blockedFace, startFaces) = true;

    // Iterate until the front is empty
    while (returnReduce(frontFaces.size(), sumOp<label>()) > 0)
    {
        // Transfer front faces across couplings
        boolList bndFaceIsFront(mesh_.nFaces() - mesh_.nInternalFaces(), false);
        forAll(frontFaces, i)
        {
            const label facei = frontFaces[i];
            const label intoOwner = frontFaceIntoOwners[i];

            if (!mesh_.isInternalFace(facei) && !intoOwner)
            {
                bndFaceIsFront[facei - mesh_.nInternalFaces()] = true;
            }
        }
        syncTools::swapBoundaryFaceList(mesh_, bndFaceIsFront);

        // Set new front faces
        forAll(bndFaceIsFront, bndFacei)
        {
            const label facei = mesh_.nInternalFaces() + bndFacei;

            if (bndFaceIsFront[bndFacei] && !blockedFace[facei])
            {
                blockedFace[facei] = true;
                frontFaces.append(facei);
                frontFaceIntoOwners.append(true);
            }
        }

        // Transfer across cells
        newFrontFaces.clear();
        newFrontFaceIntoOwners.clear();

        forAll(frontFaces, i)
        {
            const label facei = frontFaces[i];
            const label intoOwner = frontFaceIntoOwners[i];

            const label celli =
                intoOwner
              ? mesh_.faceOwner()[facei]
              : mesh_.isInternalFace(facei) ? mesh_.faceNeighbour()[facei] : -1;

            if (celli != -1)
            {
                cellIsLayer[celli] = true;

                const label oppositeFacei =
                    mesh_.cells()[celli]
                   .opposingFaceLabel(facei, mesh_.faces());
                const bool oppositeIntoOwner =
                    mesh_.faceOwner()[oppositeFacei] != celli;

                if (oppositeFacei == -1)
                {
                    FatalErrorInFunction
                        << "Cannot find face on cell " << mesh_.cells()[celli]
                        << " opposing face " << mesh_.faces()[facei]
                        << ". Mesh is not layered." << exit(FatalError);
                }
                else if (!blockedFace[oppositeFacei])
                {
                    blockedFace[oppositeFacei] = true;
                    newFrontFaces.append(oppositeFacei);
                    newFrontFaceIntoOwners.append(oppositeIntoOwner);
                }
            }
        }

        frontFaces.transfer(newFrontFaces);
        frontFaceIntoOwners.transfer(newFrontFaceIntoOwners);
    }

    // Determine whether cells on the other sides of couplings are in layers
    boolList bndCellIsLayer(mesh_.nFaces() - mesh_.nInternalFaces(), false);
    forAll(bndCellIsLayer, bndFacei)
    {
        const label facei = mesh_.nInternalFaces() + bndFacei;
        const label owni = mesh_.faceOwner()[facei];

        bndCellIsLayer[bndFacei] = cellIsLayer[owni];
    }
    syncTools::swapBoundaryFaceList(mesh_, bndCellIsLayer);

    // Block faces where one adjacent cell is in a layer and the other is not
    forAll(mesh_.faces(), facei)
    {
        const label owni = mesh_.faceOwner()[facei];
        if (mesh_.isInternalFace(facei))
        {
            const label nbri = mesh_.faceNeighbour()[facei];
            if (cellIsLayer[owni] != cellIsLayer[nbri])
            {
                blockedFace[facei] = true;
            }
        }
        else
        {
            const label bndFacei = facei - mesh_.nInternalFaces();
            if (cellIsLayer[owni] != bndCellIsLayer[bndFacei])
            {
                blockedFace[facei] = true;
            }
        }
    }
}


void Foam::functionObjects::layerAverage::calcLayers()
{
    DynamicList<label> startFaces;
    DynamicList<bool> startFaceIntoOwners;

    forAll(patchIDs_, i)
    {
        const polyPatch& pp = mesh_.boundaryMesh()[patchIDs_[i]];
        startFaces.append(identity(pp.size()) + pp.start());
        startFaceIntoOwners.append(boolList(pp.size(), true));
    }

    forAll(zoneIDs_, i)
    {
        const faceZone& fz = mesh_.faceZones()[zoneIDs_[i]];
        startFaces.append(fz);
        startFaceIntoOwners.append(fz.flipMap());
    }

    // Identify faces that separate the layers
    boolList blockedFace(mesh_.nFaces(), false);
    boolList cellIsLayer(mesh_.nCells(), false);
    walkOppositeFaces
    (
        startFaces,
        startFaceIntoOwners,
        blockedFace,
        cellIsLayer
    );

    // Do analysis for connected layers
    regionSplit rs(mesh_, blockedFace);
    nLayers_ = rs.nRegions();
    cellLayer_.transfer(rs);

    // Get rid of regions that are not layers
    label layeri0 = labelMax, layeri1 = -labelMax;
    forAll(cellLayer_, celli)
    {
        if (cellIsLayer[celli])
        {
            layeri0 = min(layeri0, cellLayer_[celli]);
            layeri1 = max(layeri1, cellLayer_[celli]);
        }
        else
        {
            cellLayer_[celli] = -1;
        }
    }
    reduce(layeri0, minOp<label>());
    reduce(layeri1, maxOp<label>());
    nLayers_ = layeri0 != labelMax ? 1 + layeri1 - layeri0 : 0;
    forAll(cellLayer_, celli)
    {
        if (cellLayer_[celli] != -1)
        {
            cellLayer_[celli] -= layeri0;
        }
    }

    // Report
    if (nLayers_ != 0)
    {
        Info<< "    Detected " << nLayers_ << " layers" << nl << endl;
    }
    else
    {
        WarningInFunction
            << "No layers detected" << endl;
    }

    // Sum number of entries per layer
    layerCount_ = sum(scalarField(mesh_.nCells(), 1));

    // Average the cell centres
    const pointField layerCentres(sum(mesh_.cellCentres())/layerCount_);

    // Sort by the direction component
    x_ = layerCentres.component(axis_);
    sortMap_ = identity(layerCentres.size());
    stableSort(sortMap_, UList<scalar>::less(x_));
    x_.map(x_, sortMap_);

    // If symmetric, keep only half of the coordinates
    if (symmetric_)
    {
        x_.setSize(nLayers_/2);
    }
}


template<>
Foam::vector
Foam::functionObjects::layerAverage::symmetricCoeff<Foam::vector>() const
{
    vector result = vector::one;
    result[axis_] = -1;
    return result;
}


template<>
Foam::symmTensor
Foam::functionObjects::layerAverage::symmetricCoeff<Foam::symmTensor>() const
{
    return sqr(symmetricCoeff<vector>());
}


template<>
Foam::tensor
Foam::functionObjects::layerAverage::symmetricCoeff<Foam::tensor>() const
{
    return symmetricCoeff<symmTensor>();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::layerAverage::layerAverage
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::layerAverage::~layerAverage()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::functionObjects::layerAverage::read(const dictionary& dict)
{
    Info<< type() << " " << name() << ":" << nl;

    patchIDs_ =
        mesh_.boundaryMesh().patchSet
        (
            dict.lookupOrDefault<wordReList>("patches", wordReList())
        ).toc();

    zoneIDs_ =
        findStrings
        (
            dict.lookupOrDefault<wordReList>("zones", wordReList()),
            mesh_.faceZones().names()
        );

    if (patchIDs_.empty() && zoneIDs_.empty())
    {
        WarningInFunction
            << "No patches or zones specified" << endl;
    }

    symmetric_ = dict.lookupOrDefault<bool>("symmetric", false);

    axis_ = axisNames_.read(dict.lookup("axis"));

    fields_ = dict.lookup<wordList>("fields");

    formatter_ = setWriter::New(dict.lookup("setFormat"), dict);

    calcLayers();

    return true;
}


Foam::wordList Foam::functionObjects::layerAverage::fields() const
{
    return fields_;
}


bool Foam::functionObjects::layerAverage::execute()
{
    return true;
}


bool Foam::functionObjects::layerAverage::write()
{
    // Create list of available fields
    wordList fieldNames;
    forAll(fields_, fieldi)
    {
        if
        (
            false
            #define FoundTypeField(Type, nullArg) \
              || foundObject<VolField<Type>>(fields_[fieldi])
            FOR_ALL_FIELD_TYPES(FoundTypeField)
            #undef FoundTypeField
        )
        {
            fieldNames.append(fields_[fieldi]);
        }
        else
        {
            cannotFindObject(fields_[fieldi]);
        }
    }

    // Collapse the fields
    #define DeclareTypeValueSets(Type, nullArg) \
        PtrList<Field<Type>> Type##ValueSets(fieldNames.size());
    FOR_ALL_FIELD_TYPES(DeclareTypeValueSets);
    #undef DeclareTypeValueSets
    forAll(fieldNames, fieldi)
    {
        #define CollapseTypeFields(Type, nullArg)                           \
            if (mesh_.foundObject<VolField<Type>>(fieldNames[fieldi]))      \
            {                                                               \
                const VolField<Type>& field =                               \
                    mesh_.lookupObject<VolField<Type>>(fieldNames[fieldi]); \
                                                                            \
                Type##ValueSets.set                                         \
                (                                                           \
                    fieldi,                                                 \
                    average(field.primitiveField())                         \
                );                                                          \
            }
        FOR_ALL_FIELD_TYPES(CollapseTypeFields);
        #undef CollapseTypeFields
    }

    // Output directory
    fileName outputPath =
        mesh_.time().globalPath()/writeFile::outputPrefix/name();
    if (mesh_.name() != fvMesh::defaultRegion)
    {
        outputPath = outputPath/mesh_.name();
    }

    // Write
    if (Pstream::master() && x_.size())
    {
        formatter_->write
        (
            outputPath/mesh_.time().timeName(),
            typeName,
            coordSet(true, axisNames_[axis_], x_),
            fieldNames
            #define TypeValueSetsParameter(Type, nullArg) , Type##ValueSets
            FOR_ALL_FIELD_TYPES(TypeValueSetsParameter)
            #undef TypeValueSetsParameter
        );
    }

    return true;
}


void Foam::functionObjects::layerAverage::updateMesh(const mapPolyMesh&)
{
    Info<< type() << " " << name() << ":" << nl;
    calcLayers();
}


void Foam::functionObjects::layerAverage::movePoints(const polyMesh&)
{
    Info<< type() << " " << name() << ":" << nl;
    calcLayers();
}


// ************************************************************************* //
