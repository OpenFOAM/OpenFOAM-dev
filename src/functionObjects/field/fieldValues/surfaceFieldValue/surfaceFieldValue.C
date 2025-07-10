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

#include "surfaceFieldValue.H"
#include "processorFvPatch.H"
#include "processorCyclicFvPatch.H"
#include "sampledSurface.H"
#include "generatedFaceZone.H"
#include "mergePoints.H"
#include "indirectPrimitivePatch.H"
#include "PatchTools.H"
#include "fvMeshStitcher.H"
#include "polyTopoChangeMap.H"
#include "polyMeshMap.H"
#include "polyDistributionMap.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace fieldValues
{
    defineTypeNameAndDebug(surfaceFieldValue, 0);
    addToRunTimeSelectionTable(fieldValue, surfaceFieldValue, dictionary);
    addToRunTimeSelectionTable(functionObject, surfaceFieldValue, dictionary);
}
}
}

const Foam::NamedEnum
<
    Foam::functionObjects::fieldValues::surfaceFieldValue::selectionTypes,
    4
> Foam::functionObjects::fieldValues::surfaceFieldValue::selectionTypeNames
{
    "faceZone",
    "patch",
    "patches",
    "sampledSurface"
};

const Foam::NamedEnum
<
    Foam::functionObjects::fieldValues::surfaceFieldValue::operationType,
    15
> Foam::functionObjects::fieldValues::surfaceFieldValue::operationTypeNames_
{
    "none",
    "sum",
    "sumMag",
    "orientedSum",
    "average",
    "areaAverage",
    "areaIntegrate",
    "min",
    "max",
    "minMag",
    "maxMag",
    "CoV",
    "UI",
    "areaNormalAverage",
    "areaNormalIntegrate"
};


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::labelList Foam::functionObjects::fieldValues::surfaceFieldValue::patchis
(
    const wordReList& patchNames
) const
{
    labelList patchis;

    forAll(patchNames, i)
    {
        const labelList patchiis =
            mesh_.boundaryMesh().findIndices(patchNames[i]);

        if (patchiis.empty())
        {
            FatalErrorInFunction
                << type() << ' ' << this->name() << ": "
                << selectionTypeNames[selectionType_]
                << "(" << patchNames[i] << "):" << nl
                << "    Unknown patch name: " << patchNames[i]
                << ". Valid patch names are: "
                << mesh_.boundaryMesh().names() << nl
                << exit(FatalError);
        }

        patchis.append(patchiis);
    }

    return patchis;
}


void Foam::functionObjects::fieldValues::surfaceFieldValue::setFaceZoneFaces()
{
    const faceZone& zone = faceZonePtr_->zone();

    // Ensure addressing is built on all processes
    mesh_.polyBFacePatches();
    mesh_.polyBFacePatchFaces();

    DynamicList<label> faceIds(zone.size());
    DynamicList<label> facePatchIds(zone.size());
    DynamicList<label> faceSigns(zone.size());

    forAll(zone, zoneFacei)
    {
        const label facei = zone[zoneFacei];
        const label faceSign = zone.flipMap()[zoneFacei] ? -1 : 1;

        if (mesh_.isInternalFace(facei))
        {
            faceIds.append(facei);
            facePatchIds.append(-1);
            faceSigns.append(faceSign);
        }
        else
        {
            const label bFacei = facei - mesh_.nInternalFaces();

            const labelUList patches = mesh_.polyBFacePatches()[bFacei];
            const labelUList patchFaces = mesh_.polyBFacePatchFaces()[bFacei];

            forAll(patches, i)
            {
                // Don't include processor patch faces twice
                const fvPatch& fvp = mesh_.boundary()[patches[i]];
                if
                (
                    isType<processorFvPatch>(fvp)
                 && refCast<const processorFvPatch>(fvp).neighbour()
                )
                {
                    continue;
                }

                faceIds.append(patchFaces[i]);
                facePatchIds.append(patches[i]);
                faceSigns.append(faceSign);
            }
        }
    }

    faceId_.transfer(faceIds);
    facePatchId_.transfer(facePatchIds);
    faceSign_.transfer(faceSigns);

    nFaces_ = returnReduce(faceId_.size(), sumOp<label>());
}


void Foam::functionObjects::fieldValues::surfaceFieldValue::setPatchesFaces()
{
    const labelList patchis = this->patchis(patchNames_);

    faceId_.clear();
    facePatchId_.clear();
    faceSign_.clear();

    forAll(patchis, i)
    {
        const label patchi = patchis[i];
        const fvPatch& fvp = mesh_.boundary()[patchi];

        faceId_.append(identityMap(fvp.size()));
        facePatchId_.append(labelList(fvp.size(), patchi));
        faceSign_.append(labelList(fvp.size(), 1));

        // If we have selected a cyclic, also include any associated processor
        // cyclic faces
        forAll(mesh_.boundary(), pcPatchj)
        {
            const fvPatch& pcFvp = mesh_.boundary()[pcPatchj];

            if
            (
                isA<processorCyclicFvPatch>(pcFvp)
             && refCast<const processorCyclicFvPatch>(pcFvp).referPatchIndex()
             == patchi
            )
            {
                faceId_.append(identityMap(pcFvp.size()));
                facePatchId_.append(labelList(pcFvp.size(), pcPatchj));
                faceSign_.append(labelList(pcFvp.size(), 1));
            }
        }
    }

    nFaces_ = returnReduce(faceId_.size(), sumOp<label>());
}


void
Foam::functionObjects::fieldValues::surfaceFieldValue::setSampledSurfaceFaces()
{
    surfacePtr_().update();

    nFaces_ = returnReduce(surfacePtr_().faces().size(), sumOp<label>());
}


void Foam::functionObjects::fieldValues::surfaceFieldValue::combineMeshGeometry
(
    faceList& faces,
    pointField& points
) const
{
    List<faceList> allFaces(Pstream::nProcs());
    List<pointField> allPoints(Pstream::nProcs());

    labelList globalFacesIs(faceId_);
    forAll(globalFacesIs, i)
    {
        if (facePatchId_[i] != -1)
        {
            label patchi = facePatchId_[i];
            globalFacesIs[i] += mesh_.boundaryMesh()[patchi].start();
        }
    }

    // Add local faces and points to the all* lists
    indirectPrimitivePatch pp
    (
        IndirectList<face>(mesh_.faces(), globalFacesIs),
        mesh_.points()
    );
    allFaces[Pstream::myProcNo()] = pp.localFaces();
    allPoints[Pstream::myProcNo()] = pp.localPoints();

    Pstream::gatherList(allFaces);
    Pstream::gatherList(allPoints);

    // Renumber and flatten
    label nFaces = 0;
    label nPoints = 0;
    forAll(allFaces, proci)
    {
        nFaces += allFaces[proci].size();
        nPoints += allPoints[proci].size();
    }

    faces.setSize(nFaces);
    points.setSize(nPoints);

    nFaces = 0;
    nPoints = 0;

    // My own data first
    {
        const faceList& fcs = allFaces[Pstream::myProcNo()];
        forAll(fcs, i)
        {
            const face& f = fcs[i];
            face& newF = faces[nFaces++];
            newF.setSize(f.size());
            forAll(f, fp)
            {
                newF[fp] = f[fp] + nPoints;
            }
        }

        const pointField& pts = allPoints[Pstream::myProcNo()];
        forAll(pts, i)
        {
            points[nPoints++] = pts[i];
        }
    }

    // Other proc data follows
    forAll(allFaces, proci)
    {
        if (proci != Pstream::myProcNo())
        {
            const faceList& fcs = allFaces[proci];
            forAll(fcs, i)
            {
                const face& f = fcs[i];
                face& newF = faces[nFaces++];
                newF.setSize(f.size());
                forAll(f, fp)
                {
                    newF[fp] = f[fp] + nPoints;
                }
            }

            const pointField& pts = allPoints[proci];
            forAll(pts, i)
            {
                points[nPoints++] = pts[i];
            }
        }
    }

    // Merge
    labelList oldToNew;
    pointField newPoints;
    const bool hasMerged = mergePoints
    (
        points,
        small,
        false,
        oldToNew,
        newPoints
    );

    if (hasMerged)
    {
        points.transfer(newPoints);
        forAll(faces, i)
        {
            inplaceRenumber(oldToNew, faces[i]);
        }
    }
}


void Foam::functionObjects::fieldValues::surfaceFieldValue::
combineSurfaceGeometry
(
    faceList& faces,
    pointField& points
) const
{
    if (selectionType_ == selectionTypes::sampledSurface)
    {
        const sampledSurface& s = surfacePtr_();

        if (Pstream::parRun())
        {
            // Dimension as fraction of mesh bounding box
            scalar mergeDim = 1e-10*mesh_.bounds().mag();

            labelList pointsMap;

            PatchTools::gatherAndMerge
            (
                mergeDim,
                primitivePatch
                (
                    SubList<face>(s.faces(), s.faces().size()),
                    s.points()
                ),
                points,
                faces,
                pointsMap
            );
        }
        else
        {
            faces = s.faces();
            points = s.points();
        }
    }
}


Foam::scalar
Foam::functionObjects::fieldValues::surfaceFieldValue::area() const
{
    if (selectionType_ == selectionTypes::sampledSurface)
    {
        return gSum(surfacePtr_().magSf());
    }
    else
    {
        return gSum(filterField(mesh_.magSf()));
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::fieldValues::surfaceFieldValue::writeFileHeader
(
    const label i
)
{
    if (operation_ != operationType::none)
    {
        writeCommented(file(), "Selection");
        file()
            << setw(1) << ':' << setw(1) << ' '
            << selectionTypeNames[selectionType_] << "("
            << selectionName_.c_str() << ")" << endl;

        writeHeaderValue(file(), "Faces", nFaces_);

        writeHeaderValue(file(), "Area", area_);

        writeCommented(file(), "Time");

        if (writeNFaces_) file() << tab << "Faces";
        if (writeArea_) file() << tab << "Area";

        forAll(fields_, fieldi)
        {
            file()
                << tab << operationTypeNames_[operation_]
                << "(" << fields_[fieldi] << ")";
        }

        file() << endl;
    }
}


bool Foam::functionObjects::fieldValues::surfaceFieldValue::processValues
(
    const Field<scalar>& values,
    const scalarField& signs,
    const scalarField& weights,
    const vectorField& Sf,
    scalar& result
) const
{
    switch (operation_)
    {
        default:
        {
            // Fall through to same-type operations
            return processValuesTypeType
            (
                values,
                signs,
                weights,
                Sf,
                result
            );
        }
    }
}


bool Foam::functionObjects::fieldValues::surfaceFieldValue::processValues
(
    const Field<vector>& values,
    const scalarField& signs,
    const scalarField& weights,
    const vectorField& Sf,
    scalar& result
) const
{
    switch (operation_)
    {
        case operationType::areaNormalAverage:
        {
            result = gSum(weights*values & Sf)/gSum(mag(weights*Sf));
            return true;
        }
        case operationType::areaNormalIntegrate:
        {
            result = gSum(weights*values & Sf);
            return true;
        }
        default:
        {
            // No fall through. Different types.
            return false;
        }
    }
}


void Foam::functionObjects::fieldValues::surfaceFieldValue::moveMesh()
{
    switch (selectionType_)
    {
        case selectionTypes::faceZone:
            if (faceZonePtr_->zone().moveUpdate()) setFaceZoneFaces();
            break;
        case selectionTypes::patch:
        case selectionTypes::patches:
            break;
        case selectionTypes::sampledSurface:
            surfacePtr_->expire();
            setSampledSurfaceFaces();
            break;
    }

    area_ = nFaces_ ? area() : NaN;
}


void Foam::functionObjects::fieldValues::surfaceFieldValue::changeMesh()
{
    switch (selectionType_)
    {
        case selectionTypes::faceZone:
            if (!faceZonePtr_->zone().moveUpdate()) setFaceZoneFaces();
            break;
        case selectionTypes::patch:
        case selectionTypes::patches:
            setPatchesFaces();
            break;
        case selectionTypes::sampledSurface:
            break;
    }

    moveMesh();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fieldValues::surfaceFieldValue::surfaceFieldValue
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldValue(name, runTime, dict, typeName),
    surfaceWriterPtr_(nullptr),
    selectionType_(selectionTypeNames.select(dict)),
    selectionName_(string::null),
    operation_(operationTypeNames_.read(dict.lookup("operation"))),
    weightFieldNames_(),
    nFaces_(0),
    area_(NaN),
    writeNFaces_(dict.lookupOrDefault("writeNumberOfFaces", false)),
    writeArea_(dict.lookupOrDefault("writeArea", false)),
    faceZonePtr_(nullptr),
    patchNames_(),
    faceId_(),
    facePatchId_(),
    faceSign_(),
    surfacePtr_(nullptr)
{
    read(dict);
}

Foam::functionObjects::fieldValues::surfaceFieldValue::surfaceFieldValue
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    fieldValue(name, obr, dict, typeName),
    surfaceWriterPtr_(nullptr),
    selectionType_(selectionTypeNames.select(dict)),
    selectionName_(string::null),
    operation_(operationTypeNames_.read(dict.lookup("operation"))),
    weightFieldNames_(),
    nFaces_(0),
    area_(NaN),
    writeNFaces_(dict.lookupOrDefault("writeNumberOfFaces", false)),
    writeArea_(dict.lookupOrDefault("writeArea", false)),
    faceZonePtr_(nullptr),
    patchNames_(),
    faceId_(),
    facePatchId_(),
    faceSign_(),
    surfacePtr_(nullptr)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::fieldValues::surfaceFieldValue::~surfaceFieldValue()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::fieldValues::surfaceFieldValue::read
(
    const dictionary& dict
)
{
    fieldValue::read(dict);

    const word selection(selectionTypeNames[selectionType_]);

    switch (selectionType_)
    {
        case selectionTypes::faceZone:
        {
            faceZonePtr_.reset(new generatedFaceZone(mesh_, dict));

            break;
        }
        case selectionTypes::patch:
        {
            const wordRe patchName = dict.lookup<wordRe>(selection);

            patchNames_ = wordReList(1, patchName);

            selectionName_ =
                patchName.isPattern()
              ? '"' + patchName + '"'
              : string(patchName);

            break;
        }
        case selectionTypes::patches:
        {
            patchNames_ = dict.lookup<wordReList>(selection);

            const labelList patchis = this->patchis(patchNames_);

            selectionName_.clear();

            forAll(patchis, i)
            {
                const fvPatch& fvp = mesh_.boundary()[patchis[i]];
                selectionName_.append((i ? " " : "") + fvp.name());
            }

            break;
        }
        case selectionTypes::sampledSurface:
        {
            surfacePtr_ =
                sampledSurface::New
                (
                    name(),
                    mesh_,
                    dict.subDict
                    (
                        selectionTypeNames[selectionTypes::sampledSurface]
                    )
                );

            selectionName_ = surfacePtr_().name();

            break;
        }
    }

    if (dict.found("weightFields"))
    {
        dict.lookup("weightFields") >> weightFieldNames_;
    }
    else if (dict.found("weightField"))
    {
        weightFieldNames_.setSize(1);

        dict.lookup("weightField") >> weightFieldNames_[0];
    }

    if (writeFields_)
    {
        const word surfaceFormat(dict.lookup("surfaceFormat"));

        surfaceWriterPtr_.reset
        (
            surfaceWriter::New(surfaceFormat, dict).ptr()
        );
    }

    changeMesh();

    // Report configuration
    Info<< type() << ' ' << name() << " read:" << nl;
    Info<< "    number of faces = " << nFaces_ << nl;
    if (nFaces_)
    {
        Info<< "    area = " << area_ << nl;
    }
    Info<< "    operation = " << operationTypeNames_[operation_] << nl;
    if (weightFieldNames_.size() == 1)
    {
        Info<< "    weight field = " << weightFieldNames_[0] << nl;
    }
    if (weightFieldNames_.size() > 1)
    {
        Info<< "    weight fields =";
        forAll(weightFieldNames_, i) Info<< ' ' << weightFieldNames_[i];
        Info<< nl;
    }
    Info<< endl;

    return true;
}


bool Foam::functionObjects::fieldValues::surfaceFieldValue::write()
{
    if (nFaces_ == 0)
    {
        return false;
    }

    // Look to see if any fields exist. Use the flag to suppress output later.
    bool anyFields = false;
    forAll(fields_, i)
    {
        #define validFieldType(fieldType, none)                          \
            anyFields = anyFields || validField<fieldType>(fields_[i]);
        FOR_ALL_FIELD_TYPES(validFieldType);
        #undef validFieldType
    }
    if (!anyFields && fields_.size() > 1) // (error for 1 will happen below)
    {
        cannotFindObjects(fields_);
    }

    // Initialise the file, write the header, etc...
    if (anyFields && operation_ != operationType::none)
    {
        fieldValue::write();
    }

    // Update the surface
    if (selectionType_ == selectionTypes::sampledSurface)
    {
        surfacePtr_().update();
    }

    // Write the time
    if (anyFields && operation_ != operationType::none && Pstream::master())
    {
        writeTime(file());
    }

    // Write the number of faces and/or the area if necessary
    if (anyFields && operation_ != operationType::none && Pstream::master())
    {
        if (writeNFaces_)
        {
            file() << tab << nFaces_;
        }
        if (writeArea_)
        {
            file() << tab << area_;
        }
    }
    if (writeNFaces_)
    {
        Log << "    number of faces = " << nFaces_ << endl;
    }
    if (writeArea_)
    {
        Log << "    area = " << area_ << endl;
    }

    // Construct the sign and weight fields and the surface normals
    const scalarField signs
    (
        selectionType_ == selectionTypes::sampledSurface
      ? scalarField(surfacePtr_().Sf().size(), 1)
      : List<scalar>(faceSign_)
    );
    scalarField weights(signs.size(), 1);
    forAll(weightFieldNames_, i)
    {
        weights *= getFieldValues<scalar>(weightFieldNames_[i]);
    }
    const vectorField Sf
    (
        selectionType_ == selectionTypes::sampledSurface
      ? surfacePtr_().Sf()
      : (signs*filterField(mesh_.Sf()))()
    );

    // Create storage for the values
    #define DeclareValues(fieldType, nullArg) \
        PtrList<Field<fieldType>> fieldType##Values(fields_.size());
    FOR_ALL_FIELD_TYPES(DeclareValues);
    #undef DeclareValues

    // Process the fields in turn. Get the values and store if necessary.
    // Compute the operation and write into the file and the log.
    forAll(fields_, i)
    {
        const word& fieldName = fields_[i];
        bool ok = false;

        #define writeValuesFieldType(fieldType, none)                          \
        {                                                                      \
            const bool typeOk = validField<fieldType>(fieldName);              \
                                                                               \
            if (typeOk)                                                        \
            {                                                                  \
                tmp<Field<fieldType>> values =                                 \
                    getFieldValues<fieldType>(fieldName);                      \
                                                                               \
                writeValues<fieldType>                                         \
                (                                                              \
                    fieldName,                                                 \
                    values(),                                                  \
                    signs,                                                     \
                    weights,                                                   \
                    Sf                                                         \
                );                                                             \
                                                                               \
                if (writeFields_)                                              \
                {                                                              \
                    fieldType##Values.set                                      \
                    (                                                          \
                        i,                                                     \
                        getFieldValues<fieldType>(fieldName).ptr()             \
                    );                                                         \
                }                                                              \
            }                                                                  \
                                                                               \
            ok = ok || typeOk;                                                 \
        }
        FOR_ALL_FIELD_TYPES(writeValuesFieldType);
        #undef writeValuesFieldType

        if (!ok)
        {
            cannotFindObject(fieldName);
        }
    }

    // Finalise the file and the log
    if (anyFields && operation_ != operationType::none && Pstream::master())
    {
        file() << endl;
    }
    Log << endl;

    // Write a surface file with the values if specified
    if (writeFields_)
    {
        faceList faces;
        pointField points;

        if (selectionType_ == selectionTypes::sampledSurface)
        {
            combineSurfaceGeometry(faces, points);
        }
        else
        {
            combineMeshGeometry(faces, points);
        }

        if (Pstream::master())
        {
            surfaceWriterPtr_->write
            (
                baseFileDir()/name()/time_.name(),
                word(selectionTypeNames[selectionType_])
              + "("
              + (
                    selectionType_ == selectionTypes::patches
                  ? selectionName_.replace(" ", ",").c_str()
                  : selectionName_.c_str()
                )
              + ")",
                points,
                faces,
                fields_,
                false
                #define ValuesParameter(fieldType, nullArg) \
                    , fieldType##Values
                FOR_ALL_FIELD_TYPES(ValuesParameter)
                #undef ValuesParameter
            );
        }
    }

    return true;
}


void Foam::functionObjects::fieldValues::surfaceFieldValue::movePoints
(
    const polyMesh& mesh
)
{
    if (&mesh == &this->mesh())
    {
        fieldValue::movePoints(mesh);
        if (selectionType_ == selectionTypes::faceZone)
        {
            faceZonePtr_->movePoints();
        }
        moveMesh();
    }
}


void Foam::functionObjects::fieldValues::surfaceFieldValue::topoChange
(
    const polyTopoChangeMap& map
)
{
    if (&map.mesh() == &mesh())
    {
        fieldValue::topoChange(map);
        if (selectionType_ == selectionTypes::faceZone)
        {
            faceZonePtr_->topoChange(map);
        }
        changeMesh();
    }
}


void Foam::functionObjects::fieldValues::surfaceFieldValue::mapMesh
(
    const polyMeshMap& map
)
{
    if (&map.mesh() == &mesh())
    {
        fieldValue::mapMesh(map);
        if (selectionType_ == selectionTypes::faceZone)
        {
            faceZonePtr_->mapMesh(map);
        }
        changeMesh();
    }
}


void Foam::functionObjects::fieldValues::surfaceFieldValue::distribute
(
    const polyDistributionMap& map
)
{
    if (&map.mesh() == &mesh())
    {
        fieldValue::distribute(map);
        if (selectionType_ == selectionTypes::faceZone)
        {
            faceZonePtr_->distribute(map);
        }
        changeMesh();
    }
}


// ************************************************************************* //
