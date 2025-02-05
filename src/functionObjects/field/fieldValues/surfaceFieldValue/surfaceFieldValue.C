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

template<>
const char* Foam::NamedEnum
<
    Foam::functionObjects::fieldValues::surfaceFieldValue::selectionTypes,
    4
>::names[] =
{
    "faceZone",
    "patch",
    "patches",
    "sampledSurface"
};

template<>
const char* Foam::NamedEnum
<
    Foam::functionObjects::fieldValues::surfaceFieldValue::operationType,
    14
>::names[] =
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
    "areaNormalAverage",
    "areaNormalIntegrate"
};

const Foam::NamedEnum
<
    Foam::functionObjects::fieldValues::surfaceFieldValue::selectionTypes,
    4
> Foam::functionObjects::fieldValues::surfaceFieldValue::selectionTypeNames;

const Foam::NamedEnum
<
    Foam::functionObjects::fieldValues::surfaceFieldValue::operationType,
    14
> Foam::functionObjects::fieldValues::surfaceFieldValue::operationTypeNames_;


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::fieldValues::surfaceFieldValue::setFaceZoneFaces
(
    const word& zoneName
)
{
    const label zoneId = mesh_.faceZones().findIndex(zoneName);

    if (zoneId < 0)
    {
        FatalErrorInFunction
            << type() << " " << name() << ": "
            << selectionTypeNames[selectionType_]
            << "(" << zoneName << "):" << nl
            << "    Unknown face zone name: " << zoneName
            << ". Valid face zones are: " << mesh_.faceZones().toc()
            << nl << exit(FatalError);
    }

    selectionName_ = zoneName;

    // Ensure addressing is built on all processes
    mesh_.polyBFacePatches();
    mesh_.polyBFacePatchFaces();

    const faceZone& zone = mesh_.faceZones()[zoneId];

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


void Foam::functionObjects::fieldValues::surfaceFieldValue::setPatchesFaces
(
    const labelList& patchis
)
{
    selectionName_.clear();

    faceId_.clear();
    facePatchId_.clear();
    faceSign_.clear();

    forAll(patchis, i)
    {
        const label patchi = patchis[i];
        const fvPatch& fvp = mesh_.boundary()[patchi];

        selectionName_.append((i ? " " : "") + fvp.name());

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


void Foam::functionObjects::fieldValues::surfaceFieldValue::setPatchesFaces
(
    const wordReList& patchNames
)
{
    labelList patchis;

    forAll(patchNames, i)
    {
        const labelList patchiis =
            mesh_.boundaryMesh().findIndices(patchNames[i]);

        if (patchiis.empty())
        {
            FatalErrorInFunction
                << type() << " " << this->name() << ": "
                << selectionTypeNames[selectionType_]
                << "(" << patchNames[i] << "):" << nl
                << "    Unknown patch name: " << patchNames[i]
                << ". Valid patch names are: "
                << mesh_.boundaryMesh().names() << nl
                << exit(FatalError);
        }

        patchis.append(patchiis);
    }

    setPatchesFaces(patchis);
}


void Foam::functionObjects::fieldValues::surfaceFieldValue::setPatchFaces
(
    const wordRe& patchName
)
{
    setPatchesFaces(wordReList(1, patchName));

    selectionName_ =
        patchName.isPattern() ? '"' + patchName + '"' : string(patchName);
}


void
Foam::functionObjects::fieldValues::surfaceFieldValue::setSampledSurfaceFaces
(
    const dictionary& dict
)
{
    surfacePtr_ = sampledSurface::New(name(), mesh_, dict);

    surfacePtr_().update();

    selectionName_ = surfacePtr_().name();

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
    bool hasMerged = mergePoints
    (
        points,
        small,
        false,
        oldToNew,
        newPoints
    );

    if (hasMerged)
    {
        if (debug)
        {
            Pout<< "Merged from " << points.size()
                << " down to " << newPoints.size() << " points" << endl;
        }

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
Foam::functionObjects::fieldValues::surfaceFieldValue::totalArea() const
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

void Foam::functionObjects::fieldValues::surfaceFieldValue::initialise
(
    const dictionary& dict
)
{
    switch (selectionType_)
    {
        case selectionTypes::faceZone:
        {
            setFaceZoneFaces
            (
                dict.lookupBackwardsCompatible<word>({"faceZone", "name"})
            );
            break;
        }
        case selectionTypes::patch:
        {
            setPatchFaces
            (
                dict.lookupBackwardsCompatible<wordRe>({"patch", "name"})
            );
            break;
        }
        case selectionTypes::patches:
        {
            setPatchesFaces(dict.lookup<wordReList>("patches"));
            break;
        }
        case selectionTypes::sampledSurface:
        {
            setSampledSurfaceFaces(dict.subDict("sampledSurfaceDict"));
            break;
        }
        default:
        {
            FatalErrorInFunction
                << type() << " " << name() << ": "
                << selectionTypeNames[selectionType_]
                << "    Unknown selection type. Valid selection types are:"
                << selectionTypeNames.sortedToc() << nl << exit(FatalError);
        }
    }

    if (nFaces_ == 0 && (!mesh_.stitcher().stitches() || !mesh_.conformal()))
    {
        FatalErrorInFunction
            << type() << " " << name() << ": "
            << selectionTypeNames[selectionType_]
            << "(" << selectionName_.c_str() << "):" << nl
            << " selection has no faces" << exit(FatalError);
    }

    if (selectionType_ == selectionTypes::sampledSurface)
    {
        surfacePtr_().update();
    }

    totalArea_ = totalArea();

    Info<< type() << " " << name() << ":" << nl
        << "    total faces  = " << nFaces_
        << nl
        << "    total area   = " << totalArea_
        << nl;

    if (dict.readIfPresent("weightFields", weightFieldNames_))
    {
        Info<< name() << " " << operationTypeNames_[operation_]
            << " weight fields " << weightFieldNames_;
    }
    else if (dict.found("weightField"))
    {
        weightFieldNames_.setSize(1);
        dict.lookup("weightField") >> weightFieldNames_[0];

        Info<< name() << " " << operationTypeNames_[operation_]
            << " weight field " << weightFieldNames_[0];
    }

    Info<< nl << endl;

    if (writeFields_)
    {
        const word surfaceFormat(dict.lookup("surfaceFormat"));

        surfaceWriterPtr_.reset
        (
            surfaceWriter::New(surfaceFormat, dict).ptr()
        );
    }
}


void Foam::functionObjects::fieldValues::surfaceFieldValue::writeFileHeader
(
    const label i
)
{
    if (operation_ != operationType::none)
    {
        writeCommented(file(), "Selection : ");
        file()
            << selectionTypeNames[selectionType_] << "("
            << selectionName_.c_str() << ")" << endl;
        writeCommented(file(), "Faces  : ");
        file() << nFaces_ << endl;
        writeCommented(file(), "Area   : ");
        file() << totalArea_ << endl;

        writeCommented(file(), "Time");
        if (writeArea_)
        {
            file() << tab << "Area";
        }

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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fieldValues::surfaceFieldValue::surfaceFieldValue
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldValue(name, runTime, dict, typeName),
    dict_(dict),
    surfaceWriterPtr_(nullptr),
    selectionType_
    (
        selectionTypeNames.read
        (
            dict.lookupBackwardsCompatible({"select", "regionType"})
        )
    ),
    selectionName_(string::null),
    operation_(operationTypeNames_.read(dict.lookup("operation"))),
    weightFieldNames_(),
    writeArea_(dict.lookupOrDefault("writeArea", false)),
    nFaces_(0),
    faceId_(),
    facePatchId_(),
    faceSign_()
{
    read(dict_);
}

Foam::functionObjects::fieldValues::surfaceFieldValue::surfaceFieldValue
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    fieldValue(name, obr, dict, typeName),
    dict_(dict),
    surfaceWriterPtr_(nullptr),
    selectionType_
    (
        selectionTypeNames.read
        (
            dict.lookupBackwardsCompatible({"select", "regionType"})
        )
    ),
    selectionName_(string::null),
    operation_(operationTypeNames_.read(dict.lookup("operation"))),
    weightFieldNames_(),
    writeArea_(dict.lookupOrDefault("writeArea", false)),
    nFaces_(0),
    faceId_(),
    facePatchId_(),
    faceSign_()
{
    read(dict_);
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
    if (dict != dict_)
    {
        dict_ = dict;
    }

    fieldValue::read(dict);
    initialise(dict);

    return true;
}


bool Foam::functionObjects::fieldValues::surfaceFieldValue::write()
{
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

    // Write the total area if necessary
    if (writeArea_)
    {
        totalArea_ = totalArea();
        if (anyFields && operation_ != operationType::none && Pstream::master())
        {
            file() << tab << totalArea_;
        }
        Log << "    total area = " << totalArea_ << endl;
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
                  : selectionName_
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
    if (&mesh == &mesh_)
    {
        // It may be necessary to reset if the mesh moves. The total area might
        // change, as might non-conformal faces.
        initialise(dict_);
    }
}


void Foam::functionObjects::fieldValues::surfaceFieldValue::topoChange
(
    const polyTopoChangeMap& map
)
{
    if (&map.mesh() == &mesh_)
    {
        initialise(dict_);
    }
}


void Foam::functionObjects::fieldValues::surfaceFieldValue::mapMesh
(
    const polyMeshMap& map
)
{
    if (&map.mesh() == &mesh_)
    {
        initialise(dict_);
    }
}


void Foam::functionObjects::fieldValues::surfaceFieldValue::distribute
(
    const polyDistributionMap& map
)
{
    if (&map.mesh() == &mesh_)
    {
        initialise(dict_);
    }
}


// ************************************************************************* //
