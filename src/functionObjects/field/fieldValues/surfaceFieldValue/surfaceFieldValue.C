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

#include "surfaceFieldValue.H"
#include "emptyPolyPatch.H"
#include "coupledPolyPatch.H"
#include "sampledSurface.H"
#include "mergePoints.H"
#include "indirectPrimitivePatch.H"
#include "PatchTools.H"
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
    Foam::functionObjects::fieldValues::surfaceFieldValue::regionTypes,
    3
>::names[] =
{
    "faceZone",
    "patch",
    "sampledSurface"
};

template<>
const char* Foam::NamedEnum
<
    Foam::functionObjects::fieldValues::surfaceFieldValue::operationType,
    16
>::names[] =
{
    "none",
    "sum",
    "sumMag",
    "sumDirection",
    "sumDirectionBalance",
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
    Foam::functionObjects::fieldValues::surfaceFieldValue::regionTypes,
    3
> Foam::functionObjects::fieldValues::surfaceFieldValue::regionTypeNames_;

const Foam::NamedEnum
<
    Foam::functionObjects::fieldValues::surfaceFieldValue::operationType,
    16
> Foam::functionObjects::fieldValues::surfaceFieldValue::operationTypeNames_;


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::fieldValues::surfaceFieldValue::setFaceZoneFaces()
{
    label zoneId = mesh_.faceZones().findZoneID(regionName_);

    if (zoneId < 0)
    {
        FatalErrorInFunction
            << type() << " " << name() << ": "
            << regionTypeNames_[regionType_] << "(" << regionName_ << "):" << nl
            << "    Unknown face zone name: " << regionName_
            << ". Valid face zones are: " << mesh_.faceZones().names()
            << nl << exit(FatalError);
    }

    const faceZone& fZone = mesh_.faceZones()[zoneId];

    DynamicList<label> faceIds(fZone.size());
    DynamicList<label> facePatchIds(fZone.size());
    DynamicList<label> faceSigns(fZone.size());

    forAll(fZone, i)
    {
        label facei = fZone[i];

        label faceId = -1;
        label facePatchId = -1;
        if (mesh_.isInternalFace(facei))
        {
            faceId = facei;
            facePatchId = -1;
        }
        else
        {
            facePatchId = mesh_.boundaryMesh().whichPatch(facei);
            const polyPatch& pp = mesh_.boundaryMesh()[facePatchId];
            if (isA<coupledPolyPatch>(pp))
            {
                if (refCast<const coupledPolyPatch>(pp).owner())
                {
                    faceId = pp.whichFace(facei);
                }
                else
                {
                    faceId = -1;
                }
            }
            else if (!isA<emptyPolyPatch>(pp))
            {
                faceId = facei - pp.start();
            }
            else
            {
                faceId = -1;
                facePatchId = -1;
            }
        }

        if (faceId >= 0)
        {
            if (fZone.flipMap()[i])
            {
                faceSigns.append(-1);
            }
            else
            {
                faceSigns.append(1);
            }
            faceIds.append(faceId);
            facePatchIds.append(facePatchId);
        }
    }

    faceId_.transfer(faceIds);
    facePatchId_.transfer(facePatchIds);
    faceSign_.transfer(faceSigns);
    nFaces_ = returnReduce(faceId_.size(), sumOp<label>());

    if (debug)
    {
        Pout<< "Original face zone size = " << fZone.size()
            << ", new size = " << faceId_.size() << endl;
    }
}


void Foam::functionObjects::fieldValues::surfaceFieldValue::setPatchFaces()
{
    const label patchid = mesh_.boundaryMesh().findPatchID(regionName_);

    if (patchid < 0)
    {
        FatalErrorInFunction
            << type() << " " << name() << ": "
            << regionTypeNames_[regionType_] << "(" << regionName_ << "):" << nl
            << "    Unknown patch name: " << regionName_
            << ". Valid patch names are: "
            << mesh_.boundaryMesh().names() << nl
            << exit(FatalError);
    }

    const polyPatch& pp = mesh_.boundaryMesh()[patchid];

    label nFaces = pp.size();
    if (isA<emptyPolyPatch>(pp))
    {
        nFaces = 0;
    }

    faceId_.setSize(nFaces);
    facePatchId_.setSize(nFaces);
    faceSign_.setSize(nFaces);
    nFaces_ = returnReduce(faceId_.size(), sumOp<label>());

    forAll(faceId_, facei)
    {
        faceId_[facei] = facei;
        facePatchId_[facei] = patchid;
        faceSign_[facei] = 1;
    }
}


void Foam::functionObjects::fieldValues::surfaceFieldValue::sampledSurfaceFaces
(
    const dictionary& dict
)
{
    surfacePtr_ = sampledSurface::New
    (
        name(),
        mesh_,
        dict.subDict("sampledSurfaceDict")
    );
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
    if (regionType_ == regionTypes::sampledSurface)
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
    if (regionType_ == regionTypes::sampledSurface)
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
    switch (regionType_)
    {
        case regionTypes::faceZone:
        {
            dict.lookup("name") >> regionName_;
            setFaceZoneFaces();
            break;
        }
        case regionTypes::patch:
        {
            dict.lookup("name") >> regionName_;
            setPatchFaces();
            break;
        }
        case regionTypes::sampledSurface:
        {
            sampledSurfaceFaces(dict);
            regionName_ = surfacePtr_().name();
            break;
        }
        default:
        {
            FatalErrorInFunction
                << type() << " " << name() << ": "
                << regionTypeNames_[regionType_] << "(" << regionName_ << "):"
                << nl << "    Unknown region type. Valid region types are:"
                << regionTypeNames_.sortedToc() << nl << exit(FatalError);
        }
    }

    if (nFaces_ == 0)
    {
        FatalErrorInFunction
            << type() << " " << name() << ": "
            << regionTypeNames_[regionType_] << "(" << regionName_ << "):" << nl
            << "    Region has no faces" << exit(FatalError);
    }

    if (regionType_ == regionTypes::sampledSurface)
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

    if (dict.readIfPresent("scaleFactor", scaleFactor_))
    {
        Info<< "    scale factor = " << scaleFactor_ << nl;
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
        writeCommented(file(), "Region type : ");
        file() << regionTypeNames_[regionType_] << " " << regionName_ << endl;
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
        case operationType::sumDirection:
        {
            const vector n(dict_.lookup("direction"));
            result = gSum(weights*pos0(values*(Sf & n))*mag(values));
            return true;
        }
        case operationType::sumDirectionBalance:
        {
            const vector n(dict_.lookup("direction"));
            const scalarField nv(values*(Sf & n));
            result = gSum(weights*(pos0(nv)*mag(values) - neg(nv)*mag(values)));
            return true;
        }
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
            return false;
        }
    }
}


bool Foam::functionObjects::fieldValues::surfaceFieldValue::processValues
(
    const Field<vector>& values,
    const scalarField& signs,
    const scalarField& weights,
    const vectorField& Sf,
    vector& result
) const
{
    switch (operation_)
    {
        case operationType::sumDirection:
        {
            const vector n = normalised(dict_.lookup<vector>("direction"));
            const scalarField nv(n & values);
            result = gSum(weights*pos0(nv)*n*(nv));
            return true;
        }
        case operationType::sumDirectionBalance:
        {
            const vector n = normalised(dict_.lookup<vector>("direction"));
            const scalarField nv(n & values);
            result = gSum(weights*pos0(nv)*n*(nv));
            return true;
        }
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
    regionType_(regionTypeNames_.read(dict.lookup("regionType"))),
    operation_(operationTypeNames_.read(dict.lookup("operation"))),
    weightFieldNames_(),
    scaleFactor_(1),
    writeArea_(dict.lookupOrDefault("writeArea", false)),
    nFaces_(0),
    faceId_(),
    facePatchId_(),
    faceSign_()
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
    regionType_(regionTypeNames_.read(dict.lookup("regionType"))),
    operation_(operationTypeNames_.read(dict.lookup("operation"))),
    weightFieldNames_(),
    scaleFactor_(1),
    writeArea_(dict.lookupOrDefault("writeArea", false)),
    nFaces_(0),
    faceId_(),
    facePatchId_(),
    faceSign_()
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
    initialise(dict);

    return true;
}


bool Foam::functionObjects::fieldValues::surfaceFieldValue::write()
{
    if (operation_ != operationType::none)
    {
        fieldValue::write();
    }

    if (regionType_ == regionTypes::sampledSurface)
    {
        surfacePtr_().update();
    }

    if (operation_ != operationType::none && Pstream::master())
    {
        writeTime(file());
    }

    if (writeArea_)
    {
        totalArea_ = totalArea();
        if (operation_ != operationType::none && Pstream::master())
        {
            file() << tab << totalArea_;
        }
        Log << "    total area = " << totalArea_ << endl;
    }

    // Write the surface geometry
    if (writeFields_)
    {
        faceList faces;
        pointField points;

        if (regionType_ == regionTypes::sampledSurface)
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
                outputDir(),
                regionTypeNames_[regionType_] + ("_" + regionName_),
                points,
                faces
            );
        }
    }

    // Construct the sign and weight fields and the surface normals
    const scalarField signs
    (
        regionType_ == regionTypes::sampledSurface
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
        regionType_ == regionTypes::sampledSurface
      ? surfacePtr_().Sf()
      : (signs*filterField(mesh_.Sf()))()
    );

    // Process the fields
    forAll(fields_, i)
    {
        const word& fieldName = fields_[i];
        bool ok = false;

        #define writeValuesFieldType(fieldType, none)                          \
            ok =                                                               \
                ok                                                             \
             || writeValues<fieldType>                                         \
                (                                                              \
                    fieldName,                                                 \
                    signs,                                                     \
                    weights,                                                   \
                    Sf                                                         \
                );
        FOR_ALL_FIELD_TYPES(writeValuesFieldType);
        #undef writeValuesFieldType

        if (!ok)
        {
            WarningInFunction
                << "Requested field " << fieldName
                << " not found in database and not processed"
                << endl;
        }
    }

    if (operation_ != operationType::none && Pstream::master())
    {
        file() << endl;
    }

    Log << endl;

    return true;
}


// ************************************************************************* //
