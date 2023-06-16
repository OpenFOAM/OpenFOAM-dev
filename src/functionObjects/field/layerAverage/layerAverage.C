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

#include "layerAverage.H"
#include "FaceCellWave.H"
#include "layerInfo.H"
#include "regionSplit.H"
#include "syncTools.H"
#include "volFields.H"
#include "writeFile.H"
#include "polyTopoChangeMap.H"
#include "polyMeshMap.H"
#include "polyDistributionMap.H"
#include "OSspecific.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace functionObjects
    {
        defineTypeNameAndDebug(layerAverage, 0);
        addToRunTimeSelectionTable(functionObject, layerAverage, dictionary);
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::layerAverage::calcLayers()
{
    // Initialise the faces from which the layers extrude
    DynamicList<label> startFaces;
    DynamicList<layerInfo> startFacesInfo;
    forAll(patchIDs_, i)
    {
        const polyPatch& pp = mesh_.boundaryMesh()[patchIDs_[i]];
        forAll(pp, j)
        {
            startFaces.append(pp.start() + j);
            startFacesInfo.append(layerInfo(0, -1));
        }
    }
    forAll(zoneIDs_, i)
    {
        const faceZone& fz = mesh_.faceZones()[zoneIDs_[i]];
        forAll(fz, j)
        {
            startFaces.append(fz[j]);
            startFacesInfo.append(layerInfo(0, fz.flipMap()[j] ? -1 : +1));
        }
    }

    // Wave to generate layer indices
    List<layerInfo> allFaceLayerInfo(mesh_.nFaces());
    List<layerInfo> allCellLayerInfo(mesh_.nCells());
    FaceCellWave<layerInfo> wave
    (
        mesh_,
        startFaces,
        startFacesInfo,
        allFaceLayerInfo,
        allCellLayerInfo,
        mesh_.globalData().nTotalCells() + 1
    );

    // Copy indices out of the wave and determine the total number of layers
    nLayers_ = 0;
    cellLayer_ = labelList(mesh_.nCells(), -1);
    forAll(cellLayer_, celli)
    {
        if (allCellLayerInfo[celli].valid(wave.data()))
        {
            const label layeri = allCellLayerInfo[celli].cellLayer();
            nLayers_ = max(nLayers_, layeri + 1);
            cellLayer_[celli] = layeri;
        }
    }
    reduce(nLayers_, maxOp<label>());

    // Report
    if (nLayers_ != 0)
    {
        Info<< "    Detected " << nLayers_ << " layers" << nl << endl;
    }
    else
    {
        WarningInFunction<< "No layers detected" << endl;
    }

    // Write the indices for debugging
    if (debug)
    {
        tmp<volScalarField> cellLayer =
            volScalarField::New
            (
                "cellLayer",
                mesh_,
                dimensionedScalar(dimless, 0)
            );
        cellLayer.ref().primitiveFieldRef() = List<scalar>(cellLayer_);
        cellLayer.ref().write();
    }

    // Sum number of entries per layer
    layerCount_ = sum(scalarField(mesh_.nCells(), 1));

    // Average the cell centres
    layerCentre_ = sum(mesh_.cellCentres())/layerCount_;

    // If symmetric, keep only half of the coordinates
    if (symmetric_)
    {
        layerCentre_.setSize(nLayers_/2);
    }
}


template<>
Foam::vector
Foam::functionObjects::layerAverage::symmetricCoeff<Foam::vector>() const
{
    direction i = -1;

    switch (axis_)
    {
        case coordSet::axisType::X:
        case coordSet::axisType::Y:
        case coordSet::axisType::Z:
            i = label(axis_) - label(coordSet::axisType::X);
            break;
        case coordSet::axisType::XYZ:
        case coordSet::axisType::DISTANCE:
        case coordSet::axisType::DEFAULT:
            FatalErrorInFunction
                << "Symmetric layer average requested with "
                << coordSet::axisTypeNames_[axis_] << " axis. Symmetric "
                << "averaging is only possible along coordinate axes."
                << exit(FatalError);
            break;
    }

    vector result = vector::one;
    result[i] = -1;
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

    axis_ =
        coordSet::axisTypeNames_
        [
            dict.lookupOrDefault<word>
            (
                "axis",
                coordSet::axisTypeNames_[coordSet::axisType::DEFAULT]
            )
        ];

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

    // Write
    if (Pstream::master() && layerCentre_.size())
    {
        // Make output directory
        const fileName outputPath =
            time_.globalPath()
           /writeFile::outputPrefix
           /(mesh_.name() != polyMesh::defaultRegion ? mesh_.name() : word())
           /name()
           /time_.name();
        mkDir(outputPath);

        scalarField layerDistance(layerCentre_.size(), 0);
        for (label i = 1; i < layerCentre_.size(); ++ i)
        {
            layerDistance[i] =
                layerDistance[i-1] + mag(layerCentre_[i] - layerCentre_[i-1]);
        }

        formatter_->write
        (
            outputPath,
            typeName,
            coordSet
            (
                identityMap(layerCentre_.size()),
                word::null,
                layerCentre_,
                coordSet::axisTypeNames_[coordSet::axisType::DISTANCE],
                layerDistance,
                coordSet::axisTypeNames_[axis_]
            ),
            fieldNames
            #define TypeValueSetsParameter(Type, nullArg) , Type##ValueSets
            FOR_ALL_FIELD_TYPES(TypeValueSetsParameter)
            #undef TypeValueSetsParameter
        );
    }

    return true;
}


void Foam::functionObjects::layerAverage::movePoints(const polyMesh& mesh)
{
    if (&mesh == &mesh_)
    {
        Info<< type() << " " << name() << ":" << nl;
        calcLayers();
    }
}


void Foam::functionObjects::layerAverage::topoChange
(
    const polyTopoChangeMap& map
)
{
    if (&map.mesh() == &mesh_)
    {
        Info<< type() << " " << name() << ":" << nl;
        calcLayers();
    }
}


void Foam::functionObjects::layerAverage::mapMesh(const polyMeshMap& map)
{
    if (&map.mesh() == &mesh_)
    {
        Info<< type() << " " << name() << ":" << nl;
        calcLayers();
    }
}


void Foam::functionObjects::layerAverage::distribute
(
    const polyDistributionMap& map
)
{
    if (&map.mesh() == &mesh_)
    {
        Info<< type() << " " << name() << ":" << nl;
        calcLayers();
    }
}


// ************************************************************************* //
