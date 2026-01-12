/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025-2026 OpenFOAM Foundation
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

#include "cutLayerAverage.H"
#include "OSspecific.H"
#include "cellCutPlot.H"
#include "volPointInterpolation.H"
#include "writeFile.H"
#include "polyTopoChangeMap.H"
#include "polyMeshMap.H"
#include "polyDistributionMap.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(cutLayerAverage, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        cutLayerAverage,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::fileName Foam::functionObjects::cutLayerAverage::outputPath() const
{
    return
        time_.globalPath()
       /writeFile::outputPrefix
       /(mesh_.name() != polyMesh::defaultRegion ? mesh_.name() : word())
       /name()
       /time_.name();
}


void Foam::functionObjects::cutLayerAverage::calcWeights()
{
    const pointField& points = mesh_.points();

    // Calculate or lookup the point coordinates
    tmp<pointScalarField> tdistance =
        distanceName_ == word::null
      ? tmp<pointScalarField>(nullptr)
      : mesh_.foundObject<volScalarField>(distanceName_)
      ? volPointInterpolation::New(mesh_).interpolate
        (
            mesh_.lookupObject<volScalarField>(distanceName_)
        )
      : tmp<pointScalarField>
        (
            mesh_.lookupObject<pointScalarField>(distanceName_)
        );
    tmp<scalarField> tpointXs =
        distanceName_ == word::null
      ? points & direction_
      : tmp<scalarField>(tdistance().primitiveField());
    const scalarField& pointXs = tpointXs();

    // Calculate the cut coordinates and the weights
    const scalarField cutXs
    (
        cellCutPlot::calcCutXs
        (
            mesh(),
            zone_,
            pointXs,
            interpolate_,
            interpolate_ ? nLayers_ : nLayers_ + 1,
            nOptimiseIter_,
            debug,
            name(),
            formatter_()
        )
    );
    weights_.reset
    (
        new List<cellCutPlot::weight>
        (
            cellCutPlot::calcWeights
            (
                mesh(),
                zone_,
                pointXs,
                cutXs,
                interpolate_
            )
        )
    );

    // Calculate plot coordinates and widths
    if (interpolate_)
    {
        layerDistances_.reset(new scalarField(cutXs));
    }
    else
    {
        const SubField<scalar> distance0s(cutXs, nLayers_);
        const SubField<scalar> distance1s(cutXs, nLayers_, 1);
        layerDistances_.reset(((distance0s + distance1s)/2).ptr());
        layerThicknesses_.reset((distance1s - distance0s).ptr());
    }

    // Calculate the layer positions
    layerPositions_.reset
    (
        cutPlot::applyWeights<vector>
        (
            nLayers_,
            weights_,
            mesh().cellCentres()
        ).ptr()
    );

    if (debug)
    {
        cellCutPlot::writeLayers
        (
            mesh(),
            cellCutPlot::calcWeights
            (
                mesh(),
                zone_,
                pointXs,
                cutXs,
                interpolate_,
                false
            ),
            name()
        );
    }
}


void Foam::functionObjects::cutLayerAverage::clear()
{
    weights_.clear();
    layerDistances_.clear();
    layerThicknesses_.clear();
    layerPositions_.clear();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::cutLayerAverage::cutLayerAverage
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    zone_(mesh())
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::cutLayerAverage::~cutLayerAverage()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::cutLayerAverage::read(const dictionary& dict)
{
    zone_.read(dict);

    const bool haveDirection = dict.found("direction");
    const bool haveDistance = dict.found("distance");
    if (haveDirection == haveDistance)
    {
        FatalIOErrorInFunction(dict)
            << "keywords direction and distance both "
            << (haveDirection ? "" : "un") << "defined in "
            << "dictionary " << dict.name()
            << exit(FatalIOError);
    }
    else if (haveDirection)
    {
        direction_ = normalised(dict.lookup<vector>("direction"));
        distanceName_ = word::null;
    }
    else if (haveDistance)
    {
        direction_ = vector::nan;
        distanceName_ = dict.lookup<word>("distance");
    }

    nLayers_ = dict.lookup<label>("nPoints");

    interpolate_ = dict.lookupOrDefault<bool>("interpolate", false);

    fields_ = dict.lookup<wordList>("fields");

    axis_ =
        coordSet::axisTypeNames_
        [
            dict.lookupOrDefault<word>
            (
                "axis",
                coordSet::axisTypeNames_[coordSet::axisType::DEFAULT]
            )
        ];

    formatter_ = setWriter::New(dict.lookup("setFormat"), dict);

    nOptimiseIter_ = dict.lookupOrDefault("nOptimiseIter", 2);

    clear();

    return true;
}


Foam::wordList Foam::functionObjects::cutLayerAverage::fields() const
{
    wordList result(fields_);

    if (distanceName_ != word::null)
    {
        result.append(distanceName_);
    }

    return result;
}


bool Foam::functionObjects::cutLayerAverage::execute()
{
    return true;
}


bool Foam::functionObjects::cutLayerAverage::write()
{
    if (!weights_.valid())
    {
        calcWeights();
    }

    const bool writeThickness =
        !interpolate_
     && (
            axis_ == coordSet::axisType::DEFAULT
         || axis_ == coordSet::axisType::DISTANCE
        );

    // Create list of field names
    wordList fieldNames;
    if (writeThickness)
    {
        fieldNames.append("thickness");
    }
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

    // Calculate the field values
    #define DeclareTypeFieldValues(Type, nullArg) \
        PtrList<Field<Type>> Type##FieldValues(fieldNames.size());
    FOR_ALL_FIELD_TYPES(DeclareTypeFieldValues);
    #undef DeclareTypeFieldValues
    if (writeThickness)
    {
        scalarFieldValues.set(0, new scalarField(layerThicknesses_));
    }
    for (label fieldi = writeThickness; fieldi < fieldNames.size(); ++ fieldi)
    {
        #define CollapseTypeFields(Type, GeoField)                             \
            if (mesh_.foundObject<GeoField<Type>>(fieldNames[fieldi]))        \
            {                                                                 \
                const GeoField<Type>& field =                                 \
                    mesh_.lookupObject<GeoField<Type>>(fieldNames[fieldi]);   \
                                                                              \
                Type##FieldValues.set                                         \
                (                                                             \
                    fieldi,                                                   \
                    cutPlot::applyWeights<Type>(nLayers_, weights_, field)    \
                );                                                            \
            }
        FOR_ALL_FIELD_TYPES(CollapseTypeFields, VolField);
        FOR_ALL_FIELD_TYPES(CollapseTypeFields, VolInternalField);
        #undef CollapseTypeFields
    }

    // Write
    if (Pstream::master() && layerPositions_->size())
    {
        mkDir(outputPath());

        formatter_->write
        (
            outputPath(),
            typeName,
            coordSet
            (
                identityMap(layerPositions_->size()),
                word::null,
                layerPositions_,
                coordSet::axisTypeNames_[coordSet::axisType::DISTANCE],
                layerDistances_,
                coordSet::axisTypeNames_[axis_]
            ),
            fieldNames
            #define TypeFieldValuesParameter(Type, nullArg) , Type##FieldValues
            FOR_ALL_FIELD_TYPES(TypeFieldValuesParameter)
            #undef TypeFieldValuesParameter
        );
    }

    return true;
}


void Foam::functionObjects::cutLayerAverage::movePoints
(
    const polyMesh& mesh
)
{
    if (&mesh == &mesh_)
    {
        zone_.movePoints();
        clear();
    }
}


void Foam::functionObjects::cutLayerAverage::topoChange
(
    const polyTopoChangeMap& map
)
{
    if (&map.mesh() == &mesh_)
    {
        zone_.topoChange(map);
        clear();
    }
}


void Foam::functionObjects::cutLayerAverage::mapMesh
(
    const polyMeshMap& map
)
{
    if (&map.mesh() == &mesh_)
    {
        zone_.mapMesh(map);
        clear();
    }
}


void Foam::functionObjects::cutLayerAverage::distribute
(
    const polyDistributionMap& map
)
{
    if (&map.mesh() == &mesh_)
    {
        zone_.distribute(map);
        clear();
    }
}


// ************************************************************************* //
