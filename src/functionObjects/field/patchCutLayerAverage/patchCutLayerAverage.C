/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2026 OpenFOAM Foundation
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

#include "patchCutLayerAverage.H"
#include "OSspecific.H"
#include "patchCutPlot.H"
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
    defineTypeNameAndDebug(patchCutLayerAverage, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        patchCutLayerAverage,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::fileName Foam::functionObjects::patchCutLayerAverage::outputPath() const
{
    return
        time_.globalPath()
       /writeFile::outputPrefix
       /(mesh_.name() != polyMesh::defaultRegion ? mesh_.name() : word())
       /name()
       /time_.name();
}


void Foam::functionObjects::patchCutLayerAverage::calcWeights()
{
    const polyPatch& pp = mesh_.boundaryMesh()[patchName_];
    const pointField& points = pp.localPoints();

    // Calculate or lookup the point coordinates
    tmp<scalarField> tpointXs =
        distanceName_ == word::null
      ? points & direction_
      : (
            mesh_.foundObject<volScalarField>(distanceName_)
          ? volPointInterpolation::New(mesh_).interpolate
            (
                mesh_.lookupObject<volScalarField>(distanceName_)
            )
          : tmp<pointScalarField>
            (
                mesh_.lookupObject<pointScalarField>(distanceName_)
            )
        )().boundaryField()[pp.index()].patchInternalField();
    const scalarField& pointXs = tpointXs();

    // Calculate the cut coordinates and the weights
    const scalarField cutXs
    (
        patchCutPlot::calcCutXs
        (
            pp,
            pointXs,
            interpolate_,
            interpolate_ ? nLayers_ : nLayers_ + 1,
            nOptimiseIter_,
            debug,
            name(),
            mesh(),
            formatter_()
        )
    );
    weights_.reset
    (
        new List<patchCutPlot::weight>
        (
            patchCutPlot::calcWeights
            (
                pp,
                pointXs,
                cutXs,
                interpolate_
            )
        )
    );

    // Calculate layer coordinates and widths
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
        cutPlot::applyWeights
        (
            nLayers_,
            weights_,
            pointField(pp.faceCentres())
        ).ptr()
    );

    // Write the layers as a tensor field for debugging
    if (debug)
    {
        patchCutPlot::writeLayers
        (
            pp,
            patchCutPlot::calcWeights
            (
                pp,
                pointXs,
                cutXs,
                interpolate_,
                false
            ),
            name(),
            mesh()
        );
    }
}


void Foam::functionObjects::patchCutLayerAverage::clear()
{
    weights_.clear();
    layerDistances_.clear();
    layerThicknesses_.clear();
    layerPositions_.clear();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::patchCutLayerAverage::patchCutLayerAverage
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

Foam::functionObjects::patchCutLayerAverage::~patchCutLayerAverage()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::patchCutLayerAverage::read(const dictionary& dict)
{
    patchName_ = dict.lookup<word>("patch");

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


Foam::wordList Foam::functionObjects::patchCutLayerAverage::fields() const
{
    wordList result(fields_);

    if (distanceName_ != word::null)
    {
        result.append(distanceName_);
    }

    return result;
}


bool Foam::functionObjects::patchCutLayerAverage::execute()
{
    return true;
}


bool Foam::functionObjects::patchCutLayerAverage::write()
{
    if (!weights_.valid())
    {
        calcWeights();
    }

    const polyPatch& pp = mesh_.boundaryMesh()[patchName_];

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
        #define CollapseTypeFields(Type, nullArg)                             \
            if (mesh_.foundObject<VolField<Type>>(fieldNames[fieldi]))        \
            {                                                                 \
                const VolField<Type>& field =                                 \
                    mesh_.lookupObject<VolField<Type>>(fieldNames[fieldi]);   \
                                                                              \
                Type##FieldValues.set                                         \
                (                                                             \
                    fieldi,                                                   \
                    cutPlot::applyWeights                                     \
                    (                                                         \
                        nLayers_,                                             \
                        weights_,                                             \
                        field.boundaryField()[pp.index()]                     \
                    )                                                         \
                );                                                            \
            }
        FOR_ALL_FIELD_TYPES(CollapseTypeFields);
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


void Foam::functionObjects::patchCutLayerAverage::movePoints
(
    const polyMesh& mesh
)
{
    if (&mesh == &mesh_)
    {
        clear();
    }
}


void Foam::functionObjects::patchCutLayerAverage::topoChange
(
    const polyTopoChangeMap& map
)
{
    if (&map.mesh() == &mesh_)
    {
        clear();
    }
}


void Foam::functionObjects::patchCutLayerAverage::mapMesh
(
    const polyMeshMap& map
)
{
    if (&map.mesh() == &mesh_)
    {
        clear();
    }
}


void Foam::functionObjects::patchCutLayerAverage::distribute
(
    const polyDistributionMap& map
)
{
    if (&map.mesh() == &mesh_)
    {
        clear();
    }
}


// ************************************************************************* //
