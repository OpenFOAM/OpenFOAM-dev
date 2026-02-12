/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "sectionalForceGraph.H"
#include "fvMesh.H"
#include "OSspecific.H"
#include "patchCutPlot.H"
#include "setWriter.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(sectionalForceGraph, 0);
    addToRunTimeSelectionTable(functionObject, sectionalForceGraph, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::sectionalForceGraph::clear()
{
    sectionalForcesBase::clear();

    distancesPtr_.clear();
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::vector Foam::functionObjects::sectionalForceGraph::normal() const
{
    return normal_;
}


Foam::point Foam::functionObjects::sectionalForceGraph::origin() const
{
    return origin_;
}


Foam::tmp<Foam::scalarField>
Foam::functionObjects::sectionalForceGraph::distances() const
{
    if (!distancesPtr_.valid())
    {
        distancesPtr_.set
        (
            patchCutPlot::calcCutXs
            (
                patch(),
                patchPointDistances(),
                false,
                nPoints_,
                nOptimiseIter_,
                debug,
                name(),
                mesh(),
                formatter_()
            ).ptr()
        );
    }

    return distancesPtr_();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::sectionalForceGraph::sectionalForceGraph
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    sectionalForcesBase(name, runTime, dict),
    origin_(vector::uniform(NaN)),
    distancesPtr_(nullptr),
    nPoints_(-1),
    formatter_(nullptr),
    nOptimiseIter_(-1)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::sectionalForceGraph::~sectionalForceGraph()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::sectionalForceGraph::read(const dictionary& dict)
{
    sectionalForcesBase::read(dict);

    normal_ = normalised(dict.lookup<vector>("normal"));
    origin_ = dict.lookupBackwardsCompatible<vector>({"origin", "CofR"});

    nPoints_ = dict.lookup<label>("nPoints");

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

    nOptimiseIter_ = dict.lookupOrDefault<label>("nOptimiseIter", 2);

    return true;
}


bool Foam::functionObjects::sectionalForceGraph::write()
{
    // Get the coordinates
    tmp<scalarField> tdistances = this->distances();
    const scalarField& distances = tdistances();

    // Construct result fields
    wordList fieldNames({"force", "moment"});

    #define DeclareTypeFieldValues(Type, nullArg) \
        PtrList<Field<Type>> Type##FieldValues(2);
    FOR_ALL_FIELD_TYPES(DeclareTypeFieldValues);
    #undef DeclareTypeFieldValues

    vectorFieldValues.set(0, new vectorField(distances.size(), vector::zero));
    vectorFieldValues.set(1, new vectorField(distances.size(), vector::zero));

    vectorField& force = vectorFieldValues[0];
    vectorField& moment = vectorFieldValues[1];

    // Add the fluid contribution
    sectionalForcesBase::addFluid(force, moment);

    // Write out the graph
    if (Pstream::master())
    {
        mkDir(outputPath());

        formatter_->write
        (
            outputPath(),
            typeName,
            coordSet
            (
                identityMap(distances.size()),
                word::null,
                origin() + distances*normal_,
                coordSet::axisTypeNames_[coordSet::axisType::DISTANCE],
                distances,
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


// ************************************************************************* //
