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

#include "rigidBodySectionalForceGraph.H"
#include "motionSolver_fvMeshMover.H"
#include "OSspecific.H"
#include "patchCutPlot.H"
#include "setWriter.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(rigidBodySectionalForceGraph, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        rigidBodySectionalForceGraph,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::rigidBodySectionalForceGraph::clear()
{
    sectionalForcesBase::clear();

    distancesPtr_.clear();
}


Foam::tmp<Foam::scalarField>
Foam::functionObjects::rigidBodySectionalForceGraph::distances() const
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

Foam::functionObjects::rigidBodySectionalForceGraph::
rigidBodySectionalForceGraph
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    rigidBodySectionalForcesBase(name, runTime, dict),
    distancesPtr_(nullptr),
    nPoints_(-1),
    formatter_(nullptr),
    nOptimiseIter_(-1)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::rigidBodySectionalForceGraph::
~rigidBodySectionalForceGraph()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::rigidBodySectionalForceGraph::read
(
    const dictionary& dict
)
{
    rigidBodySectionalForcesBase::read(dict);

    nPoints_ = dict.lookup<label>("nPoints");

    formatter_ = setWriter::New(dict.lookup("setFormat"), dict);

    nOptimiseIter_ = dict.lookupOrDefault<label>("nOptimiseIter", 2);

    return true;
}


bool Foam::functionObjects::rigidBodySectionalForceGraph::write()
{
    // Get the coordinates
    tmp<scalarField> tdistances = this->distances();
    const scalarField& distances = tdistances();

    // Construct result fields
    wordList fieldNames
    ({
        "fluidForce",
        "bodyForce",
        "totalForce",
        "fluidMoment",
        "bodyMoment",
        "totalMoment"
    });

    #define DeclareTypeFieldValues(Type, nullArg) \
        PtrList<Field<Type>> Type##FieldValues(fieldNames.size());
    FOR_ALL_FIELD_TYPES(DeclareTypeFieldValues);
    #undef DeclareTypeFieldValues

    forAll(vectorFieldValues, i)
    {
        vectorFieldValues.set
        (
            i,
            new vectorField(distances.size(), vector::zero)
        );
    }

    vectorField& fluidForce = vectorFieldValues[0];
    vectorField& bodyForce = vectorFieldValues[1];
    vectorField& totalForce = vectorFieldValues[2];
    vectorField& fluidMoment = vectorFieldValues[3];
    vectorField& bodyMoment = vectorFieldValues[4];
    vectorField& totalMoment = vectorFieldValues[5];

    // Calculate the fluid contribution
    rigidBodySectionalForcesBase::addFluid(fluidForce, fluidMoment);

    // Calculate the body contribution
    rigidBodySectionalForcesBase::addBody(bodyForce, bodyMoment);

    // Sum to create total sectional forces and moments
    totalForce = fluidForce + bodyForce;
    totalMoment = fluidMoment + bodyMoment;

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
                true,
                axisName(),
                distances + localOrigin()[axisi()],
                coordSet::axisTypeNames_[coordSet::axisType::DISTANCE]
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
