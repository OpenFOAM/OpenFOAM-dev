/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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

#include "linearInterpolationWeights.H"
#include "addToRunTimeSelectionTable.H"
#include "ListOps.H"
#include "Pair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(linearInterpolationWeights, 0);
addToRunTimeSelectionTable
(
    interpolationWeights,
    linearInterpolationWeights,
    word
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::Pair<Foam::scalar> linearInterpolationWeights::integrationWeights
(
    const label i,
    const scalar t
) const
{
    // t is in range samples_[i] .. samples_[i+1]

    scalar s = (t-samples_[i])/(samples_[i+1]-samples_[i]);

    if (s < -small || s > 1+small)
    {
        FatalErrorInFunction
            << "Value " << t << " outside range " << samples_[i]
            << " .. " << samples_[i+1]
            << exit(FatalError);
    }

    scalar d = samples_[i+1]-t;

    return Pair<scalar>(d*0.5*(1-s), d*0.5*(1+s));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

linearInterpolationWeights::linearInterpolationWeights
(
    const scalarField& samples
)
:
    interpolationWeights(samples),
    index_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool linearInterpolationWeights::valueWeights
(
    const scalar t,
    labelList& indices,
    scalarField& weights
) const
{
    bool indexChanged = false;

    // Check if current timeIndex is still valid
    if
    (
        index_ >= 0
     && index_ < samples_.size()
     && (
            samples_[index_] <= t
         && (index_ == samples_.size()-1 || t <= samples_[index_+1])
        )
    )
    {
        // index_ still at correct slot
    }
    else
    {
        // search for correct index
        index_ = findLower(samples_, t);
        indexChanged = true;
    }


    if (index_ == -1)
    {
        // Use first element only
        indices.setSize(1);
        weights.setSize(1);
        indices[0] = 0;
        weights[0] = 1.0;
    }
    else if (index_ == samples_.size()-1)
    {
        // Use last element only
        indices.setSize(1);
        weights.setSize(1);
        indices[0] = samples_.size()-1;
        weights[0] = 1.0;
    }
    else
    {
        // Interpolate
        indices.setSize(2);
        weights.setSize(2);

        indices[0] = index_;
        indices[1] = index_+1;

        scalar t0 = samples_[indices[0]];
        scalar t1 = samples_[indices[1]];
        scalar deltaT = t1-t0;

        weights[0] = (t1-t)/deltaT;
        weights[1] = 1.0-weights[0];
    }

    return indexChanged;
}


bool linearInterpolationWeights::integrationWeights
(
    const scalar t1,
    const scalar t2,
    labelList& indices,
    scalarField& weights
) const
{
    if (t2 < t1-vSmall)
    {
        FatalErrorInFunction
            << "Integration should be in positive direction."
            <<  " t1:" << t1 << " t2:" << t2
            << exit(FatalError);
    }

    // Currently no fancy logic on cached index like in value

    //- Find lower or equal index
    label i1 = findLower(samples_, t1, 0, lessEqOp<scalar>());
    //- Find lower index
    label i2 = findLower(samples_, t2);

    // For now just fail if any outside table
    if (i1 == -1 || i2 == samples_.size()-1)
    {
        FatalErrorInFunction
            << "Integrating outside table " << samples_[0] << ".."
            << samples_.last() << " not implemented."
            << " t1:" << t1 << " t2:" << t2 << exit(FatalError);
    }

    label nIndices = i2-i1+2;


    // Determine if indices already correct
    bool anyChanged = false;

    if (nIndices != indices.size())
    {
        anyChanged = true;
    }
    else
    {
        // Closer check

        label index = i1;
        forAll(indices, i)
        {
            if (indices[i] != index)
            {
                anyChanged = true;
                break;
            }
            index++;
        }
    }

    indices.setSize(nIndices);
    weights.setSize(nIndices);
    weights = 0.0;

    // Sum from i1+1 to i2+1
    for (label i = i1+1; i <= i2; i++)
    {
        scalar d = samples_[i+1]-samples_[i];
        indices[i-i1] = i;
        weights[i-i1] += 0.5*d;
        indices[i+1-i1] = i+1;
        weights[i+1-i1] += 0.5*d;
    }

    // Add from i1 to t1
    {
        Pair<scalar> i1Tot1 = integrationWeights(i1, t1);
        indices[0] = i1;
        weights[0] += i1Tot1.first();
        indices[1] = i1+1;
        weights[1] += i1Tot1.second();
    }

    // Subtract from t2 to i2+1
    {
        Pair<scalar> wghts = integrationWeights(i2, t2);
        indices[i2-i1] = i2;
        weights[i2-i1] += -wghts.first();
        indices[i2-i1+1] = i2+1;
        weights[i2-i1+1] += -wghts.second();
    }

    return anyChanged;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
