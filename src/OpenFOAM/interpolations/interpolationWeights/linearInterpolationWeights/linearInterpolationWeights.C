/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2019 OpenFOAM Foundation
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
    // Check if current index is still valid
    bool changed = false;
    if
    (
        (
            index_ == -1
         && t <= samples_.first()
        )
     || (
            index_ >= 0
         && index_ < samples_.size() - 1
         && t >= samples_[index_]
         && t <= samples_[index_ + 1]
        )
     || (
            index_ == samples_.size() - 1
         && t >= samples_.last()
        )
    )
    {
        // The index is still in the correct slot
    }
    else
    {
        // The index is no longer in the correct slot, so search for a new one
        index_ = findLower(samples_, t);
        changed = true;
    }

    // Calculate the number of indices and weights and resize the result
    const label n = index_ == -1 || index_ == samples_.size() - 1 ? 1 : 2;
    indices.resize(n);
    weights.resize(n);

    // Compute the value
    if (index_ == -1)
    {
        // Use the first value
        indices[0] = 0;
        weights[0] = 1;
    }
    else if (index_ == samples_.size() - 1)
    {
        // Use the last value
        indices[0] = samples_.size() - 1;
        weights[0] = 1;
    }
    else
    {
        // Interpolate within the interval
        const scalar f =
            (t - samples_[index_])/(samples_[index_ + 1] - samples_[index_]);
        indices[0] = index_;
        weights[0] = 1 - f;
        indices[1] = index_ + 1;
        weights[1] = f;
    }

    return changed;
}


bool linearInterpolationWeights::integrationWeights
(
    scalar t1,
    scalar t2,
    labelList& indices,
    scalarField& weights
) const
{
    // If the arguments are in descending order, then swap them and set the
    // weights' sign negative
    label sign = +1;
    if (t1 > t2)
    {
        Swap(t1, t2);
        sign = -1;
    }

    //- Search for lower indices
    //  Note: currently there is no caching of this search like in valueWeights
    const label i1 = findLower(samples_, t1);
    const label i2 = findLower(samples_, t2);
    const label iClip1 = min(max(i1, 0), samples_.size() - 2);
    const label iClip2 = min(max(i2, 0), samples_.size() - 2);
    const label n = max(i2 - i1 + (i1 == iClip1) + (i2 == iClip2), 1);

    // Check if anything changed
    bool changed = false;
    if (indices.size() != n)
    {
        changed = true;
    }
    else
    {
        forAll(indices, indexi)
        {
            if (indices[indexi] == indexi + max(i1, 0))
            {
                changed = true;
                break;
            }
        }
    }

    // Resize the result arrays
    indices.resize(n);
    indices = -1;
    weights.resize(n);
    weights = 0;

    // Add out of bounds interval below the table
    if (i1 == -1)
    {
        indices[0] = 0;
        weights[0] += sign*(samples_[0] - t1);
    }
    if (i2 == -1)
    {
        indices[0] = 0;
        weights[0] -= sign*(samples_[0] - t2);
    }

    // Add partial interval from t1 to i1 + 1
    if (i1 == iClip1)
    {
        const scalar f = (t1 - samples_[i1])/(samples_[i1 + 1] - samples_[i1]);
        const scalar d = samples_[i1 + 1] - t1;
        indices[0] = i1;
        weights[0] += sign*(1 - f)*d/2;
        indices[1] = i1 + 1;
        weights[1] += sign*(1 + f)*d/2;
    }

    // Sum whole intervals from i1 + 1 to i2
    if (i1 != i2) for (label i = i1 + 1; i <= iClip2; i ++)
    {
        const scalar d = samples_[i + 1] - samples_[i];
        indices[i - iClip1] = i;
        weights[i - iClip1] += sign*d/2;
        indices[i - iClip1 + 1] = i + 1;
        weights[i - iClip1 + 1] += sign*d/2;
    }

    // Subtract partial interval from t2 to i2 + 1
    if (i2 == iClip2)
    {
        const scalar f = (t2 - samples_[i2])/(samples_[i2 + 1] - samples_[i2]);
        const scalar d = samples_[i2 + 1] - t2;
        indices[n - 2] = i2;
        weights[n - 2] -= sign*(1 - f)*d/2;
        indices[n - 1] = i2 + 1;
        weights[n - 1] -= sign*(1 + f)*d/2;
    }

    // Add out of bounds interval above the table
    if (i1 == samples_.size() - 1)
    {
        indices[n - 1] = samples_.size() - 1;
        weights[n - 1] -= sign*(t1 - samples_.last());
    }
    if (i2 == samples_.size() - 1)
    {
        indices[n - 1] = samples_.size() - 1;
        weights[n - 1] += sign*(t2 - samples_.last());
    }

    return changed;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
