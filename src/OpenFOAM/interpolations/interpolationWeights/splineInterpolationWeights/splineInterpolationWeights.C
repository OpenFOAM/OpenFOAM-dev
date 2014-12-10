/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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

#include "splineInterpolationWeights.H"
#include "addToRunTimeSelectionTable.H"
#include "ListOps.H"
#include "linearInterpolationWeights.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(splineInterpolationWeights, 0);
addToRunTimeSelectionTable
(
    interpolationWeights,
    splineInterpolationWeights,
    word
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

splineInterpolationWeights::splineInterpolationWeights
(
    const scalarField& samples,
    const bool checkEqualDistance
)
:
    interpolationWeights(samples),
    index_(-1)
{
    if (checkEqualDistance && samples_.size() > 2)
    {
        const scalar interval = samples_[1]-samples[0];
        for (label i = 2; i < samples_.size(); i++)
        {
            scalar d = samples_[i]-samples[i-1];

            if (mag(d-interval) > SMALL)
            {
                WarningIn
                (
                    "splineInterpolationWeights::splineInterpolationWeights"
                    "(const scalarField&)"
                )   << "Spline interpolation only valid for constant intervals."
                    << nl
                    << "Interval 0-1 : " << interval << nl
                    << "Interval " << i-1 << '-' << i << " : "
                    << d << endl;
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool splineInterpolationWeights::valueWeights
(
    const scalar t,
    labelList& indices,
    scalarField& weights
) const
{
    bool indexChanged = false;

    // linear interpolation
    if (samples_.size() <= 2)
    {
        return linearInterpolationWeights(samples_).valueWeights
        (
            t,
            indices,
            weights
        );
    }

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


    // Clamp if outside table
    if (index_ == -1)
    {
        indices.setSize(1);
        weights.setSize(1);

        indices[0] = 0;
        weights[0] = 1;
        return indexChanged;
    }
    else if (index_ == samples_.size()-1)
    {
        indices.setSize(1);
        weights.setSize(1);

        indices[0] = samples_.size()-1;
        weights[0] = 1;
        return indexChanged;
    }



    label lo = index_;
    label hi = index_+1;

    // weighting
    scalar mu = (t - samples_[lo])/(samples_[hi] - samples_[lo]);

    scalar w0 = 0.5*(mu*(-1+mu*(2-mu)));            // coeff of lo-1
    scalar w1 = 0.5*(2+mu*(mu*(-5 + mu*(3))));      // coeff of lo
    scalar w2 = 0.5*(mu*(1 + mu*(4 + mu*(-3))));    // coeff of hi
    scalar w3 = 0.5*(mu*mu*(-1 + mu));              // coeff of hi+1

    if (lo > 0)
    {
        if (hi < samples_.size()-1)
        {
            // Four points available
            indices.setSize(4);
            weights.setSize(4);

            indices[0] = lo-1;
            indices[1] = lo;
            indices[2] = hi;
            indices[3] = hi+1;

            weights[0] = w0;
            weights[1] = w1;
            weights[2] = w2;
            weights[3] = w3;
        }
        else
        {
            // No y3 available. Extrapolate: y3=3*y2-y1
            indices.setSize(3);
            weights.setSize(3);

            indices[0] = lo-1;
            indices[1] = lo;
            indices[2] = hi;

            weights[0] = w0;
            weights[1] = w1 - w3;
            weights[2] = w2 + 2*w3;
        }
    }
    else
    {
        // No y0 available. Extrapolate: y0=2*y1-y2;
        if (hi < samples_.size()-1)
        {
            indices.setSize(3);
            weights.setSize(3);

            indices[0] = lo;
            indices[1] = hi;
            indices[2] = hi+1;

            weights[0] = w1 + 2*w0;
            weights[1] = w2 - w0;
            weights[2] = w3;
        }
        else
        {
            indices.setSize(2);
            weights.setSize(2);

            indices[0] = lo;
            indices[1] = hi;

            weights[0] = w1 + 2*w0 - w3;
            weights[1] = w2 - w0 + 2*w3;
        }
    }

    return indexChanged;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
