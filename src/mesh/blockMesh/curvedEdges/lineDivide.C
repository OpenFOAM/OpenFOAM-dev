/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "lineDivide.H"
#include "curvedEdge.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    //- Calculate the geometric expension factor from the expansion ratio
    inline scalar calcGexp(const scalar expRatio, const label nDiv)
    {
        return nDiv > 1 ? pow(expRatio, 1.0/(nDiv - 1)) : 0.0;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::lineDivide::lineDivide
(
    const curvedEdge& cedge,
    const label nDiv,
    const gradingDescriptors& gd
)
:
    points_(nDiv + 1),
    divisions_(nDiv + 1)
{
    divisions_[0]    = 0.0;
    divisions_[nDiv] = 1.0;

    scalar secStart = divisions_[0];
    label secnStart = 1;

    // Check that there are more divisions than sections
    if (nDiv >= gd.size())
    {
        forAll(gd, sectioni)
        {
            scalar blockFrac = gd[sectioni].blockFraction();
            scalar nDivFrac = gd[sectioni].nDivFraction();
            scalar expRatio = gd[sectioni].expansionRatio();

            label secnDiv = label(nDivFrac*nDiv + 0.5);
            if (sectioni == gd.size() - 1)
            {
                secnDiv = nDiv - secnStart + 1;
            }
            label secnEnd = secnStart + secnDiv;

            // Calculate the spacing
            if (expRatio == 1.0)
            {
                for (label i = secnStart; i < secnEnd; i++)
                {
                    divisions_[i] =
                        secStart
                      + blockFrac*scalar(i - secnStart + 1)/secnDiv;
                }
            }
            else
            {
                // Calculate geometric expansion factor from the expansion ratio
                const scalar expFact = calcGexp(expRatio, secnDiv);

                for (label i = secnStart; i < secnEnd; i++)
                {
                    divisions_[i] =
                        secStart
                      + blockFrac*(1.0 - pow(expFact, i - secnStart + 1))
                    /(1.0 - pow(expFact, secnDiv));
                }
            }

            secStart = divisions_[secnEnd - 1];
            secnStart = secnEnd;
        }
    }
    // Otherwise mesh uniformly
    else
    {
        for (label i=1; i < nDiv; i++)
        {
            divisions_[i] = scalar(i)/nDiv;
        }
    }


    // Calculate the points
    for (label i = 0; i <= nDiv; i++)
    {
        points_[i] = cedge.position(divisions_[i]);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::pointField& Foam::lineDivide::points() const
{
    return points_;
}


const Foam::scalarList& Foam::lineDivide::lambdaDivisions() const
{
    return divisions_;
}


// ************************************************************************* //
