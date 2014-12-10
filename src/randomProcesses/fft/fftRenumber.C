/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Description
    Multi-dimensional renumbering used in the Numerical Recipes
   fft routine. This version is recursive, so works in n-d :
   determined by the length of array nn

\*---------------------------------------------------------------------------*/

#include "fftRenumber.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// recursively evaluate the indexing necessary to do the folding of
// the fft data. We recurse until we have the indexing ncessary for
// the folding in all directions.

void fftRenumberRecurse
(
    List<complex>& data,
    List<complex>& renumData,
    const labelList& nn,
    label nnprod,
    label ii,
    label l1,
    label l2
)
{
    if (ii == nn.size())
    {
        // we've worked out the renumbering scheme. Now copy
        // the components across

        data[l1] = complex(renumData[l2].Re(),renumData[l2].Im());
    }
    else
    {
        // do another level of folding. First work out the
        // multiplicative value of the index

        nnprod /= nn[ii];
        label i_1(0);

        for (label i=0; i<nn[ii]; i++)
        {
            // now evaluate the indices (both from array 1 and to
            // array 2). These get multiplied by nnprod to (cumulatively)
            // find the real position in the list corresponding to
            // this set of indices.

            if (i<nn[ii]/2)
            {
                i_1 = i + nn[ii]/2;
            }
            else
            {
                i_1 = i - nn[ii]/2;
            }


            // go to the next level of recursion.

            fftRenumberRecurse
            (
                data,
                renumData,
                nn,
                nnprod,
                ii+1,
                l1+i*nnprod,
                l2+i_1*nnprod
            );
        }
    }
}


// fftRenumber : fold the n-d data array to get the fft components in
// the right places.

#include "fftRenumber.H"

void fftRenumber
(
    List<complex>& data,
    const labelList& nn
)
{
    List<complex> renumData(data);

    label nnprod(1);
    forAll(nn, i)
    {
        nnprod *= nn[i];
    }

    label ii(0), l1(0), l2(0);

    fftRenumberRecurse
    (
        data,
        renumData,
        nn,
        nnprod,
        ii,
        l1,
        l2
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
