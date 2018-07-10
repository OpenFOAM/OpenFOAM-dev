/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Equality operator for cellShape class

\*---------------------------------------------------------------------------*/

#include "cellShape.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool Foam::operator==(const cellShape& a, const cellShape& b)
{
    // Basic rule: we assume that the sequence of labels in each list
    // will be circular in the same order (but not necessarily in the
    // same direction). The limitation of this method is that with 3D
    // topologies I cannot guarantee that a congruent but not
    // identical cellShape (i.e. one sharing the same points but in a
    // different order) will necessarily be matched.

    const labelList& labsA = a;
    const labelList& labsB = b;

    // Trivial reject: faces are different size
    label sizeA = labsA.size();
    label sizeB = labsB.size();
    if (sizeA != sizeB)
    {
        return false;
    }


    // First we look for the occurrence of the first label in A, in B

    label Bptr = -1;
    label firstA = labsA[0];
    forAll(labsB, i)
    {
        if (labsB[i] == firstA)
        {
            Bptr = i;                // Denotes 'found match' at element 'i'
            break;
        }
    }

    // If no match was found, exit false
    if (Bptr < 0)
    {
        return false;
    }

    // Now we must look for the direction, if any
    label secondA = labsA[1];
    label dir = 0;

    // Check whether at top of list
    Bptr++;
    if (Bptr == labsB.size())
    {
        Bptr = 0;
    }

    // Test whether upward label matches second A label
    if (labsB[Bptr] == secondA)
    {
        // Yes - direction is 'up'
        dir = 1;
    }
    else
    {
        // No - so look downwards, checking whether at bottom of list
        Bptr -= 2;
        if (Bptr < 0)
        {
            // Case (1) Bptr=-1
            if (Bptr == -1)
            {
                Bptr = labsB.size() - 1;
            }

            // Case (2) Bptr = -2
            else
            {
                Bptr = labsB.size() - 2;
            }
        }

        // Test whether downward label matches second A label
        if (labsB[Bptr] == secondA)
        {
            // Yes - direction is 'down'
            dir = -1;
        }
    }

    // Check whether a match was made at all, and exit false if not
    if (dir == 0)
    {
        return false;
    }

    // Decrement size by 2 to account for first searches
    sizeA -= 2;

    // We now have both direction of search and next element
    // to search, so we can continue search until no more points.
    label Aptr = 1;
    if (dir > 0)
    {
        while (sizeA--)
        {
            Aptr++;
            if (Aptr >= labsA.size())
            {
                Aptr = 0;
            }

            Bptr++;
            if (Bptr >= labsB.size())
            {
                Bptr = 0;
            }
            if (labsA[Aptr] != labsB[Bptr])
            {
                return false;
            }
        }
    }
    else
    {
        while (sizeA--)
        {
            Aptr++;
            if (Aptr >= labsA.size())
            {
                Aptr = 0;
            }

            Bptr--;
            if (Bptr < 0)
            {
                Bptr = labsB.size() - 1;
            }
            if (labsA[Aptr] != labsB[Bptr])
            {
                return false;
            }
        }
    }

    // They must be equal
    return true;
}


// ************************************************************************* //
