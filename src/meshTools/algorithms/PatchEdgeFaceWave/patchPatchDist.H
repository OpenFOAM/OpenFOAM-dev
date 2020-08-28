/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

Class
    Foam::patchPatchDist

Description
    Like wallDist but calculates on a patch the distance to nearest neighbouring
    patches. Uses PatchEdgeFaceWave to do actual calculation.

SourceFiles
    patchPatchDist.C

\*---------------------------------------------------------------------------*/

#ifndef patchPatchDist_H
#define patchPatchDist_H

#include "scalarField.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//class polyMesh;
class polyPatch;

/*---------------------------------------------------------------------------*\
                       Class patchPatchDist Declaration
\*---------------------------------------------------------------------------*/

class patchPatchDist
:
    public scalarField
{
    // Private Member Data

        //- Patch to operate on
        const polyPatch& patch_;

        //- Patches to determine the distance to
        const labelHashSet nbrPatchIDs_;

        //- Number of unset faces.
        label nUnset_;

public:

    // Constructors

        //- Construct from patch and neighbour patches.
        patchPatchDist
        (
            const polyPatch& pp,
            const labelHashSet& nbrPatchIDs
        );


    //- Destructor
    virtual ~patchPatchDist();


    // Member Functions

        const scalarField& y() const
        {
            return *this;
        }

        label nUnset() const
        {
            return nUnset_;
        }

        //- Correct for mesh geom/topo changes
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
