/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
    Foam::patchWave

Description
    Takes a set of patches to start MeshWave from. After construction holds
    distance at cells and distance at patches. Is e.g. used by wallDist to
    construct volScalarField with correct distance to wall.

SourceFiles
    patchWave.C

\*---------------------------------------------------------------------------*/

#ifndef patchWave_H
#define patchWave_H

#include "cellDistFuncs.H"
#include "FieldField.H"
#include "MeshWave.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;
class wallPoint;

/*---------------------------------------------------------------------------*\
                          Class patchWave Declaration
\*---------------------------------------------------------------------------*/

class patchWave
:
    public cellDistFuncs
{
    // Private Data

        //- Current patch subset (stored as patchIDs)
        labelHashSet patchIDs_;

        //- Do accurate distance calculation for near-wall cells.
        bool correctWalls_;

        //- Number of cells/faces unset after MeshWave has finished
        label nUnset_;

        //- Distance at cell centres
        scalarField distance_;

        //- Distance at patch faces
        FieldField<Field, scalar> patchDistance_;


    // Private Member Functions

        //- Set initial set of changed faces (= all faces of patches in
        //  patchIDs). Set changedFaces to labels of changed faces,
        //  changedInfo to face centres.
        void setChangedFaces
        (
            const labelHashSet& patchIDs,
            labelList& changedFaces,
            List<wallPoint>& changedInfo
        ) const;

        //- Copy MeshWave cell values. Return number of illegal/unset
        //  cells.
        label getValues(const MeshWave<wallPoint>&);


public:

    // Constructors

        //- Construct from mesh and patches to initialise to 0 and flag
        //  whether or not to correct wall.
        //  Calculate for all cells. correctWalls : correct wall (face&point)
        //  cells for correct distance, searching neighbours.
        patchWave
        (
            const polyMesh& mesh,
            const labelHashSet& patchIDs,
            bool correctWalls = true
        );

    //- Destructor
    virtual ~patchWave();


    // Member Functions

        //- Correct for mesh geom/topo changes
        virtual void correct();


        label nUnset() const
        {
            return nUnset_;
        }

        const scalarField& distance() const
        {
            return distance_;
        }

        //- Non const access so we can 'transfer' contents for efficiency
        scalarField& distance()
        {
            return distance_;
        }

        const FieldField<Field, scalar>& patchDistance() const
        {
            return patchDistance_;
        }

        FieldField<Field, scalar>& patchDistance()
        {
            return patchDistance_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
