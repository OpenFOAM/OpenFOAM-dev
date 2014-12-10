/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::patchDataWave

Description
    Takes a set of patches to start MeshWave from.

    Holds after construction distance at cells and distance at patches
    (like patchWave), but also additional transported data.
    It is used, for example, in the y+ calculation.

See Also
   The patchWave class.

SourceFiles
    patchDataWave.C

\*---------------------------------------------------------------------------*/

#ifndef patchDataWave_H
#define patchDataWave_H

#include "cellDistFuncs.H"
#include "FieldField.H"
#include "UPtrList.H"
#include "MeshWave.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;
class wallPoint;

/*---------------------------------------------------------------------------*\
                        Class patchDataWave Declaration
\*---------------------------------------------------------------------------*/

template<class TransferType>
class patchDataWave
:
    public cellDistFuncs
{

private:

    typedef typename TransferType::dataType Type;


    // Private Member Data

        //- Current patch subset (stored as patchIDs)
        labelHashSet patchIDs_;

        //- Reference to initial extra data at patch faces
        const UPtrList<Field<Type> >& initialPatchValuePtrs_;

        //- Do accurate distance calculation for near-wall cells.
        bool correctWalls_;

        //
        // After construction:
        //

        //- Number of cells/faces unset after MeshWave has finished
        label nUnset_;

        //- Distance at cell centres
        scalarField distance_;

        //- Distance at patch faces
        FieldField<Field, scalar> patchDistance_;

        //- Extra data at cell centres
        Field<Type> cellData_;

        //- Extra data at patch faces
        FieldField<Field, Type> patchData_;


    // Private Member Functions

        //- Set initial set of changed faces
        void setChangedFaces
        (
            const labelHashSet& patchIDs,
            labelList&,
            List<TransferType>&
        ) const;

        //- Copy MeshWave values into *this
        label getValues(const MeshWave<TransferType>&);


public:

    // Constructors

        //- Construct from mesh, information on patches to initialize and flag
        //  whether or not to correct wall.
        //  Calculate for all cells. correctWalls : correct wall (face&point)
        //  cells for correct distance, searching neighbours.
        patchDataWave
        (
            const polyMesh& mesh,
            const labelHashSet& patchIDs,
            const UPtrList<Field<Type> >& initialPatchValuePtrs,
            bool correctWalls = true
        );


    //- Destructor
    virtual ~patchDataWave();


    // Member Functions

        //- Correct for mesh geom/topo changes
        virtual void correct();


        const scalarField& distance() const
        {
            return distance_;
        }

        //- Non const access so we can 'transfer' contents for efficiency.
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

        const Field<Type>& cellData() const
        {
            return cellData_;
        }

        Field<Type>& cellData()
        {
            return cellData_;
        }

        const FieldField<Field, Type>& patchData() const
        {
            return patchData_;
        }

        FieldField<Field, Type>& patchData()
        {
            return patchData_;
        }

        label nUnset() const
        {
            return nUnset_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "patchDataWave.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
