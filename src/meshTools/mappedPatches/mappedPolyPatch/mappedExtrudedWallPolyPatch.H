/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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
    Foam::mappedExtrudedWallPolyPatch

Description
    Wall patch which holds a mapping engine to map values from another patch,
    and which offsets the geometry to account for the extrusion thickness

SourceFiles
    mappedExtrudedWallPolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef mappedExtrudedWallPolyPatch_H
#define mappedExtrudedWallPolyPatch_H

#include "mappedWallPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMesh;

/*---------------------------------------------------------------------------*\
                     Class mappedExtrudedWallPolyPatch Declaration
\*---------------------------------------------------------------------------*/

class mappedExtrudedWallPolyPatch
:
    public mappedWallPolyPatch
{
private:

    // Private Member Data

        //- The bottom patch on the other side of the extruded region. Empty if
        //  this is not the extruded region.
        word bottomPatch_;

        //- The face areas on the bottom patch
        mutable autoPtr<vectorField> bottomFaceAreasPtr_;

        //- The face centres on the bottom patch
        mutable autoPtr<pointField> bottomFaceCentresPtr_;

        //- The local points on the bottom patch
        mutable autoPtr<pointField> bottomLocalPointsPtr_;


protected:

    // Protected Member Functions

        //- Initialise the calculation of the patch geometry
        virtual void initCalcGeometry(PstreamBuffers&);

        //- Calculate the patch geometry
        virtual void calcGeometry(PstreamBuffers&);

        //- Initialise the patches for moving points
        virtual void initMovePoints(PstreamBuffers&, const pointField&);

        //- Correct patches after moving points
        virtual void movePoints(PstreamBuffers&, const pointField&);

        //- Initialise the update of the patch topology
        virtual void initTopoChange(PstreamBuffers&);

        //- Update of the patch topology
        virtual void topoChange(PstreamBuffers&);

        //- Get the face-areas for the patch
        virtual tmp<vectorField> patchFaceAreas() const;

        //- Get the face-centres for the patch
        virtual tmp<pointField> patchFaceCentres() const;

        //- Get the local points for the patch
        virtual tmp<pointField> patchLocalPoints() const;


public:

    //- Runtime type information
    TypeName("mappedExtrudedWall");


    // Constructors

        //- Construct from components
        mappedExtrudedWallPolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType
        );

        //- Construct from dictionary
        mappedExtrudedWallPolyPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType
        );

        //- Construct as copy, resetting the boundary mesh
        mappedExtrudedWallPolyPatch
        (
            const mappedExtrudedWallPolyPatch&,
            const polyBoundaryMesh&
        );

        //- Construct given the original patch and resetting the
        //  face list and boundary mesh information
        mappedExtrudedWallPolyPatch
        (
            const mappedExtrudedWallPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        );

        //- Construct and return a clone, resetting the boundary mesh
        virtual autoPtr<polyPatch> clone(const polyBoundaryMesh& bm) const
        {
            return autoPtr<polyPatch>
            (
                new mappedExtrudedWallPolyPatch(*this, bm)
            );
        }

        //- Construct and return a clone, resetting the face list
        //  and boundary mesh
        virtual autoPtr<polyPatch> clone
        (
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        ) const
        {
            return autoPtr<polyPatch>
            (
                new mappedExtrudedWallPolyPatch
                (
                    *this,
                    bm,
                    index,
                    newSize,
                    newStart
                )
            );
        }


    //- Destructor
    virtual ~mappedExtrudedWallPolyPatch();


    // Member Functions

        //- Write the polyPatch data as a dictionary
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
