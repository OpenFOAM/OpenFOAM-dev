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
    Foam::regionCoupledWallPolyPatch

Description
    regionCoupledWall

SourceFiles
    regionCoupledWallPolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef regionCoupledWallPolyPatch_H
#define regionCoupledWallPolyPatch_H

#include "wallPolyPatch.H"
#include "regionCoupledBase.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMesh;

/*---------------------------------------------------------------------------*\
                      Class regionCoupledWallPolyPatch Declaration
\*---------------------------------------------------------------------------*/

class regionCoupledWallPolyPatch
:
    public wallPolyPatch,
    public regionCoupledBase
{

protected:

        //- Initialise the calculation of the patch geometry
        virtual void initCalcGeometry(PstreamBuffers&);

        //- Calculate the patch geometry
        virtual void calcGeometry(PstreamBuffers&);

        //- Initialise the patches for moving points
        virtual void initMovePoints(PstreamBuffers&, const pointField&);

        //- Correct patches after moving points
        virtual void movePoints(PstreamBuffers&, const pointField&);

        //- Initialise the update of the patch topology
        virtual void initUpdateMesh(PstreamBuffers&);

        //- Update of the patch topology
        virtual void updateMesh(PstreamBuffers&);


public:

    //- Runtime type information
    TypeName("regionCoupledWall");


    // Constructors

        //- Construct from components
        regionCoupledWallPolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType
        );


        //- Construct from dictionary
        regionCoupledWallPolyPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType
        );

        //- Construct as copy, resetting the boundary mesh
        regionCoupledWallPolyPatch
        (
            const regionCoupledWallPolyPatch&,
            const polyBoundaryMesh&
        );

        //- Construct given the original patch and resetting the
        //  face list and boundary mesh information
        regionCoupledWallPolyPatch
        (
            const regionCoupledWallPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        );

        //- Construct given the original patch and a map
        regionCoupledWallPolyPatch
        (
            const regionCoupledWallPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const labelUList& mapAddressing,
            const label newStart
        );

        //- Construct and return a clone, resetting the boundary mesh
        virtual autoPtr<wallPolyPatch::polyPatch> clone
        (
            const polyBoundaryMesh& bm
        ) const
        {
            return autoPtr<wallPolyPatch::polyPatch>
            (
                new regionCoupledWallPolyPatch(*this, bm)
            );
        }

        //- Construct and return a clone, resetting the face list
        //  and boundary mesh
        virtual autoPtr<wallPolyPatch::polyPatch> clone
        (
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        ) const
        {
            return autoPtr<wallPolyPatch::polyPatch>
            (
                new regionCoupledWallPolyPatch
                (
                    *this,
                    bm,
                    index,
                    newSize,
                    newStart
                )
            );
        }

        //- Construct and return a clone, resetting the face list
        //  and boundary mesh
        virtual autoPtr<polyPatch> clone
        (
            const polyBoundaryMesh& bm,
            const label index,
            const labelUList& mapAddressing,
            const label newStart
        ) const
        {
            return autoPtr<polyPatch>
            (
                new regionCoupledWallPolyPatch
                (
                    *this,
                    bm,
                    index,
                    mapAddressing,
                    newStart
                )
            );
        }


    //- Destructor
    virtual ~regionCoupledWallPolyPatch();


    // Member Functions

        //- Return the type
        virtual const word& regionCoupleType() const
        {
            return type();
        }

        //- Write the polyPatch data as a dictionary
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
