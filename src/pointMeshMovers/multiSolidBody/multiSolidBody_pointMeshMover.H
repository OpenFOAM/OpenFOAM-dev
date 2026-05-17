/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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
    Foam::pointMeshMovers::multiSolidBody

Description
    Solid-body motion of the mesh parts specified by a run-time selectable
    motion functions.

SourceFiles
    multiSolidBody_pointMeshMover.C

\*---------------------------------------------------------------------------*/

#ifndef multiSolidBody_pointMeshMover_H
#define multiSolidBody_pointMeshMover_H

#include "displacementPoints0.H"
#include "solidBodyMotionFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace pointMeshMovers
{

/*---------------------------------------------------------------------------*\
                     Class multiSolidBody Declaration
\*---------------------------------------------------------------------------*/

class multiSolidBody
:
    public displacementPoints0
{
    // Private Data

        //- The motion control function
        PtrList<solidBodyMotionFunction> SBMFs_;

        //- Specified cellZones
        labelList zoneIndices_;

        //- Points to move per cellZone
        labelListList zonePoints_;

        //- Cache the current transforms for topoChange
        mutable List<septernion> transforms_;


    // Private Member Functions

        //- Re-calculate the zone point indices from the cell-zones
        void updateZonePointIndices();


public:

    //- Runtime type information
    TypeName("multiSolidBody");


    // Constructors

        //- Construct from mesh and dictionary
        multiSolidBody(const polyMesh&, const dictionary&);

        //- Disallow default bitwise copy construction
        multiSolidBody(const multiSolidBody&) = delete;


    //- Destructor
    ~multiSolidBody();


    // Member Functions

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> newPoints();

        //- This is a solid body motion
        virtual bool solidBodyMotion() const
        {
            return true;
        }

        //- Update local data for topology changes
        virtual void topoChange(const polyTopoChangeMap&);

        //- Update corresponding to the given distribution map
        virtual void distribute(const polyDistributionMap&);

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const multiSolidBody&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace pointMeshMovers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
