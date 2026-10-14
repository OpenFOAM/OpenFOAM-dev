/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2026 OpenFOAM Foundation
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
    Foam::pointMeshMover::displacementPoints0

Description
    Abstract base class for displacement pointMesh movers.

SourceFiles
    displacementPoints0.C

\*---------------------------------------------------------------------------*/

#ifndef displacementPoints0_H
#define displacementPoints0_H

#include "pointMeshMover.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace pointMeshMovers
{

/*---------------------------------------------------------------------------*\
                      Class displacementPoints0 Declaration
\*---------------------------------------------------------------------------*/

class displacementPoints0
:
    public pointMeshMover
{

protected:

    // Protected data

        //- Starting points
        pointVectorField points0_;


public:

    //- Runtime type information
    TypeName("displacementPoints0");


    // Constructors

        //- Construct from mesh and dictionary
        displacementPoints0
        (
            const polyMesh&,
            const dictionary&,
            const word& type
        );

        //- Disallow default bitwise copy construction
        displacementPoints0(const displacementPoints0&) = delete;


    //- Destructor
    virtual ~displacementPoints0();


    // Member Functions

        static pointVectorField readPoints0(const polyMesh& mesh);

        //- Return reference to the reference field
        pointField& points0()
        {
            return points0_;
        }

        //- Return reference to the reference field
        const pointField& points0() const
        {
            return points0_;
        }

        //- Update local data for geometry changes
        virtual void movePoints(const pointField&);

        //-  Update local data for topology changes
        virtual void topoChange(const polyTopoChangeMap&);

        //- Update from another mesh using the given map
        //  Resets points0 to the points of the new mesh
        virtual void mapMesh(const polyMeshMap&) = 0;

        //- Update corresponding to the given distribution map
        virtual void distribute(const polyDistributionMap&);

        //- Write points0 if the mesh topology changed
        virtual bool write() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const displacementPoints0&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace pointMeshMovers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
