/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2026 OpenFOAM Foundation
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
    Foam::pointMeshMovers::displacement

Description
    Abstract base class for displacement pointMesh movers.

    The boundary displacement is set as a boundary condition
    on the pointDisplacement pointVectorField.

SourceFiles
    displacement.C

\*---------------------------------------------------------------------------*/

#ifndef displacement_pointMeshMover_H
#define displacement_pointMeshMover_H

#include "displacementPoints0.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace pointMeshMovers
{

/*---------------------------------------------------------------------------*\
                  Class displacement Declaration
\*---------------------------------------------------------------------------*/

class displacement
:
    public displacementPoints0
{
protected:

    // Protected data

        //- Point motion field
        mutable pointVectorField pointDisplacement_;

        //- Return the points corresponding to the current displacement
        tmp<pointField> points() const;


public:

    //- Runtime type information
    TypeName("displacement");


    // Constructors

        //- Construct from mesh and dictionary
        displacement
        (
            const polyMesh&,
            const dictionary&,
            const word& type
        );

        //- Disallow default bitwise copy construction
        displacement(const displacement&) = delete;


    //- Destructor
    virtual ~displacement();


    // Member Functions

        //- Return reference to the point motion displacement field
        pointVectorField& pointDisplacement()
        {
            return pointDisplacement_;
        }

        //- Return const reference to the point motion displacement field
        const pointVectorField& pointDisplacement() const
        {
            return pointDisplacement_;
        }

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> newPoints() = 0;

        //- Update from another mesh using the given map (not implemented)
        virtual void mapMesh(const polyMeshMap&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const displacement&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace pointMeshMovers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
