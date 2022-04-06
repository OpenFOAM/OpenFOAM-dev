/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2022 OpenFOAM Foundation
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
    Foam::velocityMotionSolver

Description
    Virtual base class for velocity motion solver

    The boundary displacement is set as a boundary condition
    on the pointMotionU pointVectorField.

SourceFiles
    velocityMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef velocityMotionSolver_H
#define velocityMotionSolver_H

#include "motionSolver.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyTopoChangeMap;

/*---------------------------------------------------------------------------*\
                    Class velocityMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class velocityMotionSolver
:
    public motionSolver
{
protected:

    // Protected data

        //- Point motion field
        mutable pointVectorField pointMotionU_;


public:

    //- Runtime type information
    TypeName("velocityMotionSolver");


    // Constructors

        //- Construct from mesh and dictionary
        velocityMotionSolver
        (
            const word& name,
            const polyMesh&,
            const dictionary&,
            const word& type
        );

        //- Disallow default bitwise copy construction
        velocityMotionSolver
        (
            const velocityMotionSolver&
        );


    //- Destructor
    virtual ~velocityMotionSolver();


    // Member Functions

        //- Return reference to the point motion velocity field
        pointVectorField& pointMotionU()
        {
            return pointMotionU_;
        }

        //- Return const reference to the point motion velocity field
        const pointVectorField& pointMotionU() const
        {
            return pointMotionU_;
        }

        //- Update local data for geometry changes
        virtual void movePoints(const pointField&);

        //- Update corresponding to the given map
        virtual void topoChange(const polyTopoChangeMap&);

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&);

        //- Update corresponding to the given distribution map
        virtual void distribute(const polyDistributionMap&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const velocityMotionSolver&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
