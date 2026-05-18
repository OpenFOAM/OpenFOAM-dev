/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2026 OpenFOAM Foundation
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
    Foam::pointMeshMovers::functionalRigidBody

Description
    Multiple rigid body mesh motion in which movement is specified as a
    solidBodyMotionFunction for each body and optionally for the exterior of the
    region.

    Applies distance weighted SLERP interpolation between the septernions
    returned by each of the body movement functions.

    The mesh patches for each body must be specified and optionally pointZones
    that move with each body.  The motion of the mesh around each body may be
    constrained by the optional \c innerDistance and \c outerDistance, the mesh
    within the \c innerDistance moves with the body without distortion and the
    mesh beyond the \c outerDistance does not move with body but may move by the
    proximity of other bodies or by the optionally specified exterior motion.

    Example of the motion of a single body in a stationary region:
    \verbatim
    mover
    {
        type            pointMeshMover;

        libs            ("libpointMeshMovers.so");

        pointMeshMover
        {
            type            functionalRigidBody;

            exterior
            {
                innerDistance   0.05;
                outerDistance   0.35;
            }

            patches         (floatingObject);
            innerDistance   0.05;
            outerDistance   0.35;

            function
            {
                type            oscillatingRotatingMotion;
                origin          (0.5 0.45 0.35);
                amplitude       (0 10 0);
                omega           2;
            }
        }
    }
    \endverbatim

    Example of the motion of a two bodies in a moving region all of which
    rotate around the same axis with an oscillatory manner:
    \verbatim
    mover
    {
        type            pointMeshMover;

        libs            ("libpointMeshMovers.so");

        pointMeshMover
        {
            type            functionalRigidBody;

            exterior
            {
                innerDistance   0.05;
                outerDistance   0.35;

                function
                {
                    type            oscillatingRotatingMotion;
                    origin          (0.5 0.45 0.35);
                    amplitude       (0 1 0);
                    omega           1;
                }
            }

            bodies
            {
                floatingObject1
                {
                    patches         (floatingObject1);
                    innerDistance   0.05;
                    outerDistance   0.3;

                    function
                    {
                        type            oscillatingRotatingMotion;
                        origin          (0.5 0.45 0.35);
                        amplitude       (0 5 0);
                        omega           2;
                    }
                }

                floatingObject2
                {
                    patches         (floatingObject2);
                    innerDistance   0.05;
                    outerDistance   0.3;

                    function
                    {
                        type            oscillatingRotatingMotion;
                        origin          (0.5 0.45 0.35);
                        amplitude       (0 5 0);
                        omega           -2;
                    }
                }
            }
        }
    }
    \endverbatim

See also
    Foam::multiRigidBodyMeshMotion

SourceFiles
    functionalRigidBody.C

\*---------------------------------------------------------------------------*/

#ifndef functionalRigidBody_pointMeshMover_H
#define functionalRigidBody_pointMeshMover_H

#include "multiRigidBody_pointMeshMover.H"
#include "solidBodyMotionFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace pointMeshMovers
{

/*---------------------------------------------------------------------------*\
                Class functionalRigidBody Declaration
\*---------------------------------------------------------------------------*/

class functionalRigidBody
:
    public multiRigidBody
{
    // Private Data

        //- Body motion control functions
        PtrList<solidBodyMotionFunction> SBMFs_;


    // Private Member Functions

        //- Return the list of transformation for each body from time 0
        virtual List<septernion> transforms0() const;

        //- Move the bodies
        virtual void moveBodies();


public:

    //- Runtime type information
    TypeName("functionalRigidBody");


    // Constructors

        //- Construct from polyMesh and dictionary
        functionalRigidBody(const polyMesh&, const dictionary& dict);

        //- Disallow default bitwise copy construction
        functionalRigidBody(const functionalRigidBody&) = delete;


    //- Destructor
    ~functionalRigidBody();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const functionalRigidBody&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace pointMeshMovers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
