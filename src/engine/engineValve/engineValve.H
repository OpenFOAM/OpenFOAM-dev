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
    Foam::engineValve

Description
    Foam::engineValve

SourceFiles
    engineValve.C

\*---------------------------------------------------------------------------*/

#ifndef engineValve_H
#define engineValve_H

#include "word.H"
#include "coordinateSystem.H"
#include "polyPatchID.H"
#include "Table.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;
class engineTime;

/*---------------------------------------------------------------------------*\
                         Class engineValve Declaration
\*---------------------------------------------------------------------------*/

class engineValve
{
    // Private Data

        //- Name of valve
        word name_;

        //- Reference to engine mesh
        const polyMesh& mesh_;

        //- Reference to engine database
        const engineTime& engineDB_;

        //- Coordinate system
        autoPtr<coordinateSystem> csPtr_;


        // Patch and zone names

            //- Valve bottom patch
            polyPatchID bottomPatch_;

            //- Valve poppet patch
            polyPatchID poppetPatch_;

            //- Valve stem patch
            polyPatchID stemPatch_;

            //- Valve curtain manifold patch
            polyPatchID curtainInPortPatch_;

            //- Valve curtain cylinder patch
            polyPatchID curtainInCylinderPatch_;

            //- Valve detach in cylinder patch
            polyPatchID detachInCylinderPatch_;

            //- Valve detach in port patch
            polyPatchID detachInPortPatch_;

            //- Faces to detach
            labelList detachFaces_;


        // Valve lift data

            //- Valve lift profile
            Function1s::Table<scalar> liftProfile_;

            //- Lift curve start angle
            scalar liftProfileStart_;

            //- Lift curve end angle
            scalar liftProfileEnd_;

            //- Minimum valve lift.  On this lift the valve is considered closed
            const scalar minLift_;


            // Valve layering data

                //- Min top layer thickness
                const scalar minTopLayer_;

                //- Max top layer thickness
                const scalar maxTopLayer_;

                //- Min bottom layer thickness
                const scalar minBottomLayer_;

                //- Max bottom layer thickness
                const scalar maxBottomLayer_;


            //- Valve diameter
            const scalar diameter_;


    // Private Member Functions

        //- Adjust crank angle to drop within the limits of the lift profile
        scalar adjustCrankAngle(const scalar theta) const;


public:

    // Constructors

        //- Construct from components
        engineValve
        (
            const word& name,
            const polyMesh& mesh,
            const autoPtr<coordinateSystem>& valveCS,
            const word& bottomPatchName,
            const word& poppetPatchName,
            const word& stemPatchName,
            const word& curtainInPortPatchName,
            const word& curtainInCylinderPatchName,
            const word& detachInCylinderPatchName,
            const word& detachInPortPatchName,
            const labelList& detachFaces,
            const Function1s::Table<scalar>& liftProfile,
            const scalar minLift,
            const scalar minTopLayer,
            const scalar maxTopLayer,
            const scalar minBottomLayer,
            const scalar maxBottomLayer,
            const scalar diameter

        );

        //- Construct from dictionary
        engineValve
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        engineValve(const engineValve&) = delete;


    // Member Functions

        //- Return name
        const word& name() const
        {
            return name_;
        }

        //- Return coordinate system
        const coordinateSystem& cs() const
        {
            return csPtr_();
        }

        //- Return lift profile
        const Function1s::Table<scalar>& liftProfile() const
        {
            return liftProfile_;
        }

        //- Return valve diameter
        scalar diameter() const
        {
            return diameter_;
        }


        // Valve patches

            //- Return ID of bottom patch
            const polyPatchID& bottomPatchID() const
            {
                return bottomPatch_;
            }

            //- Return ID of poppet patch
            const polyPatchID& poppetPatchID() const
            {
                return poppetPatch_;
            }

            //- Return ID of stem patch
            const polyPatchID& stemPatchID() const
            {
                return stemPatch_;
            }

            //- Return ID of curtain in cylinder patch
            const polyPatchID& curtainInCylinderPatchID() const
            {
                return curtainInCylinderPatch_;
            }

            //- Return ID of curtain in port patch
            const polyPatchID& curtainInPortPatchID() const
            {
                return curtainInPortPatch_;
            }


            //- Return ID of detach in cylinder patch
            const polyPatchID& detachInCylinderPatchID() const
            {
                return detachInCylinderPatch_;
            }

            //- Return ID of detach in port patch
            const polyPatchID& detachInPortPatchID() const
            {
                return detachInPortPatch_;
            }

            //- Return face labels of detach curtain
            const labelList& detachFaces() const
            {
                return detachFaces_;
            }


        // Valve layering thickness

            scalar minTopLayer() const
            {
                return minTopLayer_;
            }

            scalar maxTopLayer() const
            {
                return maxTopLayer_;
            }

            scalar minBottomLayer() const
            {
                return minBottomLayer_;
            }

            scalar maxBottomLayer() const
            {
                return maxBottomLayer_;
            }


        // Valve position and velocity

            //- Return valve lift given crank angle in degrees
            scalar lift(const scalar theta) const;

            //- Is the valve open?
            bool isOpen() const;

            //- Return current lift
            scalar curLift() const;

            //- Return valve velocity for current time-step
            scalar curVelocity() const;

            //- Return list of active patch labels for the valve head
            //  (stem is excluded)
            labelList movingPatchIDs() const;


        //- Write dictionary
        void writeDict(Ostream&) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const engineValve&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
