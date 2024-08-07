/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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
    Foam::fvMeshMovers::multiValveEngine

Description
    A mesh mover using explicit node translation based on scaled distance
    functions per moving object. The mover supports any number of valves
    together with piston motion and following features:

    - Piston motion: Function1 of user-time, may be set to
      crankConnectingRodMotion for standard crank and connecting rod motion.

    - Valve motion: Function1, may be set to table if the valve lift date is
      provided in the form of a table.

    - Smooth mesh motion between a moving object and other patches.

    - linerPatches: the set of patches corresponding to the cylinder liner
      Used by createEngineZones

    - slidingPatches: a set of patches along which mesh is allowed
      to deform. For example, on the cylinder liner, it is desired to
      slide mesh nodes while piston is moving.

    - frozenZones: list of pointZones the points of which are frozen,
      i.e. do not move with respect to any moving object.

    - Run-time clearance estimation based on patch-to-patch distances printed.

    - Supports cellSet and cellZone definitions to restrict mesh motion.

    - Supports domains with nonConformalCoupling (NCC) interfaces,
      enabling e.g. nodes to slide along with the interface.

    - Closing the valve can be achieved by meshToMesh mapping onto a new
      grid with closed valve geometry at user given time.

    - Mesh motion can be controlled per moving object by setting:

        - patches: list of patches defining the object.

        - motion: a Function1 which returns the object position
          as a function of time.

        - movingZones: list of pointZones the points of which move with the
          object.

        - frozenZones: list of pointZones the points of which are frozen,
          i.e. do not move with respect to this moving object.

        - maxMotionDistance: a distance away from the moving object
          after nodes are not allowed to move. (Default inf.)

        - movingFrozenLayerThickness: thickness of layer in which points move
          with the moving object. (Default 0)

        - staticFrozenLayerThickness: thickness of layer in which points
          are fixed with respect to static patches (e.g. walls). (Default 0)

        - cosineScaling: a switch whether nodal translation is weighted by
          its distance from the moving object. The objective is to yield less
          deformation near the moving object and sustain e.g. boundary layer.
          (Default no, i.e. linear weighting)

        - travelInterval: part of the stroke travelled after
          which the cached motion scaling weights are recalculated

        For valve object only:

            - minLift: a minimum valve lift value after considered closed.


    Some of the above parameters are highlighted in a given schematic
    piston-valve configuration w.r.t entries used to control piston motion.
    Furthermore, an example dictionary entries are provided below.
    \verbatim
                      |             |         |             |
                      |             |         |             |
                      |             |    S    |             |
                      |             |    T    |             |
                      |             |    E    |             |
                      |             |    M    |             |
                     /              |         |              \
                    /               |         |               \
                   /                |         |                \
     _____________/                 |         |                 \_____________
    |        :                      |         |                      :        |
    |        :      /```````````````           ```````````````\      :        |
    |        :     /                VALVE HEAD                 \     :        |
    | L      :    /_____________________________________________\    :        |
    | I      :                         /\                            :        |
    | N      :                         || staticFrozenLayerThickness :        |
    | E      : NCC (optional)          \/ (w.r.t. piston motion)     :        |
    | R      :                      ``````````                       :        |
    |        :                                                       :        |
    |        :                                                       :        |
    |........:.......................................................:........|
    |        :                         /\                            :        |
    |        :                         || movingFrozenLayerThickness :        |
    |________:_________________________\/____________________________:________|
                                       PISTON
    \endverbatim

    \verbatim
    mover
    {
        type                multiValveEngine;
        libs                ("libfvMeshMoversMultiValveEngine.so");

        frozenZones         (frozenZone1 frozenZone2);

        slidingPatches
        (
            liner
            valveStem
            "nonCouple.*"
        );

        linerPatches        (liner);

        piston
        {
            patches             (piston);
            axis                (0 0 1);

            motion
            {
                type                crankConnectingRodMotion;

                conRodLength        1e3;
                stroke              1.0;
            }

            // Move the points in the piston bowl with the piston
            movingZones         (pistonBowl);

            // Freeze the points in the cylinder head
            frozenZones         (cylinderHead);

            // Optional
            maxMotionDistance    1e30;
            movingFrozenLayerThickness  0;
            staticFrozenLayerThickness  0;

            travelInterval      0.1;

            cosineScaling       yes;
        }

        valves
        {
            iv
            {
                patches     (valveHead);
                axis        (0 0 1);

                // Optional
                maxMotionDistance   1e30;
                movingFrozenLayerThickness  0;
                staticFrozenLayerThickness  0;

                travelInterval      0.01;

                cosineScaling       yes;

                minLift     0.001;

                motion
                {
                    type    table;
                    values
                    (
                        (0      0)
                        (480    0.1)
                        (720    0)
                    );
                    // For multi-cycle simulations, use repeat
                    outOfBounds     repeat;
                    interpolationScheme linear;
                }
            }
        }
    }
    \endverbatim

    Note:
      The implementation utilises pointDist objects for distance computation,
      resulting distance fields do not propagate through NCC interfaces.  Hence,
      there should be no horizontal NCC interface separating piston from
      cylinder head as it would result in potentially ill defined mesh
      deformation. Due to same feature, in a schematic case setup above, valve
      motion affects only cells between NCC patches even though no cellSet is
      explicitly defined.

SourceFiles
    multiValveEngine.C

\*---------------------------------------------------------------------------*/

#ifndef multiValveEngine_H
#define multiValveEngine_H

#include "fvMeshMover.H"
#include "Function1.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class pointDist;

namespace fvMeshMovers
{

/*---------------------------------------------------------------------------*\
                      Class multiValveEngine Declaration
\*---------------------------------------------------------------------------*/

class multiValveEngine
:
    public fvMeshMover
{
public:

    class movingObject
    {
        friend class multiValveEngine;

    protected:

        // Protected member data

            //- Reference to engine mesh mover
            const multiValveEngine& meshMover_;


    public:

            //- Name of the object
            const word name;

            //- Axis
            const vector axis;


    protected:

        // Protected member data

            //- Piston motion function
            autoPtr<Function1<scalar>> motion_;

            //- Object patchNames
            wordReList patchNames_;

            //- Object patchSet
            labelHashSet patchSet_;

            //- Patches which must not deform during mesh motion
            labelHashSet staticPatchSet_;

            //- Distance from object beyond which the mesh does not deform
            scalar maxMotionDistance_;

            //- Distance from object within mesh nodes move
            //  along with the patch
            scalar movingFrozenLayerThickness_;

            //- Distance from static wall patch at domain perimeter
            //  within mesh nodes do not move
            scalar staticFrozenLayerThickness_;

            wordReList movingPointZones_;

            wordReList frozenPointZones_;

            //- Points to move when cell zone is supplied
            labelList pointIDs_;

            //- Interpolation scale (1 at moving patches, 0 at far-field)
            pointScalarField scale_;

            //- Switch to use cosine-based scaling
            //  for the boundary point motion
            Switch cosine_;

            //- Update the scale_ field when the travel exceeds travelInterval
            const scalar travelInterval_;

            //- Count of the scale_ field updates
            mutable label executionCount_;

            //- Position at last scale update
            mutable scalar position0_;


        // Protected Member Functions

            //- Scale the mesh point deformation with distance functions
            //  w.r.t. moving patches, static patches and a frozen zones.
            void calcScale
            (
                const pointMesh& pMesh,
                const scalarField& pDistMoving,
                const scalarField& pDistStatic,
                scalar dMoving,
                scalar dMax,
                scalar dStatic
            );

            void transformPoints
            (
                pointField& newPoints,
                const vector& translationVector
            );

            //- Generate staticPatchSet_ based on patch entries
            void createStaticPatchSet();

            //- Patch-set construction
            void initPatchSets();

            labelHashSet movingPointZones() const;

            labelHashSet staticPointZones() const;


    public:

        //- Object patchSet
        const labelHashSet& patchSet;

        // Constructors

            //- Construct from dictionary
            movingObject
            (
                const word& name,
                const multiValveEngine& engine,
                const dictionary& dict
            );

        //- Destructor
        virtual ~movingObject()
        {}


        // Member Functions

            //- Update from another mesh using the given map
            virtual void mapMesh(const polyMeshMap&);
    };


    class pistonObject
    :
        public movingObject
    {
        // Private member data

            //- Bore
            scalar bore_;

            //- Clearance
            scalar clearance_;


        // Private member functions

            //- Return the bore
            //  calculated from the lateral extent of the piston
            scalar bore() const;

            //- Calculate clearance estimate based on minimum distance
            //  between piston and any other patch excluding slidingPatches.
            void correctClearance(pointDist&);


    public:

        //- Name of the piston bowl pointZone
        static word pistonBowlName;

        // Constructors

            //- Construct from dictionary
            pistonObject
            (
                const word& name,
                const multiValveEngine& engine,
                const dictionary& dict
            );


        // Member Functions

            //- Return the piston position for the given CA theta
            scalar position(const scalar theta) const;

            //- Return the current piston position
            scalar position() const;

            //- Return piston displacement for current time-step
            scalar displacement() const;

            //- Return piston position for current time-step
            scalar speed() const;

            //- Return clearance estimate
            scalar clearance() const;

            //- update points due to piston motion
            void updatePoints(pointField&);

            //- Update from another mesh using the given map
            virtual void mapMesh(const polyMeshMap&);
    };


    class valveObject
    :
        public movingObject
    {
        // Private member data

            //- Minimum valve lift.
            //  On this lift the valve is considered closed
            const scalar minLift_;


        // Private member functions


    public:

        // Constructors

            //- Construct from dictionary
            valveObject
            (
                const word& name,
                const multiValveEngine& engine,
                const dictionary& dict
            );

            //- Dummy clone function for PtrList
            autoPtr<valveObject> clone() const
            {
                NotImplemented;
                return autoPtr<valveObject>(nullptr);
            }


        // Member Functions

            // Valve position and velocity

                //- Return valve position for the given time
                scalar lift(const scalar theta) const;

                //- Return current valve position
                scalar lift() const;

                //- Return current valve speed
                scalar speed() const;

                //- Return valve displacement for current time-step
                scalar displacement() const;

                //- Is the valve open?
                bool isOpen() const;


            //- update points due to valve motion
            void updatePoints(pointField&);
    };


    class valveList
    :
        public PtrList<valveObject>
    {
    public:

        // Constructors

            //- Construct from Istream
            valveList
            (
                const multiValveEngine& engine,
                const dictionary& dict
            );
    };


    friend class movingObject;
    friend class pistonObject;
    friend class valveObject;


private:

    // Private member data

        //- User-defined liner patches
        wordReList linerPatches_;

        //- User-defined liner patch set
        labelHashSet linerPatchSet_;

        //- User-defined patches which the mesh can slide along
        wordReList slidingPatches_;

        //- User-defined patch set which the mesh can slide along
        labelHashSet slidingPatchSet_;

        //- Piston object
        pistonObject piston_;

        //- Container for all valves
        valveList valves_;

        //- Static patch set
        labelHashSet staticPatchSet_;

        wordReList frozenPointZones_;


    // Private member functions

        //- Lookup and return the liner patch set
        labelHashSet findLinerPatchSet() const;

        //- Lookup and return the sliding patch set
        labelHashSet findSlidingPatchSet();

        //- Find and return the static patch set
        labelHashSet findStaticPatchSet();


public:

    //- Runtime type information
    TypeName("multiValveEngine");

    //- Name of the cylinder head pointZone
    static word cylinderHeadName;

    // Constant access to member data

        //- User-defined liner patches
        const labelHashSet& linerPatchSet;

        //- User-defined patches which the mesh can slide along
        const labelHashSet& slidingPatchSet;

        //- Piston object
        const pistonObject& piston;

        //- Container for all valves
        const valveList& valves;

        //- Static patch set
        const labelHashSet& staticPatchSet;


    // Constructors

        //- Construct from fvMesh
        multiValveEngine(fvMesh& mesh, const dictionary& dict);

        //- Disallow default bitwise copy construction
        multiValveEngine(const multiValveEngine&) = delete;


    //- Destructor
    ~multiValveEngine();


    // Member Functions

        //- Return current user-time, CAD, s or ...
        scalar userTime() const;

        //- Return current user-time-step, CAD, s, ...
        scalar userDeltaT() const;

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();

        //- Update corresponding to the given map
        virtual void topoChange(const polyTopoChangeMap&);

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&);

        //- Update corresponding to the given distribution map
        virtual void distribute(const polyDistributionMap&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const multiValveEngine&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace fvMeshMovers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
