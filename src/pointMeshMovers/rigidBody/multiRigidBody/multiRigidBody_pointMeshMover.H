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
    Foam::pointMeshMovers::multiRigidBody

Description
    Abstract base-class for multiple rigid body mesh motion.

    Applies distance weighted SLERP interpolation between the septernions
    returned by each of the body movement functions.

SourceFiles
    multiRigidBody.C

\*---------------------------------------------------------------------------*/

#ifndef multiRigidBody_pointMeshMover_H
#define multiRigidBody_pointMeshMover_H

#include "displacementPoints0.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace pointMeshMovers
{

/*---------------------------------------------------------------------------*\
                    Class multiRigidBody Declaration
\*---------------------------------------------------------------------------*/

class multiRigidBody
:
    public displacementPoints0
{
protected:

    //- Class containing the patches and point motion weighting for each body
    class bodyMesh
    {
        const polyMesh& mesh_;

        //- Name of the body
        const word name_;

        //- List of mesh patch names associated with this body
        const wordReList patches_;

        //- Set of mesh patches associated with this body
        labelHashSet patchSet_;

        //- Names of the pointZones which move with the body
        wordReList pointZones_;

        //- pointZones which move with the body
        labelHashSet pointZoneSet_;

        //- Inner morphing distance (limit of solid-body region)
        const scalar di_;

        //- Outer morphing distance (limit of linear interpolation region)
        const scalar do_;

        //- Current interpolation weight
        //  (1 at patches and within di_, 0 at do_ and beyond)
        pointScalarField weight_;

        void calcWeights(const pointField& points0);


    public:

        //- ID of the body in the RBD::rigidBodyMotion
        label bodyIndex;

        friend class multiRigidBody;

        bodyMesh
        (
            const polyMesh& mesh,
            const word& name,
            const dictionary& dict
        );

        bodyMesh
        (
            const polyMesh& mesh,
            const word& name
        );

        template<class Type>
        inline Type weight(const Type& pDist) const;

        inline const word& name()
        {
            return name_;
        }

        inline const wordReList& patches() const
        {
            return patches_;
        }
    };


    // Protected Data

        //- List of the bodyMeshes containing the patches and point motion
        //  weighting for each body
        PtrList<bodyMesh> bodyMeshes_;


private:

    // Private Member Functions

        //- Return the list of weights for each body to pointi
        //   w is used as the work-space for the weights and returned
        List<scalar>& weights(const label pointi, List<scalar>& w) const;

        //- Return the list of transformation for each body from time 0
        virtual List<septernion> transforms0() const = 0;

        //- Move the bodies
        virtual void moveBodies() = 0;


public:

    //- Runtime type information
    TypeName("multiRigidBody");


    // Constructors

        //- Construct from polyMesh and dictionary
        multiRigidBody(const polyMesh&, const dictionary& dict);

        //- Disallow default bitwise copy construction
        multiRigidBody(const multiRigidBody&) = delete;


    //- Destructor
    ~multiRigidBody();


    // Member Functions

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> newPoints();

        //-  Update local data for topology changes
        virtual void topoChange(const polyTopoChangeMap&);

        //- Update from another mesh using the given map
        //  Resets points0 to the points of the new mesh
        virtual void mapMesh(const polyMeshMap&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const multiRigidBody&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace pointMeshMovers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
