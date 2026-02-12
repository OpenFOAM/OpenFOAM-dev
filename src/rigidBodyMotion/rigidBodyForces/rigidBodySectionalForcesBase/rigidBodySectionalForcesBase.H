/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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
    Foam::functionObjects::rigidBodySectionalForcesBase

Description
    Base class for rigid-body sectional forces function objects

SourceFiles
    rigidBodySectionalForcesBase.C

\*---------------------------------------------------------------------------*/

#ifndef rigidBodySectionalForcesBase_functionObject_H
#define rigidBodySectionalForcesBase_functionObject_H

#include "sectionalForcesBase.H"
#include "coordSet.H"
#include "rigidBodyMotion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class rigidBodySectionalForcesBase Declaration
\*---------------------------------------------------------------------------*/

class rigidBodySectionalForcesBase
:
    public sectionalForcesBase
{
    // Private Data

        //- Name of the rigid body
        word body_;

        //- The axis on which to plot the results
        coordSet::axisType axis_;

        //- Local origin of the plot
        point localOrigin_;

        //- The axis on which to plot the results
        label axisi_;


protected:

    // Protected Member Functions

        //- Return the motion system
        const RBD::rigidBodyMotion& motion() const;

        //- Return the index of the body
        label bodyIndex() const;

        //- Return the body
        const RBD::rigidBody& body() const;

        //- Return the axis on which to plot the result
        label axisi() const;

        //- Return the name of the axis on which to plot the result
        const word& axisName() const;

        //- Return the vector axis on which to plot the result. Note that this
        //  is in the local coordinate system of the body. The global axis is
        //  given by the 'normal' method.
        vector axis() const;

        //- Return the local origin of the plot
        const point& localOrigin() const;

        //- Return the normal to the section planes
        virtual vector normal() const;

        //- Return the origin of the plot
        virtual point origin() const;

        //- Calculate the forces and moments at the cuts
        void addFluid
        (
            vectorField& force,
            vectorField& moment
        ) const;

        //- Calculate the forces and moments at the cuts
        void addBody
        (
            vectorField& force,
            vectorField& moment
        ) const;


public:

    //- Runtime type information
    TypeName("rigidBodySectionalForcesBase");


    // Constructors

        //- Construct from Time and dictionary
        rigidBodySectionalForcesBase
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        rigidBodySectionalForcesBase
        (
            const rigidBodySectionalForcesBase&
        ) = delete;


    //- Destructor
    virtual ~rigidBodySectionalForcesBase();


    // Member Functions

        //- Read the rigidBodySectionalForcesBase data
        virtual bool read(const dictionary&);

        //- Update for mesh point-motion
        virtual void movePoints(const polyMesh&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const rigidBodySectionalForcesBase&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
