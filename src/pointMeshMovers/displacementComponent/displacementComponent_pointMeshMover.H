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
    Foam::pointMeshMovers::displacementComponent

Description
    Abstract base class for component displacement pointMesh movers

    The boundary displacement is set as a boundary condition
    on the pointDisplacementX pointScalarField.

SourceFiles
    displacementComponent.C

\*---------------------------------------------------------------------------*/

#ifndef displacementComponent_pointMeshMover_H
#define displacementComponent_pointMeshMover_H

#include "pointMeshMover.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace pointMeshMovers
{

/*---------------------------------------------------------------------------*\
              Class displacementComponent Declaration
\*---------------------------------------------------------------------------*/

class displacementComponent
:
    public pointMeshMover
{
protected:

    // Protected data

        //- The component name to solve for
        word cmptName_;

        //- The component to solve for
        direction cmpt_;

        //- Reference point field for this component
        scalarField points0_;

        //- Point motion field
        mutable pointScalarField pointDisplacement_;


private:

    // Private Member Functions

        //- Return the component corresponding to the given component name
        direction cmpt(const word& cmptName) const;


public:

    //- Runtime type information
    TypeName("displacementComponent");


    // Constructors

        //- Construct from polyMesh and dictionary and type
        displacementComponent
        (
            const polyMesh&,
            const dictionary&,
            const word& type
        );

        //- Disallow default bitwise copy construction
        displacementComponent
        (
            const displacementComponent&
        ) = delete;


    //- Destructor
    virtual ~displacementComponent();


    // Member Functions

        //- Return reference to the reference field
        scalarField& points0()
        {
            return points0_;
        }

        //- Return reference to the reference field
        const scalarField& points0() const
        {
            return points0_;
        }

        //- Update local data for geometry changes
        virtual void movePoints(const pointField&);

        //-  Update local data for topology changes
        virtual void topoChange(const polyTopoChangeMap&);

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&);

        //- Update corresponding to the given distribution map
        virtual void distribute(const polyDistributionMap&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const displacementComponent&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace pointMeshMovers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
