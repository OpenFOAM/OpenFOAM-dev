/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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
    Foam::componentDisplacementMotionSolver

Description
    Virtual base class for displacement motion solver

    The boundary displacement is set as a boundary condition
    on the pointDisplacementX pointScalarField.

SourceFiles
    componentDisplacementMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef componentDisplacementMotionSolver_H
#define componentDisplacementMotionSolver_H

#include "motionSolver.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                   Class componentDisplacementMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class componentDisplacementMotionSolver
:
    public motionSolver
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

        //- Disallow default bitwise copy construct
        componentDisplacementMotionSolver
        (
            const componentDisplacementMotionSolver&
        );

        //- Disallow default bitwise assignment
        void operator=(const componentDisplacementMotionSolver&);


public:

    //- Runtime type information
    TypeName("componentDisplacementMotionSolver");


    // Constructors

        //- Construct from polyMesh and dictionary and type
        componentDisplacementMotionSolver
        (
            const polyMesh&,
            const IOdictionary&,
            const word& type
        );


    //- Destructor
    virtual ~componentDisplacementMotionSolver();


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
        virtual void updateMesh(const mapPolyMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
