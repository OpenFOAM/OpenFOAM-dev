/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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
    Foam::displacementSBRStressFvMotionSolver

Description
    Mesh motion solver for an fvMesh.  Based on solving the cell-centre
    solid-body rotation stress equations for the motion displacement.

SourceFiles
    displacementSBRStressFvMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef displacementSBRStressFvMotionSolver_H
#define displacementSBRStressFvMotionSolver_H

#include "displacementMotionSolver.H"
#include "fvMotionSolverCore.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward class declarations
class motionDiffusivity;

/*---------------------------------------------------------------------------*\
             Class displacementSBRStressFvMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class displacementSBRStressFvMotionSolver
:
    public displacementMotionSolver,
    public fvMotionSolverCore
{
    // Private data

        //- Cell-centre motion field
        mutable volVectorField cellDisplacement_;

        //- Diffusivity used to control the motion
        autoPtr<motionDiffusivity> diffusivityPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        displacementSBRStressFvMotionSolver
        (
            const displacementSBRStressFvMotionSolver&
        );

        //- Disallow default bitwise assignment
        void operator=(const displacementSBRStressFvMotionSolver&);


public:

    //- Runtime type information
    TypeName("displacementSBRStress");


    // Constructors

        //- Construct from polyMesh and IOdictionary
        displacementSBRStressFvMotionSolver
        (
            const polyMesh&,
            const IOdictionary&
        );


    //- Destructor
    ~displacementSBRStressFvMotionSolver();


    // Member Functions

        //- Return reference to the cell motion displacement field
        volVectorField& cellDisplacement()
        {
            return cellDisplacement_;
        }

        //- Return const reference to the cell motion displacement field
        const volVectorField& cellDisplacement() const
        {
            return cellDisplacement_;
        }

        //- Return diffusivity
        motionDiffusivity& diffusivity()
        {
            return diffusivityPtr_();
        }

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> curPoints() const;

        //- Solve for motion
        virtual void solve();

        //- Update topology
        virtual void updateMesh(const mapPolyMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
