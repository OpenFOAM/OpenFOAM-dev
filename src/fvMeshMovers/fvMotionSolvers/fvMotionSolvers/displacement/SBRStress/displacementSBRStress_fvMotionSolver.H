/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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
    Foam::fvMotionSolvers::displacementSBRStress

Description
    Mesh motion solver for an fvMesh.  Based on solving the cell-centre
    solid-body rotation stress equations for the motion displacement.

SourceFiles
    displacementSBRStress_fvMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef displacementSBRStress_fvMotionSolver_H
#define displacementSBRStress_fvMotionSolver_H

#include "displacement_pointMeshMover.H"
#include "fvMotionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward class declarations
class motionDiffusivity;

namespace fvMotionSolvers
{

/*---------------------------------------------------------------------------*\
                 Class displacementSBRStress Declaration
\*---------------------------------------------------------------------------*/

class displacementSBRStress
:
    public fvMotionSolver,
    public pointMeshMovers::displacement
{
    // Private Data

        //- Cell-centre motion field
        mutable volVectorField cellDisplacement_;

        //- Diffusivity type and data stream
        ITstream diffusivityType_;

        //- Diffusivity used to control the motion
        autoPtr<motionDiffusivity> diffusivityPtr_;


public:

    //- Runtime type information
    TypeName("displacementSBRStress");


    // Constructors

        //- Construct from fvMesh and dictionary
        displacementSBRStress(fvMesh&, const dictionary&);

        //- Construct from polyMesh and dictionary
        displacementSBRStress(const polyMesh&, const dictionary&);

        //- Disallow default bitwise copy construction
        displacementSBRStress(const displacementSBRStress&) = delete;


    //- Destructor
    ~displacementSBRStress();


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
        virtual tmp<pointField> newPoints();

        //- Update corresponding to the given map
        virtual void topoChange(const polyTopoChangeMap&);

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const displacementSBRStress&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvMotionSolvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
