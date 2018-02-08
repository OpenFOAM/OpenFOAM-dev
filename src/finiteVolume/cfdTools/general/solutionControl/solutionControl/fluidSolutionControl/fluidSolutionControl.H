/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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
    Foam::fluidSolutionControl

Description
    Fluid solution control class. Holds flags specific to the solution of
    fluids.

SourceFiles
    fluidSolutionControlI.H
    fluidSolutionControl.C

\*---------------------------------------------------------------------------*/

#ifndef fluidSolutionControl_H
#define fluidSolutionControl_H

#include "nonOrthogonalSolutionControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class fluidSolutionControl Declaration
\*---------------------------------------------------------------------------*/

class fluidSolutionControl
:
    public nonOrthogonalSolutionControl
{
protected:

    // Protected data

        //- Flag to indicate to solve for the flow
        bool solveFlow_;

        //- Flag to indicate to solve for momentum
        bool momentumPredictor_;

        //- Flag to indicate to solve using transonic algorithm
        bool transonic_;

        //- Flag to indicate to relax pressure using the "consistent" approach
        //  of simplec
        bool consistent_;


public:

    // Static Data Members

        //- Run-time type information
        TypeName("fluidSolutionControl");


    // Constructors

        //- Construct from a mesh and the name of the algorithm
        fluidSolutionControl(fvMesh& mesh, const word& algorithmName);


    //- Destructor
    virtual ~fluidSolutionControl();


    // Member Functions

        // IO

            //- Read controls
            virtual bool read();

        // Access

            //- Flag to indicate to solve for the flow
            inline bool solveFlow() const;

            //- Flag to indicate to solve for momentum
            inline bool momentumPredictor() const;

            //- Flag to indicate to solve using transonic algorithm
            inline bool transonic() const;

            //- Flag to indicate to relax pressure using the "consistent"
            //  approach of simplec
            inline bool consistent() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fluidSolutionControlI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
