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
    Foam::pimpleLoop

Description
    PIMPLE loop class. Implements the logic which controls the PIMPLE loop
    generically for a given corrector convergence control. Can therefore be
    used be either single- or multi-region control classes.

SourceFiles
    pimpleLoopI.H
    pimpleLoop.C

\*---------------------------------------------------------------------------*/

#ifndef pimpleLoop_H
#define pimpleLoop_H

#include "correctorConvergenceControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class pimpleLoop Declaration
\*---------------------------------------------------------------------------*/

class pimpleLoop
{
protected:

    // Protected data

        //- Reference to the solution control
        const solutionControl& control_;

        //- Maximum number of PIMPLE correctors
        label nCorrPIMPLE_;

        //- Current PIMPLE corrector
        label corrPIMPLE_;

        //- Flag to indicate convergence of the PIMPLE loop
        bool converged_;


    // Protected Member Functions

        //- Read controls
        virtual bool read();


public:

    // Static data members

        //- Run-time type information
        TypeName("pimpleLoop");


    // Constructors

        //- Construct from a solution control
        pimpleLoop(const solutionControl& control);


    //- Destructor
    virtual ~pimpleLoop();


    // Member Functions

        // Access

            //- Maximum number of PIMPLE correctors
            inline label nCorrPIMPLE() const;

            //- Flag to indicate the first PIMPLE iteration
            inline bool firstPIMPLEIter() const;

            //- Flag to indicate the last PIMPLE iteration
            inline bool finalPIMPLEIter() const;

            //- Flag to indicate the first iteration
            inline bool firstIter() const;

            //- Flag to indicate the last iteration
            inline bool finalIter() const;


        // Evolution

            //- PIMPLE loop
            bool loop(correctorConvergenceControl& convergence);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pimpleLoopI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
