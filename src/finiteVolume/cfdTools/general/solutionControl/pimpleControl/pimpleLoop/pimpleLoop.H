/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2023 OpenFOAM Foundation
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
    Pimple loop class. Implements the logic which controls the pimple loop
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

        //- Maximum number of correctors
        label nCorr_;

        //- Current corrector
        label corr_;

        //- Flag to indicate convergence of the pimple loop
        bool converged_;


    // Protected Member Functions

        //- Read controls
        virtual bool read();


public:

    // Static Data Members

        //- Run-time type information
        TypeName("pimpleLoop");


    // Constructors

        //- Construct from a solution control
        pimpleLoop(const solutionControl& control);


    //- Destructor
    virtual ~pimpleLoop();


    // Member Functions

        // Access

            //- Maximum number of correctors
            inline label nCorr() const;

            //- Flag to indicate any iteration
            inline bool anyIter() const;

            //- Flag to indicate the first iteration
            inline bool firstIter() const;

            //- Flag to indicate the last iteration
            inline bool finalIter() const;


        // Evolution

            //- Pimple loop
            bool loop(correctorConvergenceControl& convergence);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pimpleLoopI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
