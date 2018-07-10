/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
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
    Foam::pimpleControl

Description
    Pimple control class. Provides time-loop control methods which exit the
    simulation once convergence criteria have been reached. Also provides
    Pimple-loop control methods which exit the iteration once corrector
    convergence criteria have been met. Example usage:

    \verbatim
    pimpleControl pimple(mesh);

    while (pimple.run(runTime))
    {
        // pre-time-increment operations ...

        runTime ++;
        Info<< "Time = " << runTime.timeName() << nl << endl;

        while (pimple.loop())
        {
            // solve ...
        }

        // post-solve operations ...
    }
    \endverbatim

SourceFiles
    pimpleControlI.H
    pimpleControl.C

\*---------------------------------------------------------------------------*/

#ifndef pimpleControl_H
#define pimpleControl_H

#include "pimpleNoLoopControl.H"
#include "pimpleLoop.H"
#include "singleRegionConvergenceControl.H"
#include "singleRegionCorrectorConvergenceControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class pimpleControl Declaration
\*---------------------------------------------------------------------------*/

class pimpleControl
:
    public pimpleNoLoopControl,
    public pimpleLoop
{
public:

    // Static data members

        //- Run-time type information
        TypeName("pimpleControl");


    // Constructors

        //- Construct from a mesh and the name of the algorithm
        pimpleControl(fvMesh& mesh, const word& algorithmName="PIMPLE");


    //- Destructor
    virtual ~pimpleControl();


    // Member Functions

        // IO

            //- Read controls
            virtual bool read();

        // Access

            //- Flag to indicate whether to solve the turbulence
            inline bool turbCorr() const;

        // Evolution

            //- Pimple loop
            bool loop();

            //- Time run loop
            bool run(Time& time);

            //- Time loop loop
            bool loop(Time& time);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pimpleControlI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
