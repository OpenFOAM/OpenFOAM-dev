/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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
    Foam::pimpleSingleRegionControl

Description
    Pimple single region control class. Provides time-loop control methods
    which exit the simulation once convergence criteria have been reached. Also
    provides Pimple-loop control methods which exit the iteration once corrector
    convergence criteria have been met. Example usage:

    \verbatim
    pimpleSingleRegionControl pimple(mesh);

    while (pimple.run(runTime))
    {
        // pre-time-increment operations ...

        runTime ++;
        Info<< "Time = " << runTime.userTimeName() << nl << endl;

        while (pimple.loop())
        {
            // solve ...
        }

        // post-solve operations ...
    }
    \endverbatim

SourceFiles
    pimpleSingleRegionControlI.H
    pimpleSingleRegionControl.C

\*---------------------------------------------------------------------------*/

#ifndef pimpleSingleRegionControl_H
#define pimpleSingleRegionControl_H

#include "pimpleLoop.H"
#include "pimpleNoLoopControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class pimpleSingleRegionControl Declaration
\*---------------------------------------------------------------------------*/

class pimpleSingleRegionControl
:
    public pimpleLoop
{
    pimpleNoLoopControl& pimple_;


protected:

    // Protected Member Functions

        //- Read controls
        virtual bool read();


public:

    // Static Data Members

        //- Run-time type information
        TypeName("pimpleSingleRegionControl");


    // Constructors

        //- Construct from single region pimple
        pimpleSingleRegionControl(pimpleNoLoopControl& pimple);

        //- Disallow default bitwise copy construction
        pimpleSingleRegionControl(const pimpleSingleRegionControl&) = delete;


    //- Destructor
    virtual ~pimpleSingleRegionControl();


    // Member Functions

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

#include "pimpleSingleRegionControlI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
