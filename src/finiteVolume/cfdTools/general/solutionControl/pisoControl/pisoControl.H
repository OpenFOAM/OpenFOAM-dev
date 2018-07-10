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
    Foam::pisoControl

Description
    PISO control class. Provides time-loop and PISO-loop control methods. No
    convergence checking is done.

SourceFiles
    pisoControlI.H
    pisoControl.C

\*---------------------------------------------------------------------------*/

#ifndef pisoControl_H
#define pisoControl_H

#include "fluidSolutionControl.H"

#define PISO_CONTROL

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class pisoControl Declaration
\*---------------------------------------------------------------------------*/

class pisoControl
:
    public fluidSolutionControl
{
protected:

    // Protected data

        //- Maximum number of PISO correctors
        label nCorrPISO_;

        //- Current PISO corrector
        label corrPISO_;


public:

    // Static data members

        //- Run-time type information
        TypeName("pisoControl");


    // Constructors

        //- Construct from a mesh and the name of the algorithm
        pisoControl(fvMesh& mesh, const word& algorithmName="PISO");


    //- Destructor
    virtual ~pisoControl();


    // Member Functions

        // IO

            //- Read controls
            virtual bool read();

        // Access

            //- Maximum number of PISO correctors
            inline label nCorrPISO() const;

            //- Flag to indicate the first PISO iteration
            inline bool firstPISOIter() const;

            //- Flag to indicate the last PISO iteration
            inline bool finalPISOIter() const;

            //- Flag to indicate the last inner iteration (last PISO and last
            //  non-orthogonal)
            inline bool finalInnerIter() const;


        // Evolution

            //- PISO loop
            bool correct();

            //- Time run loop
            bool run(Time& time);

            //- Time loop loop
            bool loop(Time& time);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pisoControlI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
