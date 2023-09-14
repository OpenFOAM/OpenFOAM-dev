/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2023 OpenFOAM Foundation
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
    Foam::displacedPhaseInterface

Description
    Class to represent an interface between phases which has been displaced to
    some extent by a third phase.

SourceFiles
    displacedPhaseInterface.C

\*---------------------------------------------------------------------------*/

#ifndef displacedPhaseInterface_H
#define displacedPhaseInterface_H

#include "phaseInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class displacedPhaseInterface Declaration
\*---------------------------------------------------------------------------*/

class displacedPhaseInterface
:
    virtual public phaseInterface
{
private:

    // Private Data

        //- Displacing phase
        const phaseModel& displacing_;


public:

    //- Runtime type information
    TypeName("displacedPhaseInterface");


    // Constructors

        //- Construct from phases
        displacedPhaseInterface
        (
            const phaseModel& phase1,
            const phaseModel& phase2,
            const phaseModel& displacing
        );

        //- Construct from fluid and name
        displacedPhaseInterface
        (
            const phaseSystem& fluid,
            const word& name
        );


    //- Destructor
    virtual ~displacedPhaseInterface();


    // Static Member Functions

        //- Return the separator that delimits this interface's name
        static word separator()
        {
            return "displacedBy";
        }


    // Member Functions

        //- Name
        virtual word name() const;


        // Access

            //- Displacing phase
            const phaseModel& displacing() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
