/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2019 OpenFOAM Foundation
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
    Foam::adaptiveLinear

Description
    Produces a linear ramp which adapts its gradient to changes in
    endTime and deltaT to always arrive at the relaxationEnd value at the end of
    the run

SourceFiles
    adaptiveLinear.C

\*---------------------------------------------------------------------------*/

#ifndef adaptiveLinear_H
#define adaptiveLinear_H

#include "relaxationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class adaptiveLinear Declaration
\*---------------------------------------------------------------------------*/

class adaptiveLinear
:
    public relaxationModel
{

private:

    // Private Data

        //- Relaxation coefficient at the start of the iteration sequence.
        scalar relaxationStart_;

        //- Relaxation coefficient at the end of the iteration sequence.
        scalar relaxationEnd_;

        //- Store the time when the last request was made for relaxation,
        //  prevents multiple calls to relaxation in a timestep from
        //  incrementing the value
        scalar lastTimeValue_;

        //- Current relaxation value
        scalar relaxation_;


public:

    //- Runtime type information
    TypeName("adaptiveLinear");

    // Constructors

        //- Construct from components
        adaptiveLinear
        (
            const dictionary& relaxationDict,
            const Time& runTime
        );


    //- Destructor
    virtual ~adaptiveLinear()
    {}


    // Member Functions

        //- Return the current relaxation coefficient
        virtual scalar relaxation();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
