/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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
    Foam::rampHoldFall

Description
    Piecewise linear function with a ramp from a start value to a plateaux
    value, holding at this, then a linear fall to an end value.

SourceFiles
    rampHoldFall.C

\*---------------------------------------------------------------------------*/

#ifndef rampHoldFall_H
#define rampHoldFall_H

#include "relaxationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class rampHoldFall Declaration
\*---------------------------------------------------------------------------*/

class rampHoldFall
:
    public relaxationModel
{

private:

    // Private data

        //- Relaxation coefficient at the start of the ramp
        scalar rampStartRelaxation_;

        //- Relaxation coefficient for the hold portion
        scalar holdRelaxation_;

        //- Relaxation coefficient at the end of the fall
        scalar fallEndRelaxation_;

        //- Fraction through the run where the ramp ends and the hold starts
        scalar rampEndFraction_;

        //- Fraction through the run where the hold ends and the fall starts
        scalar fallStartFraction_;

        //- Gradient in ramp portion, normalised against time
        scalar rampGradient_;

        //- Gradient in fall portion, normalised against time
        scalar fallGradient_;

public:

    //- Runtime type information
    TypeName("rampHoldFall");

    // Constructors

        //- Construct from components
        rampHoldFall
        (
            const dictionary& relaxationDict,
            const Time& runTime
        );


    //- Destructor
    virtual ~rampHoldFall()
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
