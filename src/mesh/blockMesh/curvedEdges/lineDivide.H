/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::lineDivide

Description
    Divides a line into segments

SourceFiles
    lineDivide.C

\*---------------------------------------------------------------------------*/

#ifndef lineDivide_H
#define lineDivide_H

#include "pointField.H"
#include "scalarList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class curvedEdge;

/*---------------------------------------------------------------------------*\
                         Class lineDivide Declaration
\*---------------------------------------------------------------------------*/

class lineDivide
{
    // Private data

        pointField points_;

        scalarList divisions_;

public:

    // Constructors

        //- Construct from components with discretization and expansion ratio
        lineDivide
        (
            const curvedEdge&,
            const label ndiv,
            const scalar xratio = 1.0
        );


    // Member Functions

        //- Return the points
        const pointField& points() const;

        //- Return the list of lambda values
        const scalarList& lambdaDivisions() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
