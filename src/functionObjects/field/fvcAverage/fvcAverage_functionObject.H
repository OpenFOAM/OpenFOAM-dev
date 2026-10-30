/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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
    Foam::functionObjects::fvcAverage

Description
    Calculates a cell-centred volField by averaging the face values of a
    surfaceField using fvc::average.

See also
    Foam::functionObjects::reconstruct

SourceFiles
    fvcAverage.C

\*---------------------------------------------------------------------------*/

#ifndef fvcAverage_functionObject_H
#define fvcAverage_functionObject_H

#include "fieldExpression.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class fvcAverage Declaration
\*---------------------------------------------------------------------------*/

class fvcAverage
:
    public fieldExpression
{
    // Private Member Functions

        //- Calculate the magnitude of the field and register the result
        template<class Type>
        inline bool calcFvcAverage();

        //- Calculate the fvcAverageion field and return true if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("fvcAverage");


    // Constructors

        //- Construct from Time and dictionary
        fvcAverage
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~fvcAverage();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
