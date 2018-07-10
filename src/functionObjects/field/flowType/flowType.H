/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
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
    Foam::functionObjects::flowType

Description
    Calculates and writes the flowType of a velocity field.

    The flow type parameter is obtained according to the following equation:
    \verbatim
                 |D| - |Omega|
        lambda = -------------
                 |D| + |Omega|

        -1 = rotational flow
         0 = simple shear flow
         1 = planar extensional flow
    \endverbatim

See also
    Foam::functionObjects::fieldExpression
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    flowType.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_flowType_H
#define functionObjects_flowType_H

#include "fieldExpression.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class flowType Declaration
\*---------------------------------------------------------------------------*/

class flowType
:
    public fieldExpression
{
    // Private Member Functions

        //- Calculate the flowType field and return true if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("flowType");


    // Constructors

        //- Construct from Time and dictionary
        flowType
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~flowType();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
