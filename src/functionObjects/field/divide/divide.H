/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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
    Foam::functionObjects::divide

Description
    Divide a list of fields.

    The numerator can be any volume or surface field and the denominators can
    be scalar volume or surface fields.

    Example of function object specification:
    \verbatim
    Ttot
    {
        type            divide;
        libs            ("libfieldFunctionObjects.so");
        fields          (rho p);
        result          psi;
        executeControl  writeTime;
        writeControl    writeTime;
    }
    \endverbatim

See also
    Foam::functionObjects::fieldsExpression
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    divide.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_divide_H
#define functionObjects_divide_H

#include "fieldsExpression.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                             Class divide Declaration
\*---------------------------------------------------------------------------*/

class divide
:
    public fieldsExpression
{
    // Private Member Functions

        //- Multiply the list of fields and return true if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("divide");


    // Constructors

        //- Construct from Time and dictionary
        divide
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~divide();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
