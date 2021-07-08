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
    Foam::functionObjects::multiply

Description
    Multiply a list of fields.

    The operation can be applied to any volume or surface fields provided that
    the resulting multiplication does not produce a tensor of rank greater than
    two.

    Example of function object specification:
    \verbatim
    Ttot
    {
        type            multiply;
        libs            ("libfieldFunctionObjects.so");
        fields          (rho U);
        result          rhoU;
        executeControl  writeTime;
        writeControl    writeTime;
    }
    \endverbatim

See also
    Foam::functionObjects::fieldsExpression
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    multiply.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_multiply_H
#define functionObjects_multiply_H

#include "fieldsExpression.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                             Class multiply Declaration
\*---------------------------------------------------------------------------*/

class multiply
:
    public fieldsExpression
{
    // Private Member Functions

        //- Multiply the list of fields and return true if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("multiply");


    // Constructors

        //- Construct from Time and dictionary
        multiply
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~multiply();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
