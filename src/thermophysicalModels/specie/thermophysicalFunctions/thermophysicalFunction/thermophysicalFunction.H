/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::thermophysicalFunction

Description
    Abstract base class for thermo-physical functions

\*---------------------------------------------------------------------------*/

#ifndef thermophysicalFunction_H
#define thermophysicalFunction_H

#include "scalar.H"
#include "IOstreams.H"
#include "typeInfo.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
class thermophysicalFunction;
Ostream& operator<<(Ostream& os, const thermophysicalFunction& l);

/*---------------------------------------------------------------------------*\
                   Class thermophysicalFunction Declaration
\*---------------------------------------------------------------------------*/

class thermophysicalFunction
{

public:

    //- Runtime type information
    TypeName("thermophysicalFunction");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            thermophysicalFunction,
            dictionary,
            (const dictionary& dict),
            (dict)
        );


    // Constructors

        //- Construct null
        thermophysicalFunction()
        {}

        //- Return pointer to new thermophysicalFunction created from dict
        static autoPtr<thermophysicalFunction> New
        (
            const dictionary& dict,
            const word& name
        );


    //- Destructor
    virtual ~thermophysicalFunction()
    {}


    // Member Functions

        //- Evaluate the function and return the result
        virtual scalar f(scalar p, scalar T) const = 0;

        //- Write the function coefficients
        virtual void write(Ostream& os) const = 0;

        //- Write the function coefficients within the property sub-dictionary
        virtual void write(Ostream& os, const word& name) const;


    // Ostream Operator

        friend Ostream& operator<<
        (
            Ostream& os,
            const thermophysicalFunction& f
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
