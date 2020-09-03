/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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
    Foam::thermophysicalFunctions::table2D

Description
    Tabulated property function that linearly interpolates between
    the table2D values.

Usage
    \table
        Property    | Description
        Tlow        | Lower temperature limit of the table2D
        Thigh       | Upper temperature limit of the table2D
        values      | Property values, assumed uniformly distributed
    \endtable

    Example for the density of water between 280 and 350K
    \verbatim
    rho
    {
        type    table2D;

        pLow    1000;
        pHigh   1e6;

        Tlow    280;
        Thigh   350;

        values
        (
            999.87
            991.53
            973.7
        );
    }
    \endverbatim

\*---------------------------------------------------------------------------*/

#ifndef table2DThermophysicalFunction_H
#define table2DThermophysicalFunction_H

#include "thermophysicalFunction.H"
#include "scalarMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace thermophysicalFunctions
{

/*---------------------------------------------------------------------------*\
                            Class table2D Declaration
\*---------------------------------------------------------------------------*/

class table2D
:
    public thermophysicalFunction
{
    // Private member data

        //- Name of dictionary from which this function is instantiated
        word name_;

        //- Lowest pressure in the table
        scalar pLow_;

        //- Highest pressure in the table
        scalar pHigh_;

        //- Lowest temperature in the table
        scalar Tlow_;

        //- Highest temperature in the table
        scalar Thigh_;

        //- Table values
        scalarRectangularMatrix values_;

        //- pressure increment derived from pLow_, pHigh_ and values_.m()
        scalar deltap_;

        //- Temperature increment derived from Tlow_, Thigh_ and values_.n()
        scalar deltaT_;


    // Private member functions

        inline void checkRange
        (
            scalar p,
            scalar ndp,
            label ip,
            scalar T,
            scalar ndT,
            label iT
        ) const;


public:

    //- Runtime type information
    TypeName("table2D");


    // Constructors

        //- Construct from entry name and dictionary
        table2D(const word& name, const dictionary& dict);

        //- Construct from dictionary
        table2D(const dictionary& dict);


    // Member Functions

        //- Return the non-uniform table of values
        const scalarRectangularMatrix& values() const
        {
            return values_;
        }

        //- Evaluate the function and return the result
        scalar f(scalar p, scalar T) const;

        //- Evaluate the derivative of the function w.r.t. p
        scalar dfdp(scalar p, scalar T) const;

        //- Evaluate the derivative of the function w.r.t. T
        scalar dfdT(scalar p, scalar T) const;

        //- Write the function coefficients
        void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace thermophysicalFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
