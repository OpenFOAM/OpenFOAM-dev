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
    Foam::Function2s::UniformTable

Description
    Tabulated property function that linearly interpolates between
    the uniformTable values.

Usage
    \table
        Property    | Description
        Tlow        | Lower temperature limit of the uniformTable
        Thigh       | Upper temperature limit of the uniformTable
        values      | Property values, assumed uniformly distributed
    \endtable

    Example for the density of water between 280 and 350K
    \verbatim
    rho
    {
        type    uniformTable;

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

#ifndef UniformTable2_H
#define UniformTable2_H

#include "Function2.H"
#include "RectangularMatrix.H"
#include "Pair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function2s
{

/*---------------------------------------------------------------------------*\
                            Class UniformTable Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class UniformTable
:
    public FieldFunction2<Type, UniformTable<Type>>
{
    // Private member data

        //- Lowest values in the table
        Pair<scalar> low_;

        //- Highest values in the table
        Pair<scalar> high_;

        //- Table values
        RectangularMatrix<Type> values_;

        //- x increment derived from low_, high_ and values_.m()
        scalar deltax_;

        //- y increment derived from low_, high_ and values_.n()
        scalar deltay_;


    // Private member functions

        inline void checkRange
        (
            scalar x,
            scalar ndx,
            label ix,
            scalar y,
            scalar ndy,
            label iy
        ) const;


public:

    //- Runtime type information
    TypeName("uniformTable");


    // Constructors

        //- Construct from name and dictionary
        UniformTable(const word& name, const dictionary& dict);


    // Member Functions

        //- Return the non-uniform table of values
        const RectangularMatrix<Type>& values() const
        {
            return values_;
        }

        //- Evaluate the function and return the result
        virtual Type value(scalar x, scalar y) const;

        //- Evaluate the derivative of the function w.r.t. p
        Type dfdp(scalar p, scalar T) const;

        //- Evaluate the derivative of the function w.r.t. T
        Type dfdT(scalar p, scalar T) const;

        //- Write the function coefficients
        void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function2s
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "UniformTable2.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
