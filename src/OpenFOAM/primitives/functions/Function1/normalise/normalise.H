/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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
    Foam::Function1s::normalise

Description
    Function1 which scales a given 'value' function so that its integral within
    specified bounds is equal to one

Usage
    Example of use within a scale function to scale a given profile to a
    desired total integral value:
    \verbatim
        <name>
        {
            type    scale;

            scale
            {
                type    normalise;

                bounds  (0 1.5);

                value
                {
                    type    table;
                    values
                    (
                        (0    0.1272)
                        (0.25 8.7011)
                        (0.5  8.8882)
                        (0.75 6.7608)
                        (1    4.9547)
                        (1.25 3.9213)
                        (1.5  0     )
                    );
                }
            }

            value   constant 6e-6;
        }
    \endverbatim

    Where:
    \table
        Property | Description               | Data type       | Default
        bounds   | Bounds of the integration | Pair<scalar>    |
        value    | Function to normalise     | Function1<Type> |
    \endtable

SourceFiles
    normalise.C

\*---------------------------------------------------------------------------*/

#ifndef normalise_H
#define normalise_H

#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1s
{

/*---------------------------------------------------------------------------*\
                         Class ScaleIntegral Declaration
\*---------------------------------------------------------------------------*/

class normalise
:
    public FieldFunction1<scalar, normalise>
{
    // Private Data

        //- Bounds over which to normalise
        const Pair<scalar> bounds_;

        //- Value function
        const autoPtr<Function1<scalar>> value_;

        //- Computed scale
        const scalar scale_;


public:

    // Runtime type information
    TypeName("normalise");


    // Constructors

        //- Construct from name and dictionary
        normalise
        (
            const word& name,
            const unitConversions& units,
            const dictionary& dict
        );

        //- Copy constructor
        normalise(const normalise& se);


    //- Destructor
    virtual ~normalise();


    // Member Functions

        //- Return value
        virtual inline scalar value(const scalar x) const;

        //- Integrate between two values
        virtual inline scalar integral(const scalar x1, const scalar x2) const;

        //- Write data to dictionary stream
        virtual void write(Ostream& os, const unitConversions& units) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const normalise&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1s
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "normaliseI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
