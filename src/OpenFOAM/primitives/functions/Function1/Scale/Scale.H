/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2021 OpenFOAM Foundation
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
    Foam::Function1s::Scale

Description
    Function1 which scales a given 'value' function by a 'scale' scalar function
    and scales the 'x' argument of the 'value' and 'scale' functions by the
    optional 'xScale' scalar function.

    This is particularly useful to ramp a time-varying value by one of the
    monotonic ramp functions.

    Usage for a vector:
    \verbatim
        <name>
        {
            type      scale;

            scale
            {
                type        linearRamp;

                start       0;
                duration    10;
            }

            value
            {
                type        sine;

                frequency   10;
                amplitude   1;
                scale       (1 0.1 0);
                level       (10 1 0);
            }
        }
    \endverbatim

    Simplified usage to scale by a constant factor, e.g. 2:
    \verbatim
        <name>
        {
            type      scale;

            scale     2;

            value
            {
                type        sine;

                frequency   10;
                amplitude   1;
                scale       (1 0.1 0);
                level       (10 1 0);
            }
        }
    \endverbatim

    Usage including the optional 'xScale' function:
    \verbatim
        <name>
        {
            type      scale;

            xScale    0.5;
            scale     2;

            value
            {
                type        sine;

                frequency   10;
                amplitude   1;
                scale       (1 0.1 0);
                level       (10 1 0);
            }
        }
    \endverbatim

    Where:
    \table
        Property | Description                                    | Required
        value    | Function of type Function1<Type>               | yes
        scale    | Scaling function of type Function1<scalar>     | yes
        xScale   | 'x' scaling function of type Function1<scalar> | no
    \endtable

See also
    Foam::Function1s::ramp
    Foam::Function1s::reverseRamp

SourceFiles
    Scale.C

\*---------------------------------------------------------------------------*/

#ifndef Scale_H
#define Scale_H

#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1s
{

/*---------------------------------------------------------------------------*\
                            Class Scale Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Scale
:
    public FieldFunction1<Type, Scale<Type>>
{
    // Private Data

        //- Scalar scaling function
        autoPtr<Function1<scalar>> scale_;

        //- Argument scaling function
        autoPtr<Function1<scalar>> xScale_;

        //- Value function
        autoPtr<Function1<Type>> value_;

        //- Is this function integrable?
        bool integrableScale_;

        //- Is this function integrable?
        bool integrableValue_;


    // Private Member Functions

        //- Read the coefficients from the given dictionary
        void read(const dictionary& dict);


public:

    // Runtime type information
    TypeName("scale");


    // Constructors

        //- Construct from name and functions
        Scale
        (
            const word& name,
            const Function1<scalar>& scale,
            const Function1<scalar>& xScale,
            const Function1<Type>& value
        );

        //- Construct from name and dictionary
        Scale
        (
            const word& name,
            const dictionary& dict
        );

        //- Copy constructor
        Scale(const Scale<Type>& se);


    //- Destructor
    virtual ~Scale();


    // Member Functions

        //- Return value
        virtual inline Type value(const scalar x) const;

        //- Integrate between two values
        virtual inline Type integral(const scalar x1, const scalar x2) const;

        //- Write data to dictionary stream
        virtual void write(Ostream& os) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const Scale<Type>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1s
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ScaleI.H"

#ifdef NoRepository
    #include "Scale.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
