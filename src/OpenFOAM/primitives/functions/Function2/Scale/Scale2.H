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
    Foam::Function2s::Scale

Description
    Function2 which scales a given 'value' function by a 'scale' scalar function
    and scales the 'x' and 'y' arguments of the 'value' and 'scale' functions
    by the optional 'xScale' and 'yScale' scalar functions.

See also
    Foam::Function2s::ramp
    Foam::Function2s::reverseRamp

SourceFiles
    Scale.C

\*---------------------------------------------------------------------------*/

#ifndef Scale2_H
#define Scale2_H

#include "Function2.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function2s
{

/*---------------------------------------------------------------------------*\
                            Class Scale Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Scale
:
    public FieldFunction2<Type, Scale<Type>>
{
    // Private Data

        //- Scalar scaling function
        autoPtr<Function2<scalar>> scale_;

        //- Argument scaling function
        autoPtr<Function1<scalar>> xScale_;

        //- Argument scaling function
        autoPtr<Function1<scalar>> yScale_;

        //- Value function
        autoPtr<Function2<Type>> value_;


    // Private Member Functions

        //- Read the coefficients from the given dictionary
        void read(const dictionary& dict);


public:

    // Runtime type information
    TypeName("scale");


    // Constructors

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
        virtual inline Type value(const scalar x, const scalar y) const;

        //- Write data to dictionary stream
        virtual void write(Ostream& os) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const Scale<Type>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function2s
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "Scale2I.H"

#ifdef NoRepository
    #include "Scale2.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
