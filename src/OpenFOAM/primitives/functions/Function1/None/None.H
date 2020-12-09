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
    Foam::Function1s::None

Description
    Undefined Function1 which returns an error when called.

\*---------------------------------------------------------------------------*/

#ifndef None_H
#define None_H

#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1s
{

/*---------------------------------------------------------------------------*\
                            Class None Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class None
:
    public FieldFunction1<Type, None<Type>>
{
    // Private member data

        //- Name of dictionary from which this function is instantiated
        fileName dictName_;


public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from name and dictionary
        None(const word& name, const dictionary& dict);


    //- Destructor
    virtual ~None();


    // Member Functions

        //- Generates an error
        virtual Type value(const scalar x) const;

        //- Generates an error
        virtual Type integral(const scalar x1, const scalar x2) const;

        //- Write in dictionary format
        virtual void write(Ostream& os) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const None<Type>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1s
} // End namespace Foam

#ifdef NoRepository
    #include "None.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
