/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2021 OpenFOAM Foundation
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
    Foam::Function2s::Coded

Description
    Constructs a dynamically compiled function of two variables.

See also
    Foam::dynamicCode
    Foam::functionEntries::codeStream

SourceFiles
    CodedFunction2I.H
    CodedFunction2.C

\*---------------------------------------------------------------------------*/

#ifndef CodedFunction2_H
#define CodedFunction2_H

#include "Function2.H"
#include "codedBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function2s
{

/*---------------------------------------------------------------------------*\
                            Class coded Declaration
\*---------------------------------------------------------------------------*/

class coded
{
public:

    ClassNameNoDebug("Function2");
};


/*---------------------------------------------------------------------------*\
                            Class Coded Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Coded
:
    public Function2<Type>,
    public codedBase
{
    // Private Data

        //- The dynamically generated Function2 pointer
        mutable autoPtr<Function2<Type>> redirectFunction2Ptr_;


    // Private Member Functions

        //- Get the keywords associated with source code
        virtual wordList codeKeys() const;

        //- Adapt the context for the current object
        virtual void prepare(dynamicCode&, const dynamicCodeContext&) const;

        //- Clear the ptr to the redirected object
        virtual void clearRedirect() const;

        //- Compile, link and return the now coded Function2
        autoPtr<Function2<Type>> compileNew();


public:

    // Runtime type information
    TypeName("coded");


    // Constructors

        //- Construct from name and dictionary
        Coded
        (
            const word& name,
            const dictionary& dict
        );

        //- Copy constructor
        Coded(const Coded<Type>& cf1);

        //- Construct and return a clone
        virtual tmp<Function2<Type>> clone() const;


    //- Destructor
    virtual ~Coded();


    // Member Functions

        //- Return value as a function of two scalars
        virtual inline Type value(const scalar x, const scalar y) const;

        //- Return value as a function of two scalar fields
        virtual tmp<Field<Type>> value
        (
            const scalarField& x,
            const scalarField& y
        ) const;

        //- Write data to dictionary stream
        virtual void write(Ostream& os) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const Coded<Type>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function2s
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "CodedFunction2I.H"

#ifdef NoRepository
    #include "CodedFunction2.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
