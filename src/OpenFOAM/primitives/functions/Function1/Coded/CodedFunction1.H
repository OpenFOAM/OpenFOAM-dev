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
    Foam::Function1s::Coded

Description
    Constructs on-the-fly a Function1 which is then used to evaluate.

Usage
    Example:
    \verbatim
    <patchName>
    {
        type            codedFixedValue;
        value           uniform 0;
        name            rampedFixedValue;   // name of generated BC

        code
        #{
            operator==(min(10, 0.1*this->db().time().value()));
        #};

        // codeInclude
        //#{
        //    #include "fvCFD.H"
        //#};

        // codeOptions
        //#{
        //    -I$(LIB_SRC)/finiteVolume/lnInclude
        //#};
    }
    \endverbatim

    A special form is if the 'code' section is not supplied. In this case
    the code is read from a (runTimeModifiable!) dictionary system/codeDict
    which would have a corresponding entry:

    \verbatim
    <patchName>
    {
        code
        #{
            operator==(min(10, 0.1*this->db().time().value()));
        #};
    }
    \endverbatim

See also
    Foam::dynamicCode
    Foam::functionEntries::codeStream

SourceFiles
    CodedI.H
    Coded.C

\*---------------------------------------------------------------------------*/

#ifndef CodedFunction1_H
#define CodedFunction1_H

#include "Function1.H"
#include "codedBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1s
{

/*---------------------------------------------------------------------------*\
                           Class Coded Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Coded
:
    public Function1<Type>,
    public codedBase
{
    // Private static data

        //- The keywords associated with source code
        static const wordList codeKeys_;


    // Private Data

        //- Dictionary contents for the Function1
        const dictionary dict_;

        //- Name of the dynamically generated Function1
        const word name_;

        //- The dynamically generated Function1 pointer
        mutable autoPtr<Function1<Type>> redirectFunction1Ptr_;


    // Private Member Functions

        //- Adapt the context for the current object
        virtual void prepare(dynamicCode&, const dynamicCodeContext&) const;

        //- Return a description (type + name) for the output
        virtual string description() const;

        //- Clear the ptr to the redirected object
        virtual void clearRedirect() const;

        //- Get the dictionary to initialize the codeContext
        virtual const dictionary& codeDict() const;

        //- Get the keywords associated with source code
        virtual const wordList& codeKeys() const;


public:

    // Static Data Members

        //- Name of the C code template to be used
        static const word codeTemplateC;

        //- Name of the H code template to be used
        static const word codeTemplateH;


    // Runtime type information
    TypeName("coded");


    // Constructors

        //- Construct from entry name and dictionary
        Coded
        (
            const word& entryName,
            const dictionary& dict
        );

        //- Copy constructor
        Coded(const Coded<Type>& cf1);

        //- Construct and return a clone
        virtual tmp<Function1<Type>> clone() const;


    //- Destructor
    virtual ~Coded();


    // Member Functions

        //- Return value
        virtual inline Type value(const scalar x) const;

        //- Return value as a function of (scalar) independent variable
        virtual tmp<Field<Type>> value(const scalarField& x) const;

        //- Integrate between two values
        virtual Type integrate(const scalar x1, const scalar x2) const;

        //- Integrate between two (scalar) values
        virtual tmp<Field<Type>> integrate
        (
            const scalarField& x1,
            const scalarField& x2
        ) const;

        //- Write in dictionary format
        virtual void writeData(Ostream& os) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const Coded<Type>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1s
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "CodedFunction1I.H"

#ifdef NoRepository
    #include "CodedFunction1.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
