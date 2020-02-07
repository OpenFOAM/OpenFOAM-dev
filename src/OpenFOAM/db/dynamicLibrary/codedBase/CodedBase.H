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
    Foam::CodedBase

Description
    Templated wrapper class to reduce code duplication and simplify maintenance
    of coded classes.

See also
    Foam::codeBase

SourceFiles
    CodedBase.C

\*---------------------------------------------------------------------------*/

#ifndef CodedBase_H
#define CodedBase_H

#include "codedBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class CodedType>
class CodedBase
:
    public codedBase
{
    // Private static data

        //- The keywords associated with source code
        static const wordList codeKeys_;

        //- Name of the dynamically generated CodedType
        const word codeName_;

        //- Dictionary contents for the CodedType
        const dictionary dict_;


    // Private Member Functions

        //- Get the keywords associated with source code
        virtual const wordList& codeKeys() const;

        //- Return a description (type + name) for the output
        virtual string description() const;


public:

    // Static Data Members

        //- Name of the C code template to be used
        static const word codeTemplateC;

        //- Name of the H code template to be used
        static const word codeTemplateH;


    // Constructors

        //- Construct null
        CodedBase()
        {}

        //- Construct from dictionary
        CodedBase(const dictionary& dict);

        //- Copy constructor
        CodedBase(const CodedBase<CodedType>& cb);


    //- Destructor
    virtual ~CodedBase();


    // Member Functions

        //- Name of the dynamically generated CodedType
        virtual const word& codeName() const;

        //- Get the dictionary to initialize the codeContext
        virtual const dictionary& codeDict() const;

        void writeCode(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CodedBase.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
