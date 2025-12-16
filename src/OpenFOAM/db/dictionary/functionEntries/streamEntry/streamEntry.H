/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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
    Foam::functionEntries::streamEntry

Description
    Compiles and executes C++ OpenFOAM code string expressions

    inserting the result into the dictionary or dictionary entry.  The
    functionality of \c #stream is equivalent to \c #codeStream but with a
    simpler syntax, relying on separate \c #codeInclude entries to provide
    required include files.  For more details see
    Foam::functionEntries::codeStream.

Usage
    Example to find the root of f(x) = x + B*x/sqrt(1 + sqr(x)) - A
    using Newton-Raphson iteration

    \verbatim
        A       1;
        B       1;

        x       #stream
        #{
            // Lookup the coefficients
            const scalar A = $A;
            const scalar B = $B;

            // Initial guess for x
            scalar x = 0;

            scalar x0 = x;

            do
            {
                // Store the previous iteration x for the convergence check
                x0 = x;

                // Temporary sub-function evaluations
                const scalar f1 = 1 + sqr(x);
                const scalar f2 = sqrt(f1);

                // Evaluate the function
                const scalar f = x + B*x/f2 - A;

                // Evaluate the derivative
                const scalar df = 1 + B/(f1*f2);

                // Update x
                x = x0 - f/df;

                // Test for convergence
            } while (mag(x - x0) > 1e-6);

            os << x;
        #};
    \endverbatim

See also
    Foam::functionEntries::codeIncludeEntry
    Foam::functionEntries::codeStream
    Foam::functionEntries::calc

SourceFiles
    streamEntry.C

\*---------------------------------------------------------------------------*/

#ifndef streamEntry_H
#define streamEntry_H

#include "functionEntry.H"
#include "OTstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{

// Forward declaration of friend classes
class codeBlockEntry;

/*---------------------------------------------------------------------------*\
                          Class streamEntry Declaration
\*---------------------------------------------------------------------------*/

class streamEntry
:
    public functionEntry
{
protected:

    // Protected Member Functions

        //- Return the code string
        static string codeString
        (
            const label index,
            const dictionary& codeDict,
            Istream&,
            const string& startString = string::null,
            const string& endString = string::null
        );

        //- Perform the calculation and return the result in an OTstream
        static OTstream resultStream
        (
            const dictionary& dict,
            Istream& is,
            const string& startString = string::null,
            const string& endString = string::null
        );


public:

    // Related types

        //- Declare friendship with the streamEntry class
        friend class codeBlockEntry;


    //- Runtime type information
    FunctionTypeName("#stream");


    // Constructors

        //- Construct from line number, dictionary and Istream
        streamEntry
        (
            const label lineNumber,
            const dictionary& parentDict,
            Istream& is
        );

        //- Copy construct
        streamEntry(const streamEntry&) = default;

        //- Clone
        virtual autoPtr<entry> clone(const dictionary&) const
        {
            return autoPtr<entry>(new streamEntry(*this));
        }


    // Member Functions

        //- Expand the functionEntry into the contextDict
        virtual bool execute(dictionary& contextDict, Istream&);

        //- Expand the functionEntry into the contextEntry
        static bool execute
        (
            const dictionary& contextDict,
            primitiveEntry& contextEntry,
            Istream&
        );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const streamEntry&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionEntries
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
