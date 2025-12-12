/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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
    Foam::functionEntries::calcEntry

Description
    Compiles and executes code string expressions,
    returning the result to the dictionary entry.

    \c \#calc reads the following code string to generate the library source
    code stored in the local \c dynamicCode directory with a subdirectory name
    corresponding to the SHA1 of the code.  The code is then compiled into a
    dynamically loaded library libcodeStream_<SHA1>.so stored in the \c
    dynamicCode/platforms/$WM_OPTIONS/lib directory using 'wmake libso'.  The
    resulting library is loaded in executed with arguments
    \code
        (const dictionary& dict, Ostream& os)
    \endcode
    where the dictionary is the current dictionary.  The code writes results to
    the current entry via the \c Ostream \c os.

    The verbatim string format \c \#{ ... \c \#} can optionally be used to allow
    string/word comparison without the need to escape the quotes and multi-line
    input without the need to escape the newlines.

    Dictionary entries constructed \c \#calc can conveniently access and
    use typed variables.  This means calculations involving vectors and tensors
    and list etc. are possible.  To access a variable and construct it as a
    given type within a \c \#calc entry, put the type immediately after the $
    symbol inside angled brackets <>. So, $<vector>var or $<vector>{var}
    substitutes a variable named var as a vector.

    The variable values are no longer embedded into the code but looked-up at
    run-time making the code corresponding to each \c \#calc independent of the
    values in the dictionary and of each other.  Hence the \c \#calc code does
    not need to be recompiled when the values in the dictionary are changed,
    only if the code is changed.

Usage
    Simple variable multiplication
    \verbatim
    a 1.1;
    b 3.2;
    c #calc "$a*$b";
    \endverbatim

    Special care is required for calc entries that include a division since
    "/" is also used as the scoping operator to identify keywords in
    sub-dictionaries. For example, "$a/b" expects a keyword "b" within a
    sub-dictionary named "a". A division can be correctly executed by using a
    space between a variables and "/", e.g.
    \verbatim
    c #calc "$a / $b";
    \endverbatim

    or bracketing the variable, e.g.
    \verbatim
    c #calc "($a)/$b";
    \endverbatim

    // Example using the string type
    \verbatim
    s   "field";
    fieldName #calc "$<string>s + \"Name\" ";
    \endverbatim

    Additional include files for the \c \#calc code compilation can be specified
    using the \c \#calcInclude entry, e.g. if functions from transform.H are
    used
    \verbatim
    angleOfAttack   5; // degs

    angle           #calc "-degToRad($angleOfAttack)";

    #calcInclude    "transform.H"
    liftDir         #calc "transform(Ry($angle), vector(0, 0, 1))";
    dragDir         #calc "transform(Ry($angle), vector(1, 0, 0))";
    \endverbatim

    Calculate the magnitude of the velocity and turbulent kinetic energy
    where the velocity is looked-up from testCalc2
    \verbatim
    magU        #calc "mag($<vector>testCalc2!U)";
    k           #calc "1.5*magSqr(0.05*$<vector>{${FOAM_CASE}/testCalc2!U})";
    \endverbatim

    Example to generate a single block blockMeshDict for use with snappyHexMesh
    with no redundant information
    \verbatim
        min         (-2.5 -1.2 -3.0);   // Minimum coordinates of the block
        max         (2.5 1.2 3.0);      // Maximum coordinates of the block
        nCellsByL   33.3333;            // Number of cells per unit length

        // Calculate the number of cells in each block direction
        nCells      #calc
    "Vector<label>($nCellsByL*($<vector>max - $<vector>min) + vector::one/2)";

        // Generate the vertices using a boundBox
        #calcInclude "boundBox.H"
        vertices #calc "boundBox($<vector>min, $<vector>max).points()";

        blocks
        (
            hex #calc "identityMap(8)" $nCells simpleGrading (1 1 1)
        );

        defaultPatch
        {
            type patch;
        }

        boundary
        ();
    \endverbatim

See also
    Foam::functionEntries::calcIncludeEntry
    Foam::functionEntries::codeStream

SourceFiles
    calcEntry.C

\*---------------------------------------------------------------------------*/

#ifndef calcEntry_H
#define calcEntry_H

#include "functionEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{

// Forward declaration of friend classes
class codeBlockEntry;

/*---------------------------------------------------------------------------*\
                          Class calcEntry Declaration
\*---------------------------------------------------------------------------*/

class calcEntry
:
    public functionEntry
{
    // Private Member Functions

        //- Add the code entry to the codeDict
        static string codeString
        (
            const label index,
            const dictionary& codeDict,
            Istream&
        );

        //- Perform the calculation and return the resulting string
        static string calc
        (
            const dictionary& dict,
            Istream& is
        );


public:

    // Related types

        //- Declare friendship with the calcEntry class
        friend class codeBlockEntry;


    //- Runtime type information
    FunctionTypeName("#calc");


    // Constructors

        //- Copy construct
        calcEntry(const calcEntry&) = default;

        //- Clone
        virtual autoPtr<entry> clone(const dictionary&) const
        {
            return autoPtr<entry>(new calcEntry(*this));
        }


    // Member Functions

        //- Expand the functionEntry into the contextDict
        virtual bool execute(dictionary& contextDict, Istream&)
        {
            NotImplemented;
            return false;
        }

        //- Expand the functionEntry into the contextEntry
        static bool execute
        (
            const dictionary& contextDict,
            primitiveEntry& contextEntry,
            Istream&
        );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const calcEntry&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionEntries
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
