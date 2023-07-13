/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    multi-line input without the need to escape the newlines but in this case
    the result must be written to the \c Ostream \c os explicitly.

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

    If the code string is delimited by \#{ ... \c \#} multiple lines
    and multiple code statements can be used to generate the entry using
    'os << ...;'. This is equivalent to \#codeStream but with a more
    compact syntax.

    \verbatim
    #calcInclude    "transform.H"
    maxAngle        30;
    nAngles         7;
    Us              #calc
    #{
        const vector U($<vector>testCalc2!U);
        const int nAngles = $nAngles;
        const scalar angleStep = ($<scalar>maxAngle)/(nAngles - 1);
        List<vector> Us(nAngles);
        for(int i=0; i<nAngles; i++)
        {
            const scalar angle = degToRad(i*angleStep);
            Us[i] = transform(Ry(angle), U);
        }
        os << Us;
    #};
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

class dlLibraryTable;

namespace functionEntries
{

/*---------------------------------------------------------------------------*\
                          Class calcEntry Declaration
\*---------------------------------------------------------------------------*/

class calcEntry
:
    public functionEntry
{
    // Private Member Functions

        //- Perform the calculation and return the resulting string
        static string calc
        (
            const dictionary& dict,
            Istream& is
        );


public:

    //- Runtime type information
    ClassName("calc");


    // Constructors

        //- Disallow default bitwise copy construction
        calcEntry(const calcEntry&) = delete;


    // Member Functions

        //- Execute the functionEntry in a sub-dict context
        static bool execute(dictionary& dict, Istream&);

        //- Execute the functionEntry in a primitiveEntry context
        static bool execute
        (
            const dictionary& dict,
            primitiveEntry&,
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
