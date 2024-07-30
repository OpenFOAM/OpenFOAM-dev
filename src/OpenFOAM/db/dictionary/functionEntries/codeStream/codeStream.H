/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
    Foam::functionEntries::codeStream

Description
    Dictionary entry that contains C++ OpenFOAM code that is compiled to
    generate the entry itself.

    \c \#codeStream reads three entries: \c code, \c codeInclude (optional),
    \c codeOptions (optional) to generate the library source code stored in the
    local \c dynamicCode directory with a subdirectory name corresponding to the
    SHA1 of the code.  The code is then compiled into a dynamically loaded
    library libcodeStream_<SHA1>.so stored in the
    \c dynamicCode/platforms/$WM_OPTIONS/lib directory using 'wmake libso'.
    The resulting library is loaded in executed with arguments
    \code
        (const dictionary& dict, Ostream& os)
    \endcode
    where the dictionary is the current dictionary.  The code writes results to
    the current entry via the \c Ostream \c os.

    The verbatim string format \c \#{ ... \c \#} is used to allow multi-line
    input without the need to escape the newlines.

    Dictionary entries constructed with \c \#codeStream can conveniently access
    and use typed variables.  This means calculations involving vectors and
    tensors and list etc. are possible.  To access a variable and construct it
    as a given type within a \c \#codeStream entry, put the type immediately
    after the $ symbol inside angled brackets <>. So, $<vector>var or
    $<vector>{var} substitutes a variable named var as a vector.

    The variable values are no longer embedded into the code but looked-up at
    run-time making the code corresponding to each \c \#codeStream independent
    of the values in the dictionary and of each other.  Hence the \c
    \#codeStream code does not need to be recompiled when the values in the
    dictionary are changed, only if the code is changed.

Usage
    Example to set the internal field of a field:
    \verbatim
        internalField  #codeStream
        {
            code
            #{
                const IOdictionary& d = static_cast<const IOdictionary&>(dict);
                const fvMesh& mesh = refCast<const fvMesh>(d.db());
                scalarField fld(mesh.nCells(), 12.34);
                writeEntry(os, "", fld);
            #};

            //- Optional:
            codeInclude
            #{
                #include "volFields.H"
            #};

            //- Optional:
            codeOptions
            #{
                -I$(LIB_SRC)/finiteVolume/lnInclude
            #};
        };
    \endverbatim

    Example to rotate a list of points around an axis by a given angle
    \verbatim
        points  ((3 0 0) (2 1 1) (1 2 2) (0 3 3));
        rotation
        {
            axis    (0 1 1);
            angle   45;
        }

        #codeStream
        {
            codeInclude
            #{
                #include "pointField.H"
                #include "transform.H"
            #};

            code
            #{
                const pointField points($<List<point>>points);
                const vector axis = $<vector>!rotation/axis;
                const scalar angle = degToRad($!rotation/angle);
                os << "pointsRotated" << nl
                   << (Ra(axis, angle) & points)() << ";";
            #};
        };
    \endverbatim

    Example Compute the centre and trianglation of a polygon
    \verbatim
       polygon  ((0 0 0) (1 0 0) (2 1 0) (0 2 0) (-1 1 0));

       #codeStream
       {
           codeInclude
           #{
               #include "polygonTriangulate.H"
           #};

           code
           #{
               const List<point> polygon($<List<point>>polygon);
               writeEntry(os, "polygonCentre", face::centre(polygon));

               polygonTriangulate triEngine;
               triEngine.triangulate(polygon);
               os << "polygonTris" << ' ' << triEngine.triPoints() << ";";
           #};
        };
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
        vertices    #codeStream
        {
            codeInclude
            #{
                #include "boundBox.H"
            #};

            code
            #{
                os << boundBox($<vector>min, $<vector>max).points();
            #};
        };

        blocks
        (
            hex (0 1 2 3 4 5 6 7) $nCells simpleGrading (1 1 1)
        );

        defaultPatch
        {
            type patch;
        }

        boundary
        ();
    \endverbatim

See also
    Foam::functionEntries::calcEntry

SourceFiles
    codeStream.C

\*---------------------------------------------------------------------------*/

#ifndef codeStream_H
#define codeStream_H

#include "functionEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
class dlLibraryTable;

namespace functionEntries
{

// Forward declaration of friend classes
class calcEntry;

/*---------------------------------------------------------------------------*\
                         Class codeStream Declaration
\*---------------------------------------------------------------------------*/

class codeStream
:
    public functionEntry
{
    //- Interpreter function type
    typedef void (*streamingFunctionType)(Ostream&, const dictionary&);


    // Private Member Functions

        //- Helper function: parent (of parent etc.) of dictionary up to the top
        static const dictionary& topDict(const dictionary&);

        //- Return true for master only reading of global dictionary
        static bool masterOnlyRead(const dictionary& dict);

        //- Construct, compile, load and return streaming function
        static streamingFunctionType getFunction
        (
            const dictionary& contextDict,
            const dictionary& codeDict
        );

        //- Compile and execute the code and return the resulting string
        static string run
        (
            const dictionary& contextDict,
            Istream& is
        );


public:

    // Static Data Members

        //- Keywords associated with source code
        static const wordList codeKeys;

        //- Name of the dictionary variables in the source code
        static const wordList codeDictVars;

        //- Name of the C code template to be used
        static const word codeTemplateC;


    // Related types

        //- Declare friendship with the calcEntry class
        friend class calcEntry;


    //- Runtime type information
    ClassName("codeStream");


    // Constructors

        //- Disallow default bitwise copy construction
        codeStream(const codeStream&) = delete;


    // Member Functions

        //- Execute the functionEntry in a sub-dict context
        static bool execute(dictionary& contextDict, Istream&);

        //- Execute the functionEntry in a primitiveEntry context
        static bool execute
        (
            const dictionary& contextDict,
            primitiveEntry&,
            Istream&
        );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const codeStream&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionEntries
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
