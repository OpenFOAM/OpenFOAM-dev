/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

%{

#undef yyFlexLexer

 /* ------------------------------------------------------------------------ *\
   ------ local definitions
 \* ------------------------------------------------------------------------ */

#include "STLsurfaceFormatCore.H"

using namespace Foam;

// Dummy yyFlexLexer::yylex() to keep the linker happy. It is not called
//! \cond dummy
int yyFlexLexer::yylex()
{
    FatalErrorIn("yyFlexLexer::yylex()")
        << "Should not have called this function"
        << abort(FatalError);
    return 0;
}
//! \endcond

// Dummy yywrap to keep yylex happy at compile time.
// It is called by yylex but is not used as the mechanism to change file.
// See <<EOF>>
//! \cond dummy
#if YY_FLEX_SUBMINOR_VERSION < 34
extern "C" int yywrap()
#else
int yyFlexLexer::yywrap()
#endif
{
    return 1;
}
//! \endcond


//- A lexer for parsing STL ASCII files.
//  Returns DynamicList(s) of points and facets (zoneIds).
//  The facets are within a solid/endsolid grouping
class STLASCIILexer
:
    public yyFlexLexer
{
    // Private data

        bool  sorted_;
        label groupID_;      // current solid group
        label lineNo_;
        word  startError_;

        DynamicList<point> points_;
        DynamicList<label> facets_;
        DynamicList<word>  names_;
        DynamicList<label> sizes_;
        HashTable<label>   lookup_;

public:

    // Constructors

        //- From input stream and the approximate number of vertices in the STL
        STLASCIILexer(istream* is, const label approxNpoints);


    // Member Functions

        //- The lexer function itself
        int lex();

    // Access

        //- Do all the solid groups appear in order
        bool sorted() const
        {
            return sorted_;
        }

        //- A list of points corresponding to a pointField
        DynamicList<point>& points()
        {
            return points_;
        }

        //- A list of facet IDs (group IDs)
        //  corresponds to the number of triangles
        DynamicList<label>& facets()
        {
            return facets_;
        }

        //- Names
        DynamicList<word>& names()
        {
            return names_;
        }

        //- Sizes
        DynamicList<label>& sizes()
        {
            return sizes_;
        }
};


STLASCIILexer::STLASCIILexer(istream* is, const label approxNpoints)
:
    yyFlexLexer(is),
    sorted_(true),
    groupID_(-1),
    lineNo_(1),
    points_(approxNpoints),
    facets_(approxNpoints)
{}


 /* ------------------------------------------------------------------------ *\
   ------ cppLexer::yylex()
 \* ------------------------------------------------------------------------ */

#define YY_DECL int STLASCIILexer::lex()

%}

one_space             [ \t\f\r]
space                 {one_space}*
some_space            {one_space}+

alpha                 [_A-Za-z]
digit                 [0-9]

integer               {digit}+
signedInteger         [-+]?{integer}

word                  ([[:alnum:]]|[[:punct:]])*
string                {word}({some_space}{word})*

exponent_part         [eE][-+]?{digit}+
fractional_constant   [-+]?(({digit}*"."{digit}+)|({digit}+"."?))

floatNum              (({fractional_constant}{exponent_part}?)|({digit}+{exponent_part}))

x                     {floatNum}
y                     {floatNum}
z                     {floatNum}

solid                 {space}("solid"|"SOLID"){space}
color                 {space}("color"|"COLOR"){some_space}{floatNum}{some_space}{floatNum}{some_space}{floatNum}{space}
facet                 {space}("facet"|"FACET"){space}
normal                {space}("normal"|"NORMAL"){space}
point                 {space}{x}{some_space}{y}{some_space}{z}{space}
outerloop             {space}("outer"{some_space}"loop")|("OUTER"{some_space}"LOOP"){space}
vertex                {space}("vertex"|"VERTEX"){space}
endloop               {space}("endloop"|"ENDLOOP"){space}
endfacet              {space}("endfacet"|"ENDFACET"){space}
endsolid              {space}("endsolid"|"ENDSOLID")({some_space}{word})*


 /* ------------------------------------------------------------------------ *\
                      -----  Exclusive start states -----
 \* ------------------------------------------------------------------------ */

%option stack

%x readSolidName
%x readFacet
%x readNormal
%x readVertices
%x readVertex
%x stlError

%%

%{
    // End of read character pointer returned by strtof
    // char* endPtr;

    STLpoint normal;
    STLpoint vertex;
    label cmpt = 0;   // component index used for reading vertex

    static const char* stateNames[7] =
    {
        "reading solid",
        "reading solid name",
        "reading facet",
        "reading normal",
        "reading vertices",
        "reading vertex",
        "error"
    };

    static const char* stateExpects[7] =
    {
        "'solid', 'color', 'facet' or 'endsolid'",
        "<string>",
        "'normal', 'outer loop' or 'endfacet'",
        "<float> <float> <float>",
        "'vertex' or 'endloop'",
        "<float> <float> <float>",
        ""
    };
%}


 /* ------------------------------------------------------------------------ *\
                            ------ Start Lexing ------
 \* ------------------------------------------------------------------------ */

 /*                      ------ Reading control header ------                */

{solid} {
        BEGIN(readSolidName);
    }

<readSolidName>{string} {
        word name(Foam::string::validate<word>(YYText()));

        HashTable<label>::const_iterator fnd = lookup_.find(name);
        if (fnd != lookup_.end())
        {
            if (groupID_ != fnd())
            {
                // group appeared out of order
                sorted_ = false;
            }
            groupID_ = fnd();
        }
        else
        {
            groupID_ = sizes_.size();
            lookup_.insert(name, groupID_);
            names_.append(name);
            sizes_.append(0);
        }
        BEGIN(INITIAL);
    }

<readSolidName>{space}\n {
        word name("solid");

        HashTable<label>::const_iterator fnd = lookup_.find(name);
        if (fnd != lookup_.end())
        {
            if (groupID_ != fnd())
            {
                // group appeared out of order
                sorted_ = false;
            }
            groupID_ = fnd();
        }
        else
        {
            groupID_ = sizes_.size();
            lookup_.insert(name, groupID_);
            names_.append(name);
            sizes_.append(0);
        }

        lineNo_++;
        BEGIN(INITIAL);
    }

{color} {
    }

{facet} {
        BEGIN(readFacet);
    }

<readFacet>{normal} {
        BEGIN(readNormal);
    }

<readNormal>{point} {
        /*
         skip reading normals:
         normal.x() = strtof(YYText(), &endPtr);
         normal.y() = strtof(endPtr, &endPtr);
         normal.z() = strtof(endPtr, &endPtr);
         normals_.append(normal);
         */
        BEGIN(readFacet);
    }

<readFacet>{outerloop} {
        BEGIN(readVertices);
    }

<readVertices>{vertex} {
        BEGIN(readVertex);
    }

<readVertex>{space}{signedInteger}{space} {
        vertex[cmpt++] = atol(YYText());

        if (cmpt == 3)
        {
            cmpt = 0;
            points_.append(vertex);
            BEGIN(readVertices);
        }
    }

<readVertex>{space}{floatNum}{space} {
        vertex[cmpt++] = atof(YYText());

        if (cmpt == 3)
        {
            cmpt = 0;
            points_.append(vertex);
            BEGIN(readVertices);
        }
    }

<readVertices>{endloop} {
        BEGIN(readFacet);
    }

<readFacet>{endfacet} {
        facets_.append(groupID_);
        sizes_[groupID_]++;
        BEGIN(INITIAL);
    }

{endsolid} {
    }


 /* ------------------ Ignore remaining space and \n s. -------------------- */

<*>{space} {}
<*>\n      { lineNo_++; }


 /* ------------------- Any other characters are errors -------------------- */

<*>. {
        startError_ = YYText();
        yy_push_state(stlError);
    }


 /* ---------------------------- Error handler ----------------------------- */

<stlError>.* {
        yy_pop_state();
        FatalErrorIn
        (
            "fileFormats::STLsurfaceFormatCore::readASCII(const fileName&)"
        )   << "while " << stateNames[YY_START] << " on line " << lineNo_ << nl
            << "    expected " << stateExpects[YY_START]
            << " but found '" << startError_.c_str() << YYText() << "'"
            << exit(FatalError);
    }


 /*  ------------------------ On EOF terminate ----------------------------  */

<<EOF>> {
            yyterminate();
    }
%%


//
// member function
//
bool Foam::fileFormats::STLsurfaceFormatCore::readASCII
(
    istream& is,
    const off_t dataFileSize
)
{
    // Create the lexer with the approximate number of vertices in the STL
    // from the file size
    STLASCIILexer lexer(&is, dataFileSize/400);
    while (lexer.lex() != 0) {}

    sorted_ = lexer.sorted();

    // transfer to normal lists
    points_.transfer(lexer.points());
    zoneIds_.transfer(lexer.facets());
    names_.transfer(lexer.names());
    sizes_.transfer(lexer.sizes());

    return true;
}

 /* ------------------------------------------------------------------------ *\
    ------ End of STLfileFormatASCII.L
 \* ------------------------------------------------------------------------ */
