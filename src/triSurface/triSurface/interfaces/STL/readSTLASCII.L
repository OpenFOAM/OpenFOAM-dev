/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "IFstream.H"
#include "triSurface.H"
#include "floatVector.H"
#include "OSspecific.H"
#include "mergePoints.H"
//#include "memInfo.H"

using namespace Foam;

// Dummy yyFlexLexer::yylex() to keep the linker happy. It is not called
//! \cond dummy
int yyFlexLexer::yylex()
{
    FatalErrorInFunction
        << "Should not have called this function"
        << abort(FatalError);
    return 0;
}
//! \endcond

// Dummy yywrap to keep yylex happy at compile time.
// It is called by yylex but is not used as the mechanism to change file.
// See <<EOF>>
//! \cond dummy
#if YY_FLEX_MINOR_VERSION < 6 && YY_FLEX_SUBMINOR_VERSION < 34
extern "C" int yywrap()
#else
int yyFlexLexer::yywrap()
#endif
{
    return 1;
}
//! \endcond


class STLLexer
:
    public yyFlexLexer
{
    // Private data

        label nTriangles_;
        short region_;       // current region
        short maxRegion_;    // max region

        label lineNo_;
        word startError_;

        DynamicList<floatVector> STLpoints_;
        // DynamicList<floatVector > STLnormals_;
        DynamicList<label> STLlabels_;
        HashTable<label, word> STLsolidNames_;


public:

    // Constructors

        //- From input stream and the approximate number of vertices in the STL
        STLLexer(istream* is, const label approxNpoints);


    // Member Functions

        //- The lexer function itself
        int lex();

    // Access

        label nTriangles() const
        {
            return nTriangles_;
        }

        DynamicList<floatVector>& STLpoints()
        {
            return STLpoints_;
        }

        // DynamicList<floatVector>& STLnormals()
        //{
        //    return STLnormals_;
        //}

        DynamicList<label>& STLlabels()
        {
            return STLlabels_;
        }

        const HashTable<label, word>& STLsolidNames() const
        {
            return STLsolidNames_;
        }
};


STLLexer::STLLexer(istream* is, const label approxNpoints)
:
    yyFlexLexer(is),

    nTriangles_(0),
    region_(-1),
    maxRegion_(0),

    lineNo_(1),

    STLpoints_(approxNpoints),
    // STLnormals_(approxNpoints),
    STLlabels_(approxNpoints)
{}


 /* ------------------------------------------------------------------------ *\
   ------ cppLexer::yylex()
 \* ------------------------------------------------------------------------ */

#define YY_DECL int STLLexer::lex()

%}

one_space             [ \t\f\r]
space                 {one_space}*
some_space            {one_space}+

alpha                 [_A-Za-z]
digit                 [0-9]

integer               {digit}+
signedInteger         [-+]?{digit}+

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
%x stlerror

%%

%{
    // End of read character pointer returned by strtof
    // char* endPtr;

    floatVector normal;
    floatVector vertex;
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
        word solidName(Foam::string::validate<word>(YYText()));
        if (STLsolidNames_.found(solidName))
        {
            region_ = STLsolidNames_[solidName];
        }
        else
        {
            region_ = maxRegion_++;
            STLsolidNames_.insert(solidName, region_);
        }
        BEGIN(INITIAL);
    }

<readSolidName>{space}\n {
        word solidName("solid");
        if (STLsolidNames_.found(solidName))
        {
            region_ = STLsolidNames_[solidName];
        }
        else
        {
            region_ = maxRegion_++;
            STLsolidNames_.insert(solidName, region_);
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
        normal.x() = strtof(YYText(), &endPtr);
        normal.y() = strtof(endPtr, &endPtr);
        normal.z() = strtof(endPtr, &endPtr);
        STLnormals_.append(normal);
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
            STLpoints_.append(vertex);
            BEGIN(readVertices);
        }
    }

<readVertex>{space}{floatNum}{space} {
        vertex[cmpt++] = atof(YYText());

        if (cmpt == 3)
        {
            cmpt = 0;
            STLpoints_.append(vertex);
            BEGIN(readVertices);
        }
    }

<readVertices>{endloop} {
        BEGIN(readFacet);
    }

<readFacet>{endfacet} {
        nTriangles_++;
        STLlabels_.append(region_);
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
        yy_push_state(stlerror);
    }


 /* ---------------------------- Error handler ----------------------------- */

<stlerror>.* {
        yy_pop_state();
        FatalErrorInFunction
            << "while " << stateNames[YY_START] << " on line " << lineNo_ << nl
            << "    expected " << stateExpects[YY_START]
            << " but found '" << startError_.c_str() << YYText() << "'"
            << exit(FatalError);
    }


 /*  ------------------------ On EOF terminate ----------------------------  */

<<EOF>> {
            yyterminate();
    }
%%


#include <fstream>

bool triSurface::readSTLASCII(const fileName& STLfileName)
{
    IFstream STLstream(STLfileName);

    if (!STLstream)
    {
        FatalErrorInFunction
            << "file " << STLfileName << " not found"
            << exit(FatalError);
    }

    // memInfo memStat;

    // memStat.update();
    // Pout<< "At start:" << memStat.rss() << endl;

    // Create the lexer obtaining the approximate number of vertices in the STL
    // from the file size
    STLLexer lexer(&STLstream.stdStream(), Foam::fileSize(STLfileName)/400);
    while (lexer.lex() != 0)
    {}

    // memStat.update();
    // Pout<< "After lexing:" << memStat.rss() << endl;

    DynamicList<floatVector>& STLpoints = lexer.STLpoints();
    DynamicList<label>& STLlabels = lexer.STLlabels();

    /*
    DynamicList<floatVector>& STLnormals = lexer.STLnormals();

    if (STLpoints.size() != 3*STLnormals.size())
    {
        FatalErrorInFunction
            << "in file " << STLfileName << endl
            << "Problem: read " << STLnormals.size() << " normals"
            << " but " << STLpoints.size() << " points"
            << exit(FatalError);
    }
    */

    labelList pointMap;
    label nUniquePoints = mergePoints
    (
        STLpoints,
        100*SMALL,              // merge distance
        false,                  // verbose
        pointMap
    );

    // memStat.update();
    // Pout<< "After merging:" << memStat.rss() << endl;


    pointField& sp = storedPoints();

    setSize(lexer.nTriangles());
    sp.setSize(nUniquePoints);
    forAll(STLpoints, pointi)
    {
        const floatVector& pt = STLpoints[pointi];
        sp[pointMap[pointi]] = vector
        (
            scalar(pt.x()),
            scalar(pt.y()),
            scalar(pt.z())
        );
    }

    // Assign triangles
    label pointi = 0;
    forAll(*this, i)
    {
        operator[](i)[0] = pointMap[pointi++];
        operator[](i)[1] = pointMap[pointi++];
        operator[](i)[2] = pointMap[pointi++];
        operator[](i).region() = STLlabels[i];
    }

    // memStat.update();
    // Pout<< "After assigning:" << memStat.rss() << endl;


    STLpoints.clear();
    STLlabels.clear();


    // Convert solidNames into regionNames
    patches_.setSize(lexer.STLsolidNames().size());

    for
    (
        HashTable<label, word>::const_iterator iter =
            lexer.STLsolidNames().begin();
        iter != lexer.STLsolidNames().end();
        ++iter
    )
    {
        patches_[iter()].name() = iter.key();
    }

    // Fill in the missing information in the patches
    setDefaultPatches();

    // memStat.update();
    // Pout<< "After patchifying:" << memStat.rss() << endl;

    return true;
}


 /* ------------------------------------------------------------------------ *\
    ------ End of STLToFoam.L
 \* ------------------------------------------------------------------------ */
