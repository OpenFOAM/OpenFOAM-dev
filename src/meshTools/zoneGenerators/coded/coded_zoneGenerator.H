/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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
    Foam::zoneGenerators::coded

Description
    Coded zoneGenerator

    The code provided should set \c pointIndices and/or \c cellIndices
    and/or \c faceIndices and optionally \c faceMap.  These index lists
    are automatically converted into the corresponding zones in the \c zoneSet
    which is returned by the \c zoneGenerator.

    The entries are:
    \plaintable
       code        | Code to set pointIndices, cellIndices, faceIndices
       codeInclude | Include files
       codeOptions | Include paths; inserted into EXE_INC in Make/options
       codeLibs    | Link line; inserted into LIB_LIBS in Make/options
    \endplaintable

    Example of a zoneGenerators::coded which creates a \c cellZone of the cells
    in which the value of the turbulent kinetic energy is within the highest
    10%:
    \verbatim
    cellZone
    {
        type        coded;

        // Regenerate cellZone for every evaluation
        regenerate  yes;

        code
        #{
            // Lookup field from the objectRegistry, e.g. k
            const volScalarField& k = mesh().lookupObject<volScalarField>("k");

            // Calculate the minimum and maximum k values
            const scalar mink = min(k).value();
            const scalar maxk = max(k).value();

            // Fraction of the range of k
            const scalar f = 0.1;

            // Lowest k in selection
            const scalar lowk = f*mink + (1.0 - f)*maxk;

            cellIndices = selectIndices
            (
                k,
                [&](const scalar k)
                {
                    return k > lowk;
                }
            );
        #};
    }
    \endverbatim

See also
    Foam::zoneGenerator
    Foam::codedBase

SourceFiles
    coded_zoneGenerator.C

\*---------------------------------------------------------------------------*/

#ifndef coded_zoneGenerator_H
#define coded_zoneGenerator_H

#include "zoneGenerator.H"
#include "codedBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace zoneGenerators
{

/*---------------------------------------------------------------------------*\
                         Class coded Declaration
\*---------------------------------------------------------------------------*/

class coded
:
    public zoneGenerator,
    public codedBase
{
    // Private Data

        //- Keywords associated with source code
        static const wordList codeKeys;

        //- Name of the dictionary variables in the source code
        static const wordList codeDictVars;

        // Name of the code options template file
        static const word codeOptions;

        // Names of the source files to compile
        static const wordList compileFiles;

        // Names of the source files to copy (header files)
        static const wordList copyFiles;

        //- The dynamically generated zoneGenerator pointer
        mutable autoPtr<zoneGenerator> redirectZoneGeneratorPtr_;


public:

    //- Runtime type information
    TypeName("coded");


    // Constructors

        //- Construct from name, polyMesh and dictionary
        coded
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        coded(const coded&) = delete;


    //- Destructor
    virtual ~coded();


    // Member Functions

        //- Generate and return the zoneSet
        virtual zoneSet generate() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const coded&) = delete;
};


//- Helper struct to enable the zone type for the corresponding list of indices
//  which allows for the generation of empty zones
struct enabledLabelList
:
    public labelList
{
    bool enabled_;

    enabledLabelList()
    :
        enabled_(false)
    {}

    void operator=(const labelList& indices)
    {
        labelList::operator=(indices);
        enabled_ = true;
    }

    bool enabled() const
    {
        return enabled_ || size();
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace zoneGenerators
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
