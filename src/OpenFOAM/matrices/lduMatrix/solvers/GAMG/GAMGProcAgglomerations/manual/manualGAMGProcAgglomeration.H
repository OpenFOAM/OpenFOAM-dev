/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2025 OpenFOAM Foundation
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
    Foam::manualGAMGProcAgglomeration

Description
    Manual processor agglomerator

    In the GAMG control dictionary:

Usage
    Example of manual processor agglomeration
    \verbatim
    processorAgglomeration
    {
        agglomerator    manual;

        // List of level + list of processor clusters
        //
        // Each cluster of processor is specified as a list of IDs which are
        // merged onto the first element in the list.

        processorAgglomeration
        (
            (
                3           // at level 3
                (
                    (0 1)   // coarse 0 from 0,1 (and moved onto 0)
                    (3 2)   // coarse 1 from 2,3 (and moved onto 3)
                )
            )
            (
                6           // at level6
                (
                    (0 1)   // coarse 0 from 0,1 (and moved onto 0)
                )
            )
        );
    }
    \endverbatim

SourceFiles
    manualGAMGProcAgglomeration.C

\*---------------------------------------------------------------------------*/

#ifndef manualGAMGProcAgglomeration_H
#define manualGAMGProcAgglomeration_H

#include "GAMGProcAgglomeration.H"
#include "DynamicList.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class GAMGAgglomeration;

/*---------------------------------------------------------------------------*\
                 Class manualGAMGProcAgglomeration Declaration
\*---------------------------------------------------------------------------*/

class manualGAMGProcAgglomeration
:
    public GAMGProcAgglomeration
{
    // Private Data

        //- Per level the agglomeration map
        const List<Tuple2<label, List<labelList>>> procAgglomMaps_;

        //- Any allocated communicators
        DynamicList<label> comms_;


public:

    //- Runtime type information
    TypeName("manual");


    // Constructors

        //- Construct given agglomerator and controls
        manualGAMGProcAgglomeration
        (
            GAMGAgglomeration& agglom,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        manualGAMGProcAgglomeration
        (
            const manualGAMGProcAgglomeration&
        ) = delete;


    //- Destructor
    virtual ~manualGAMGProcAgglomeration();


    // Member Functions

       //- Modify agglomeration. Return true if modified
        virtual bool agglomerate();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const manualGAMGProcAgglomeration&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
