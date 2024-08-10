/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2024 OpenFOAM Foundation
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
    Foam::sampledSets::cellSet

Description
    Samples at the cell-centres of a given cell set

Usage
    \table
        Property    | Description                            | Req'd? | Default
        set         | Name of the cell set                   | yes    |
    \endtable

    Example specification:
    \verbatim
    {
        type        cellSet;
        set         c0;
    }
    \endverbatim

SourceFiles
    cellSet.C

\*---------------------------------------------------------------------------*/

#ifndef cellSet_sampledSet_H
#define cellSet_sampledSet_H

#include "sampledSet.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace sampledSets
{

/*---------------------------------------------------------------------------*\
                      Class cellSet Declaration
\*---------------------------------------------------------------------------*/

class cellSet
:
    public sampledSet
{
    // Private Data

        //- Name of the cell set
        const word setName_;


    // Private Member Functions

        //- Generate the samples
        virtual void genSamples();


public:

    //- Runtime type information
    TypeName("cellSet");


    // Constructors

        //- Construct from dictionary
        cellSet
        (
            const word& name,
            const polyMesh& mesh,
            const meshSearch& searchEngine,
            const dictionary& dict
        );


    //- Destructor
    virtual ~cellSet();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace sampledSets
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
