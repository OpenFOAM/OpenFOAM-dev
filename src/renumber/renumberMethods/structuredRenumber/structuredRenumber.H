/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2024 OpenFOAM Foundation
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
    Foam::structuredRenumber

Description
    Renumbering according to mesh layers.
    depthFirst = true:
        first column gets ids 0..nLayer-1,
        second nLayers..2*nLayers-1 etc.
    depthFirst = false:
        first layer gets ids 0,1,2 etc.

SourceFiles
    structuredRenumber.C

\*---------------------------------------------------------------------------*/

#ifndef structuredRenumber_H
#define structuredRenumber_H

#include "renumberMethod.H"
#include "topoDistanceData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class structuredRenumber Declaration
\*---------------------------------------------------------------------------*/

class structuredRenumber
:
    public renumberMethod
{
    // Private constructors

        //- Construct given the renumber dictionary
        structuredRenumber
        (
            const dictionary& renumberDict,
            const dictionary& methodDict
        );


public:

    // Public classes

        //- Less function class that can be used for sorting according to
        //  column and layer
        class layerLess
        {
            const Switch depthFirst_;
            const labelList& order_;
            const List<topoDistanceData>& distance_;

        public:

            layerLess
            (
                const Switch depthFirst,
                const labelList& order,
                const List<topoDistanceData>& distance
            )
            :
                depthFirst_(depthFirst),
                order_(order),
                distance_(distance)
            {}

            bool operator()(const label a, const label b);
        };


    // Private Data

        const wordReList patches_;

        const label nLayers_;

        const Switch depthFirst_;

        const autoPtr<renumberMethod> method_;

        const Switch reverse_;


public:

    //- Runtime type information
    TypeName("structured");


    // Constructors

        //- Construct given the renumber dictionary
        structuredRenumber(const dictionary& renumberDict);

        //- Disallow default bitwise copy construction
        structuredRenumber(const structuredRenumber&) = delete;


    //- Destructor
    virtual ~structuredRenumber()
    {}


    // Member Functions

        //- Return the order in which cells need to be visited, i.e.
        //  from ordered back to original cell label.
        //  This is only defined for geometric renumberMethods.
        virtual labelList renumber(const pointField&) const
        {
            NotImplemented;
            return labelList(0);
        }

        //- Return the order in which cells need to be visited, i.e.
        //  from ordered back to original cell label.
        //  Use the mesh connectivity (if needed)
        virtual labelList renumber
        (
            const polyMesh& mesh,
            const pointField& cc
        ) const;

        //- Return the order in which cells need to be visited, i.e.
        //  from ordered back to original cell label.
        //  The connectivity is equal to mesh.cellCells() except
        //  - the connections are across coupled patches
        virtual labelList renumber
        (
            const labelListList& cellCells,
            const pointField& cc
        ) const
        {
            NotImplemented;
            return labelList(0);
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const structuredRenumber&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
