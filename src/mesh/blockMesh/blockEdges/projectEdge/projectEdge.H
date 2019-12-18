/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2019 OpenFOAM Foundation
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
    Foam::projectEdge

Description
    Defines the edge from the projection onto a surface (single surface)
    or intersection of two surfaces.

SourceFiles
    projectEdge.C

\*---------------------------------------------------------------------------*/

#ifndef projectEdge_H
#define projectEdge_H

#include "blockEdge.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class pointConstraint;

/*---------------------------------------------------------------------------*\
                         Class projectEdge Declaration
\*---------------------------------------------------------------------------*/

class projectEdge
:
    public blockEdge
{
    // Private Data

        const searchableSurfaces& geometry_;

        //- The indices of surfaces onto which the points are projected
        labelList surfaces_;


    // Private Member Functions

        //- Single point find nearest
        void findNearest(const point&, point& near, pointConstraint&) const;


public:

    //- Runtime type information
    TypeName("project");


    // Constructors

        //- Construct from Istream setting pointsList
        projectEdge
        (
            const dictionary& dict,
            const label index,
            const searchableSurfaces& geometry,
            const pointField& points,
            Istream&
        );

        //- Disallow default bitwise copy construction
        projectEdge(const projectEdge&) = delete;


    //- Destructor
    virtual ~projectEdge()
    {}


    // Member Functions

        //- Return the point positions corresponding to the curve parameters
        //  0 <= lambda <= 1
        virtual point position(const scalar) const;

        //- Return the point positions corresponding to the curve parameters
        //  0 <= lambda <= 1
        virtual tmp<pointField> position(const scalarList&) const;

        //- Return the length of the curve
        virtual scalar length() const
        {
            NotImplemented;
            return 1;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const projectEdge&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
