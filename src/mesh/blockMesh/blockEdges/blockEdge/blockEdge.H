/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::blockEdge

Description
    Define a curved edge that is parameterized for 0<lambda<1
    between the start and end point.

SourceFiles
    blockEdge.C

\*---------------------------------------------------------------------------*/

#ifndef blockEdge_H
#define blockEdge_H

#include "searchableSurfaces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class blockEdge;

Ostream& operator<<(Ostream&, const blockEdge&);


/*---------------------------------------------------------------------------*\
                          Class blockEdge Declaration
\*---------------------------------------------------------------------------*/

class blockEdge
{
protected:

    // Protected data

        const pointField& points_;

        const label start_;
        const label end_;


    // Protected Member Functions

        //- Return a complete point field by appending the start/end points
        //  to the given list
        static pointField appendEndPoints
        (
            const pointField&,
            const label start,
            const label end,
            const pointField& otherKnots
        );


public:

    //- Runtime type information
    TypeName("blockEdge");

    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            blockEdge,
            Istream,
            (
                const dictionary& dict,
                const label index,
                const searchableSurfaces& geometry,
                const pointField& points,
                Istream& is
            ),
            (dict, index, geometry, points, is)
        );


    // Constructors

        //- Construct from components
        blockEdge
        (
            const pointField& points,
            const label start,
            const label end
        );

        //- Construct from Istream setting pointsList
        blockEdge
        (
            const dictionary& dict,
            const label index,
            const pointField&,
            Istream&
        );

        //- Clone function
        virtual autoPtr<blockEdge> clone() const;

        //- New function which constructs and returns pointer to a blockEdge
        static autoPtr<blockEdge> New
        (
            const dictionary& dict,
            const label index,
            const searchableSurfaces& geometry,
            const pointField&,
            Istream&
        );

        //- Class used for the read-construction of
        //  PtrLists of blockEdge
        class iNew
        {
            const dictionary& dict_;
            const searchableSurfaces& geometry_;
            const pointField& points_;
            mutable label index_;

        public:

            iNew
            (
                const dictionary& dict,
                const searchableSurfaces& geometry,
                const pointField& points
            )
            :
                dict_(dict),
                geometry_(geometry),
                points_(points),
                index_(0)
            {}

            autoPtr<blockEdge> operator()(Istream& is) const
            {
                return blockEdge::New(dict_, index_++, geometry_, points_, is);
            }
        };


    //- Destructor
    virtual ~blockEdge()
    {}


    // Member Functions

        //- Return label of start point
        inline label start() const;

        //- Return label of end point
        inline label end() const;

        //- Compare the given start and end points with this curve
        //  Return:
        //  -  0: different
        //  - +1: identical
        //  - -1: same edge, but different orientation
        inline int compare(const blockEdge&) const;

        //- Compare the given start and end points with this curve
        //  Return:
        //  -  0: different
        //  - +1: identical
        //  - -1: same edge, but different orientation
        inline int compare(const edge&) const;

        //- Compare the given start and end points with this curve
        //  Return:
        //  -  0: different
        //  - +1: identical
        //  - -1: same edge, but different orientation
        inline int compare(const label start, const label end) const;

        //- Return the point position corresponding to the curve parameter
        //  0 <= lambda <= 1
        virtual point position(const scalar) const = 0;

        //- Return the point positions corresponding to the curve parameters
        //  0 <= lambda <= 1
        virtual tmp<pointField> position(const scalarList&) const;

        //- Return the length of the curve
        virtual scalar length() const = 0;

        //- Write edge with variable backsubstitution
        void write(Ostream&, const dictionary&) const;


    // Ostream operator

        friend Ostream& operator<<(Ostream&, const blockEdge&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "blockEdgeI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
