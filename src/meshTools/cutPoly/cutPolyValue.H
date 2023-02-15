/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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

Namespace
    Foam::cutPoly

Description
    Functions and classes for extracting values from cut edges, faces and cells

SourceFiles
    cutPolyValueI.H
    cutPolyValueTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef cutPolyValue_H
#define cutPolyValue_H

#include "cutPoly.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace cutPoly
{

//- Get the local coordinate within an edge, given end point values and an
//  iso-value
inline scalar edgeCutLambda
(
    const edge& e,
    const scalarField& pAlphas,
    const scalar isoAlpha
);

//- Linearly interpolate a value from the end points to the cut point of an
//  edge, given a local coordinate within the edge
template<class Type>
Type edgeCutValue
(
    const edge& e,
    const scalar lambda,
    const Field<Type>& pPsis
);

//- Linearly interpolate a value from the end points to the cut point of an
//  edge, given end point values and an iso-value
template<class Type>
Type edgeCutValue
(
    const edge& e,
    const scalarField& pAlphas,
    const scalar isoAlpha,
    const Field<Type>& pPsis
);


/*---------------------------------------------------------------------------*\
                       Class FaceCutValues Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class FaceCutValues
{
    // Private Data

        //- The face
        const face& f_;

        //- The face's cuts
        const List<labelPair>& fCuts_;

        //- The point values that we want to iterate over
        const Field<Type>& pPsis_;

        //- The point values that define the iso-surface
        const scalarField& pAlphas_;

        //- The value that defines the iso-surface
        const scalar isoAlpha_;

        //- Do we want to iterate over the cut below or above the iso surface?
        const bool below_;


public:

    // Public Typedefs

        //- The value type
        typedef Type value_type;


    // Public Classes

        //- Forward iterator
        class const_iterator
        {
            // Private Data

                //- Reference to the face-values
                const FaceCutValues<Type>& fValues_;

                //- Index of the sub-face and sub-face-point
                label i_, j_;

                //- Cache of values on cut edges for this sub-face segment
                mutable Tuple2<label, Pair<Type>> iAndCutPsis_;


            // Private Member Functions

                //- Get the values on the cut edges for this sub-face segment
                const Pair<Type>& cutPsis(const label i) const;

                //- Return the size of the given sub-face
                label size(const label i) const;

                //- Return the value for the given sub-face and sub-face-point
                const Type psi(const label i, const label j) const;


        public:

            // Constructors

                //- Construct from components
                const_iterator
                (
                    const FaceCutValues<Type>& fValues,
                    const label i,
                    const label j
                );


            // Member Functions

                //- Get the next value around the sub-face
                Type next() const;


            // Member Operators

                //- Equality comparison
                bool operator==(const const_iterator& it) const;

                //- Inequality comparison
                bool operator!=(const const_iterator& it) const;

                //- Dereference
                Type operator*() const;

                //- Increment
                inline const_iterator& operator++();

                //- Increment
                inline const_iterator operator++(int);
        };


    // Constructors

        //- Construct from components
        FaceCutValues
        (
            const face& f,
            const List<labelPair>& fCuts,
            const Field<Type>& pPsis,
            const scalarField& pAlphas,
            const scalar isoAlpha,
            const bool below
        );


    // Member Functions

        //- Get the beginning of the iteration over the values
        const_iterator begin() const;

        //- Get the end of the iteration over the values
        const_iterator end() const;
};


/*---------------------------------------------------------------------------*\
                       Class CellCutValues Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class CellCutValues
{
    // Private Data

        //- The cell
        const cell& c_;

        //- The cell's edge addressing
        const cellEdgeAddressing& cAddr_;

        //- The cell's cuts
        const labelListList& cCuts_;

        //- The faces
        const faceList& fs_;

        //- The point values that we want to iterate over
        const Field<Type>& pPsis_;

        //- The point values that define the iso-surface
        const scalarField& pAlphas_;

        //- The value that defines the iso-surface
        const scalar isoAlpha_;


public:

    // Public Typedefs

        //- The value type
        typedef Type value_type;


    // Public Classes

        //- Forward iterator
        class const_iterator
        {
            // Private Data

                //- Reference to the cell-values
                const CellCutValues<Type>& cValues_;

                //- Index of the loop and loop-point
                label i_, j_;


            // Private Member Functions

                //- The values at this loop-point and the next
                Pair<Type> psis_;

                //- Return the value for the given loop and loop-point
                Type psi(const label i, const label j);


        public:

            // Constructors

                //- Construct from components
                const_iterator
                (
                    const CellCutValues<Type>& cValues,
                    const label i,
                    const label j
                );


            // Member Functions

                //- Get the next value around the loop
                const Type& next() const;


            // Member Operators

                //- Equality comparison
                bool operator==(const const_iterator& it) const;

                //- Inequality comparison
                bool operator!=(const const_iterator& it) const;

                //- Dereference
                const Type& operator*() const;

                //- Increment
                inline const_iterator& operator++();

                //- Increment
                inline const_iterator operator++(int);
        };


    // Constructors

        //- Construct from components
        CellCutValues
        (
            const cell& c,
            const cellEdgeAddressing& cAddr,
            const labelListList& cCuts,
            const faceList& fs,
            const Field<Type>& pPsis,
            const scalarField& pAlphas,
            const scalar isoAlpha
        );


    // Member Functions

        //- Get the beginning of the iteration over the values
        const_iterator begin() const;

        //- Get the end of the iteration over the values
        const_iterator end() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace cutPoly
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cutPolyValueI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "cutPolyValueTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
