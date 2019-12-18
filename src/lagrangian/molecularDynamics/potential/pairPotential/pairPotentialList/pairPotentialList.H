/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::pairPotentialList

Description

SourceFiles
    pairPotentialListI.H
    pairPotentialList.C

\*---------------------------------------------------------------------------*/

#ifndef pairPotentialList_H
#define pairPotentialList_H

#include "PtrList.H"
#include "word.H"
#include "pairPotential.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class pairPotentialList Declaration
\*---------------------------------------------------------------------------*/

class pairPotentialList
:
    public PtrList<pairPotential>
{
    // Private Data

        label nIds_;

        scalar rCutMax_;

        scalar rCutMaxSqr_;

        autoPtr<pairPotential> electrostaticPotential_;


    // Private Member Functions

        inline label pairPotentialIndex
        (
            const label a,
            const label b
        ) const;

        void readPairPotentialDict
        (
            const List<word>& idList,
            const dictionary& pairPotentialDict,
            const polyMesh& mesh
        );


public:

    // Constructors

        pairPotentialList();

        //- Construct from idList and potential dictionaries
        pairPotentialList
        (
            const List<word>& idList,
            const dictionary& pairPotentialDict,
            const polyMesh& mesh
        );

        //- Disallow default bitwise assignment
        void operator=(const pairPotentialList&) = delete;


    //- Destructor
    ~pairPotentialList();


    // Member Functions

        void buildPotentials
        (
            const List<word>& idList,
            const dictionary& pairPotentialDict,
            const polyMesh& mesh
        );


        // Access

            inline scalar rCutMax() const;

            inline scalar rCutMaxSqr() const;

            bool rCutMaxSqr(scalar rIJMagSqr) const;

            const pairPotential& pairPotentialFunction
            (
                const label a,
                const label b
            ) const;

            // Return true if rIJ is within rCut for this pair.
            bool rCutSqr
            (
                const label a,
                const label b,
                const scalar rIJMagSqr
            ) const;

            scalar rMin(const label a, const label b) const;

            scalar dr(const label a, const label b) const;

            scalar rCutSqr(const label a, const label b) const;

            scalar rCut(const label a, const label b) const;

            scalar force
            (
                const label a,
                const label b,
                const scalar rIJMag
            ) const;

            scalar energy
            (
                const label a,
                const label b,
                const scalar rIJMag
            ) const;

            inline const pairPotential& electrostatic() const;


    // Member Operators

        //- Disallow default bitwise copy construction
        pairPotentialList(const pairPotentialList&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pairPotentialListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
