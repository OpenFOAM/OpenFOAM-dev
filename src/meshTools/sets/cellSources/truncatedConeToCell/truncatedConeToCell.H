/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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
    Foam::truncatedConeToCell

Description
    A topoSetSource to select cells with centres inside a truncated cone.

SourceFiles
    truncatedConeToCell.C

\*---------------------------------------------------------------------------*/

#ifndef truncatedConeToCell_H
#define truncatedConeToCell_H

#include "topoSetSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class truncatedConeToCell Declaration
\*---------------------------------------------------------------------------*/

class truncatedConeToCell
:
    public topoSetSource
{
    // Private Data

        //- First point on cylinder axis
        vector point1_;

        //- Second point on cylinder axis
        vector point2_;

        //- Radius 1
        scalar radius1_;

        //- Radius 2
        scalar radius2_;


    // Private Member Functions

        void combine(topoSet& set, const bool add) const;


public:

    //- Runtime type information
    TypeName("truncatedConeToCell");


    // Constructors

        //- Construct from components
        truncatedConeToCell
        (
            const polyMesh& mesh,
            const vector& point1,
            const vector& point2,
            const scalar radius1,
            const scalar radius2
        );

        //- Construct from dictionary
        truncatedConeToCell
        (
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~truncatedConeToCell();


    // Member Functions

        virtual sourceType setType() const
        {
            return CELLSETSOURCE;
        }

        virtual void applyToSet
        (
            const topoSetSource::setAction action,
            topoSet&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
