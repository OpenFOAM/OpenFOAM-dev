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
    Foam::nearestToPoint

Description
    A topoSetSource to select points nearest to points.

SourceFiles
    nearestToPoint.C

\*---------------------------------------------------------------------------*/

#ifndef nearestToPoint_H
#define nearestToPoint_H

#include "topoSetSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class nearestToPoint Declaration
\*---------------------------------------------------------------------------*/

class nearestToPoint
:
    public topoSetSource
{

    // Private Data

        //- Add usage string
        static addToUsageTable usage_;

        //- Points to select nearest to
        pointField points_;


    // Private Member Functions

        void combine(topoSet& set, const bool add) const;


public:

    //- Runtime type information
    TypeName("nearestToPoint");

    // Constructors

        //- Construct from components
        nearestToPoint
        (
            const polyMesh& mesh,
            const pointField& points
        );

        //- Construct from dictionary
        nearestToPoint
        (
            const polyMesh& mesh,
            const dictionary& dict
        );

        //- Construct from Istream
        nearestToPoint
        (
            const polyMesh& mesh,
            Istream&
        );


    //- Destructor
    virtual ~nearestToPoint();


    // Member Functions

        virtual sourceType setType() const
        {
            return POINTSETSOURCE;
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
