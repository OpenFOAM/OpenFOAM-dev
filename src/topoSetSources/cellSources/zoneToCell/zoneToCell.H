/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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
    Foam::zoneToCell

Description
    A topoSetSource to select faces based on cellZone.

SourceFiles
    zoneToCell.C

\*---------------------------------------------------------------------------*/

#ifndef zoneToCell_H
#define zoneToCell_H

#include "topoSetSource.H"
#include "wordRe.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class zoneToCell Declaration
\*---------------------------------------------------------------------------*/

class zoneToCell
:
    public topoSetSource
{
    // Private Data

        //- Name/regular expression of cellZone
        wordRe zoneName_;


    // Private Member Functions

        void combine(topoSet& set, const bool add) const;


public:

    //- Runtime type information
    TypeName("zoneToCell");


    // Constructors

        //- Construct from components
        zoneToCell
        (
            const polyMesh& mesh,
            const word& zoneName
        );

        //- Construct from dictionary
        zoneToCell
        (
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~zoneToCell();


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
