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
    Foam::shapeToCell

Description
    A topoSetSource to select cells based on cell shape.

    Handles all ones from cellModeller and splitHex with 10 degrees
    feature angle.

SourceFiles
    shapeToCell.C

\*---------------------------------------------------------------------------*/

#ifndef shapeToCell_H
#define shapeToCell_H

#include "topoSetSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class shapeToCell Declaration
\*---------------------------------------------------------------------------*/

class shapeToCell
:
    public topoSetSource
{

    // Private Data

        //- Add usage string
        static addToUsageTable usage_;

        //- Name of cell type
        word type_;

    // Private Member Functions

        //- Depending on cell type add to or delete from cellSet.
        void combine(topoSet& set, const bool add) const;


public:

    //- Runtime type information
    TypeName("shapeToCell");


    // Static data

        //- Cos of feature angle for polyHedral to be splitHex
        static scalar featureCos;

    // Constructors

        //- Construct from components
        shapeToCell
        (
            const polyMesh& mesh,
            const word& type
        );

        //- Construct from dictionary
        shapeToCell
        (
            const polyMesh& mesh,
            const dictionary& dict
        );

        //- Construct from Istream
        shapeToCell
        (
            const polyMesh& mesh,
            Istream&
        );


    //- Destructor
    virtual ~shapeToCell();


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
