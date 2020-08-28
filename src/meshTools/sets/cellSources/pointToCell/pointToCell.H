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
    Foam::pointToCell

Description
    A topoSetSource to select cells based on usage of points.

SourceFiles
    pointToCell.C

\*---------------------------------------------------------------------------*/

#ifndef pointToCell_H
#define pointToCell_H

#include "topoSetSource.H"
#include "NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class pointToCell Declaration
\*---------------------------------------------------------------------------*/

class pointToCell
:
    public topoSetSource
{
public:
        //- Enumeration defining the valid options
        enum pointAction
        {
            ANY,    // Cells using any point in set
            EDGE    // Cells using an edge with both points in set
            // ALL   // Possible extension: cells whose all points are in set
        };

private:

        //- Add usage string
        static addToUsageTable usage_;

        static const NamedEnum<pointAction, 2> pointActionNames_;

        //- Name of set to use
        word setName_;

        //- Option
        pointAction option_;


    // Private Member Functions

        //- Depending on point-to-cell option add to or delete from cellSet.
        void combine(topoSet& set, const bool add) const;


public:

    //- Runtime type information
    TypeName("pointToCell");

    // Constructors

        //- Construct from components
        pointToCell
        (
            const polyMesh& mesh,
            const word& setName,
            const pointAction option
        );

        //- Construct from dictionary
        pointToCell
        (
            const polyMesh& mesh,
            const dictionary& dict
        );

        //- Construct from Istream
        pointToCell
        (
            const polyMesh& mesh,
            Istream&
        );


    //- Destructor
    virtual ~pointToCell();


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
