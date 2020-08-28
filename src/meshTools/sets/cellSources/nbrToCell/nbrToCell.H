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
    Foam::nbrToCell

Description
    A topoSetSource to select cells based on number of neighbouring cells
    (i.e. number of internal or coupled faces)

SourceFiles
    nbrToCell.C

\*---------------------------------------------------------------------------*/

#ifndef nbrToCell_H
#define nbrToCell_H

#include "topoSetSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class nbrToCell Declaration
\*---------------------------------------------------------------------------*/

class nbrToCell
:
    public topoSetSource
{

    // Private Data

        //- Add usage string
        static addToUsageTable usage_;

        //- Number of internal faces on cell
        label minNbrs_;


    // Private Member Functions

        void combine(topoSet& set, const bool add) const;


public:

    //- Runtime type information
    TypeName("nbrToCell");

    // Constructors

        //- Construct from components
        nbrToCell
        (
            const polyMesh& mesh,
            const label minNbrs
        );

        //- Construct from dictionary
        nbrToCell
        (
            const polyMesh& mesh,
            const dictionary& dict
        );

        //- Construct from Istream
        nbrToCell
        (
            const polyMesh& mesh,
            Istream&
        );


    //- Destructor
    virtual ~nbrToCell();


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
