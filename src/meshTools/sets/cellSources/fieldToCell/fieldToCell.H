/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
    Foam::fieldToCell

Description
    A topoSetSource to select cells based on field values.

SourceFiles
    fieldToCell.C

\*---------------------------------------------------------------------------*/

#ifndef fieldToCell_H
#define fieldToCell_H

#include "topoSetSource.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class fieldToCell Declaration
\*---------------------------------------------------------------------------*/

class fieldToCell
:
    public topoSetSource
{
    // Private Data

        //- Name of volScalarField
        word fieldName_;

        //- Min to subset field values with
        scalar min_;

        //- Max to subset field values with
        scalar max_;


    // Private Member Functions

        //- Depending on field values add to or delete from cellSet.
        void applyToSet
        (
            const topoSetSource::setAction action,
            const scalarField& field,
            topoSet& set
        ) const;


public:

    //- Runtime type information
    TypeName("fieldToCell");


    // Constructors

        //- Construct from components
        fieldToCell
        (
            const polyMesh& mesh,
            const word& fieldName,
            const scalar min,
            const scalar max
        );

        //- Construct from dictionary
        fieldToCell
        (
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~fieldToCell();


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
