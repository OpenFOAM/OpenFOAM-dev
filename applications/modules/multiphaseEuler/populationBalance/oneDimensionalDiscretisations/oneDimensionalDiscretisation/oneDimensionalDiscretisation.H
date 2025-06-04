/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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
    Foam::oneDimensionalDiscretisation

Description
    Base class for ways in which to generate coordinates in one-dimensional
    space. Used to space the representative spherical diameters of the groups
    in a population balance model.

SourceFiles
    oneDimensionalDiscretisation.C

\*---------------------------------------------------------------------------*/

#ifndef oneDimensionalDiscretisation_H
#define oneDimensionalDiscretisation_H

#include "runTimeSelectionTables.H"
#include "dimensionedScalar.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class oneDimensionalDiscretisation Declaration
\*---------------------------------------------------------------------------*/

class oneDimensionalDiscretisation
{
private:

    // Private Member Data

        //- Name
        const word name_;

        //- Dimensions
        const dimensionSet dims_;

        //- Coordinates
        const scalarField coordinates_;


public:

    //- Runtime type information
    TypeName("oneDimensionalDiscretisation");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            oneDimensionalDiscretisation,
            dictionary,
            (
                const word& name,
                const dimensionSet& dims,
                const label n,
                const dictionary& dict
            ),
            (name, dims, n, dict)
        );


    // Constructors

        //- Construct from components
        oneDimensionalDiscretisation
        (
            const word& name,
            const dimensionSet& dims,
            const tmp<scalarField>& coordinates
        );


    //- Selector
    static autoPtr<oneDimensionalDiscretisation> New
    (
        const word& name,
        const dimensionSet& dims,
        const label n,
        const dictionary& dict
    );


    //- Destructor
    virtual ~oneDimensionalDiscretisation();


    // Member Operators

        //- Return the coordinates
        inline const scalarField& coordinates() const;

        //- Return the coordinates as a list of dimensioned scalars
        PtrList<dimensionedScalar> dimensionedCoordinates() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "oneDimensionalDiscretisationI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
