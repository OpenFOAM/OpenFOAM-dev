/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2025 OpenFOAM Foundation
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
    Foam::populationBalance::shapeModels::spherical

Description
    Spherical shape model

SourceFiles
    spherical.C

\*---------------------------------------------------------------------------*/

#ifndef spherical_H
#define spherical_H

#include "shapeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalance
{
namespace shapeModels
{

/*---------------------------------------------------------------------------*\
                          Class spherical Declaration
\*---------------------------------------------------------------------------*/

class spherical
:
    public shapeModel
{
public:

    //- Runtime type information
    TypeName("spherical");


    // Constructors

        //- Construct from a dictionary
        spherical
        (
            const dictionary& dict,
            const populationBalanceModel& popBal
        );


    //- Destructor
    virtual ~spherical();


    // Member Functions

        //- Return the representative surface area of a group
        virtual tmp<volScalarField> a(const label i) const;

        //- Return the representative diameter of a group
        virtual tmp<volScalarField> d(const label i) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace shapeModels
} // End namespace populationBalance
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
