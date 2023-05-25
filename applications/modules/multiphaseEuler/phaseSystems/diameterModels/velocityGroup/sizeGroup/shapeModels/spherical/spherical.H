/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2023 OpenFOAM Foundation
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
    Foam::diameterModels::shapeModels::spherical

Description
    Uniform spherical shape.

SourceFiles
    spherical.C

\*---------------------------------------------------------------------------*/

#ifndef spherical_H
#define spherical_H

#include "shapeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
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

        //- Construct from dictionary and sizeGroup
        spherical
        (
            const dictionary& dict,
            const sizeGroup& group
        );


    //- Destructor
    virtual ~spherical();


    // Member Functions

        //- Return representative surface area of the sizeGroup
        virtual const tmp<volScalarField> a() const;

        //- Return representative diameter of the sizeGroup
        virtual const tmp<volScalarField> d() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace shapeModels
} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
