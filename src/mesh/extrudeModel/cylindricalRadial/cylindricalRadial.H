/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2023 OpenFOAM Foundation
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
    Foam::extrudeModels::cylindricalRadial

Description
    Extrudes by transforming points in the cylindrical radial direction

    Generates layers at radii specified by a Foam::Function1
    of the layer index.

\*---------------------------------------------------------------------------*/

#ifndef cylindricalRadial_H
#define cylindricalRadial_H

#include "extrudeModel.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace extrudeModels
{

/*---------------------------------------------------------------------------*\
                      Class cylindricalRadial Declaration
\*---------------------------------------------------------------------------*/

class cylindricalRadial
:
    public extrudeModel
{
    // Private Data

        //- Point on axis
        const point axisPt_;

        //- Normalised axis
        vector axis_;

        //- Radial distribution of layers
        autoPtr<Function1<scalar>> R_;


public:

    //- Runtime type information
    TypeName("cylindricalRadial");

    // Constructors

        //- Construct from dictionary
        cylindricalRadial(const dictionary& dict);


    //-Destructor
    virtual ~cylindricalRadial();


    // Member Operators

        //- Return the new point corresponding to the surfacePoint
        //  on the specified layer.
        //  The surfaceNormal is not used.
        point operator()
        (
            const point& surfacePoint,
            const vector& surfaceNormal,
            const label layer
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace extrudeModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
