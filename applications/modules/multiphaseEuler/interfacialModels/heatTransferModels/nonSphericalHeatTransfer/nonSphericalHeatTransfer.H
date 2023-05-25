/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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
    Foam::heatTransferModels::nonSphericalHeatTransfer

Description
    Model which applies an analytical solution for heat transfer from the
    surface of a sphere to the matter within the sphere, modified by a
    non-spherical factor.

SourceFiles
    nonSphericalHeatTransfer.C

\*---------------------------------------------------------------------------*/

#ifndef nonSphericalHeatTransfer_H
#define nonSphericalHeatTransfer_H

#include "sphericalHeatTransfer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace heatTransferModels
{

/*---------------------------------------------------------------------------*\
                    Class nonSphericalHeatTransfer Declaration
\*---------------------------------------------------------------------------*/

class nonSphericalHeatTransfer
:
    public sphericalHeatTransfer
{
    // Private Data

        //- Non-spherical factor
        const dimensionedScalar factor_;


public:

    //- Runtime type information
    TypeName("nonSpherical");


    // Constructors

        //- Construct from a dictionary and an interface
        nonSphericalHeatTransfer
        (
            const dictionary& dict,
            const phaseInterface& interface,
            const bool registerObject
        );


    //- Destructor
    virtual ~nonSphericalHeatTransfer();


    // Member Functions

        //- The heat transfer function K used in the enthalpy equation
        tmp<volScalarField> K(const scalar residualAlpha) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace heatTransferModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
