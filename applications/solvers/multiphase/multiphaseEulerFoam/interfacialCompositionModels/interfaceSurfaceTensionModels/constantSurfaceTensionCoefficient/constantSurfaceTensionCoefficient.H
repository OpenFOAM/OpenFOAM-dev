/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2022 OpenFOAM Foundation
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
    Foam::interfaceSurfaceTensionModels::constantSurfaceTensionCoefficient

Description
    Constant value surface tension model.

SourceFiles
    constantSurfaceTensionCoefficient.C

\*---------------------------------------------------------------------------*/

#ifndef constantSurfaceTensionCoefficient_H
#define constantSurfaceTensionCoefficient_H

#include "interfaceSurfaceTensionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace interfaceSurfaceTensionModels
{

/*---------------------------------------------------------------------------*\
              Class constantSurfaceTensionCoefficient Declaration
\*---------------------------------------------------------------------------*/

class constantSurfaceTensionCoefficient
:
    public interfaceSurfaceTensionModel
{
    // Private Data

        //- Constant surface tension value
        const dimensionedScalar sigma_;


public:

    //- Runtime type information
    TypeName("constant");


    // Constructors

        //- Construct from a dictionary and an interface
        constantSurfaceTensionCoefficient
        (
            const dictionary& dict,
            const phaseInterface& interface
        );


    //- Destructor
    virtual ~constantSurfaceTensionCoefficient();


    // Member Functions

        //- Surface tension
        virtual tmp<volScalarField> sigma() const;

        //- Surface tension for a patch
        virtual tmp<scalarField> sigma(const label patchi) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace interfaceSurfaceTensionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
