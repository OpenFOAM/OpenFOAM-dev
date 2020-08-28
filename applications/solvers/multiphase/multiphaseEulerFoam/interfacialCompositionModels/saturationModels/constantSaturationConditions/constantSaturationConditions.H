/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2020 OpenFOAM Foundation
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
    Foam::saturationModels::constantSaturationConditions

Description
    Constant saturation pressure and temperature.

SourceFiles
    constantSaturationConditions.C

\*---------------------------------------------------------------------------*/

#ifndef constantSaturationConditions_H
#define constantSaturationConditions_H

#include "saturationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace saturationModels
{

/*---------------------------------------------------------------------------*\
                Class constantSaturationConditions Declaration
\*---------------------------------------------------------------------------*/

class constantSaturationConditions
:
    public saturationModel
{
protected:

    // Private Data

        //- Constant saturation pressure
        dimensionedScalar pSat_;

        //- Constant saturation temperature
        dimensionedScalar Tsat_;


public:

    //- Runtime type information
    TypeName("constant");

    // Constructors

        //- Construct from a dictionary
        constantSaturationConditions
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~constantSaturationConditions();


    // Member Functions

        //- Saturation pressure
        virtual tmp<volScalarField> pSat(const volScalarField& T) const;

        //- Saturation pressure derivative w.r.t. temperature
        virtual tmp<volScalarField> pSatPrime(const volScalarField& T) const;

        //- Natural log of the saturation pressure
        virtual tmp<volScalarField> lnPSat(const volScalarField& T) const;

        //- Saturation temperature
        virtual tmp<volScalarField> Tsat(const volScalarField& p) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace saturationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
