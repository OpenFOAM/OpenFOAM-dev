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
    Foam::diameterModels::fixedInterfacialArea

Description
    fixedInterfacialArea dispersed-phase diameter model.
    The interfacial are is set by providing phase surface area divided by phase
    volume, AvbyAlpha, etiher as a constant value or as a field.

SourceFiles
    fixedInterfacialAreaDiameter.C

\*---------------------------------------------------------------------------*/

#ifndef fixedInterfacialAreaDiameter_H
#define fixedInterfacialAreaDiameter_H

#include "diameterModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{

/*---------------------------------------------------------------------------*\
                    Class fixedInterfacialArea Declaration
\*---------------------------------------------------------------------------*/

class fixedInterfacialArea
:
    public diameterModel
{
    // Private Data

        //- Uniform fixed area by volume of the stationary phase for calculation
        //  of the interfacial area
        dimensionedScalar AvbyAlpha_;

        //- Optional area by volume field of the stationary phase for
        //  calculation of interfacial area
        autoPtr<volScalarField> AvbyAlphaFieldPtr_;


public:

    //- Runtime type information
    TypeName("fixedInterfacialArea");


    // Constructors

        //- Construct from dictionary and phase
        fixedInterfacialArea
        (
            const dictionary& diameterProperties,
            const phaseModel& phase
        );


    //- Destructor
    virtual ~fixedInterfacialArea();


    // Member Functions

        //- Get the diameter field
        virtual tmp<volScalarField> d() const;

        //- Get the surface area per unit volume field
        virtual tmp<volScalarField> Av() const;

        //- Read diameterProperties dictionary
        virtual bool read(const dictionary& diameterProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
