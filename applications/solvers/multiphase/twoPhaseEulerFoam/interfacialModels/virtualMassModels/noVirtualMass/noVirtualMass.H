/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 OpenFOAM Foundation
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
    Foam::virtualMassModels::noVirtualMass

Description

SourceFiles
    noVirtualMass.C

\*---------------------------------------------------------------------------*/

#ifndef noVirtualMass_H
#define noVirtualMass_H

#include "virtualMassModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace virtualMassModels
{

/*---------------------------------------------------------------------------*\
                         Class noVirtualMass Declaration
\*---------------------------------------------------------------------------*/

class noVirtualMass
:
    public virtualMassModel
{
public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from a dictionary and a phase pair
        noVirtualMass
        (
            const dictionary& dict,
            const phasePair& pair,
            const bool registerObject
        );


    //- Destructor
    virtual ~noVirtualMass();


    // Member Functions

        //- Virtual mass coefficient
        virtual tmp<volScalarField> Cvm() const;

        //- The virtual mass function K used in the momentum equation
        virtual tmp<volScalarField> K() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace virtualMassModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
