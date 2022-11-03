/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2022 OpenFOAM Foundation
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
    Foam::virtualMassModels::Lamb

Description
    Virtual mass model of Lamb.

    Reference:
    \verbatim
        Lamb, H. (1993).
        Hydrodynamics.
        Cambridge university press.
    \endverbatim

SourceFiles
    Lamb.C

\*---------------------------------------------------------------------------*/

#ifndef Lamb_H
#define Lamb_H

#include "dispersedVirtualMassModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class aspectRatioModel;

namespace virtualMassModels
{

/*---------------------------------------------------------------------------*\
                            Class Lamb Declaration
\*---------------------------------------------------------------------------*/

class Lamb
:
    public dispersedVirtualMassModel
{
private:

    //- Private Data

        //- The aspect ratio model
        autoPtr<aspectRatioModel> aspectRatio_;


public:

    //- Runtime type information
    TypeName("Lamb");


    // Constructors

        //- Construct from a dictionary and an interface
        Lamb
        (
            const dictionary& dict,
            const phaseInterface& interface,
            const bool registerObject
        );


    //- Destructor
    virtual ~Lamb();


    // Member Functions

        //- Virtual mass coefficient
        virtual tmp<volScalarField> Cvm() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace virtualMassModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
