/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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
    Foam::psiThermoCombustion

Description
    Compressibility-based thermo model wrapper for combustion models

SourceFiles
    psiThermoCombustion.C

\*---------------------------------------------------------------------------*/

#ifndef psiThermoCombustion_H
#define psiThermoCombustion_H

#include "autoPtr.H"
#include "psiCombustionModel.H"
#include "psiReactionThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                    class psiThermoCombustion Declaration
\*---------------------------------------------------------------------------*/

class psiThermoCombustion
:
    public psiCombustionModel
{
    // Private Member Functions

        //- Construct as copy (not implemented)
        psiThermoCombustion(const psiThermoCombustion&);

        //- Disallow default bitwise assignment
        void operator=(const psiThermoCombustion&);


protected:

    // Protected data

        //- Pointer to chemistry model
        autoPtr<psiReactionThermo> thermoPtr_;


public:

    // Constructors

        //- Construct from components and thermo
        psiThermoCombustion(const word& modelType, const fvMesh& mesh);


    //- Destructor
    virtual ~psiThermoCombustion();


    // Member Functions

        //- Return access to the thermo package
        virtual psiReactionThermo& thermo();

        //- Return const access to the thermo package
        virtual const psiReactionThermo& thermo() const;

        //- Return const access to the density field
        virtual tmp<volScalarField> rho() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
