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
    Foam::rhoThermoCombustion

Description
    Density-based thermo model wrapper for combustion models

SourceFiles
    rhoThermoCombustion.C

\*---------------------------------------------------------------------------*/

#ifndef rhoThermoCombustion_H
#define rhoThermoCombustion_H

#include "autoPtr.H"
#include "rhoCombustionModel.H"
#include "rhoChemistryModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                    class rhoThermoCombustion Declaration
\*---------------------------------------------------------------------------*/

class rhoThermoCombustion
:
    public rhoCombustionModel
{
    // Private Member Functions

        //- Construct as copy (not implemented)
        rhoThermoCombustion(const rhoThermoCombustion&);

        //- Disallow default bitwise assignment
        void operator=(const rhoThermoCombustion&);


protected:

    // Protected data

        //- Pointer to thermo model
        autoPtr<rhoReactionThermo> thermoPtr_;


public:

    // Constructors

        //- Construct from components and thermo
        rhoThermoCombustion(const word& modelType, const fvMesh& mesh);


    //- Destructor
    virtual ~rhoThermoCombustion();


    // Member Functions

        //- Return access to the thermo package
        virtual rhoReactionThermo& thermo();

        //- Return const access to the thermo package
        virtual const rhoReactionThermo& thermo() const;

        //- Return const access to the density field
        virtual tmp<volScalarField> rho() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
