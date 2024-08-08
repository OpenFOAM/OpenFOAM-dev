/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
    Foam::radiationModels::absorptionEmissionModels::binary

Description
    Radiation coefficient based on two absorption models

SourceFiles
    binary.C

\*---------------------------------------------------------------------------*/

#ifndef binary_H
#define binary_H

#include "absorptionEmissionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{
namespace absorptionEmissionModels
{

/*---------------------------------------------------------------------------*\
                           Class binary Declaration
\*---------------------------------------------------------------------------*/

class binary
:
    public absorptionEmissionModel
{
    // Private Data

        //- First absorption model
        autoPtr<absorptionEmissionModel> model1_;

        //- Second absorption model
        autoPtr<absorptionEmissionModel> model2_;


public:

    //- Runtime type information
    TypeName("binary");


    // Constructors

        //- Construct from components
        binary(const dictionary& dict, const fvMesh& mesh);


    //- Destructor
    virtual ~binary();


    // Member Functions

        // Access

            // Absorption coefficient

                //- Absorption coefficient for continuous phase
                virtual tmp<volScalarField> aCont(const label bandI = 0) const;

                //- Absorption coefficient for dispersed phase
                virtual tmp<volScalarField> aDisp(const label bandI = 0) const;


            // Emission coefficient

                //- Emission coefficient for continuous phase
                virtual tmp<volScalarField> eCont(const label bandI = 0) const;

                //- Emission coefficient for dispersed phase
                virtual tmp<volScalarField> eDisp(const label bandI = 0) const;


            // Emission contribution

                //- Emission contribution for continuous phase
                virtual tmp<volScalarField> ECont(const label bandI = 0) const;

                //- Emission contribution for continuous phase
                virtual tmp<volScalarField> EDisp(const label bandI = 0) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace absorptionEmissionModels
} // End namespace radiationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
