/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::radiationModels::absorptionEmissionModels::constant

Description
    Constant radiation absorption and emission coefficients for continuous
    phase

SourceFiles
    constantAbsorptionEmission.C

\*---------------------------------------------------------------------------*/

#ifndef constantAbsorptionEmission_H
#define constantAbsorptionEmission_H

#include "absorptionEmissionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{
namespace absorptionEmissionModels
{

/*---------------------------------------------------------------------------*\
                          Class constant Declaration
\*---------------------------------------------------------------------------*/

class constant
:
    public absorptionEmissionModel
{
    // Private Data

        //- Absorption model dictionary
        dictionary coeffsDict_;

        //- Absorption coefficient / [1/m]
        dimensionedScalar a_;

        //- Emission coefficient / [1/m]
        dimensionedScalar e_;

        //- Emission contribution / [kg/m/s^3]
        dimensionedScalar E_;


public:

    //- Runtime type information
    TypeName("constant");


    // Constructors

        //- Construct from components
        constant(const dictionary& dict, const fvMesh& mesh);


    //- Destructor
    virtual ~constant();


    // Member Functions

        // Access

            // Absorption coefficient

                //- Absorption coefficient for continuous phase
                tmp<volScalarField> aCont(const label bandI = 0) const;


            // Emission coefficient

                //- Emission coefficient for continuous phase
                tmp<volScalarField> eCont(const label bandI = 0) const;


            // Emission contribution

                //- Emission contribution for continuous phase
                tmp<volScalarField> ECont(const label bandI = 0) const;


    // Member Functions

        inline bool isGrey() const
        {
            return true;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace absorptionEmissionModels
} // End namespace radiationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
