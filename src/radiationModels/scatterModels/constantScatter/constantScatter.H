/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::radiationModels::scatterModels::constant

Description
    Constant radiation scatter coefficient

SourceFiles
    scatterModel.C

\*---------------------------------------------------------------------------*/

#ifndef constantScatter_H
#define constantScatter_H

#include "scatterModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{
namespace scatterModels
{

/*---------------------------------------------------------------------------*\
                          Class constant Declaration
\*---------------------------------------------------------------------------*/

class constant
:
    public scatterModel
{
    // Private Data

        //- Coefficients dictionary
        dictionary coeffsDict_;

        //- Scattering coefficient / [1/m]
        dimensionedScalar sigma_;

        //- Linear-anisotropic phase function coefficient / []
        //  -1 < C < 1
        //  - = backward scattering
        //  0 = isotropic scattering (reasonable default value)
        //  + = forward scattering
        dimensionedScalar C_;


public:

    //- Runtime type information
    TypeName("constant");


    // Constructors

        //- Construct from components
        constant(const dictionary& dict, const fvMesh& mesh);


    //- Destructor
    virtual ~constant();


    // Member Functions

        //- Return scatter coefficient
        tmp<volScalarField> sigmaEff() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace scatterModels
} // End namespace radiationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
