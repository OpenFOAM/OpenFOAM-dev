/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024-2025 OpenFOAM Foundation
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
    Foam::fv::homogeneousCondensation

Description
    Model for the homogeneous nucleation of liquid droplets out of a gaseous
    mixture

Usage
    Example usage:
    \verbatim
    homogeneousCondensation
    {
        type            homogeneousCondensation;
        libs            ("libmultiphaseEulerFvModels.so");

        // Phases between which the transfer occurs. The first phase is the
        // gas, and the second is the condensed liquid.
        phases          (gas water);

        // The specie that is condensing
        specie          H2O;

        // Linearise the latent heat contribution into the energy equation?
        energySemiImplicit no;

        // Saturation curve for the specie in the gaseous phase
        pSat            ArdenBuck;
    }
    \endverbatim

SourceFiles
    homogeneousCondensation.C

\*---------------------------------------------------------------------------*/

#ifndef homogeneousCondensation_H
#define homogeneousCondensation_H

#include "homogeneousNucleation.H"
#include "saturationPressureModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                   Class homogeneousCondensation Declaration
\*---------------------------------------------------------------------------*/

class homogeneousCondensation
:
    public homogeneousNucleation
{
private:

    // Private Data

        //- Saturation pressure model
        autoPtr<saturationPressureModel> saturationModel_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs(const dictionary& dict);

        //- Calculate and return the diameter of nucleates and the phase change
        //  rate per quantity of the solution phase
        virtual Pair<tmp<volScalarField::Internal>>
            dAndMDotByAlphaSolution() const;


public:

    //- Runtime type information
    TypeName("homogeneousCondensation");


    // Constructors

        //- Construct from explicit source name and mesh
        homogeneousCondensation
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Member Functions

        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
