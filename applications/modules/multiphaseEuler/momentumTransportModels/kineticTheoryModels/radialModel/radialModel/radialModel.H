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
    Foam::kineticTheoryModels::radialModel

SourceFiles
    radialModel.C

\*---------------------------------------------------------------------------*/

#ifndef radialModel_H
#define radialModel_H

#include "dictionary.H"
#include "volFields.H"
#include "dimensionedTypes.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{

/*---------------------------------------------------------------------------*\
                         Class radialModel Declaration
\*---------------------------------------------------------------------------*/

class radialModel
{
    // Private member functions

        //- Read the coefficients from coeffDict
        virtual bool readCoeffs(const dictionary& coeffDict)
        {
            return true;
        }


public:

    //- Runtime type information
    TypeName("radialModel");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        radialModel,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );


    // Constructors

        //- Construct from the coefficients dictionary
        radialModel(const dictionary& coeffDict);

        //- Disallow default bitwise copy construction
        radialModel(const radialModel&) = delete;


    // Selectors

        static autoPtr<radialModel> New
        (
            const dictionary& dict
        );


    //- Destructor
    virtual ~radialModel();


    // Member Functions

        //- Read the coefficients from the optional ".*Coeffs" sub-dictionary
        //  of the given dictionary
        bool read(const dictionary& dict);

        //- Radial distribution function
        virtual tmp<volScalarField> g0
        (
            const volScalarField& alpha,
            const dimensionedScalar& alphaMinFriction,
            const dimensionedScalar& alphaMax
        ) const = 0;

        //- Derivative of the radial distribution function
        virtual tmp<volScalarField> g0prime
        (
            const volScalarField& alpha,
            const dimensionedScalar& alphaMinFriction,
            const dimensionedScalar& alphaMax
        ) const = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const radialModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace kineticTheoryModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
