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
    Foam::kineticTheoryModels::viscosityModel

Description

SourceFiles
    viscosityModel.C

\*---------------------------------------------------------------------------*/

#ifndef kineticTheoryViscosityModel_H
#define kineticTheoryViscosityModel_H

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
                       Class viscosityModel Declaration
\*---------------------------------------------------------------------------*/

class viscosityModel
{
    // Private member functions

        //- Read the coefficients from coeffDict
        virtual bool readCoeffs(const dictionary& coeffDict)
        {
            return true;
        }


public:

    //- Runtime type information
    TypeName("viscosityModel");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        viscosityModel,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );


    // Constructors

        //- Construct from the coefficients dictionary
        viscosityModel(const dictionary& coeffDict);

        //- Disallow default bitwise copy construction
        viscosityModel(const viscosityModel&) = delete;


    // Selectors

        static autoPtr<viscosityModel> New
        (
            const dictionary& dict
        );


    //- Destructor
    virtual ~viscosityModel();


    // Member Functions

        //- Read the coefficients from the optional ".*Coeffs" sub-dictionary
        //  of the given dictionary
        bool read(const dictionary& dict);

        virtual tmp<volScalarField> nu
        (
            const volScalarField& alpha1,
            const volScalarField& Theta,
            const volScalarField& g0,
            const volScalarField& rho1,
            const volScalarField& da,
            const dimensionedScalar& e
        ) const = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const viscosityModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace kineticTheoryModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
