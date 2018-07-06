/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::radiation::cloudScatter

Description
    Cloud radiation scatter model

SourceFiles
    scatterModel.C

\*---------------------------------------------------------------------------*/

#ifndef cloudScatter_H
#define cloudScatter_H

#include "scatterModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                         Class cloudScatter Declaration
\*---------------------------------------------------------------------------*/

class cloudScatter
:
    public scatterModel
{
    // Private data

        //- Absorption model dictionary
        dictionary coeffsDict_;

        //- Cloud names
        wordList cloudNames_;


public:

    //- Runtime type information
    TypeName("cloudScatter");


    // Constructors

        //- Construct from components
        cloudScatter
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~cloudScatter();


    // Member Operators

        // Access

            //- Return scatter coefficient
            tmp<volScalarField> sigmaEff() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
