/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2025 OpenFOAM Foundation
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
    Foam::saturationModels::function1Temperature

Description
    Saturation vapour temperature in terms of the vapour pressure (in Pa). The
    saturation temperature in Kelvins is specified as a Foam::Function1 type,
    to enable use of, e.g. constant, polynomial, table values.

    Examples:

    \verbatim
        type function1;
        function polynomial
        (
            308.0422
            0.0015096
            -1.61589e-8
            1.114106e-13
            -4.52216e-19
            1.05192e-24
            -1.2953e-30
            6.5365e-37
        );
    \endverbatim

    \verbatim
        type function1;
        function
        {
            type                table;
            file                "filename.csv";
            format              csv;
            nHeaderLine         1;
            columns             (0 1);
            separator           ",";
            mergeSeparators     no;
            outOfBounds         clamp;
            interpolationScheme linear;
        }
    \endverbatim

SourceFiles
    function1Temperature.C

\*---------------------------------------------------------------------------*/

#ifndef function1Temperature_H
#define function1Temperature_H

#include "saturationTemperatureModel.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace saturationModels
{

/*---------------------------------------------------------------------------*\
                          Class function1Temperature Declaration
\*---------------------------------------------------------------------------*/

class function1Temperature
:
    public saturationTemperatureModel
{
    // Private Data

        //- Saturation temperature as a function of pressure
        autoPtr<Function1<scalar>> function_;


public:

    //- Runtime type information
    TypeName("function1");


    // Constructors

        //- Construct from a dictionary and an interface
        function1Temperature
        (
            const dictionary& dict
        );


    //- Destructor
    virtual ~function1Temperature();


    // Member Functions

        //- Saturation temperature for scalarField
        virtual tmp<scalarField> Tsat(const scalarField& p) const;

        //- Saturation temperature derivative w.r.t. pressure for scalarField
        virtual tmp<scalarField> TsatPrime(const scalarField& p) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace saturationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
