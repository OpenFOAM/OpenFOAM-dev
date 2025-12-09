/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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
    Foam::saturationModels::function1Pressure

Description
    Saturation vapour pressure in terms of temperature in Kelvins. The
    saturation pressure in Pa is specified as a Foam::Function1 type, to enable
    use of, e.g. constant, polynomial, table values.

    Examples:

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
    function1Pressure.C

\*---------------------------------------------------------------------------*/

#ifndef function1Pressure_H
#define function1Pressure_H

#include "saturationPressureModel.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace saturationModels
{

/*---------------------------------------------------------------------------*\
                          Class function1Pressure Declaration
\*---------------------------------------------------------------------------*/

class function1Pressure
:
    public saturationPressureModel
{
    // Private Data

        //- Saturation pressure as a function of temperature
        autoPtr<Function1<scalar>> function_;


public:

    //- Runtime type information
    TypeName("function1");


    // Constructors

        //- Construct from a dictionary and an interface
        function1Pressure
        (
            const dictionary& dict
        );


    //- Destructor
    virtual ~function1Pressure();


    // Member Functions

        //- Saturation pressure for scalarField
        virtual tmp<scalarField> pSat(const scalarField& T) const;

        //- Saturation pressure derivative w.r.t. temperature for scalarField
        virtual tmp<scalarField> pSatPrime(const scalarField& T) const;

        //- Natural log of the saturation pressure for scalarField
        virtual tmp<scalarField> lnPSat(const scalarField& T) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace saturationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
