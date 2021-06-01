/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2021 OpenFOAM Foundation
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
    Foam::laminarModels::generalisedNewtonianViscosityModels::strainRateFunction

Description
    Run-time selected strain-rate function generalised Newtonian viscosity model

    Example linear function of strain-rate:
    \verbatim
        generalisedNewtonianModel strainRateFunction;

        function polynomial ((0 0.1) (1 1.3));
    \endverbatim

See also
    Foam::generalisedNewtonianViscosityModel
    Foam::Function1

SourceFiles
    strainRateFunction.C

\*---------------------------------------------------------------------------*/

#ifndef strainRateFunction_H
#define strainRateFunction_H

#include "generalisedNewtonianViscosityModel.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarModels
{
namespace generalisedNewtonianViscosityModels
{

/*---------------------------------------------------------------------------*\
                     Class strainRateFunction Declaration
\*---------------------------------------------------------------------------*/

class strainRateFunction
:
    public generalisedNewtonianViscosityModel
{
    // Private Data

        //- Strain-rate function
        autoPtr<Function1<scalar>> strainRateFunction_;


public:

    //- Runtime type information
    TypeName("strainRateFunction");


    // Constructors

        //- Construct from components
        strainRateFunction(const dictionary& viscosityProperties);


    //- Destructor
    virtual ~strainRateFunction()
    {}


    // Member Functions

        //- Read transportProperties dictionary
        virtual bool read(const dictionary& viscosityProperties);

        //- Return the laminar viscosity
        virtual tmp<volScalarField> nu
        (
            const volScalarField& nu0,
            const volScalarField& strainRate
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace generalisedNewtonianViscosityModels
} // End namespace laminarModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
