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
    Foam::laminarModels::generalisedNewtonianViscosityModels::Casson

Description
    Casson generalised Newtonian viscosity model

    References:
    \verbatim
        Casson, N. (1959).
        Rheology of disperse systems.
        In Proceedings of a Conference Organised by the
        British Society of Rheology.
        Pergamon Press, New York.

        Fournier, R. L. (2011).
        Basic transport phenomena in biomedical engineering.
        CRC Press.
    \endverbatim

    Example specification for blood:
    \verbatim
         viscosityModel Casson;

         m       3.934986e-6;
         tau0    2.9032e-6;
         nuMax   13.3333e-6;
         nuMin   3.9047e-6;
    \endverbatim

SourceFiles
    Casson.C

\*---------------------------------------------------------------------------*/

#ifndef Casson_H
#define Casson_H

#include "strainRateViscosityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarModels
{
namespace generalisedNewtonianViscosityModels
{

/*---------------------------------------------------------------------------*\
                           Class Casson Declaration
\*---------------------------------------------------------------------------*/

class Casson
:
    public strainRateViscosityModel
{
    // Private Data

        dimensionedScalar m_;
        dimensionedScalar tau0_;
        dimensionedScalar nuMin_;
        dimensionedScalar nuMax_;


public:

    //- Runtime type information
    TypeName("Casson");


    // Constructors

        //- Construct from components
        Casson
        (
            const dictionary& viscosityProperties,
            const Foam::viscosity& viscosity,
            const volVectorField& U
        );


    //- Destructor
    virtual ~Casson()
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
