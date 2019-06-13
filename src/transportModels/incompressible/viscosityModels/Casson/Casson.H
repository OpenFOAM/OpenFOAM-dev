/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2019 OpenFOAM Foundation
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
    Foam::viscosityModels::Casson

Description
    An incompressible Casson non-Newtonian viscosity model.

    References:
    \verbatim
        Casson, N. (1959).
        Rheology of disperse systems.
        In Proceedings of a Conference Organized by the
        British Society of Rheology.
        Pergamon Press, New York.

        Fournier, R. L. (2011).
        Basic transport phenomena in biomedical engineering.
        CRC Press.
    \endverbatim

    Example specification for blood:
    \verbatim
         transportModel  Casson;

         CassonCoeffs
         {
             m       3.934986e-6;
             tau0    2.9032e-6;
             nuMax   13.3333e-6;
             nuMin   3.9047e-6;
         }
    \endverbatim

SourceFiles
    Casson.C

\*---------------------------------------------------------------------------*/

#ifndef Casson_H
#define Casson_H

#include "viscosityModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{


/*---------------------------------------------------------------------------*\
                           Class Casson Declaration
\*---------------------------------------------------------------------------*/

class Casson
:
    public viscosityModel
{
    // Private Data

        dictionary CassonCoeffs_;

        dimensionedScalar m_;
        dimensionedScalar tau0_;
        dimensionedScalar nuMin_;
        dimensionedScalar nuMax_;

        volScalarField nu_;


    // Private Member Functions

        //- Calculate and return the laminar viscosity
        tmp<volScalarField> calcNu() const;


public:

    //- Runtime type information
    TypeName("Casson");


    // Constructors

        //- Construct from components
        Casson
        (
            const word& name,
            const dictionary& viscosityProperties,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~Casson()
    {}


    // Member Functions

        //- Return the laminar viscosity
        virtual tmp<volScalarField> nu() const
        {
            return nu_;
        }

        //- Return the laminar viscosity for patch
        virtual tmp<scalarField> nu(const label patchi) const
        {
            return nu_.boundaryField()[patchi];
        }

        //- Correct the laminar viscosity
        virtual void correct()
        {
            nu_ = calcNu();
        }

        //- Read transportProperties dictionary
        virtual bool read(const dictionary& viscosityProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace viscosityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
