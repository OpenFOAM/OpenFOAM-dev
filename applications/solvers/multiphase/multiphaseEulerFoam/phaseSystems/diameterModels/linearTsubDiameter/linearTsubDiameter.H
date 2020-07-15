/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2020 OpenFOAM Foundation
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
    Foam::linearTsub

Description
    Vapour bubble diameter model for modelling of condensation of vapour
    bubbles. Calculates bubble diameter as a function of liquid phase
    subcooling.

    Reference:
    \verbatim
        Anglart, H., Nylund, O., Kurul, N., & Podowski, M. Z. (1997).
        CFD prediction of flow and phase distribution in fuel assemblies with
        spacers.
        Nuclear Engineering and Design, 177(1-3), 215-228.
    \endverbatim

SourceFiles
    linearTsub.C

\*---------------------------------------------------------------------------*/

#ifndef linearTsubDiameter_H
#define linearTsubDiameter_H

#include "sphericalDiameter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{

/*---------------------------------------------------------------------------*\
                           Class linearTsub Declaration
\*---------------------------------------------------------------------------*/

class linearTsub
:
    public spherical
{
    // Private Data

        //- Name of the liquid phase that is used to determine subcooling
        //  temperature
        word liquidPhaseName_;

        //- Reference diameter for low subcooling temperature
        dimensionedScalar d2_;

        //- Subcooling temperature where low subcooling diameter is reached
        dimensionedScalar Tsub2_;

        //- Reference diameter for high subcooling temperature
        dimensionedScalar d1_;

        //- Subcooling temperature where high subcooling diameter is reached
        dimensionedScalar Tsub1_;

        //- Diameter field
        volScalarField& d_;


protected:

    // Protected Member Functions

        //- Get the diameter field
        virtual tmp<volScalarField> calcD() const;


public:

    //- Runtime type information
    TypeName("linearTsub");


    // Constructors

        //- Construct from components
        linearTsub
        (
            const dictionary& diameterProperties,
            const phaseModel& phase
        );


    //- Destructor
    virtual ~linearTsub();


    // Member Functions

        //- Correct the diameter field
        virtual void correct();

        //- Read phaseProperties dictionary
        virtual bool read(const dictionary& phaseProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
