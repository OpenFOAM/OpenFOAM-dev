/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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
    Foam::liftModels::SaffmanMei

Description
    Lift model of Saffman (1965) as extended by Mei (1992). Applicable to
    spherical particles.

    References:
    \verbatim
        Saffman, P. G. T. (1965).
        The lift on a small sphere in a slow shear flow.
        Journal of fluid mechanics, 22(2), 385-400.

        Mei, R. (1992).
        An approximate expression for the shear lift force on a spherical
        particle at finite Reynolds number.
        International Journal of Multiphase Flow, 18(1), 145-147.
    \endverbatim

SourceFiles
    SaffmanMei.C

\*---------------------------------------------------------------------------*/

#ifndef SaffmanMei_H
#define SaffmanMei_H

#include "dispersedLiftModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace liftModels
{

/*---------------------------------------------------------------------------*\
                         Class SaffmanMei Declaration
\*---------------------------------------------------------------------------*/

class SaffmanMei
:
    public dispersedLiftModel
{
    // Private Data

        //- Residual Reynold's number
        const dimensionedScalar residualRe_;


public:

    //- Runtime type information
    TypeName("SaffmanMei");


    // Constructors

        //- Construct from a dictionary and an interface
        SaffmanMei
        (
            const dictionary& dict,
            const phaseInterface& interface
        );


    //- Destructor
    virtual ~SaffmanMei();


    // Member Functions

        //- Lift coefficient
        virtual tmp<volScalarField> Cl() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace liftModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
