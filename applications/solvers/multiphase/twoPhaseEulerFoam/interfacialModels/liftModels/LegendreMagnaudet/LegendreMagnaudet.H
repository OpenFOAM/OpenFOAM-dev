/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2019 OpenFOAM Foundation
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
    Foam::liftModels::LegendreMagnaudet

Description
    Lift model of Legendre and Magnaudet.

    References:
    \verbatim
        "Implementation and Comparison of Correlations for interfacial Forces
        in a Gas-Liquid System within an Euler-Euler Framework"
        Otromke, M.
        PhD Thesis
        April 2013
    \endverbatim

    \verbatim
        "The lift force on a spherical bubble in a viscous linear shear flow"
        Legendre, D., Magnaudet, J.,
        Journal of Fluid Mechanics
        Volume 368, August 1998, pp. 81-126
    \endverbatim

SourceFiles
    LegendreMagnaudet.C

\*---------------------------------------------------------------------------*/

#ifndef LegendreMagnaudet_H
#define LegendreMagnaudet_H

#include "liftModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace liftModels
{

/*---------------------------------------------------------------------------*\
                           Class LegendreMagnaudet Declaration
\*---------------------------------------------------------------------------*/

class LegendreMagnaudet
:
    public liftModel
{
    // Private Data

        //- Residual Reynold's number
        const dimensionedScalar residualRe_;


public:

    //- Runtime type information
    TypeName("LegendreMagnaudet");


    // Constructors

        //- Construct from a dictionary and a phase pair
        LegendreMagnaudet
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~LegendreMagnaudet();


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
