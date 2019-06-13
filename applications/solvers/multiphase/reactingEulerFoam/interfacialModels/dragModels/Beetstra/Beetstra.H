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
    Foam::dragModels::Beetstra

Description
    Drag model of Beetstra et al. for monodisperse gas-particle flows obtained
    with direct numerical simulations with the Lattice-Boltzmann method and
    accounting for the effect of particle ensembles.

    Reference:
    \verbatim
        Beetstra, R., van der Hoef, M. A., & Kuipers, J. a. M. (2007).
        Drag force of intermediate Reynolds number flow past mono- and
        bidisperse arrays of spheres.
        AIChE Journal, 53(2), 489–501.
    \endverbatim

SourceFiles
    Beetstra.C

\*---------------------------------------------------------------------------*/

#ifndef Beetstra_H
#define Beetstra_H

#include "dragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace dragModels
{

/*---------------------------------------------------------------------------*\
                           Class Beetstra Declaration
\*---------------------------------------------------------------------------*/

class Beetstra
:
    public dragModel
{
    // Private Data

        //- Residual Reynolds Number
        const dimensionedScalar residualRe_;


public:

    //- Runtime type information
    TypeName("Beetstra");


    // Constructors

        //- Construct from a dictionary and a phase pair
        Beetstra
        (
            const dictionary& dict,
            const phasePair& pair,
            const bool registerObject
        );


    //- Destructor
    virtual ~Beetstra();


    // Member Functions

        //- Drag coefficient
        virtual tmp<volScalarField> CdRe() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace dragModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
