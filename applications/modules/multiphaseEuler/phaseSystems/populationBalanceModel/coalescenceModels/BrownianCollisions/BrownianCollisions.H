/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2023 OpenFOAM Foundation
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
    Foam::diameterModels::coalescenceModels::BrownianCollisions

Description
    Model describing coagulation due to Brownian motion. Utilises collisional
    diameters and the Cunningham slip correction. The slip correction
    coefficient is implemented in the following form:

    \f[
        C_{c_i} = 1 + \lambda [A_1 + A_2 \exp(-A_3 d_i/\lambda)]/d_i\,.
    \f]

    The coefficients default to the values proposed by Davis (1945). The mean
    free path is computed by

    \f[
        \lambda = \frac{kT}{\sqrt{2} \pi p \sigma^{2}}\,.
    \f]

    \vartable
        A_1       | Coefficient [-]
        A_2       | Coefficient [-]
        A_3       | Coefficient [-]
        \sigma    | Lennard-Jones parameter [m]
    \endvartable

    Reference:
    \verbatim
        Davies, C. N. (1945).
        Definitive equations for the fluid resistance of spheres.
        Proceedings of the Physical Society, 57(4), 259.
    \endverbatim

Usage
    \table
        Property    | Description                | Required    | Default value
        A1          | Coefficient A1             | no          | 2.514
        A2          | Coefficient A2             | no          | 0.8
        A3          | Coefficient A2             | no          | 0.55
        sigma       | Lennard-Jones parameter    | yes         | none
    \endtable

SourceFiles
    BrownianCollisions.C

\*---------------------------------------------------------------------------*/

#ifndef BrownianCollisions_H
#define BrownianCollisions_H

#include "coalescenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace coalescenceModels
{

/*---------------------------------------------------------------------------*\
                     Class BrownianCollisions Declaration
\*---------------------------------------------------------------------------*/

class BrownianCollisions
:
    public coalescenceModel
{
    // Private Data

        //- Cunningham slip correction coefficient A1
        const dimensionedScalar A1_;

        //- Cunningham slip correction coefficient A2
        const dimensionedScalar A2_;

        //- Cunningham slip correction coefficient A3
        const dimensionedScalar A3_;

        //- Lennard-Jones sigma parameter
        const dimensionedScalar sigma_;

        //- Mean free path
        volScalarField lambda_;


public:

    //- Runtime type information
    TypeName("BrownianCollisions");

    // Constructor

        BrownianCollisions
        (
            const populationBalanceModel& popBal,
            const dictionary& dict
        );


    //- Destructor
    virtual ~BrownianCollisions()
    {}


    // Member Functions

        //- Precompute diameter independent expressions
        virtual void precompute();

        //- Add to coalescenceRate
        virtual void addToCoalescenceRate
        (
            volScalarField& coalescenceRate,
            const label i,
            const label j
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace coalescenceModels
} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
