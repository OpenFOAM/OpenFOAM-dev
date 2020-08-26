/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2020 OpenFOAM Foundation
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
    Foam::dragModels::aerosolDrag

Description
    Stokes drag with Cunningham slip correction. The slip correction coefficient
    is implemented in the following form:

    \f[
        C_c = 1 + \lambda [A_1 + A_2 \exp(-A_3 d/\lambda)]/d\,.
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
    aerosolDrag.C

\*---------------------------------------------------------------------------*/

#ifndef aerosolDrag_H
#define aerosolDrag_H

#include "dragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{

/*---------------------------------------------------------------------------*\
                         Class aerosolDrag Declaration
\*---------------------------------------------------------------------------*/

class aerosolDrag
:
    public dragModel
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


public:

    //- Runtime type information
    TypeName("aerosolDrag");


    // Constructors

        //- Construct from a dictionary and a phase pair
        aerosolDrag
        (
            const dictionary& dict,
            const phasePair& pair,
            const bool registerObject
        );


    //- Destructor
    virtual ~aerosolDrag();


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
