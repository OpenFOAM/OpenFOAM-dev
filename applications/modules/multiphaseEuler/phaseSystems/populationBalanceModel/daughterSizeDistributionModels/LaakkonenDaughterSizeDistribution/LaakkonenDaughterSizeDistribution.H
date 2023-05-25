/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2023 OpenFOAM Foundation
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
    Foam::diameterModels::daughterSizeDistributionModels::
    LaakkonenDaughterSizeDistribution

Description
    Daughter size distribution model of Laakkonen et al. (2007). Note that the
    diameters in the original expression were substituted by bubble volumes
    giving

    \f[
        (1 + C_4)(2 + C_4)(3 + C_4)(4 + C_4)
        \times \left(\frac{1}{3}\right) \left(\frac{1}{v_j}\right)
        \left(\frac{v_i}{v_j}\right)^{2} \left(1 - \frac{v_i}{v_j}\right)^{C_4}
    \f]

    where

    \vartable
        v_i         |  Volume of daughter bubble i [m3]
        v_j         |  Volume of mother bubble j [m3]
    \endvartable

    The total number of daughter bubbles generated depends on C4 and evaluates
    to

    \f[
        \frac{4}{3} + \frac{C_4}{3}
    \f]

    References:
    \verbatim
        Laakkonen, M., Moilanen, P., Alopaeus, V., & Aittamaa, J. (2007).
        Modelling local bubble size distributions in agitated vessels.
        Chemical Engineering Science, 62, 721–740.
    \endverbatim

Usage
    \table
        Property     | Description             | Required    | Default value
        C4           | Coefficient C4          | no          | 18.25
    \endtable

SourceFiles
    LaakkonenBinaryBreakup.C

\*---------------------------------------------------------------------------*/

#ifndef LaakkonenDaughterSizeDistribution_H
#define LaakkonenDaughterSizeDistribution_H

#include "daughterSizeDistributionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace daughterSizeDistributionModels
{

/*---------------------------------------------------------------------------*\
              Class LaakkonenDaughterSizeDistribution Declaration
\*---------------------------------------------------------------------------*/

class LaakkonenDaughterSizeDistribution
:
    public daughterSizeDistributionModel
{
    // Private Data

        //- Optional coefficient C4, defaults to 18.25
        dimensionedScalar C4_;


public:

    //- Runtime type information
    TypeName("Laakkonen");


    // Constructor

        LaakkonenDaughterSizeDistribution
        (
            const breakupModel& breakup,
            const dictionary& dict
        );


    //- Destructor
    virtual ~LaakkonenDaughterSizeDistribution();


    // Member Functions

        //- Calculate antiderivative
        dimensionedScalar antiderivative
        (
            const dimensionedScalar& xk,
            const dimensionedScalar& v,
            const dimensionedScalar& bndr,
            const dimensionedScalar range
        ) const;

        //- Calculate contribution to sizeGroup i due to breakup in sizeGroup k
        virtual dimensionedScalar calcNik(const label i, const label k) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace daughterSizeDistributionModels
} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
