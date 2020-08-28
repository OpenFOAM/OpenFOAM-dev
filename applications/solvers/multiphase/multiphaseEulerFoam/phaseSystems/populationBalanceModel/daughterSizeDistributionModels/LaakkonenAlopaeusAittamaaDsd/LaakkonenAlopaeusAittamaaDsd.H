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
    Foam::diameterModels::daughterSizeDistributionModels::
    LaakkonenAlopaeusAittamaaDsd

Description
    Daughter size distribution model used by Laakkonen et al. (2006). Note that
    the diameters in the original expression was substituted by bubble volumes
    giving

    \f[
        \left(9 + \frac{33}{2}C_4 + 9C_4^2 + \frac{3}{2}C_4^3\right)
        \times \left(\frac{1}{3}\right) \left(\frac{1}{v_j}\right)
        \left(\frac{v_i}{v_j}\right)^{2} \left(1 - \frac{v_i}{v_j}\right)^{C_4}
    \f]

    where

    \vartable
        v_i         |  Volume of daughter bubble i [m]
        v_j         |  Volume of mother bubble j [m]
    \endvartable

    References:
    \verbatim
        Laakkonen, M., Alopaeus, V., & Aittamaa, J. (2006).
        Validation of bubble breakage, coalescence and mass transfer models for
        gas-liquid dispersion in agitated vessel.
        Chemical engineering science, 61(1), 218-228.
        Eq. 3, p. 220.
    \endverbatim

Usage
    \table
        Property     | Description             | Required    | Default value
        C4           | Coefficient C4          | no          | 4.3
    \endtable

SourceFiles
    LaakkonenAlopaeusAittamaaBinaryBreakup.C

\*---------------------------------------------------------------------------*/

#ifndef LaakkonenAlopaeusAittamaaDsd_H
#define LaakkonenAlopaeusAittamaaDsd_H

#include "daughterSizeDistributionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace daughterSizeDistributionModels
{

/*---------------------------------------------------------------------------*\
                Class LaakkonenAlopaeusAittamaaDsd Declaration
\*---------------------------------------------------------------------------*/

class LaakkonenAlopaeusAittamaaDsd
:
    public daughterSizeDistributionModel
{
    // Private Data

        //- Optional coefficient C4, defaults to 4.3
        dimensionedScalar C4_;


public:

    //- Runtime type information
    TypeName("LaakkonenAlopaeusAittamaa");


    // Constructor

        LaakkonenAlopaeusAittamaaDsd
        (
            const breakupModel& breakup,
            const dictionary& dict
        );


    //- Destructor
    virtual ~LaakkonenAlopaeusAittamaaDsd();


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
