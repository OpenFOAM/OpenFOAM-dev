/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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
    Foam::diameterModels::coalescenceModels::CoulaloglouTavlaridesCoalescence

Description
    Model of Coulaloglou and Tavlarides (1977). The coalescence rate is
    calculated by

        \f[
            C_1 (v_i^{2/3} + v_j^{2/3}) (v_i^{2/9} + v_j^{2/9})^{1/2}
            \frac{\epsilon_c^{1/3}}{1 + \alpha_d}
            \text{exp}
            \left[
              - C_2 \frac{\mu_c \rho_c}{\sigma^2}
                \frac{\epsilon_c}{(1 + \alpha_d)^{3}}
                \left(
                    \frac{v_i^{1/3} v_j^{1/3}}{v_i^{1/3} + v_j^{1/3}}
                \right)^{4}
            \right]
        \f]

    where

    \vartable
        \sigma      |  Surface tension [N/m]
        v_i         |  Volume of droplet i [m3]
        v_j         |  Volume of droplet j [m3]
        \epsilon_c  |  Turbulent dissipation rate of continuous phase [m2/s3]
        \alpha_d    |  Total void fraction of disperse phase [-]
        \mu_c       |  Molecular dynamic viscosity of liquid phase [Pa s]
        \rho_c      |  Density of continous phase [kg/m3]
    \endvartable

    References:
    \verbatim
        "Description of interaction processes in agitated liquid-liquid
         dispersions"
        Coulaloglou, C.A., Tavlarides, L.L.
        Chemical Engineering Science, Vol. 32, Issue 11, 1977, pp. 1289-1297
        Eq. 37, p. 1294
    \endverbatim

Usage
    \table
        Property     | Description             | Required    | Default value
        C1           | Coefficient C1          | no          | 2.8
        C2           | Coefficient C2          | no          | 1.83e9
    \endtable

SourceFiles
    CoulaloglouTavlaridesCoalescence.C

\*---------------------------------------------------------------------------*/

#ifndef CoulaloglouTavlaridesCoalescence_H
#define CoulaloglouTavlaridesCoalescence_H

#include "coalescenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace coalescenceModels
{

/*---------------------------------------------------------------------------*\
              Class CoulaloglouTavlaridesCoalescence Declaration
\*---------------------------------------------------------------------------*/

class CoulaloglouTavlaridesCoalescence
:
    public coalescenceModel
{
    // Private data

        //- Optional coefficient C1, defaults to 2.8
        dimensionedScalar C1_;

        //- Optional coefficient C2, defaults to 1.83e9
        dimensionedScalar C2_;


public:

    //- Runtime type information
    TypeName("CoulaloglouTavlarides");

    // Constructor

        CoulaloglouTavlaridesCoalescence
        (
            const populationBalanceModel& popBal,
            const dictionary& dict
        );


    //- Destructor
    virtual ~CoulaloglouTavlaridesCoalescence()
    {}


    // Member Functions

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
