/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
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
    Foam::diameterModels::coalescenceModels::Luo

Description
    Model of Luo (1993). The coalescence rate is calculated by

    \f[
        \frac{\pi}{4} (d_i + d_j)^2 u_{ij}
        \mathrm{exp}
        \left[
          - C_1
            \frac
            {[0.75(1 + \xi_{ij}^2)(1 + \xi_{ij}^3)]^{1/2}}
            {(\rho_d/\rho_c + C_{vm})^{1/2} (1 + \xi_{ij})^3}
            \left(\frac{\rho_c d_i u_{ij}^2}{\sigma}\right)^{1/2}
        \right]\,,
    \f]

    where

    \f[
        u_{ij} = \sqrt{\beta} (\epsilon_c d_i)^{1/3} \sqrt{1 + \xi_{ij}^{-2/3}}
    \f]

    is the mean approach velocity of the bubbles and

    \f[
        \xi_{ij} = d_i/d_j
    \f]

    their size ratio.

    \vartable
        d_i         |  Diameter of bubble i [m]
        d_j         |  Diameter of bubble j [m]
        u_{ij}      |  Mean approach velocity [m/s]
        \xi_{ij}    |  Bubble size ratio [-]
        \rho_d      |  Density of dispersed phase [kg/m^3]
        \rho_c      |  Density of continuous phase [kg/m^3]
        \sigma      |  Surface tension [N/m]
        C_{vm}      |  Virtual mass coefficient [-]
        C_1         |  Coefficient [-]
        \beta       |  Coefficient [-]
        \epsilon_c  |  Continuous phase turbulent dissipation rate [m2/s^3]
    \endvartable

    Reference:
    \verbatim
        Luo, H. (1993).
        Coalescence, breakup and liquid circulation in bubble column reactors.
        Dr. Ing (Doctoral dissertation, Thesis, Department of Chemical
        Engineering, The Norwegian Institute of Technology, Trondheim, Norway).
    \endverbatim

Usage
    \table
        Property     | Description             | Required    | Default value
        beta         | Coefficient beta        | no          | 2.0
        C1           | Coefficient C1          | no          | 1.0
    \endtable

SourceFiles
    Luo.C

\*---------------------------------------------------------------------------*/

#ifndef Luo_H
#define Luo_H

#include "coalescenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace coalescenceModels
{

/*---------------------------------------------------------------------------*\
                             Class Luo Declaration
\*---------------------------------------------------------------------------*/

class Luo
:
    public coalescenceModel
{
    // Private Data

        //- Coefficient beta, defaults to 2.0
        dimensionedScalar beta_;

        //- Optional coefficient C1, defaults to 1.0
        dimensionedScalar C1_;


public:

    //- Runtime type information
    TypeName("Luo");

    // Constructor

        Luo
        (
            const populationBalanceModel& popBal,
            const dictionary& dict
        );


    //- Destructor
    virtual ~Luo()
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
