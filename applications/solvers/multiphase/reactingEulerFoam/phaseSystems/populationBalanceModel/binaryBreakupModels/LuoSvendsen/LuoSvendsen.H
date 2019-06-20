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
    Foam::diameterModels::binaryBreakupModels::LuoSvendsen

Description
    Model of Luo and Svendsen (1996). The breakup rate is calculated by

    \f[
        C_4 \alpha_c \left(\frac{\epsilon_c}{d_j^2}\right)^{1/3}
        \int\limits_{\xi_{min}}^{1}
        \frac{\left(1 + \xi\right)^{2}}{\xi^{11/3}}
        \mathrm{exp}
        \left(
          - \frac{12c_f\sigma}{\beta\rho_c\epsilon_c^{2/3}d_j^{5/3}\xi^{11/3}}
        \right)
        \mathrm{d} \xi
    \f]

    where

    \f[
        c_f = \left(\frac{v_i}{v_j}\right)^{2/3}
          + \left(1 - \frac{v_i}{v_j}\right)^{2/3} - 1
    \f]

    \f[
        \xi_{min} = \frac{\lambda_{min}}{d_j}\,,
    \f]

    and

    \f[
        \lambda_{min} = C_5 \eta\,.
    \f]

    The integral in the first expression is solved by means of incomplete Gamma
    functions as given by Bannari et al. (2008):

    \f[
        \frac{3}{11 b^{8/11}}
        \left(
            \left[\Gamma(8/11, b) - \Gamma(8/11, t_{m})\right]
          + 2b^{3/11} \left[\Gamma(5/11, b) - \Gamma(5/11, t_{m})\right]
          + b^{6/11} \left[\Gamma(2/11, b) - \Gamma(2/11, t_{m})\right]
        \right)
    \f]

    where

    \f[
        b = \frac{12c_f\sigma}{\beta\rho_c\epsilon_c^{2/3}d_j^{5/3}}
    \f]

    and

    \f[
        t_{min} = b \xi_{min}^{-11/3}\,.
    \f]

    Note that in the code, the upper incomplete gamma function is expressed as

    \f[
        \Gamma(a,z) = Q(a,z) \Gamma(a)
    \f]

    \vartable
        \alpha_c    |  Void fraction of continuous phase [-]
        \epsilon_c  |  Turbulent dissipation rate of continuous phase [m2/s^3]
        d_j         |  Diameter of mother bubble j [m^3]
        v_i         |  Volume of daughter bubble i [m^3]
        v_j         |  Volume of mother bubble j [m^3]
        \xi         |  Integration variable [-]
        \xi_{min}   |  Lower bound of integral [-]
        c_f         |  Increase coefficient of surface area [-]
        \sigma      |  Surface tension [N/m]
        \rho_c      |  Density of continuous phase [kg/m^3]
        \eta        |  Kolmogorov length scale [m]
        \Gamma(a,z) |  Upper incomplete gamma function
        Q(a,z)      |  Regularized upper incomplete gamma function
        \Gamma(a)   |  Gamma function
    \endvartable

    References:
    \verbatim
        Luo, H., & Svendsen, H. F. (1996).
        Theoretical model for drop and bubble breakup in turbulent dispersions.
        AIChE Journal, 42(5), 1225-1233.
        Eq. 27, p. 1229.
    \endverbatim

    \verbatim
        Bannari, R., Kerdouss, F., Selma, B., Bannari, A., & Proulx, P. (2008).
        Three-dimensional mathematical modeling of dispersed two-phase flow
        using class method of population balance in bubble columns.
        Computers & chemical engineering, 32(12), 3224-3237.
        Eq. 49, p. 3230.
    \endverbatim

Usage
    \table
        Property     | Description             | Required    | Default value
        C4           | Coefficient C4          | no          | 0.923
        beta         | Coefficient beta        | no          | 2.05
        C5           | Minimum eddy ratio      | no          | 11.4
    \endtable

SourceFiles
    LuoSvendsen.C

\*---------------------------------------------------------------------------*/

#ifndef LuoSvendsen_H
#define LuoSvendsen_H

#include "binaryBreakupModel.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace binaryBreakupModels
{

/*---------------------------------------------------------------------------*\
                         Class LuoSvendsen Declaration
\*---------------------------------------------------------------------------*/

class LuoSvendsen
:
    public binaryBreakupModel
{
private:

    // Private Data

        //- Interpolation table of Q(a,z) for a=2/11
        autoPtr<interpolationTable<scalar>> gammaUpperReg2by11_;

        //- Interpolation table of Q(a,z) for a=5/11
        autoPtr<interpolationTable<scalar>> gammaUpperReg5by11_;

        //- Interpolation table of Q(a,z) for a=8/11
        autoPtr<interpolationTable<scalar>> gammaUpperReg8by11_;

        //- Empirical constant, defaults to 0.923
        dimensionedScalar C4_;

        //- Empirical constant, defaults to 2.05
        dimensionedScalar beta_;

        //- Ratio between minimum size of eddies in the inertial subrange
        //  and Kolmogorov length scale, defaults to 11.4
        dimensionedScalar minEddyRatio_;

        //- Kolmogorov length scale
        volScalarField kolmogorovLengthScale_;


public:

    //- Runtime type information
    TypeName("LuoSvendsen");

    // Constructor

        LuoSvendsen
        (
            const populationBalanceModel& popBal,
            const dictionary& dict
        );


    //- Destructor
    virtual ~LuoSvendsen()
    {}


    // Member Functions

        //- Correct diameter independent expressions
        virtual void correct();

        //- Add to binary breakupRate
        virtual void addToBinaryBreakupRate
        (
            volScalarField& binaryBreakupRate,
            const label i,
            const label j
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace binaryBreakupModels
} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
