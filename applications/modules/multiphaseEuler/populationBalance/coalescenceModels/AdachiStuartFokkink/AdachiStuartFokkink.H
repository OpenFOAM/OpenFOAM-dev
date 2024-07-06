/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2024 OpenFOAM Foundation
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
    Foam::diameterModels::coalescenceModels::AdachiStuartFokkink

Description
    Model describing aggregation of solid particles in turbulent flows.
    Applicable when particles are smaller than the Kolmogorov length scale.
    The coalescence rate is calculated by

    \f[
        \frac{4}{3}\left(\frac{3\pi}{10}\right)^{0.5}
        \left(\frac{\epsilon}{\nu}\right)^{0.5}(d_i+d_j)^3
    \f]

    where

    \vartable
        d_i      | Diameter of particle i [m]
        d_j      | Diameter of particle j [m]
        \nu      | Kinematic viscosity of continuous phase [m2/s]
        \epsilon | Continuous phase turbulent dissipation rate [m2/s3]
    \endtable

    Reference:
    \verbatim
        Adachi, Y., Stuart, M. C., & Fokkink, R. (1994).
        Kinetics of turbulent coagulation studied by means of end-over-end
        rotation.
        Journal of colloid and interface science, 165(2), 310-317.
    \endverbatim

SourceFiles
    AdachiStuartFokkink.C

\*---------------------------------------------------------------------------*/

#ifndef AdachiStuartFokkink_H
#define AdachiStuartFokkink_H

#include "coalescenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace coalescenceModels
{

/*---------------------------------------------------------------------------*\
                     Class AdachiStuartFokkink Declaration
\*---------------------------------------------------------------------------*/

class AdachiStuartFokkink
:
    public coalescenceModel
{
public:

    //- Runtime type information
    TypeName("AdachiStuartFokkink");

    // Constructor

        AdachiStuartFokkink
        (
            const populationBalanceModel& popBal,
            const dictionary& dict
        );


    //- Destructor
    virtual ~AdachiStuartFokkink()
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
