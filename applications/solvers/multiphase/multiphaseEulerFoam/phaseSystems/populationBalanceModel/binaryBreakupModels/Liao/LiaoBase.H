/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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
    Foam::diameterModels::LiaoBase

Description
    Base class for coalescence and breakup models of Liao et al. (2015).

    Reference:
    \verbatim
        Liao, Y.; Rzehak, R.; Lucas, D.; Krepper, E. (2015).
        Baseline closure models for dispersed bubbly flow:
        Bubble coalescence and breakup.
        Chemical Engineering Science, 122, 336-349.
    \endverbatim

Usage

SourceFiles
    LiaoBase.C

\*---------------------------------------------------------------------------*/

#ifndef LiaoBase_H
#define LiaoBase_H

#include "populationBalanceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{

/*---------------------------------------------------------------------------*\
                          Class LiaoBase Declaration
\*---------------------------------------------------------------------------*/

class LiaoBase
{
protected:

    // Protected Data

        //- Reference to the populationBalanceModel
        const populationBalanceModel& populationBalance_;

        //- Kolmogorov length scale
        volScalarField kolmogorovLengthScale_;

        //- Shear strain rate
        volScalarField shearStrainRate_;

        //- Eddy strain rate
        volScalarField eddyStrainRate_;

        //- Terminal velocities
        PtrList<dimensionedScalar> uTerminal_;

        //- Drag coefficients
        PtrList<dimensionedScalar> Cd_;


public:

    // Constructor

        LiaoBase
        (
            const populationBalanceModel& popBal,
            const dictionary& dict
        );


    //- Destructor
    virtual ~LiaoBase()
    {}


    // Member Functions

        //- Precompute diameter independent expressions
        virtual void precompute();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
