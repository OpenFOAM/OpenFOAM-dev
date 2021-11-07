/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
    Foam::XiEqModels::basicSubGrid

Description
    Basic sub-grid obstacle flame-wrinkling enhancement factor model.
    Details supplied by J Puttock 2/7/06.

    <b> Sub-grid flame area generation </b>

    \f$ n = N - \hat{\dwea{\vec{U}}}.n_{s}.\hat{\dwea{\vec{U}}} \f$
    \f$ n_{r} = \sqrt{n} \f$

    where:

        \f$ \hat{\dwea{\vec{U}}} = \dwea{\vec{U}} / \vert \dwea{\vec{U}}
        \vert \f$

        \f$ b = \hat{\dwea{\vec{U}}}.B.\hat{\dwea{\vec{U}}} / n_{r} \f$

    where:

        \f$ B \f$ is the file "B".

        \f$ N \f$ is the file "N".

        \f$  n_{s} \f$ is the file "ns".

    The flame area enhancement factor \f$ \Xi_{sub} \f$ is expected to
    approach:

    \f[
        \Xi_{{sub}_{eq}} =
            1 + max(2.2 \sqrt{b}, min(0.34 \frac{\vert \dwea{\vec{U}}
            \vert}{{\vec{U}}^{'}}, 1.6)) \times min(\frac{n}{4}, 1)
    \f]


SourceFiles
    basicSubGrid.C

\*---------------------------------------------------------------------------*/

#ifndef basicXiSubXiEq_H
#define basicXiSubXiEq_H

#include "XiEqModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace XiEqModels
{

/*---------------------------------------------------------------------------*\
                        Class basicSubGrid Declaration
\*---------------------------------------------------------------------------*/

class basicSubGrid
:
    public XiEqModel
{
    // Private Data

        //- tblock
        volSymmTensorField B_;

        //- Equilibrium Xi model due to turbulence
        autoPtr<XiEqModel> XiEqModel_;


public:

    //- Runtime type information
    TypeName("basicSubGrid");


    // Constructors

        //- Construct from components
        basicSubGrid
        (
            const dictionary& XiEqProperties,
            const psiuReactionThermo& thermo,
            const compressible::RASModel& turbulence,
            const volScalarField& Su
        );

        //- Disallow default bitwise copy construction
        basicSubGrid(const basicSubGrid&);


    //- Destructor
    virtual ~basicSubGrid();


    // Member Functions

        //- Return the flame-wrinkling XiEq
        virtual tmp<volScalarField> XiEq() const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& XiEqProperties);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const basicSubGrid&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace XiEqModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
