/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::XiGModels::basicSubGrid

Description

    Basic sub-grid obstacle flame-wrinking generation rate coefficient model.
    Details supplied by J Puttock 2/7/06.

    \f$ G_{sub} \f$ denotes the generation coefficient and it is given by

    \f[
        G_{sub} = k_{1} /frac{\vert \dwea{\vec{U}} \vert}{L_{obs}}
                 \frac{/Xi_{{sub}_{eq}}-1}{/Xi_{sub}}
    \f]

    and the removal:

    \f[
        - k_{1} /frac{\vert \dwea{\vec{U}} \vert}{L_{sub}}
        \frac{\Xi_{sub}-1}{\Xi_{sub}}
    \f]

    Finally, \f$ G_{sub} \f$ is added to generation rate \f$ G_{in} \f$
    due to the turbulence.

SourceFiles
    basicSubGrid.C

\*---------------------------------------------------------------------------*/

#ifndef basicXiSubG_H
#define basicXiSubG_H

#include "XiGModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace XiGModels
{

/*---------------------------------------------------------------------------*\
                        Class basicSubGrid Declaration
\*---------------------------------------------------------------------------*/

class basicSubGrid
:
    public XiGModel
{
    // Private data

        //- Sub-grid generation rate coefficient
        scalar k1;

        //- Xi generation rate model due to turbulence
        autoPtr<XiGModel> XiGModel_;


    // Private Member Functions

        //- Disallow copy construct
        basicSubGrid(const basicSubGrid&);

        //- Disallow default bitwise assignment
        void operator=(const basicSubGrid&);


public:

    //- Runtime type information
    TypeName("basicSubGridG");


    // Constructors

        //- Construct from components
        basicSubGrid
        (
            const dictionary& XiGProperties,
            const psiuReactionThermo& thermo,
            const compressible::RASModel& turbulence,
            const volScalarField& Su
        );


    //- Destructor
    virtual ~basicSubGrid();


    // Member Functions

        //- Return the flame-wrinking generation rate
        virtual tmp<volScalarField> G() const;

        //- Return the flame diffusivity
        virtual tmp<volScalarField> Db() const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& XiGProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace XiGModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
