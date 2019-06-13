/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2019 OpenFOAM Foundation
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
    Foam::turbulentDispersionModels::LopezDeBertodano

Description
    Lopez de Bertodano (1992) turbulent dispersion model.

    \verbatim
        Lopez, D. B. M. (1993).
        Turbulent bubbly two-phase flow in a triangular duct.
        PhD Thesis, Rensselaer Polytechnic Institution.
    \endverbatim

    \verbatim
        Burns, A. D., Frank, T., Hamill, I., & Shi, J. M. (2004, May).
        The Favre averaged drag model for turbulent dispersion in Eulerian
        multi-phase flows.
        In 5th international conference on multiphase flow,
        ICMF (Vol. 4, pp. 1-17).
    \endverbatim

SourceFiles
    LopezDeBertodano.C

\*---------------------------------------------------------------------------*/

#ifndef LopezDeBertodano_H
#define LopezDeBertodano_H

#include "turbulentDispersionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace turbulentDispersionModels
{

/*---------------------------------------------------------------------------*\
              Class LopezDeBertodano Declaration
\*---------------------------------------------------------------------------*/

class LopezDeBertodano
:
    public turbulentDispersionModel
{
    // Private Data

        //- Constant turbulent dispersion coefficient
        const dimensionedScalar Ctd_;


public:

    //- Runtime type information
    TypeName("LopezDeBertodano");


    // Constructors

        //- Construct from a dictionary and a phase pair
        LopezDeBertodano
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~LopezDeBertodano();


    // Member Functions

        //- Turbulent diffusivity
        //  multiplying the gradient of the phase-fraction
        virtual tmp<volScalarField> D() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace turbulentDispersionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
