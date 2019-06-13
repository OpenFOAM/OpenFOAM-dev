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
    Foam::swarmCorrections::TomiyamaSwarm

Description
    Swarm correction of Tomiyama et al.

    Reference:
    \verbatim
        "Drag Coefficients of Bubbles. 2nd Report. Drag Coefficient for a Swarm
         of Bubbles and its Applicability to Transient Flow."
        Tomiyama, A., Kataoka, I., Fukuda, T., and Sakaguchi, T.,
        Nippon Kikai Gakkai Ronbunshu
        Volume 61, Issue 588, 1995, pp. 2810-2817
    \endverbatim

SourceFiles
    TomiyamaSwarm.C

\*---------------------------------------------------------------------------*/

#ifndef TomiyamaSwarm_H
#define TomiyamaSwarm_H

#include "swarmCorrection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace swarmCorrections
{

/*---------------------------------------------------------------------------*\
                         Class TomiyamaSwarm Declaration
\*---------------------------------------------------------------------------*/

class TomiyamaSwarm
:
    public swarmCorrection
{
    // Private Data

        //- Residual phase fraction
        const dimensionedScalar residualAlpha_;

        //- Constant exponent
        const dimensionedScalar l_;


public:

    //- Runtime type information
    TypeName("Tomiyama");


    // Constructors

        //- Construct from a dictionary and a phase pair
        TomiyamaSwarm
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~TomiyamaSwarm();


    // Member Functions

        //- Swarm correction coefficient
        virtual tmp<volScalarField> Cs() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace swarmCorrections
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
