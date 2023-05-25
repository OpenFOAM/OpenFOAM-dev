/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2023 OpenFOAM Foundation
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
    Foam::wallLubricationModels::Antal

Description
    Wall lubrication model of Antal et al.

    Reference:
    \verbatim
        Antal, S. P., Lahey Jr, R. T., & Flaherty, J. E. (1991).
        Analysis of phase distribution in fully developed laminar bubbly
        two-phase flow.
        International Journal of Multiphase Flow, 17(5), 635-652.
    \endverbatim

SourceFiles
    Antal.C

\*---------------------------------------------------------------------------*/

#ifndef Antal_H
#define Antal_H

#include "dispersedWallLubricationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace wallLubricationModels
{

/*---------------------------------------------------------------------------*\
                            Class Antal Declaration
\*---------------------------------------------------------------------------*/

class Antal
:
    public dispersedWallLubricationModel
{
    // Private Data

        //- Coefficient 1
        const dimensionedScalar Cw1_;

        //- Coefficient 2
        const dimensionedScalar Cw2_;


public:

    //- Runtime type information
    TypeName("Antal");


    // Constructors

        //- Construct from a dictionary and an interface
        Antal
        (
            const dictionary& dict,
            const phaseInterface& interface
        );


    //- Destructor
    virtual ~Antal();


    // Member Functions

        //- Return phase-intensive wall lubrication force
        tmp<volVectorField> Fi() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace wallLubricationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
