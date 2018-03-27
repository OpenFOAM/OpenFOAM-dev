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
    Foam::dragModels::AttouFerschneider

Description
    Attou and Ferschneider's Drag model for film flow through packed beds. The
    implementation follows the desciption of Gunjal and Ranade, who, in the
    reference below, formulate the model in more convenient terms.

    Reference:
    \verbatim
        "Modeling of laboratory and commercial scale hydro-processing reactors
        using CFD"
        Gunjal, P.R., Ranade, V.V.,
        Chemical Engineering Science
        Volume 62, Issues 18-20, September-October 2007, pp. 5512 - 5526
    \endverbatim

SourceFiles
    AttouFerschneider.C

\*---------------------------------------------------------------------------*/

#ifndef AttouFerschneider_H
#define AttouFerschneider_H

#include "dragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;
class phaseModel;

namespace dragModels
{

/*---------------------------------------------------------------------------*\
                           Class AttouFerschneider Declaration
\*---------------------------------------------------------------------------*/

class AttouFerschneider
:
    public dragModel
{
    // Private data

        //- Name of the gaseous phase
        const word gasName_;

        //- Name of the liquidphase
        const word liquidName_;

        //- Name of the solid phase
        const word solidName_;

        //- Ergun constant 1
        const dimensionedScalar E1_;

        //- Ergun constant 2
        const dimensionedScalar E2_;


    // Private member functions

        //- Return the momentum transfer coefficient between gas and liquid
        virtual tmp<volScalarField> KGasLiquid
        (
            const phaseModel& gas,
            const phaseModel& liquid
        ) const;

        //- Return the momentum transfer coefficient between gas and solid
        virtual tmp<volScalarField> KGasSolid
        (
            const phaseModel& gas,
            const phaseModel& solid
        ) const;

        //- Return the momentum transfer coefficient between liquid and solid
        virtual tmp<volScalarField> KLiquidSolid
        (
            const phaseModel& liquid,
            const phaseModel& solid
        ) const;


public:

    //- Runtime type information
    TypeName("AttouFerschneider");


    // Constructors

        //- Construct from a dictionary and a phase pair
        AttouFerschneider
        (
            const dictionary& dict,
            const phasePair& pair,
            const bool registerObject
        );


    //- Destructor
    virtual ~AttouFerschneider();


    // Member Functions

        //- Drag coefficient
        virtual tmp<volScalarField> CdRe() const;

        //- The drag coefficient used in the momentum equation
        virtual tmp<volScalarField> K() const;

        //- The drag coefficient used in the face-momentum equations
        virtual tmp<surfaceScalarField> Kf() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace dragModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
