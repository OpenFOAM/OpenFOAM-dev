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
    Foam::liftModels::TomiyamaLift

Description
    Lift model of Tomiyama et al.

    Reference:
    \verbatim
        Tomiyama, A., Tamai, H., Zun, I., & Hosokawa, S. (2002).
        Transverse migration of single bubbles in simple shear flows.
        Chemical Engineering Science, 57(11), 1849-1858.
    \endverbatim

    The coefficient for pow3(EoH) proposed by Tomiyama (2002) has been modified
    to make the model continuous at EoH = 10.7 while maintaining the
    lift coefficient proposed by Tomiyama (2002) when EoH >= 10.7.

SourceFiles
    TomiyamaLift.C

\*---------------------------------------------------------------------------*/

#ifndef TomiyamaLift_H
#define TomiyamaLift_H

#include "dispersedLiftModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class aspectRatioModel;

namespace liftModels
{

/*---------------------------------------------------------------------------*\
                        Class TomiyamaLift Declaration
\*---------------------------------------------------------------------------*/

class TomiyamaLift
:
    public dispersedLiftModel
{
private:

    //- Private Data

        //- The aspect ratio model
        autoPtr<aspectRatioModel> aspectRatio_;


public:

    //- Runtime type information
    TypeName("Tomiyama");


    // Constructors

        //- Construct from a dictionary and an interface
        TomiyamaLift
        (
            const dictionary& dict,
            const phaseInterface& interface
        );


    //- Destructor
    virtual ~TomiyamaLift();


    // Member Functions

        //- Lift coefficient
        virtual tmp<volScalarField> Cl() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace liftModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
