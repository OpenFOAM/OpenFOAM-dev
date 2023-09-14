/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2023 OpenFOAM Foundation
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
    Foam::diameterModels::velocityGroup

Description
    Computes the Sauter mean diameter based on a user specified size
    distribution, defined in terms of size class fractions. Intended for use
    with a population balance model to account for the evolution of a size
    distribution by means of coalescence, breakup, drift and nucleation.

Usage
    Excerpt from an exemplary phaseProperties dictionary:
    \verbatim
    diameterModel   velocityGroup;

    velocityGroupCoeffs
    {
        populationBalance    bubbles;

        shapeModel           spherical;

        sizeGroups
        (
            f1 {dSph  1e-3; value 1.0;}
            f2 {dSph  2e-3; value 0.0;}
            f3 {dSph  3e-3; value 0.0;}
            f4 {dSph  4e-3; value 0.0;}
            f5 {dSph  5e-3; value 0.0;}
            ...
        );
    }
    \endverbatim

See also
    Foam::diameterModels::sizeGroup
    Foam::diameterModels::populationBalanceModel

SourceFiles
    velocityGroup.C

\*---------------------------------------------------------------------------*/

#ifndef velocityGroup_H
#define velocityGroup_H

#include "diameterModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{

// Forward declaration of classes
class sizeGroup;

/*---------------------------------------------------------------------------*\
                        Class velocityGroup Declaration
\*---------------------------------------------------------------------------*/

class velocityGroup
:
    public diameterModel
{
    // Private Data

        //- Name of the populationBalance this velocityGroup belongs to
        word popBalName_;

        //- sizeGroups belonging to this velocityGroup
        PtrList<sizeGroup> sizeGroups_;

        //- Sauter-mean diameter of the phase
        volScalarField d_;


    // Private Member Functions

        tmp<volScalarField> dsm() const;

        tmp<volScalarField> fSum() const;

        tmp<volScalarField> N() const;

        void scale();


public:

    //- Runtime type information
    TypeName("velocityGroup");


    // Constructors

        //- Construct from dictionary and phase
        velocityGroup
        (
            const dictionary& diameterProperties,
            const phaseModel& phase
        );


    //- Destructor
    virtual ~velocityGroup();


    // Member Functions

        //- Return name of populationBalance this velocityGroup belongs to
        inline const word& popBalName() const;

        //- Return sizeGroups belonging to this velocityGroup
        inline const PtrList<sizeGroup>& sizeGroups() const;

        //- Return sizeGroups belonging to this velocityGroup
        inline PtrList<sizeGroup>& sizeGroups();

        //- Get the diameter field
        virtual tmp<volScalarField> d() const;

        //- Get the surface area per unit volume field
        virtual tmp<volScalarField> Av() const;

        //- Correct the model
        virtual void correct();

        //- Read diameterProperties dictionary
        virtual bool read(const dictionary& diameterProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "velocityGroupI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
