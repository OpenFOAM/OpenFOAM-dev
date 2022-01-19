/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2022 OpenFOAM Foundation
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
    This diameterModel is intended for use with a populationBalanceModel in
    order to simulate polydispersed bubbly or particulate flows. It can hold any
    number of sizeGroups from which the Sauter mean diameter is calculated. It
    can also be used as a diameterModel without a populationBalance and would
    then behave like a constantDiameter model. In this case, some arbitrary name
    must be entered for the populationBalance keyword.

Usage
    \table
        Property          | Description
        populationBalance | Name of the corresponding populationBalance
        sizeGroups        | List of sizeGroups
    \endtable

    Example
    \verbatim
    diameterModel   velocityGroup;
    velocityGroupCoeffs
    {
        populationBalance    bubbles;

        shapeModel           constant;

        sizeGroups
        (
            f0{dSph  1.00e-3; value 0;}
            f1{dSph  1.08e-3; value 0;}
            f2{dSph  1.16e-3; value 0.25;}
            f3{dSph  1.25e-3; value 0.5;}
            f4{dSph  1.36e-3; value 0.25;}
            f5{dSph  1.46e-3; value 0;}
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

        //- Sum of the sizeGroup volume fractions and reference field from which
        //  the sizeGroup fields are derived
        volScalarField f_;

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

        //- Return reference field for sizeGroup's
        inline const volScalarField& f() const;

        //- Return sizeGroups belonging to this velocityGroup
        inline const PtrList<sizeGroup>& sizeGroups() const;

        //- Get the diameter field
        virtual tmp<volScalarField> d() const;

        //- Get the surface area per unit volume field
        virtual tmp<volScalarField> a() const;

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
