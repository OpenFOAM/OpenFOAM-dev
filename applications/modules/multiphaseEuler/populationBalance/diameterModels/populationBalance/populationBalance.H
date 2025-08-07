/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2025 OpenFOAM Foundation
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
    Foam::diameterModels::populationBalance

Description
    This diameter model computes the diameter from multiple size group
    diameters and fractions provided by an associated population balance model.

    This model requires the name of the population balance model to be
    specified, and will trigger construction of said population balance model.
    The number of groups associated with the phase is also specified here. The
    sizes/diameters themselves, as well as all models relating to evolution of
    the size distribution (i.e., coalescence and breakup) are specified for
    population balance model and not for this diameter model.

Usage
    Excerpt from an exemplary phaseProperties dictionary:
    \verbatim
    diameterModel
    {
        type                populationBalance;

        populationBalance   bubbles;

        nGroups             10;
    }
    \endverbatim

See also
    Foam::populationBalanceModel

SourceFiles
    populationBalance.C

\*---------------------------------------------------------------------------*/

#ifndef populationBalance_H
#define populationBalance_H

#include "diameterModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class populationBalanceModel;

namespace diameterModels
{

/*---------------------------------------------------------------------------*\
                        Class populationBalance Declaration
\*---------------------------------------------------------------------------*/

class populationBalance
:
    public diameterModel
{
    // Private Data

        //- Name of the populationBalance this populationBalance belongs to
        word popBalName_;

        //- Pointer to the populationBalance this populationBalance belongs to
        mutable const populationBalanceModel* popBalPtr_;

        //- The number of groups associated with this phase
        const label nGroups_;

        //- Index of the first group of this phase
        mutable label iFirst_;

        //- Sauter-mean diameter of the phase
        volScalarField d_;


public:

    //- Runtime type information
    TypeName("populationBalance");


    // Constructors

        //- Construct from dictionary and phase
        populationBalance
        (
            const dictionary& diameterProperties,
            const phaseModel& phase
        );


    //- Destructor
    virtual ~populationBalance();


    // Member Functions

        //- Return name of populationBalance this populationBalance belongs to
        inline const word& popBalName() const;

        //- Return the populationBalance this populationBalance belongs to
        const populationBalanceModel& popBal() const;

        //- Return the number of groups in this phase
        inline label nGroups() const;

        //- Return the index of the first group of this phase
        label iFirst() const;

        //- Return the index of the last group of this phase
        inline label iLast() const;

        //- Get the diameter field
        virtual tmp<volScalarField> d() const;

        //- Get the surface area per unit volume field
        virtual tmp<volScalarField> Av() const;

        //- ...
        tmp<volScalarField::Internal> fSum() const;

        //- Correct the model
        virtual void correct();

        //- Read diameterProperties dictionary
        virtual bool read(const dictionary& diameterProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "populationBalanceI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
