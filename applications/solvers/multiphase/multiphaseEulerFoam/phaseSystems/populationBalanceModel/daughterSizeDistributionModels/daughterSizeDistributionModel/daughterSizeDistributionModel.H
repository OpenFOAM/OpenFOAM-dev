/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2020 OpenFOAM Foundation
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
    Foam::diameterModels::daughterSizeDistributionModel

Description
    Base class for daughter size distribution models. Calculates and stores
    the contribution to a sizeGroup i due to breakup in a larger sizeGroup k.

Usage
    \table
        Property     | Description              | Required    | Default value
        p            | daughter particle number | no          | 2
    \endtable

SourceFiles
    daughterSizeDistributionModel.C

\*---------------------------------------------------------------------------*/

#ifndef daughterSizeDistributionModel_H
#define daughterSizeDistributionModel_H

#include "dictionary.H"
#include "dimensionedScalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{

class breakupModel;

/*---------------------------------------------------------------------------*\
                Class daughterSizeDistributionModel Declaration
\*---------------------------------------------------------------------------*/

class daughterSizeDistributionModel
{
protected:

    // Protected Data

        //- Reference to the breakupModel
        const breakupModel& breakup_;

        //- Number of daughter particles, defaults to binary breakup
        dimensionedScalar particleNumber_;

        //- Contribution to sizeGroup i due to breakup in sizeGroup k
        PtrList<PtrList<dimensionedScalar>> nik_;


public:

    //- Runtime type information
    TypeName("daughterSizeDistributionModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            daughterSizeDistributionModel,
            dictionary,
            (
                const breakupModel& breakup,
                const dictionary& dict
            ),
            (breakup, dict)
        );


    // Constructor

        daughterSizeDistributionModel
        (
            const breakupModel& breakup,
            const dictionary& dict
        );


    //- Destructor
    virtual ~daughterSizeDistributionModel();


    // Selectors

        static autoPtr<daughterSizeDistributionModel> New
        (
            const breakupModel& breakup,
            const dictionary& dict
        );


    // Member Functions

        //- Return contribution to sizeGroup i due to breakup in sizeGroup k
        const dimensionedScalar& nik(const label i, const label k) const;

        //- Correct field independent expressions
        void correct();

        //- Calculate contribution to sizeGroup i due to breakup in sizeGroup k
        virtual dimensionedScalar calcNik
        (
            const label i,
            const label k
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
