/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
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
    Foam::diameterModels::breakupModel

Description
    Base class for breakup models which give a total breakup rate and a separate
    daughter size distribution function.

SourceFiles
    breakupModel.C

\*---------------------------------------------------------------------------*/

#ifndef breakupModel_H
#define breakupModel_H

#include "populationBalanceModel.H"
#include "daughterSizeDistributionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{

/*---------------------------------------------------------------------------*\
                        Class breakupModel Declaration
\*---------------------------------------------------------------------------*/

class breakupModel
{
protected:

    // Protected data

        //- Reference to the populationBalanceModel
        const populationBalanceModel& popBal_;

        //- Dictionary
        dictionary dict_;

        //- Daughter size distribution model
        autoPtr<daughterSizeDistributionModel> dsd_;


public:

    //- Runtime type information
    TypeName("breakupModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            breakupModel,
            dictionary,
            (
                const populationBalanceModel& popBal,
                const dictionary& dict
            ),
            (popBal, dict)
        );


    //- Class used for the read-construction of
    //  PtrLists of breakup models
    class iNew
    {
        const populationBalanceModel& popBal_;

    public:

        iNew(const populationBalanceModel& popBal)
        :
            popBal_(popBal)
        {}

        autoPtr<breakupModel> operator()(Istream& is) const
        {
            word type(is);
            dictionary dict(is);
            return breakupModel::New(type, popBal_, dict);
        }
    };


    // Constructor

        breakupModel
        (
            const populationBalanceModel& popBal,
            const dictionary& dict
        );

        autoPtr<breakupModel> clone() const
        {
            NotImplemented;
            return autoPtr<breakupModel>(nullptr);
        }


    // Selector

        static autoPtr<breakupModel> New
        (
            const word& type,
            const populationBalanceModel& popBal,
            const dictionary& dict
        );


    //- Destructor
    virtual ~breakupModel()
    {}


    // Member Functions

        //- Return reference to the populationBalanceModel
        const populationBalanceModel& popBal() const
        {
            return popBal_;
        }

        //- Return const-reference to daughter size distribution pointer
        const autoPtr<daughterSizeDistributionModel>& dsdPtr() const
        {
            return dsd_;
        }

        //- Return reference to daughter size distribution pointer
        autoPtr<daughterSizeDistributionModel>& dsdPtr()
        {
            return dsd_;
        }

        //- Correct diameter independent expressions
        virtual void correct();

        //- Set total breakupRate
        virtual void setBreakupRate
        (
            volScalarField& breakupRate,
            const label i
        ) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
