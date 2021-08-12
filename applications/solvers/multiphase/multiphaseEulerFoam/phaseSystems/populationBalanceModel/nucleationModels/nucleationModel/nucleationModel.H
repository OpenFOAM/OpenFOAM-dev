/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2021 OpenFOAM Foundation
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
    Foam::diameterModels::nucleationModel

Description
    Base class for nucleation models.

SourceFiles
    nucleationModel.C

\*---------------------------------------------------------------------------*/

#ifndef nucleationModel_H
#define nucleationModel_H

#include "populationBalanceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{

/*---------------------------------------------------------------------------*\
                       Class nucleationModel Declaration
\*---------------------------------------------------------------------------*/

class nucleationModel
{
protected:

    // Protected Data

        //- Reference to the populationBalanceModel
        const populationBalanceModel& popBal_;

        //- Model dictionary
        const dictionary dict_;

        //- Velocity group in which the nucleation occurs
        const velocityGroup& velGroup_;


public:

    //- Runtime type information
    TypeName("nucleationModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            nucleationModel,
            dictionary,
            (
                const populationBalanceModel& popBal,
                const dictionary& dict
            ),
            (popBal, dict)
        );


    //- Class used for the read-construction of
    //  PtrLists of nucleation models
    class iNew
    {
        const populationBalanceModel& popBal_;

    public:

        iNew(const populationBalanceModel& popBal)
        :
            popBal_(popBal)
        {}

        autoPtr<nucleationModel> operator()(Istream& is) const
        {
            word type(is);
            dictionary dict(is);
            return nucleationModel::New(type, popBal_, dict);
        }
    };


    // Constructor

        nucleationModel
        (
            const populationBalanceModel& popBal,
            const dictionary& dict
        );

        autoPtr<nucleationModel> clone() const
        {
            NotImplemented;
            return autoPtr<nucleationModel>(nullptr);
        }


    // Selector

        static autoPtr<nucleationModel> New
        (
            const word& type,
            const populationBalanceModel& popBal,
            const dictionary& dict
        );


    //- Destructor
    virtual ~nucleationModel()
    {}


    // Member Functions

        //- Return reference to the populationBalanceModel
        const populationBalanceModel& popBal() const
        {
            return popBal_;
        }

        //- Return reference to model dictionary
        const dictionary& dict() const
        {
            return dict_;
        }


        //- Precompute diameter independent expressions
        virtual void precompute();

        //- Add to nucleationRate
        virtual void addToNucleationRate
        (
            volScalarField& nucleationRate,
            const label i
        ) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
