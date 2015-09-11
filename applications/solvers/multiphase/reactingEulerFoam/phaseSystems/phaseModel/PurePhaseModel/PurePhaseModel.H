/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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
    Foam::PurePhaseModel

Description
    Class which represents pure phases, i.e. without any species. Returns an
    empty list of mass fractions.

SourceFiles
    PurePhaseModel.C

\*---------------------------------------------------------------------------*/

#ifndef PurePhaseModel_H
#define PurePhaseModel_H

#include "phaseModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class phaseModel Declaration
\*---------------------------------------------------------------------------*/

template<class BasePhaseModel>
class PurePhaseModel
:
    public BasePhaseModel
{
protected:

    // Protected data

        //- Empty mass fraction field list
        PtrList<volScalarField> Y_;


public:

    // Constructors

        PurePhaseModel
        (
            const phaseSystem& fluid,
            const word& phaseName,
            const label index
        );


    //- Destructor
    virtual ~PurePhaseModel();


    // Member Functions

        //- Return the species fraction equation
        virtual tmp<fvScalarMatrix> YiEqn(volScalarField& Yi);


        // Thermo

            //- Return the species mass fractions
            virtual const PtrList<volScalarField>& Y() const;

            //- Access the species mass fractions
            virtual PtrList<volScalarField>& Y();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "PurePhaseModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
