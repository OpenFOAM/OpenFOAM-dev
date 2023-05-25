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
    Foam::SidedInterfacialModel

Description

SourceFiles
    SidedInterfacialModel.C

\*---------------------------------------------------------------------------*/

#ifndef SidedInterfacialModel_H
#define SidedInterfacialModel_H

#include "phaseInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class SidedInterfacialModel Declaration
\*---------------------------------------------------------------------------*/

template<class ModelType>
class SidedInterfacialModel
:
    public regIOobject
{
    // Private Data

        //- The interface
        const phaseInterface interface_;

        //- Model for within phase 1
        autoPtr<ModelType> modelInPhase1_;

        //- Model for within phase 2
        autoPtr<ModelType> modelInPhase2_;


public:

    //- Runtime type information
    TypeName("SidedInterfacialModel");


    // Constructors

        //- Construct from a dictionary and an interface
        SidedInterfacialModel
        (
            const dictionary& dict,
            const phaseInterface& interface
        );

        //- Disallow default bitwise copy construction
        SidedInterfacialModel
        (
            const SidedInterfacialModel<ModelType>&
        ) = delete;


    //- Destructor
    ~SidedInterfacialModel();


    // Member Functions

        //- Access the interface
        const phaseInterface& interface() const;

        //- Does a model exist in the given phase?
        bool haveModelInThe(const phaseModel& phase) const;

        //- Access the model within the given phase
        const ModelType& modelInThe(const phaseModel& phase) const;

        //- Access the model within the given phase
        ModelType& modelInThe(const phaseModel& phase);

        //- Dummy write for regIOobject
        bool writeData(Ostream& os) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const SidedInterfacialModel<ModelType>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define defineSidedInterfacialModelTypeNameAndDebug(ModelType, DebugSwitch)    \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        SidedInterfacialModel<ModelType>,                                      \
        (                                                                      \
            word(SidedInterfacialModel<ModelType>::typeName_()) + "<"          \
          + ModelType::typeName + ">"                                          \
        ).c_str(),                                                             \
        DebugSwitch                                                            \
    );


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SidedInterfacialModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
