/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2025 OpenFOAM Foundation
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
    Foam::diameterModels::SecondaryPropertyModel

Description
    Base class for modelling evolution of secondary representative properties
    of a size class. By default, transport between size classes due to
    coalescence and breakup conserve the property. This behaviour can be
    overridden in a derived class.

See also
    Foam::diameterModels::sizeGroup
    Foam::diameterModels::populationBalanceModel

SourceFiles
    SecondaryPropertyModel.C

\*---------------------------------------------------------------------------*/

#ifndef SecondaryPropertyModel_H
#define SecondaryPropertyModel_H

#include "sizeGroup.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{

/*---------------------------------------------------------------------------*\
                   Class SecondaryPropertyModel Declaration
\*---------------------------------------------------------------------------*/

template<class ModelType>
class SecondaryPropertyModel
:
    public ModelType
{
private:

    // Private Member Functions

        //- Return the secondary property model of a given size group
        const SecondaryPropertyModel<ModelType>& correspondingModel
        (
            const sizeGroup& fi
        ) const;


public:

    //- Runtime type information
    virtual const word& type() const = 0;


    // Constructors

        //- Construct from sizeGroup
        SecondaryPropertyModel(const sizeGroup& group);

        //- Disallow default bitwise copy construction
        SecondaryPropertyModel
        (
            const SecondaryPropertyModel<ModelType>&
        ) = delete;


    //- Destructor
    virtual ~SecondaryPropertyModel();


    // Member Functions

        // Access

            //- Return the secondary property model of a given size group
            static const SecondaryPropertyModel<ModelType>& model
            (
                const sizeGroup& fi
            );

            //- Return reference to secondary property field
            virtual const volScalarField& fld() const = 0;

            //- Access to secondary property source
            virtual volScalarField::Internal& src() = 0;


        // Edit

            //- Add coalescence contribution to secondary property source
            virtual void addCoalescence
            (
                const volScalarField::Internal& Su,
                const sizeGroup& fj,
                const sizeGroup& fk
            );

            //- Add breakup contribution to secondary property source
            virtual void addBreakup
            (
                const volScalarField::Internal& Su,
                const sizeGroup& fj
            );

            //- Reset secondary property source
            virtual void reset();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const SecondaryPropertyModel<ModelType>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SecondaryPropertyModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
