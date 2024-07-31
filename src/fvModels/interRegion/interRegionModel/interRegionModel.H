/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
    Foam::fv::interRegionModel

Description
    Base class for inter-region exchange.

\*---------------------------------------------------------------------------*/

#include "fvModel.H"
#include "volFields.H"
#include "cellsToCells.H"

#ifndef interRegionModel_H
#define interRegionModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                      Class interRegionModel Declaration
\*---------------------------------------------------------------------------*/

class interRegionModel
:
    public fvModel
{
    // Private data

        //- Master or slave region
        bool master_;

        //- Name of the neighbour region to map
        word nbrRegionName_;

        //- Interpolation method
        word interpolationMethod_;

        //- Mesh to mesh interpolation object
        mutable autoPtr<cellsToCells> interpolationPtr_;


    // Private member functions

        //- Non-virtual read
        void readCoeffs(const dictionary& dict);


protected:

    // Protected member functions

        //- Get the neighbour interRegionModel
        const interRegionModel& nbrModel() const;


public:

    //- Runtime type information
    TypeName("interRegionModel");


    // Constructors

        //- Construct from dictionary
        interRegionModel
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~interRegionModel();


    // Member Functions

        // Access

            //- Return whether the master region
            inline bool master() const;

            //- Return const access to the neighbour region name
            inline const word& nbrRegionName() const;

            //- Return const access to the neighbour mesh
            inline const fvMesh& nbrMesh() const;

            //- Return const access to the interpolation engine
            const cellsToCells& interpolation() const;


        // Interpolation

            //- Interpolate field
            template<class Type>
            tmp<Field<Type>> interpolate
            (
                const Field<Type>& field
            ) const;

            //- Interpolate field
            template<class Type>
            void interpolate
            (
                const Field<Type>& field,
                Field<Type>& result
            ) const;


        // IO

            //- Read dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "interRegionModelI.H"

#ifdef NoRepository
    #include "interRegionModelTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
