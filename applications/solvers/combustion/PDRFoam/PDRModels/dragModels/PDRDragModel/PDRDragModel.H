/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::PDRDragModel

Description
    Base-class for sub-grid obstacle drag models. The available drag model is at
    \link basic.H \endlink.

SourceFiles
    PDRDragModel.C

\*---------------------------------------------------------------------------*/

#ifndef PDRDragModel_H
#define PDRDragModel_H

#include "IOdictionary.H"
#include "psiuReactionThermo.H"
#include "turbulentFluidThermoModel.H"
#include "multivariateSurfaceInterpolationScheme.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class PDRDragModel Declaration
\*---------------------------------------------------------------------------*/

class PDRDragModel
:
    public regIOobject
{

protected:

    // Protected data

        dictionary PDRDragModelCoeffs_;

        const compressible::RASModel& turbulence_;
        const volScalarField& rho_;
        const volVectorField& U_;
        const surfaceScalarField& phi_;

        Switch on_;


public:

    //- Runtime type information
    TypeName("PDRDragModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            PDRDragModel,
            dictionary,
            (
                const dictionary& PDRProperties,
                const compressible::RASModel& turbulence,
                const volScalarField& rho,
                const volVectorField& U,
                const surfaceScalarField& phi
            ),
            (
                PDRProperties,
                turbulence,
                rho,
                U,
                phi
            )
        );


    // Constructors

        //- Construct from components
        PDRDragModel
        (
            const dictionary& PDRProperties,
            const compressible::RASModel& turbulence,
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& phi
        );

        //- Disallow default bitwise copy construction
        PDRDragModel(const PDRDragModel&);


    // Selectors

        //- Return a reference to the selected Xi model
        static autoPtr<PDRDragModel> New
        (
            const dictionary& PDRProperties,
            const compressible::RASModel& turbulence,
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~PDRDragModel();


    // Member Functions

        //- Return true if the drag model is switched on
        bool on() const
        {
            return on_;
        }

        //- Return the momentum drag coefficient
        virtual tmp<volSymmTensorField> Dcu() const = 0;

        //- Return the momentum drag turbulence generation rate
        virtual tmp<volScalarField> Gk() const = 0;

        //- Inherit read from regIOobject
        using regIOobject::read;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& PDRProperties) = 0;

        virtual bool writeData(Ostream&) const
        {
            return true;
        }

        virtual void writeFields() const
        {
            NotImplemented;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const PDRDragModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
