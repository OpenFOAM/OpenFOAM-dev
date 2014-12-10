/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::regionModels::pyrolysisModels::pyrolysisModel

Description
    Base class for pyrolysis models

SourceFiles
    pyrolysisModelI.H
    pyrolysisModel.C

\*---------------------------------------------------------------------------*/

#ifndef pyrolysisModel_H
#define pyrolysisModel_H

#include "runTimeSelectionTables.H"
#include "volFieldsFwd.H"
#include "regionModel1D.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;
class Time;

namespace regionModels
{
namespace pyrolysisModels
{

/*---------------------------------------------------------------------------*\
                      Class pyrolysisModel Declaration
\*---------------------------------------------------------------------------*/

class pyrolysisModel
:
    public regionModel1D
{
private:

    // Private Member Functions

        //- Construct fields
        void constructMeshObjects();

        //- Read pyrolysis controls
        void readPyrolysisControls();

        //- Disallow default bitwise copy construct
        pyrolysisModel(const pyrolysisModel&);

        //- Disallow default bitwise assignment
        void operator=(const pyrolysisModel&);


protected:

    // Protected Member Functions

        //- Read control parameters
        virtual bool read();

        //- Read control parameters from dictionary
        virtual bool read(const dictionary& dict);


public:

    //- Runtime type information
    TypeName("pyrolysisModel");


    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            pyrolysisModel,
            mesh,
            (
                const word& modelType,
                const fvMesh& mesh,
                const word& regionType
            ),
            (modelType, mesh, regionType)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            pyrolysisModel,
            dictionary,
            (
                const word& modelType,
                const fvMesh& mesh,
                const dictionary& dict,
                const word& regionType
            ),
            (modelType, mesh, dict, regionType)
        );


    // Constructors

        //- Construct null from mesh
        pyrolysisModel
        (
            const fvMesh& mesh,
            const word& regionType
        );

        //- Construct from type name and mesh
        pyrolysisModel
        (
            const word& modelType,
            const fvMesh& mesh,
            const word& regionType
        );

        //- Construct from type name and mesh and dictionary
        pyrolysisModel
        (
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict,
            const word& regionType
        );

        //- Return clone
        autoPtr<pyrolysisModel> clone() const
        {
            notImplemented("autoPtr<pyrolysisModel> clone() const");
            return autoPtr<pyrolysisModel>(NULL);
        }


    // Selectors

        //- Return a reference to the selected pyrolysis model
        static autoPtr<pyrolysisModel> New
        (
            const fvMesh& mesh,
            const word& regionType = "pyrolysis"
        );

        //- Return a reference to a named selected pyrolysis model
        static autoPtr<pyrolysisModel> New
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word& regionType = "pyrolysis"
        );


    //- Destructor
    virtual ~pyrolysisModel();


    // Member Functions

        // Access

            // Fields

                //- Return density [kg/m3]
                virtual const volScalarField& rho() const = 0;

                //- Return const temperature [K]
                virtual const volScalarField& T() const = 0;

                //- Return specific heat capacity [J/kg/K]
                virtual const tmp<volScalarField> Cp() const = 0;

                //- Return the region absorptivity [1/m]
                virtual tmp<volScalarField> kappaRad() const = 0;

                //- Return the region thermal conductivity [W/m/k]
                virtual tmp<volScalarField> kappa() const = 0;

                //- Return the total gas mass flux to primary region [kg/m2/s]
                virtual const surfaceScalarField& phiGas() const = 0;


            // Sources

                //- External hook to add mass to the primary region
                virtual scalar addMassSources
                (
                    const label patchI,
                    const label faceI
                );


        // Helper function

            //- Mean diffusion number of the solid region
            virtual scalar solidRegionDiffNo() const;

            //- Return max diffusivity allowed in the solid
            virtual scalar maxDiff() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace pyrolysisModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
