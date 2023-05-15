/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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
    Foam::fv::filmCloudTransfer

Description
    Film<->cloud transfer model

Usage
    Example usage:
    \verbatim
    filmCloudTransfer
    {
        type    filmCloudTransfer;

        libs    ("libfilmCloudTransfer.so");
    }
    \endverbatim

SourceFiles
    filmCloudTransfer.C

\*---------------------------------------------------------------------------*/

#ifndef filmCloudTransfer_H
#define filmCloudTransfer_H

#include "fvModel.H"
#include "isothermalFilm.H"
#include "ejectionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class filmCloudTransfer Declaration
\*---------------------------------------------------------------------------*/

class filmCloudTransfer
:
    public fvModel
{
    // Private Data

        //- The film model
        const solvers::isothermalFilm& film_;

        //- Current time index (used for updating)
        mutable label curTimeIndex_;

        // Transfers from cloud

            scalarField massFromCloud_;
            vectorField momentumFromCloud_;
            scalarField pressureFromCloud_;
            scalarField energyFromCloud_;

        //- Pointer to the ejection model
        autoPtr<ejectionModel> ejection_;


    // Private Member Functions

        //- Cloud to film rate transfer function
        template<class Type>
        tmp<VolInternalField<Type>>
        inline CloudToFilmTransferRate
        (
            const Field<Type>& prop,
            const dimensionSet& dimProp
        ) const;

        //- Film to cloud transfer function
        template<class Type>
        tmp<Field<Type>>
        inline filmToCloudTransfer
        (
            const VolInternalField<Type>& prop
        ) const;


public:

    //- Runtime type information
    TypeName("filmCloudTransfer");


    // Constructors

        //- Construct from explicit source name and mesh
        filmCloudTransfer
        (
            const word& sourceName,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        filmCloudTransfer
        (
            const filmCloudTransfer&
        ) = delete;


    // Member Functions

        // Checks

            //- Return the list of fields for which the option adds source term
            //  to the transport equation
            virtual wordList addSupFields() const;


        // Correct

            //- Solve the film and update the sources
            virtual void correct();


        // Add explicit and implicit contributions to compressible equation

            //- Add explicit droplet impingement contribution to pressure field
            virtual void addSup
            (
                fvMatrix<scalar>& eqn,
                const word& fieldName
            ) const;

            //- Add explicit contribution to phase continuity
            virtual void addSup
            (
                const volScalarField& alpha,
                fvMatrix<scalar>& eqn,
                const word& fieldName
            ) const;

            //- Add explicit contribution to phase energy equation
            virtual void addSup
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                fvMatrix<scalar>& eqn,
                const word& fieldName
            ) const;

            //- Add implicit contribution to mixture momentum equation
            virtual void addSup
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const word& fieldName
            ) const;


        // Transfer from cloud

            //- Reset the fields accumulated cloud transfer fields
            void resetFromCloudFields();

            //- Transfer parcel properties from cloud to the film
            //  accumulated fields
            void parcelFromCloud
            (
                const label facei,
                const scalar mass,
                const vector& momentum,
                const scalar pressure,
                const scalar energy
            );


        // Transfer to cloud

            //- Return true if the film is ejecting to the cloud
            bool ejecting() const;

            //- Transfer the ejected mass to the cloud
            tmp<Field<scalar>> ejectedMassToCloud() const;

            //- Transfer the ejected droplet diameter to the cloud
            tmp<Field<scalar>> ejectedDiameterToCloud() const;

            //- Transfer the film delta field to the cloud
            tmp<Field<scalar>> deltaToCloud() const;

            //- Transfer the film velocity field to the cloud
            tmp<Field<vector>> UToCloud() const;

            //- Transfer the film density field to the cloud
            tmp<Field<scalar>> rhoToCloud() const;

            //- Transfer the film temperature field to the cloud
            tmp<Field<scalar>> TToCloud() const;

            //- Transfer the film heat capacity field to the cloud
            tmp<Field<scalar>> CpToCloud() const;


        // Mesh changes

            //- Update topology using the given map
            virtual void topoChange(const polyTopoChangeMap&);

            //- Update from another mesh using the given map
            virtual void mapMesh(const polyMeshMap&);

            //- Redistribute or update using the given distribution map
            virtual void distribute(const polyDistributionMap&);

            //- Update for mesh motion
            virtual bool movePoints();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const filmCloudTransfer&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
