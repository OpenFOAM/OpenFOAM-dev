/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "FitData.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "SVD.H"

// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

template<class Form, class ExtendedStencil, class Polynomial>
Foam::FitData<Form, ExtendedStencil, Polynomial>::FitData
(
    const fvMesh& mesh,
    const ExtendedStencil& stencil,
    const bool linearCorrection,
    const scalar linearLimitFactor,
    const scalar centralWeight
)
:
    MeshObject<fvMesh, Foam::MoveableMeshObject, Form>(mesh),
    stencil_(stencil),
    linearCorrection_(linearCorrection),
    linearLimitFactor_(linearLimitFactor),
    centralWeight_(centralWeight),
    #ifdef SPHERICAL_GEOMETRY
    dim_(2),
    #else
    dim_(mesh.nGeometricD()),
    #endif
    minSize_(Polynomial::nTerms(dim_))
{
    // Check input
    if (linearLimitFactor <= small || linearLimitFactor > 3)
    {
        FatalErrorInFunction
            << "linearLimitFactor requested = " << linearLimitFactor
            << " should be between zero and 3"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class FitDataType, class ExtendedStencil, class Polynomial>
void Foam::FitData<FitDataType, ExtendedStencil, Polynomial>::findFaceDirs
(
    vector& idir,        // value changed in return
    vector& jdir,        // value changed in return
    vector& kdir,        // value changed in return
    const label facei
)
{
    const fvMesh& mesh = this->mesh();

    idir = mesh.faceAreas()[facei];
    idir /= mag(idir);

    #ifndef SPHERICAL_GEOMETRY
    if (mesh.nGeometricD() <= 2) // find the normal direction
    {
        if (mesh.geometricD()[0] == -1)
        {
            kdir = vector(1, 0, 0);
        }
        else if (mesh.geometricD()[1] == -1)
        {
            kdir = vector(0, 1, 0);
        }
        else
        {
            kdir = vector(0, 0, 1);
        }
    }
    else // 3D so find a direction in the plane of the face
    {
        const face& f = mesh.faces()[facei];
        kdir = mesh.points()[f[0]] - mesh.faceCentres()[facei];
    }
    #else
    // Spherical geometry so kdir is the radial direction
    kdir = mesh.faceCentres()[facei];
    #endif

    if (mesh.nGeometricD() == 3)
    {
        // Remove the idir component from kdir and normalise
        kdir -= (idir & kdir)*idir;

        scalar magk = mag(kdir);

        if (magk < small)
        {
            FatalErrorInFunction
                << exit(FatalError);
        }
        else
        {
            kdir /= magk;
        }
    }

    jdir = kdir ^ idir;
}


template<class FitDataType, class ExtendedStencil, class Polynomial>
void Foam::FitData<FitDataType, ExtendedStencil, Polynomial>::calcFit
(
    scalarList& coeffsi,
    const List<point>& C,
    const scalar wLin,
    const label facei
)
{
    vector idir(1,0,0);
    vector jdir(0,1,0);
    vector kdir(0,0,1);
    findFaceDirs(idir, jdir, kdir, facei);

    // Setup the point weights
    scalarList wts(C.size(), scalar(1));
    wts[0] = centralWeight_;
    if (linearCorrection_)
    {
        wts[1] = centralWeight_;
    }

    // Reference point
    point p0 = this->mesh().faceCentres()[facei];

    // Info<< "Face " << facei << " at " << p0 << " stencil points at:\n"
    //     << C - p0 << endl;

    // p0 -> p vector in the face-local coordinate system
    vector d;

    // Local coordinate scaling
    scalar scale = 1;

    // Matrix of the polynomial components
    scalarRectangularMatrix B(C.size(), minSize_, scalar(0));

    forAll(C, ip)
    {
        const point& p = C[ip];

        d.x() = (p - p0)&idir;
        d.y() = (p - p0)&jdir;
        #ifndef SPHERICAL_GEOMETRY
        d.z() = (p - p0)&kdir;
        #else
        d.z() = mag(p) - mag(p0);
        #endif

        if (ip == 0)
        {
            scale = cmptMax(cmptMag((d)));
        }

        // Scale the radius vector
        d /= scale;

        Polynomial::addCoeffs
        (
            B[ip],
            d,
            wts[ip],
            dim_
        );
    }

    // Additional weighting for constant and linear terms
    for (label i = 0; i < B.m(); i++)
    {
        B(i, 0) *= wts[0];
        B(i, 1) *= wts[0];
    }

    // Set the fit
    label stencilSize = C.size();
    coeffsi.setSize(stencilSize);

    bool goodFit = false;
    for (int iIt = 0; iIt < 8 && !goodFit; iIt++)
    {
        SVD svd(B, small);
        scalarRectangularMatrix invB(svd.VSinvUt());

        scalar maxCoeff = 0;
        label maxCoeffi = 0;

        for (label i=0; i<stencilSize; i++)
        {
            coeffsi[i] = wts[0]*wts[i]*invB(0, i);
            if (mag(coeffsi[i]) > maxCoeff)
            {
                maxCoeff = mag(coeffsi[i]);
                maxCoeffi = i;
            }
        }

        if (linearCorrection_)
        {
            goodFit =
                (mag(coeffsi[0] - wLin) < linearLimitFactor_*wLin)
             && (mag(coeffsi[1] - (1 - wLin)) < linearLimitFactor_*(1 - wLin))
             && maxCoeffi <= 1;
        }
        else
        {
            // Upwind: weight on face is 1.
            goodFit =
                (mag(coeffsi[0] - 1.0) < linearLimitFactor_*1.0)
             && maxCoeffi <= 1;
        }

        // if (goodFit && iIt > 0)
        // {
            // Info<< "FitData<Polynomial>::calcFit"
            //     << "(const List<point>& C, const label facei" << nl
            //     << "Can now fit face " << facei << " iteration " << iIt
            //     << " with sum of weights " << sum(coeffsi) << nl
            //     << "    Weights " << coeffsi << nl
            //     << "    Linear weights " << wLin << " " << 1 - wLin << nl
            //     << "    sing vals " << svd.S() << endl;
        // }

        if (!goodFit) // (not good fit so increase weight in the centre and
                      //  weight for constant and linear terms)
        {
            // if (iIt == 7)
            // {
            //     WarningInFunction
            //         << "Cannot fit face " << facei << " iteration " << iIt
            //         << " with sum of weights " << sum(coeffsi) << nl
            //         << "    Weights " << coeffsi << nl
            //         << "    Linear weights " << wLin << " " << 1 - wLin << nl
            //         << "    sing vals " << svd.S() << endl;
            // }

            wts[0] *= 10;
            if (linearCorrection_)
            {
                wts[1] *= 10;
            }

            for (label j = 0; j < B.n(); j++)
            {
                B(0, j) *= 10;
                B(1, j) *= 10;
            }

            for (label i = 0; i < B.m(); i++)
            {
                B(i, 0) *= 10;
                B(i, 1) *= 10;
            }
        }
    }

    if (goodFit)
    {
        if (linearCorrection_)
        {
            // Remove the uncorrected linear coefficients
            coeffsi[0] -= wLin;
            coeffsi[1] -= 1 - wLin;
        }
        else
        {
            // Remove the uncorrected upwind coefficients
            coeffsi[0] -= 1.0;
        }
    }
    else
    {
        // if (debug)
        // {
            WarningInFunction
                << "Could not fit face " << facei
                << "    Weights = " << coeffsi
                << ", reverting to linear." << nl
                << "    Linear weights " << wLin << " " << 1 - wLin << endl;
        // }

        coeffsi = 0;
    }
}


template<class FitDataType, class ExtendedStencil, class Polynomial>
bool Foam::FitData<FitDataType, ExtendedStencil, Polynomial>::movePoints()
{
    calcFit();
    return true;
}

// ************************************************************************* //
