/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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

#include "CentredFitSnGradData.H"
#include "surfaceFields.H"
#include "SVD.H"
#include "extendedCentredCellToFaceStencil.H"

// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

template<class Polynomial>
Foam::CentredFitSnGradData<Polynomial>::CentredFitSnGradData
(
    const fvMesh& mesh,
    const extendedCentredCellToFaceStencil& stencil,
    const scalar linearLimitFactor,
    const scalar centralWeight
)
:
    FitData
    <
        CentredFitSnGradData<Polynomial>,
        extendedCentredCellToFaceStencil,
        Polynomial
    >
    (
        mesh, stencil, true, linearLimitFactor, centralWeight
    ),
    coeffs_(mesh.nFaces())
{
    if (debug)
    {
        InfoInFunction
            << "Constructing CentredFitSnGradData<Polynomial>" << endl;
    }

    calcFit();

    if (debug)
    {
        Info<< "    Finished constructing polynomialFit data" << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Polynomial>
void Foam::CentredFitSnGradData<Polynomial>::calcFit
(
    scalarList& coeffsi,
    const List<point>& C,
    const scalar wLin,
    const scalar deltaCoeff,
    const label facei
)
{
    vector idir(1,0,0);
    vector jdir(0,1,0);
    vector kdir(0,0,1);
    this->findFaceDirs(idir, jdir, kdir, facei);

    // Setup the point weights
    scalarList wts(C.size(), scalar(1));
    wts[0] = this->centralWeight();
    wts[1] = this->centralWeight();

    // Reference point
    point p0 = this->mesh().faceCentres()[facei];

    // p0 -> p vector in the face-local coordinate system
    vector d;

    // Local coordinate scaling
    scalar scale = 1;

    // Matrix of the polynomial components
    scalarRectangularMatrix B(C.size(), this->minSize(), scalar(0));

    forAll(C, ip)
    {
        const point& p = C[ip];
        const vector p0p = p - p0;

        d.x() = p0p & idir;
        d.y() = p0p & jdir;
        d.z() = p0p & kdir;

        if (ip == 0)
        {
            scale = cmptMax(cmptMag((d)));
        }

        // Scale the radius vector
        d /= scale;

        Polynomial::addCoeffs(B[ip], d, wts[ip], this->dim());
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

        for (label i=0; i<stencilSize; i++)
        {
            coeffsi[i] = wts[1]*wts[i]*invB(1, i)/scale;
        }

        goodFit =
        (
            mag(wts[0]*wts[0]*invB(0, 0) - wLin)
          < this->linearLimitFactor()*wLin)
         && (mag(wts[0]*wts[1]*invB(0, 1) - (1 - wLin)
        ) < this->linearLimitFactor()*(1 - wLin))
         && coeffsi[0] < 0 && coeffsi[1] > 0
         && mag(coeffsi[0] + deltaCoeff) < 0.5*deltaCoeff
         && mag(coeffsi[1] - deltaCoeff) < 0.5*deltaCoeff;

        if (!goodFit)
        {
            // (not good fit so increase weight in the centre and weight
            //  for constant and linear terms)

            WarningInFunction
                << "Cannot fit face " << facei << " iteration " << iIt
                << " with sum of weights " << sum(coeffsi) << nl
                << "    Weights " << coeffsi << nl
                << "    Linear weights " << wLin << " " << 1 - wLin << nl
                << "    deltaCoeff " << deltaCoeff << nl
                << "    sing vals " << svd.S() << nl
                << "Components of goodFit:\n"
                << "    wts[0]*wts[0]*invB(0, 0) = "
                << wts[0]*wts[0]*invB(0, 0) << nl
                << "    wts[0]*wts[1]*invB(0, 1) = "
                << wts[0]*wts[1]*invB(0, 1)
                << " dim = " << this->dim() << endl;

            wts[0] *= 10;
            wts[1] *= 10;

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
        // Remove the uncorrected coefficients
        coeffsi[0] += deltaCoeff;
        coeffsi[1] -= deltaCoeff;
    }
    else
    {
        WarningInFunction
            << "Could not fit face " << facei
            << "    Coefficients = " << coeffsi
            << ", reverting to uncorrected." << endl;

        coeffsi = 0;
    }
}


template<class Polynomial>
void Foam::CentredFitSnGradData<Polynomial>::calcFit()
{
    const fvMesh& mesh = this->mesh();

    // Get the cell/face centres in stencil order.
    // Centred face stencils no good for triangles or tets.
    // Need bigger stencils
    List<List<point>> stencilPoints(mesh.nFaces());
    this->stencil().collectData(mesh.C(), stencilPoints);

    // find the fit coefficients for every face in the mesh

    const surfaceScalarField& w = mesh.surfaceInterpolation::weights();
    const surfaceScalarField& dC = mesh.nonOrthDeltaCoeffs();

    for (label facei = 0; facei < mesh.nInternalFaces(); facei++)
    {
        calcFit
        (
            coeffs_[facei],
            stencilPoints[facei],
            w[facei],
            dC[facei],
            facei
        );
    }

    const surfaceScalarField::Boundary& bw = w.boundaryField();
    const surfaceScalarField::Boundary& bdC = dC.boundaryField();

    forAll(bw, patchi)
    {
        const fvsPatchScalarField& pw = bw[patchi];
        const fvsPatchScalarField& pdC = bdC[patchi];

        if (pw.coupled())
        {
            label facei = pw.patch().start();

            forAll(pw, i)
            {
                calcFit
                (
                    coeffs_[facei],
                    stencilPoints[facei],
                    pw[i],
                    pdC[i],
                    facei
                );
                facei++;
            }
        }
    }
}


// ************************************************************************* //
