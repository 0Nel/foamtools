/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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

#include "fishParticleCloud.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fishParticle::move(fishParticle::trackData& td)
{
    td.switchProcessor = false;
    td.keepParticle = true;

    const polyMesh& mesh = cloud().pMesh();
    const polyBoundaryMesh& pbMesh = mesh.boundaryMesh();

    scalar deltaT = mesh.time().deltaT().value();
    scalar tEnd = (1.0 - stepFraction())*deltaT;
    scalar dtMax = tEnd;

    while (td.keepParticle && !td.switchProcessor && tEnd > SMALL)
    {
        if (debug)
        {
            Info<< "Time = " << mesh.time().timeName()
                << " deltaT = " << deltaT
                << " tEnd = " << tEnd
                << " steptFraction() = " << stepFraction() << endl;
        }

        // set the lagrangian time-step
        scalar dt = min(dtMax, tEnd);

        // remember which cell the parcel is in
        // since this will change if a face is hit
        label celli = cell();

        dt *= trackToFace(position() + dt*U_, td);

        tEnd -= dt;
        stepFraction() = 1.0 - tEnd/deltaT;

        cellPointWeight cpw(mesh, position(), celli, face());
        scalar rhoc = td.rhoInterp().interpolate(cpw);
        vector Uc = td.UInterp().interpolate(cpw);
        scalar nuc = td.nuInterp().interpolate(cpw);

        scalar rhop = td.spc().rhop();
        scalar magUr = mag(Uc - U_);

        scalar ReFunc = 1.0;
        scalar Re = magUr*d_/nuc;

        if (Re > 0.01)
        {
            ReFunc += 0.15*pow(Re, 0.687);
        }

        scalar Dc = (24.0*nuc/d_)*ReFunc*(3.0/4.0)*(rhoc/(d_*rhop));

        U_ = (U_ + dt*(Dc*Uc + (1.0 - rhoc/rhop)*td.g()))/(1.0 + dt*Dc);

        if (onBoundary() && td.keepParticle)
        {
            if (isA<processorPolyPatch>(pbMesh[patch(face())]))
            {
                td.switchProcessor = true;
            }
        }
    }

    return td.keepParticle;
}


bool Foam::fishParticle::hitPatch
(
    const polyPatch&,
    fishParticle::trackData&,
    const label
)
{
    return false;
}


bool Foam::fishParticle::hitPatch
(
    const polyPatch&,
    int&,
    const label
)
{
    return false;
}


void Foam::fishParticle::hitProcessorPatch
(
    const processorPolyPatch&,
    fishParticle::trackData& td
)
{
    td.switchProcessor = true;
}


void Foam::fishParticle::hitProcessorPatch
(
    const processorPolyPatch&,
    int&
)
{}


void Foam::fishParticle::hitWallPatch
(
    const wallPolyPatch& wpp,
    fishParticle::trackData& td
)
{
    vector nw = wpp.faceAreas()[wpp.whichFace(face())];
    nw /= mag(nw);

    scalar Un = U_ & nw;
    vector Ut = U_ - Un*nw;

    if (Un > 0)
    {
        U_ -= (1.0 + td.spc().e())*Un*nw;
    }

    U_ -= td.spc().mu()*Ut;
}


void Foam::fishParticle::hitWallPatch
(
    const wallPolyPatch&,
    int&
)
{}


void Foam::fishParticle::hitPatch
(
    const polyPatch&,
    fishParticle::trackData& td
)
{
    td.keepParticle = false;
}


void Foam::fishParticle::hitPatch
(
    const polyPatch&,
    int&
)
{}


void Foam::fishParticle::transformProperties (const tensor& T)
{
    Particle<fishParticle>::transformProperties(T);
    U_ = transform(T, U_);
}


void Foam::fishParticle::transformProperties(const vector& separation)
{
    Particle<fishParticle>::transformProperties(separation);
}


// ************************************************************************* //
