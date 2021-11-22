{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

    surfaceScalarField phir
    (
        IOobject
        (
            "phir",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mixture.cAlpha()*mag(phiu/mesh.magSf())*mixture.nHatf()
    );

    for (int gCorr=0; gCorr<nAlphaCorr; gCorr++)
    {
        // Create the limiter to be used for all phase-fractions
        scalarField allLambda(mesh.nFaces(), 1.0);

        // Split the limiter into a surfaceScalarField
        slicedSurfaceScalarField lambda
        (
            IOobject
            (
                "lambda",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimless,
            allLambda,
            false   // Use slices for the couples
        );


        // Create the complete convection flux for alpha1
        surfaceScalarField alphaPhi1
        (
            fvc::flux
            (
                phiu,
                alpha1,
                alphaScheme
            )
          + fvc::flux
            (
                -fvc::flux(-phir, alpha2, alpharScheme),
                alpha1,
                alpharScheme
            )
          + fvc::flux
            (
                -fvc::flux(-phir, alpha3, alpharScheme),
                alpha1,
                alpharScheme
            )
        );

        // Create the complete flux for alpha2
        surfaceScalarField alphaPhi2
        (
            fvc::flux
            (
                phiu,
                alpha2,
                alphaScheme
            )
          + fvc::flux
            (
                -fvc::flux(phir, alpha1, alpharScheme),
                alpha2,
                alpharScheme
            )
          + fvc::flux
            (
                -fvc::flux(phi23, alpha3, alpharScheme),
                alpha2,
                alpharScheme
            )
        );

        // Create the bounded (upwind) flux for alpha2
        surfaceScalarField alphaPhi2BD
        (
            upwind<scalar>(mesh, phiu).flux(alpha2)
        );

        // Calculate the flux correction for alpha2
        alphaPhi2 -= alphaPhi2BD;

        // Further limit the limiter for alpha2
        if (LTS)
        {
            const volScalarField& rDeltaT =
                fv::localEulerDdt::localRDeltaT(mesh);

            MULES::limiter
            (
                allLambda,
                rDeltaT,
                geometricOneField(),
                alpha2,
                alphaPhi2BD,
                alphaPhi2,
                zeroField(),
                zeroField(),
                oneField(),
                zeroField()
            );
        }
        else
        {
            MULES::limiter
            (
                allLambda,
                1.0/runTime.deltaT().value(),
                geometricOneField(),
                alpha2,
                alphaPhi2BD,
                alphaPhi2,
                zeroField(),
                zeroField(),
                oneField(),
                zeroField()
            );
        }

        // Construct the limited fluxes
        alphaPhi2 = alphaPhi2BD + lambda*alphaPhi2;

         // Create the complete flux for alpha3
        surfaceScalarField alphaPhi3
        (
            fvc::flux
            (
                phiu,
                alpha3,
                alphaScheme
            )
          + fvc::flux
            (
                -fvc::flux(phir, alpha1, alpharScheme),
                alpha3,
                alpharScheme
            )
          + fvc::flux
            (
                fvc::flux(phi23, alpha2, alpharScheme),
                alpha3,
                alpharScheme
            )
        );

        // Create the bounded (upwind) flux for alpha3
        surfaceScalarField alphaPhi3BD
        (
            upwind<scalar>(mesh, phiu).flux(alpha3)
        );

        // Calculate the flux correction for alpha3
        alphaPhi3 -= alphaPhi3BD;

        // Further limit the limiter for alpha3
        if (LTS)
        {
            const volScalarField& rDeltaT =
                fv::localEulerDdt::localRDeltaT(mesh);

            MULES::limiter
            (
                allLambda,
                rDeltaT,
                geometricOneField(),
                alpha3,
                alphaPhi3BD,
                alphaPhi3,
                zeroField(),
                zeroField(),
                oneField(),
                zeroField()
            );
        }
        else
        {
            MULES::limiter
            (
                allLambda,
                1.0/runTime.deltaT().value(),
                geometricOneField(),
                alpha3,
                alphaPhi3BD,
                alphaPhi3,
                zeroField(),
                zeroField(),
                oneField(),
                zeroField()
            );
        }

        // Construct the limited fluxes
        alphaPhi3 = alphaPhi3BD + lambda*alphaPhi3;

        

        // Create the diffusion coefficients for alpha2<->alpha3
        volScalarField Dc23(D23*max(alpha3, scalar(0))*pos0(alpha2));
        volScalarField Dc32(D23*max(alpha2, scalar(0))*pos0(alpha3));


        // Solve for alpha2
        fvScalarMatrix alpha2Eqn
        (
            fvm::ddt(alpha2)
          + fvc::div(alphaPhi2)
          - fvm::laplacian(Dc23 + Dc32, alpha2)
        );
        alpha2Eqn.solve();
        alpha2 = max(alpha2, scalar(0.0));

        // Solve for alpha3
        fvScalarMatrix alpha3Eqn
        (
            fvm::ddt(alpha3)
          + fvc::div(alphaPhi3)
          - fvm::laplacian(Dc23 + Dc32, alpha3)
        );
        alpha3Eqn.solve();
        alpha3 = min(alpha3, 1.0 - alpha2);
        alpha3 = max(alpha3, scalar(0.0));


        // Construct the complete mass flux
        //rhoPhi =
        //      alphaPhi1*(rho1 - rho3)
        //    + (alphaPhi2 + alpha2Eqn.flux())*(rho2 - rho3)
        //    + phi*rho3;
          rhoPhi =
		      alphaPhi3 * rho3
		    + alphaPhi2 * rho2
		    + alphaPhi1 * rho1; 

        alpha1 = 1.0 - alpha2 - alpha3;
        //rho = alpha1*rho1 + alpha2*rho2 + alpha3*rho3;
    }

    Info<< "Air phase volume fraction = "
        << alpha1.weightedAverage(mesh.V()).value()
        << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
        << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
        << endl;

    Info<< "Slurry phase volume fraction = "
        << alpha2.weightedAverage(mesh.V()).value()
        << "  Min(" << alpha2.name() << ") = " << min(alpha2).value()
        << "  Max(" << alpha2.name() << ") = " << max(alpha2).value()
        << endl;
}
