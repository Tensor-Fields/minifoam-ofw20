/// Where are these files?
#include "argList.H"
#include "Time.H"
#include "IOobject.H"
#include "fvMesh.H"
#include "fvMatrix.H"
#include "scalar.H"
#include "fvm.H"

using namespace Foam;

int main(int argc, char* argv[]) /// iso C++
{
    /// Parse args
    argList args(argc, argv);

    /// Create time
    Time runTime("banana", args);

    /// Discretize the domain
    IOobject meshFile
    (
        word("region0"),
        runTime.constant(),
        runTime,
        IOobject::MUST_READ
    );

    fvMesh mesh(meshFile);

    /// Create properties
    IOdictionary waveProperties
    (
        IOobject
        (
            word("waveProperties"),
            runTime.constant(),
            mesh,
            IOobject::MUST_READ
        )
    );

    dimensionedScalar c2
    (
        "c2",
        pow(dimLength/dimTime,2),
        //1
        waveProperties.lookup<scalar>("c2")
    );

    /// Create field

    IOobject fieldFile
    (
        word("u"),
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );

    GeometricField<scalar, fvPatchField, volMesh> u
    (fieldFile, mesh);

    while(runTime.loop())
    {
        /// Discretize the PDE
        fvMatrix<scalar> uEq
        (
            fvm::d2dt2(u) == fvm::laplacian(c2, u)
        );

        /// Solve the system
        uEq.solve();

        /// Write
        if (runTime.write())
        {
            //u.write();
            mesh.write();
        }
    }

    Info<< "End.\n";
}
