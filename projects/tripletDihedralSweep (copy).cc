#include <iostream>
#include <string>
#include <time.h>
#include <sstream>
#include "ensemble.h"
#include "PDBInterface.h"
#include "typedef.h"
#include "ran.h"

int main (int argc, char* argv[])
{
	enum aminoAcid {A,R,N,D,Dh,C,Cx,Q,E,Eh,Hd,He,Hn,Hp,I,L,K,M,F,P,O,S,T,W,Y,V,G,dA,dR,dN,dD,dDh,dC,dCx,dQ,dE,dEh,dHd,dHe,dHn,dHp,dI,dL,dK,dM,dF,dP,dO,dS,dT,dW,dY,dV};

	string inputFileName = argv[1];

    PDBInterface* thePDB = new PDBInterface(inputFileName);
    ensemble* theEnsemble = thePDB->getEnsemblePointer();
    molecule* theMol = theEnsemble->getMoleculePointer(0);
    protein* prot = static_cast<protein*>(theMol);

    residue::setCutoffDistance(4.0);
    pmf::setScaleFactor(0.0);
    rotamer::setScaleFactor(0.0);
    microEnvironment::setScaleFactor(0.0);
    amberVDW::setScaleFactor(1.0);
    amberVDW::setRadiusScaleFactor(1.0);
    amberVDW::setLinearRepulsionDampeningOff();
    amberElec::setScaleFactor(0.0);
    solvation::setItsScaleFactor(0.0);

	prot->activateForRepacking(0,0);
	for (double psi = -180.0; psi < 180.0; psi++)
	{
		prot->setPsi(0, 0, psi);
		
		double energy = prot->intraEnergy();

		cout << "Res 1 " << "NoPhi" << " " << psi << " " << energy << endl;
		
	}
	prot->activateForRepacking(0,1);
	for (double phi = -180.0; phi < 180.0; phi ++)
	{
		for (double psi = -180.0; psi < 180.0; psi ++)
		{
			prot->setPhi(0, 1, phi);
			prot->setPsi(0, 1, psi);
		
			double energy = prot->intraEnergy();

			cout << "Res 2 " << phi << " " << psi << " " << energy << endl;
		}
	}
	prot->activateForRepacking(0,2);
	for (double phi = -180.0; phi < 180.0; phi++)
	{
		prot->setPhi(0, 2, phi);
		
		double energy = prot->intraEnergy();

		cout << "Res 3 " << phi << " " << "NoPsi" << " " << energy << endl;
		
	}
	return 0;
}
