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
	enum aminoAcid {A,R,N,D,Dh,C,Cx,Q,E,Eh,Hd,He,Hn,Hp,I,L,K,M,F,P,O,S,T,W,Y,V,G,dA,dR,dN,dD,dDh,dC,dCx,dQ,dE,dEh,dHd,dHe,dHn,dHp,dI,dL,dK,dM,dF,dP,dO,dS,dT,dAT,dW,dY,dV};

	string inputFileName = argv[1];

    PDBInterface* thePDB = new PDBInterface(inputFileName);
    ensemble* theEnsemble = thePDB->getEnsemblePointer();
    molecule* theMol = theEnsemble->getMoleculePointer(0);
    protein* prot = static_cast<protein*>(theMol);

    	string outFile;

    prot->silenceMessages();
    residue::setCutoffDistance(9.0);
    pmf::setScaleFactor(0.0);
    rotamer::setScaleFactor(0.0);
    microEnvironment::setScaleFactor(0.0);
    amberVDW::setScaleFactor(1.0); //VanDerWaals Energy Scaling
    amberVDW::setRadiusScaleFactor(1.0); //VanDerWaals Radius Energy
    amberVDW::setLinearRepulsionDampeningOff(); //VanDerWaals Clash Dampening
    amberElec::setScaleFactor(0.0); //Electrostatics
    solvation::setItsScaleFactor(0.0);
    double energy;
    double psi0, phi1, psi1, phi2, chi1, bestpsi0, bestphi1, bestpsi1, bestphi2, bestchi1, currentPsi0, currentPhi1, currentPsi1, currentPhi2, currentChi1;
	prot->activateForRepacking(0,0);
	prot->activateForRepacking(0,1);
	prot->activateForRepacking(0,2);


	
/*
	cout << "Res1_psi Res2_phi Res2_psi Res2_chi1 Res3_phi TotEnergy" << endl;
	for (double psi0 = -180.0; psi0 < 180.0; psi0+= 10.0)
	{
		for (double phi1 = -180.0; phi1 < 180.0; phi1+= 10.0)
		{
			for (double psi1 = -180.0; psi1 < 180.0; psi1+= 10.0)
			{
				for (double phi2 = -180.0; phi2 < 180.0; phi2+= 10.0)
				{
					prot->setPsi(0, 0, psi0);
					prot->setPhi(0, 1, phi1);
					prot->setPsi(0, 1, psi1);
					prot->setPhi(0, 2, phi2);
					double energy = prot->intraEnergy();

					if (energy < -3.2)
					{
						for (double chi1 = -180; chi1 < 180.0; chi1+= 5.0)
						{
							prot->setChi(0, 1, 0, 0, chi1);
		

							double energy = prot->intraEnergy();
							if (energy < -3.9)
							{
//							cout << "Res1_psi " << psi0 << " Res2_phi " << phi1 << " Res2_psi " << psi1 << " Res2_chi1 " << chi1 <<" " << " Res3_phi " << phi2 << " TotEnergy " << energy << endl;
							cout << psi0 << " " << phi1 << " " << psi1 << " " << chi1 <<" " << phi2 << " " << energy << endl;
							}
						}
					}
				}
			}
		}
	}

*/


///*
				double bestEnergy = 1E10;

			for (psi0 = -180.0; psi0 < 180.0; psi0+= 10.0)
			{
				for (phi2 = -180.0; phi2 < 180.0; phi2+= 10.0)
				{
					for (chi1 = -180.0; chi1 < 180.0; chi1+= 10.0)
{
					prot->setPsi(0, 0, psi0);
					prot->setPhi(0, 1, 60);
					prot->setPsi(0, 1, -130);
					prot->setPhi(0, 2, phi2);
					prot->setChi(0, 1, 0, 0, chi1);

					energy = prot->intraEnergy();

						if (energy < bestEnergy)
						{
							bestpsi0 = psi0;
							bestphi2 = phi2;
							bestchi1 = chi1;
							bestphi1 = phi1;
							bestpsi1 = psi1;
							bestEnergy = energy;
						}
}
}
}			

					prot->setPsi(0, 0, bestpsi0);
					prot->setPhi(0, 1, bestphi1);
					prot->setPsi(0, 1, bestpsi1);
					prot->setPhi(0, 2, bestphi2);
					prot->setChi(0, 1, 0, 0, bestchi1);

	outFile = argv[2];
	pdbWriter(prot, outFile);	



			cout << bestpsi0 << " " << bestphi1 << " " << bestpsi1 << " " << bestchi1 <<" " << bestphi2 << " " << bestEnergy << endl;			

//*/

return 0;





}
