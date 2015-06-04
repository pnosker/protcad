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
    amberVDW::setScaleFactor(1.0);//VanDerWaals Energy Scaling
    amberVDW::setRadiusScaleFactor(1.0);//VanDerWaals Radius Energy
    amberVDW::setLinearRepulsionDampeningOff();//VanDerWaals Clash Dampening
    amberElec::setScaleFactor(0.0); //Electrostatics
    solvation::setItsScaleFactor(0.0);

	prot->activateForRepacking(0,0);
	prot->activateForRepacking(0,1);
	prot->activateForRepacking(0,2);
	
/*
	cout << "Res1_psi Res2_phi Res2_psi Res2_chi1 Res3_phi TotEnergy" << endl;

    for (int psi0 = -180; psi0 < 180; psi0+= 10)
    {
        for (int phi1 = -180; phi1 < 180; phi1+= 10)
		{
            for (int psi1 = -180; psi1 < 180; psi1+= 10)
			{
                for (int phi2 = -180; phi2 < 180; phi2+= 10)
				{
                    double phi0d = phi0;
                    double phi1d = phi1;
                    double psi1d = psi1;
                    double phi2d = phi2;
                    double chi1d = chi1;

                    prot->setPsi(0, 0, psi0d);
                    prot->setPhi(0, 1, phi1d);
                    prot->setPsi(0, 1, psi1d);
                    prot->setPhi(0, 2, phi2d);
					double energy = prot->intraEnergy();

					if (energy < -3.2)
					{
                        for (int chi1 = -180; chi1 < 180; chi1+= 5)
						{
                            prot->setChi(0, 1, 0, 0, chi1d);
		

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


/*

					prot->setPsi(0, 0, -90);
					prot->setPhi(0, 1, -40);
					prot->setPsi(0, 1, -60);
					prot->setPhi(0, 2, 155);
					prot->setChi(0, 1, 0, 0, -100);

	outFile = argv[2];
	pdbWriter(prot, outFile);	

*/

return 0;





}
