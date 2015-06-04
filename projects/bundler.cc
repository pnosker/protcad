//*******************************************************************************************************
//*******************************************************************************************************
//***********************************                          ******************************************
//***********************************     protEvolver 2.0      ******************************************
//***********************************                          ******************************************
//*******************************************************************************************************
//***************   -Folding Selective Protein Evolution in Implicit Solvent-   ***********************
//*******************************************************************************************************

/////// Just specify infile structure, active chains and residues indexes, and it will evolve a sequence favorable for folding

//--Included files and functions-------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <time.h>
#include <sstream>
#include "ensemble.h"
#include "PDBInterface.h"

void randomizeSideChains(protein* _prot, UInt _chainIndex);
vector < UInt > getChainSequence(protein* _prot, UInt _chainIndex);

//--Program setup----------------------------------------------------------------------------------------
int main (int argc, char* argv[])
{
	//--Running parameters
	if (argc !=2)
	{
		cout << "protEvolver <inFile.pdb>" << endl;
		exit(1);
	}
	string infile = argv[1];
	enum aminoAcid {A,R,N,D,Dh,C,Cx,Q,E,Eh,Hd,He,Hn,Hp,I,L,K,M,F,P,O,S,T,W,Y,V,G,dA,dR,dN,dD,dDh,dC,dCx,dQ,dE,dEh,dHd,dHe,dHn,dHp,dI,dL,dK,dM,dF,dP,dO,dS,dT,dAT,dW,dY,dV};
	string aminoAcidString[] = {"A","R","N","D","Dh","C","Cx","Q","E","Eh","Hd", "He","Hn","Hp","I","L","K","M","F","P","O","S","T","W","Y", "V","G","dA","dR","dN","dD","dDh","dC","dCx","dQ","dE","dEh","dHd","dHe","dHn","dHp","dI","dL","dK","dM","dF","dP","dO","dS","dT","dW","dY","dV"};
	PDBInterface* thePDB = new PDBInterface(infile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* pMol = theEnsemble->getMoleculePointer(0);
	protein* bundle = static_cast<protein*>(pMol);
	bundle->silenceMessages();
	residue::setCutoffDistance(9.0);
	pmf::setScaleFactor(0.0);
	rotamer::setScaleFactor(0.0);
	microEnvironment::setScaleFactor(0.0);
	amberVDW::setScaleFactor(0.8);
	amberVDW::setRadiusScaleFactor(1.0);
	amberVDW::setLinearRepulsionDampeningOn();
	amberElec::setScaleFactor(1.0);
	solvation::setItsScaleFactor(0.0);
	srand (time(NULL));

	//--inputs for mutation
	UInt activeChains[] = {0,1,2,3};

    UInt allowedLResidues[] = {F,I,L,V,W,Y};
	UInt activeResidues[] = {2,6,9,12,13,19,20,23,26,27,30};
	UInt allowedDResidues[] = {G};

	double phi, bestEnergy, foldingEnergy, finalEnergy, posE, totalposE, aveposE, originalEnergy, sessEnergy, tempEnergy, startAvePosE, sysEnergy, bindingEnergy, sessBindEnergy, tempBindEnergy, originalBindEnergy, bestBindingEnergy;
	UInt nobetter = 0, activeResiduesSize = sizeof(activeResidues)/sizeof(activeResidues[0]), activeChainsSize = sizeof(activeChains)/sizeof(activeChains[0]);
	UInt dResidues = sizeof(allowedDResidues)/sizeof(allowedDResidues[0]), lResidues = sizeof(allowedLResidues)/sizeof(allowedLResidues[0]);
	UInt name, mutant, numResidues, plateau = (lResidues*activeResiduesSize*activeChainsSize*2), fib = 0, oldBuffer, newBuffer, bufferBump;
	vector < UInt > mutantPosition, chainSequence, sequencePosition;
	vector < vector < UInt > > proteinSequence, finalSequence;
	vector < double > Energy(2);
	UInt randres, randchain, rescount = 0;
	stringstream convert;
	string startstr, outFile;
	name = rand() % 1000000;
	convert << name, startstr = convert.str();
	string tempModel = startstr + "_temp.pdb";
	delete thePDB;

	//--Run multiple independent evolution cycles-----------------------------------------------------
	for (UInt a = 1; a < 200; a++)
	{
		PDBInterface* thePDB = new PDBInterface(infile);
		ensemble* theEnsemble = thePDB->getEnsemblePointer();
		molecule* pMol = theEnsemble->getMoleculePointer(0);
		protein* bundle = static_cast<protein*>(pMol);


		//--load in initial pdb and mutate in random starting sequence on active chains and residues
		proteinSequence.clear(), chainSequence.clear(), Energy.clear(), fib = 0, nobetter = 0, sysEnergy = 0, bindingEnergy = 0;
		for (UInt i = 0; i < activeChainsSize; i++)
		{
			for (UInt j = 0; j < activeResiduesSize; j++)
			{
			bundle->activateForRepacking(activeChains[i], activeResidues[j]);
            }
			chainSequence = getChainSequence(bundle, activeChains[i]);
			proteinSequence.push_back(chainSequence);
		}


		//--Determine next mutation position
		rescount = 0, totalposE = 0, mutantPosition.clear();

			randres = activeResidues[rand() % activeResiduesSize];

//cout << "Residue Number = " << randres << endl;

			if (activeChainsSize > 1)
			{
				randchain = activeChains[rand() % activeChainsSize];
			}
			else
			{
				randchain = activeChains[0];
			}



//cout << "Chain Number = " << randchain << endl;

		mutantPosition.push_back(randchain);
		mutantPosition.push_back(randres);
		pdbWriter(bundle, tempModel);

//cout << "Chain " << randchain << " Average Position Energy = " << aveposE << endl;
		
		startAvePosE = aveposE;

cout << "Temp PDB written" << endl;

//		//--set Energy startpoint
		Energy.clear();
		Energy = bundle->chainBindingEnergy();

		sysEnergy = 0;
		bindingEnergy = 0;
		sysEnergy = bundle->intraSoluteEnergy(true);
		bindingEnergy = Energy[1];
		bestEnergy = sysEnergy;
		bestBindingEnergy = bindingEnergy;
		originalEnergy = bestEnergy;
		originalBindEnergy = bestBindingEnergy;
		fib = 0.000000001;

cout << "Best Energy Set" << endl;
cout << " " << endl;


		delete thePDB;
		
		//--Run single evolution loop till hitting plateau
		do
		{  
			PDBInterface* thePDB = new PDBInterface(infile);
			ensemble* theEnsemble = thePDB->getEnsemblePointer();
			molecule* pMol = theEnsemble->getMoleculePointer(0);
			protein* bundle = static_cast<protein*>(pMol);

			//--Mutate current sequence, new mutant and optimize system


			oldBuffer = fib;
			fib = 0.00000001, nobetter++;
			newBuffer = fib;
			for (UInt i = 0; i < activeChainsSize; i++)
			{
				numResidues = bundle->getNumResidues(activeChains[i]);
				for (UInt j = 0; j < numResidues; j++)
				{
					bundle->activateForRepacking(activeChains[i],j);
					if (activeChains[i] == mutantPosition[0] && j == mutantPosition[1])
					{
						//--new mutant

						sequencePosition.push_back(i);
						sequencePosition.push_back(j);
						phi = bundle->getPhi(mutantPosition[0], mutantPosition[1]);
						if (phi > 0 && phi < 180)
						{
							mutant = allowedDResidues[(rand() % dResidues)];
							bundle->mutateWBC(mutantPosition[0],mutantPosition[1], mutant);
						}
						if (phi < 0 && phi > -180)
						{
							mutant = allowedLResidues[(rand() % lResidues)];
							bundle->mutateWBC(mutantPosition[0],mutantPosition[1], mutant);
						}

 //                                  bundle->optimizeRotamers();
					}
					else
					{


						bundle->mutateWBC(activeChains[i],j, proteinSequence[i][j]);


					}

                    bundle->activateAllForRepacking(activeChains[i]);

                    bundle->optimizeRotamersPN();

				}

			}



 //          bundle->rotOptSolvent(10);

			protein* tempBundle = new protein(*bundle);
			
			//--Determine next mutation position
			rescount = 0, totalposE = 0, mutantPosition.clear();

				randres = activeResidues[rand() % activeResiduesSize];


				if (activeChainsSize > 1)
				{
					randchain = activeChains[rand() % activeChainsSize];
				}
				else
				{
					randchain = activeChains[0];
				}



			mutantPosition.push_back(randchain);
			mutantPosition.push_back(randres);


                bundle->optimizeRotamersPN();
                bundle->rotOptSolvent(100);

			//--Energy test and determination of next mutant position
			Energy.clear();
			Energy = bundle->chainBindingEnergy();

			bindingEnergy = 0;
			bindingEnergy = Energy[1];
			
			sysEnergy = 0;
			sysEnergy = bundle->intraSoluteEnergy(true);

			sessEnergy = sysEnergy;
			sessBindEnergy = bindingEnergy;


cout << " " << endl;

cout << "Starting Bind Energy = " << originalBindEnergy << " | Best Bind Energy = " << bestBindingEnergy << " | Session Binding Energy = " << sessBindEnergy << " | Buffer = " << newBuffer << " | " << "Buffer Bump = " << newBuffer - oldBuffer << flush;

			if (sessBindEnergy < (bestBindingEnergy + newBuffer))
			{

cout << " New Binding Energy = " << sessBindEnergy << endl;

				if (sessBindEnergy < bestBindingEnergy)

				{

					bestBindingEnergy = sessBindEnergy;
					pdbWriter(tempBundle, tempModel);
				}
				fib = 0, nobetter--, proteinSequence[sequencePosition[0]][sequencePosition[1]] = mutant, foldingEnergy = Energy[1];
			}
			sequencePosition.clear();
			delete thePDB;
			delete tempBundle;
		}while (nobetter < plateau);

		//--Print final energy and write a pdb file----------------------------------------------------
		PDBInterface* theModelPDB = new PDBInterface(tempModel);
		ensemble* theModelEnsemble = theModelPDB->getEnsemblePointer();
		molecule* modelMol = theModelEnsemble->getMoleculePointer(0);
		protein* model = static_cast<protein*>(modelMol);
		name = rand() % 1000000;
		stringstream convert;
		string countstr;
		convert << name, countstr = convert.str();
		outFile = countstr + ".evo.pdb";
		pdbWriter(model, outFile);
		finalSequence.clear(), chainSequence.clear(), Energy.clear(), sysEnergy = 0;
		for (UInt i = 0; i < activeChainsSize; i++)
		{
			chainSequence = getChainSequence(model, activeChains[i]);
			finalSequence.push_back(chainSequence);
		}
		Energy = model->chainBindingEnergy();
		bindingEnergy = Energy[1];		
		finalEnergy = bindingEnergy;
		cout << name << " | Total Energy: " << finalEnergy << " | Binding Energy: " << bestBindingEnergy << endl;
		for (UInt i = 0; i < activeChainsSize; i++)
		{
			for (UInt j = 0; j < finalSequence[i].size(); j++)
			{
				cout << aminoAcidString[finalSequence[i][j]] << " ";
			}
		}
		cout << endl;
		delete theModelPDB;
	}
	return 0;
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////
//////// functions //////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

void randomizeSideChains(protein* _prot, UInt _chainIndex)
{	
	UInt allowedRotsSize, randrot, restype, numResidues;
	UIntVec allowedRots;
	numResidues = _prot->getNumResidues(_chainIndex);
	for (UInt j = 0; j < numResidues; j++)
	{
		restype = _prot->getTypeFromResNum(_chainIndex, j);
		allowedRots = _prot->getAllowedRotamers(_chainIndex, j, restype, 0);
		allowedRotsSize = allowedRots.size();
		if (allowedRotsSize > 2)
		{
			randrot = rand() % allowedRotsSize;
			_prot->setRotamerWBC(_chainIndex, j, 0, allowedRots[randrot]);
		}
	}
	return;
}

vector < UInt > getChainSequence(protein* _prot, UInt _chainIndex)
{	
	UInt restype, numResidues;
	vector < UInt > sequence;
	numResidues = _prot->getNumResidues(_chainIndex);
	for (UInt j = 0; j < numResidues; j++)
	{
		restype = _prot->getTypeFromResNum(_chainIndex, j);
		sequence.push_back(restype);
	}
	return sequence;
}
