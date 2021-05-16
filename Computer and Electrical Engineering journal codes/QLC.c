/*
*************************************************************
*	About: The code implements Quantum Linear Cryptanalysis	*
*			Last Rounds-Attack is implemented				*
*	Usage: Run in command prompt "QLC.exe"					*
*	Author: Ashwini Kumar Malviya							*
*	Email: ashwinixar@gmail.com								*
*************************************************************
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "QuEST.h"

int k[5] = { 0xE3, 0x2B, 0xF5, 0x9A, 0x75 };

int p_table[8] = { 0,4,6,5,1,3,2,7 };

int s_table[16] = { 0xE,0x4,0xD,0x1,
					0x2,0xF,0xB,0x8,
					0x3,0xA,0x6,0xC,
					0x5,0x9,0x0,0x7 };

int s_box_inv(int a)
{
	int temp = 0;
	for (int i = 0; i < 2; i++)
	{
		int idx;
		int d = (a & (0x0F << (i * 4))) >> (i * 4);
		for (int i = 0; i < 16; i++)
			if (s_table[i] == d) idx = i;
		temp = temp | (idx) << (i * 4);
	}
	return temp;
}

int p_box(int a)
{
	int temp = 0;
	for (int i = 0; i < 8; i++)
	{
		if (a & (0x1 << i))
			temp = temp | (0x1 << p_table[i]);
	}
	return temp;
}

int bit_product(int a, int mask)
{
	a = a & mask;
	mask = 0;
	while (a > 0)
	{
		mask = mask ^ (a & 0x1);
		a >>= 1;
	}
	return mask;
}

//Stores encryption oracle evaluations of indexes "i" at corresponding entry i.e. oracle_E[i] = E(i)
//E : (F_2)^m -> (F_2)^n, E is a boolean encryption function which takes "m" bit string as input and produces "n" bit output
//Here m = n = 8
int oracle_E[256] = { 	0x3C,0xBD,0x22,0x5,0x15,0xEE,0xF0,0x68,0x64,0x72,0x88,0xBE,0xCC,0xD2,0xA0,0x5C,
						0x31,0xB9,0x3E,0x5E,0x9D,0xE3,0xFB,0x11,0xE1,0x9E,0x7B,0xB3,0xC1,0xFD,0xAB,0x51,
						0xF,0x82,0x37,0x57,0x9A,0xC5,0x26,0xAF,0x8F,0x97,0xC6,0xA2,0x46,0xFA,0x4F,0xD6,
						0x49,0x9,0x7C,0xB8,0x2C,0x81,0xA1,0xE8,0xDC,0x0,0xB4,0xCB,0x23,0x10,0xD3,0x43,
						0x4,0xEC,0x5A,0xF7,0xA7,0x70,0x28,0x66,0x6F,0xCA,0x86,0x1C,0x48,0x3A,0x44,0xD8,
						0x17,0xDA,0x4E,0xAC,0xDE,0x69,0x63,0xC0,0xD,0x2E,0x8C,0x83,0x7,0x4D,0x2A,0x7A,
						0x3D,0xB1,0x2F,0x6,0x16,0xEB,0xFE,0x65,0x62,0x7F,0x85,0xBB,0xCD,0xDF,0xAE,0x5D,
						0x94,0xF4,0x47,0xA9,0xD7,0xE2,0x12,0xC3,0xA,0x27,0x89,0x75,0x58,0x4A,0x38,0x98,
						0x1B,0xD1,0xC2,0xB6,0x52,0x6D,0x6E,0xE6,0x32,0xF5,0xBF,0x8E,0xB,0xA5,0x21,0x71,
						0x1E,0xDD,0xC4,0xB5,0x54,0x6C,0x60,0xE5,0x34,0xF8,0xB2,0x80,0xE,0xA8,0x2D,0x7D,
						0x2,0x84,0x33,0x53,0x99,0xC8,0x25,0x1A,0xEA,0x93,0x77,0xA4,0x45,0xF9,0x42,0xD5,
						0x13,0xD9,0x4B,0xAD,0xDB,0x61,0x6B,0xCE,0x1,0x2B,0x8D,0x8B,0x3,0x41,0x29,0x79,
						0x92,0xF2,0x96,0xAA,0x36,0xEF,0x1F,0xC7,0xFF,0x56,0x8A,0x76,0x55,0x9F,0x35,0x95,
						0xF1,0xE9,0x24,0x8,0x18,0x73,0x5B,0xB0,0xBC,0x74,0xE0,0x19,0x9B,0xD4,0x91,0x3B,
						0x4C,0xC,0xCF,0x67,0x5F,0xE4,0x14,0x87,0x3F,0xF6,0x6A,0x78,0x20,0xA6,0xD0,0x40,
						0x39,0xBA,0x30,0x50,0x9C,0xE7,0xF3,0x1D,0xED,0x90,0x7E,0xB7,0xC9,0xFC,0xA3,0x59 };


//Using Grover's search to find last round best subkey
int find_best_subkey()
{
	//Oracle evaluation results for different subkeys
	//For each subkey (represented by index of the array), corresponding array element represents a number of "x" for which ..
	// .. "bit_product(x, 0x10) XOR bit_product(Sbox_Inverse(E(x) ^ subkey) XOR 0x0E) = 0" holds
	//Each index in array represents a 4-bit "subkey"
	//The oracle elements is found by evaluating encryption function classically (can't be done as quantum simulation)
	int oracle_subkeys[16];
	for (int i = 0; i < 16; i++)
		oracle_subkeys[i] = 0;
	for (int i = 0; i < 16; i++) //Loop for each subkey
	{
		for (int j = 0; j < 256; j += 2) //Loop for each data (ciphertext) last round backward operation
		{
			//Performing last round operation in reverse
			int c = oracle_E[j] ^ i; //Applying subkey "i" on ciphertext
			c = s_box_inv(c);
			if ((bit_product(j, 0x10) ^ bit_product(c, 0x0E)) == 0)
				oracle_subkeys[i]++;
		}
	}
	//Find a subkey with the biggest bias
	int _to_search = 0;
	int _max_bias = 0;
	for(int i = 0; i < 16; i++)
	{
		if(abs(128 - oracle_subkeys[i]) > _max_bias)
		{
			_to_search = i;
			_max_bias = abs(128 - oracle_subkeys[i]);
		}
	}

	QuESTEnv env = createQuESTEnv();

    Qureg qubits = createQureg(4, env);
    initZeroState(qubits);

    //Unitary matrix creation corresponding to the "oracle_subkeys"
	ComplexMatrixN e = createComplexMatrixN(4);
    for(int i = 0; i < 16; i++)
    {
    	if(i == _to_search) e.real[i][i] = -1; //subkey that holds for most of the "x"
    	else e.real[i][i] = 1;
    }

    //Apply Hadamard on all the n-qubits to create superposition
    //Here each element in superposition is considered to represent "x" 
    for(int i = 0; i < 4; i++)
        hadamard(qubits, i);
    
    int targs[] = { 0,1,2,3 };
    for(int gi = 0; gi < 3; gi++)
    {
    	//Marking
    	multiQubitUnitary(qubits, targs, 4, e);
    	
        //Diffusion
        for(int i = 0; i < 4; i++)
            hadamard(qubits, i);
        for(int i = 0; i < 4; i++)
            pauliX(qubits, i);
        multiControlledPhaseFlip(qubits, targs, 4);
        for(int i = 0; i < 4; i++)
            pauliX(qubits, i);
        for(int i = 0; i < 4; i++)
            hadamard(qubits, i);
    }

    qreal prob;
    int subkey = 0;
    for(int i = 0; i < 4; i++)
    {
    	int outcome = measureWithStats(qubits, i, &prob);
    	if(outcome) subkey ^= (outcome << i);
    }

    destroyQureg(qubits, env);
    destroyComplexMatrixN(e);
    destroyQuESTEnv(env);

    return subkey;
}

//Using Grover's search to find last round complete key
int find_complete_key(int P, int C, int subkey)
{
	//Create oracle corresponding to remaining first 4 key bits
	//Each index in array represents the first 4-bit "subkey"
	//Oracle stores the entry "1" corresponding to the correct key and 0, otherwise
	int oracle_rem_bits[16];
	for (int i = 0; i < 16; i++) //Loop for each subkey
	{
			int temp_C = C;
			temp_C = temp_C ^ ((i << 4) ^ (subkey)); //Applying complete key ("i" concatenated with subkey) on ciphertext
			for (int r = 3; r >= 0; r--)
			{
				if (r != 3) temp_C = p_box(temp_C);
				temp_C = s_box_inv(temp_C);
				temp_C = temp_C ^ k[r];
			}
			if(P == temp_C) oracle_rem_bits[i] = 1;
			else oracle_rem_bits[i] = 0;
	}

	QuESTEnv env = createQuESTEnv();

    Qureg qubits = createQureg(4, env);
    initZeroState(qubits);

    //Unitary matrix creation corresponding to the "oracle_rem_bits"
	ComplexMatrixN e = createComplexMatrixN(4);
    for(int i = 0; i < 16; i++)
    {
    	if(oracle_rem_bits[i] == 1) e.real[i][i] = -1;
    	else e.real[i][i] = 1;
    }

    //Apply Hadamard on all the n-qubits to create superposition
    //Here each element in superposition is considered to represent "x" 
    for(int i = 0; i < 4; i++)
        hadamard(qubits, i);
    
    int targs[] = { 0,1,2,3 };
    for(int gi = 0; gi < 3; gi++)
    {
    	//Marking
    	multiQubitUnitary(qubits, targs, 4, e);
    	
        //Diffusion
        for(int i = 0; i < 4; i++)
            hadamard(qubits, i);
        for(int i = 0; i < 4; i++)
            pauliX(qubits, i);
        multiControlledPhaseFlip(qubits, targs, 4);
        for(int i = 0; i < 4; i++)
            pauliX(qubits, i);
        for(int i = 0; i < 4; i++)
            hadamard(qubits, i);
    }

    qreal prob;
    int f_subkey = 0;
    for(int i = 0; i < 4; i++)
    {
    	int outcome = measureWithStats(qubits, i, &prob);
    	if(outcome) f_subkey ^= (outcome << i);
    }
    /*
    qreal max_prob = 0;
    for(int i = 0; i < 16; i++)
    {
    	prob = getProbAmp(qubits, i);
    	if(prob > max_prob)
    	{
    		max_prob = prob;
    		f_subkey = i;
    	}
    }*/

    destroyQureg(qubits, env);
    destroyComplexMatrixN(e);
    destroyQuESTEnv(env);

    return ((f_subkey << 4) ^ (subkey));
}

int main (int narg, char *varg[])
{
    int subkey = find_best_subkey();
    printf("\nThe best subkey is %02X", subkey);

    int complete_key = find_complete_key(56, oracle_E[56], subkey);
    printf("\nThe complete key is %02X\n", complete_key);

    return 0;
}