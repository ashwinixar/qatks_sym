/*
*****************************************************
*	About: The code implements Algorithm 3		 	*
*	Usage: Run in command prompt "A3.exe"			*
*	Author: Ashwini Kumar Malviya					*
*	Email: ashwinixar@gmail.com						*
*****************************************************
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "QuEST.h"

typedef unsigned __int8 block;

//S-Box of the 8-bit custom SPN cipher
int s_table[16] = { 0xE,0x4,0xD,0x1,
					0x2,0xF,0xB,0x8,
					0x3,0xA,0x6,0xC,
					0x5,0x9,0x0,0x7 };

//Performs inverse S-Box
block s_box_inv(block a)
{
	block temp = 0;
	for (int i = 0; i < 2; i++)
	{
		int idx;
		block d = (a & (0x0F << (i * 4))) >> (i * 4);
		for (int i = 0; i < 16; i++)
			if (s_table[i] == d) idx = i;
		temp = temp | (idx) << (i * 4);
	}
	return temp;
}

//Stores encryption of indexes "i" at corresponding entry i.e. oracle_e[i] = E(i)
int oracle_e[256] = { 	0x3C,0xBD,0x22,0x5,0x15,0xEE,0xF0,0x68,0x64,0x72,0x88,0xBE,0xCC,0xD2,0xA0,0x5C,
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

double count(int y)
{
	QuESTEnv env = createQuESTEnv();

    //Value of Pi
    const double pi = 3.14;

    const int p = 8, P = (int)pow(2, p); //Precision qubits for theta of Grover's iteration
    const int n = 8, N = 256;
    Qureg qubits = createQureg((p + n), env);
    initZeroState(qubits);
    //reportQuregParams(qubits);
    //reportQuESTEnv(env);

    //Unitary matrix creation for the key "y"
	ComplexMatrixN e = createComplexMatrixN(n);
    for(int i = 0; i < 256; i++)
    {
    	block c1 = oracle_e[i] ^ y, c2 = oracle_e[i ^ 0x40] ^ y;
		c1 = s_box_inv(c1);
		c2 = s_box_inv(c2);
		if ((c1 ^ c2) == 0x0A) e.real[i][i] = -1;
    	else e.real[i][i] = 1;
    }

    //Apply Hadamard on all the (p + n)-qubits to create superposition
    //Here each element in last n-qubit superposition is considered to represent "x" 
    for(int i = 0; i < (p + n); i++)
        hadamard(qubits, i);
    
    //Target qubits on which grover's search must be performed to count the number of correct result
    int targs[n];
    for(int i = 0; i < n; i++)
    	targs[i] = p + i;

    //Applying controlled grover's operator
    int iterations = 1;
    for(int gi = p - 1; gi >= 0; gi--)
    {
    	int ctrls[] = { gi };
    	for(int j = 0; j < iterations; j++)
    	{
    		//Marking
	    	multiControlledMultiQubitUnitary(qubits, ctrls, 1, targs, n, e);
	    	
	        //Diffusion
	        for(int i = p; i < (n + p); i++)
	            hadamard(qubits, i);
	        for(int i = p; i < (n + p); i++)
	            pauliX(qubits, i);
	        multiControlledPhaseFlip(qubits, targs, n);
	        for(int i = p; i < (n + p); i++)
	            pauliX(qubits, i);
	        for(int i = p; i < (n + p); i++)
	            hadamard(qubits, i);
    	}
    	iterations *= 2;
    }
	
    //Inverse QFT on the first p-qubits
    for(int i = 0; i < p; i++)
    {
    	for(int j = 0; j < i; j++)
    	{
    		qreal angle = -2.0 * pi / (pow(2, i - j + 1));
    		controlledPhaseShift(qubits, j, i, angle);
    	}
    	hadamard(qubits, i);
    }
    
    //Measures the qubits to get the value of first p-qubits
    qreal prob = 0, max_prob = 0;
    int p_val = 0;
    for(int i = 0; i < P; i++)
    {
    	prob = getProbAmp(qubits, i);
    	if(prob > max_prob)
    	{
    		max_prob = prob;
    		p_val = i;
    	}
    }
    if(p_val > (P / 2)) p_val = P - p_val;
    //printf("\nThe measured value of (f_tilde) is %d", p_val);
    //printf("\nThe approximate count of the correct answers is: %f\n", (N * pow(sin(p_val * pi / P), 2)));

    destroyQureg(qubits, env);
    destroyComplexMatrixN(e);
    destroyQuESTEnv(env);

    return (N * pow(sin(p_val * pi / P), 2));
}

int main (int narg, char *varg[])
{
	//Oracle evaluation results of e(x, j) for all pairs "j" corresponding to a "x"
	//It is "+1" when a fixed Pair "j" is a right pair of "x" and "+0", otherwise
	//Each index in array is a 4-bit subkey "x" and corresponding entry contains the number of right pairs
	//Delta_in is 0x40 and Delta_out is 0x0A
	int epairs_count[16] = { 0,0,0,4,0,6,0,0,0,0,0,4,0,4,2,0 };

	int m = 2; //(refer algo in original paper)
	int c = 4; //A random constant value
	int s_n = 4; //Subkey size in bits
	//int n = 8; //Block size of the cipher

	int y = 0x1; // A Random threshold key to begin with
	int R_y = round(count(y));

	printf("\nStart Count is %d for subkey 0x%02X\n", R_y, y);

	QuESTEnv env = createQuESTEnv();

    Qureg qubits = createQureg(s_n, env);
    //reportQuregParams(qubits);
    //reportQuESTEnv(env);

	for(int times = 0; times <= 2 * m * c; times++)
	{
		//Set all the qubits to 0
		initZeroState(qubits);

		//Creation of unitary corresponding to oracle O1 (according to paper) using epairs_count
		ComplexMatrixN e = createComplexMatrixN(s_n);
		for(int i = 0; i < 16; i++)
	    {
	    	if(epairs_count[i] > R_y) e.real[i][i] = -1;
	    	else e.real[i][i] = 1;
	    }

	    //Searching for a key having approx count more than current key y
	    for(int i = 0; i < s_n; i++)
        	hadamard(qubits, i);
	    int targs[] = { 0,1,2,3 };
	    for(int gi = 0; gi < 3; gi++)
	    {
	    	//Marking
	    	multiQubitUnitary(qubits, targs, s_n, e);
	    	
	        //Diffusion
	        for(int i = 0; i < s_n; i++)
	            hadamard(qubits, i);
	        for(int i = 0; i < s_n; i++)
	            pauliX(qubits, i);
	        multiControlledPhaseFlip(qubits, targs, s_n);
	        for(int i = 0; i < s_n; i++)
	            pauliX(qubits, i);
	        for(int i = 0; i < s_n; i++)
	            hadamard(qubits, i);
	    }
	    qreal prob;
	    int y_prime = 0;
	    for(int i = 0; i < s_n; i++)
	    {
	    	int outcome = measureWithStats(qubits, i, &prob);
	    	if(outcome) y_prime ^= (outcome << i);
	    }

	    //Getting count of y_prime and updating y accordingly
	    int temp = round(count(y_prime));
	    printf("Count is %d for subkey 0x%02X\n", temp, y_prime);
	    if(temp > R_y)
	    {
	    	y = y_prime;
	    	R_y = temp;
	    }
	    printf("New threshold subkey is 0x%02X\n\n", y);

	    destroyComplexMatrixN(e);
	}
	printf("\nThe final subkey with highest count is 0x%02X\n", y);

	destroyQureg(qubits, env);
    destroyQuESTEnv(env);

    return 0;
}