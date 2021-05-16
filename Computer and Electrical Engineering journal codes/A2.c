/*
*****************************************************
*	About: The code implements Algorithm 2		 	*
*	Usage: Run in command prompt "A2.exe"			*
*	Author: Ashwini Kumar Malviya					*
*	Email: ashwinixar@gmail.com						*
*****************************************************
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "QuEST.h"

typedef unsigned __int8 block;

//STARTS - Encryption Oracle Implementation

block k[5] = { 0xE3, 0x2B, 0xF5, 0x9A, 0x75 };

int s_table[16] = { 0xE,0x4,0xD,0x1,
					0x2,0xF,0xB,0x8,
					0x3,0xA,0x6,0xC,
					0x5,0x9,0x0,0x7 };

int p_table[8] = { 0,4,6,5,1,3,2,7 };

block s_box(block a)
{
	block temp = 0;
	for (int i = 0; i < 2; i++)
	{
		int idx = (a & (0x0F << (i * 4))) >> (i * 4);
		temp = temp | (s_table[idx]) << (i * 4);
	}
	return temp;
}

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

block p_box(block a)
{
	block temp = 0;
	for (int i = 0; i < 8; i++)
	{
		if (a & (0x1 << i))
			temp = temp | (0x1 << p_table[i]);
	}
	return temp;
}

block oracle_enc(block p)
{
	block c = p;
	for (int r = 0; r < 4; r++)
	{
		c = c ^ k[r];
		c = s_box(c);
		if (r != 3) c = p_box(c);
		else c = c ^ k[r + 1];
	}
	return c;
}

//ENDS - Encryption Oracle Implementation

//Oracle evaluation results of "E(x) XOR E(x XOR Delta_in) = Delta_fin"
//Delta_in is 0x40 and Delta_fin is { 0x0E, 0x0B }
//Each index in array is "x" and corresponding entry is oracle evaluation result
int oracle_e[256] = {	0x38,0x51,0x78,0xF2,0xB2,0x9E,0xD8,0xE,0xB,0xB8,0xE,0xA2,0x84,0xE8,0xE4,0x84,
						0x26,0x63,0x70,0xF2,0x43,0x8A,0x98,0xD1,0xEC,0xB0,0xF7,0x30,0xC6,0xB0,0x81,0x2B,
						0x32,0x33,0x18,0x51,0x8C,0x2E,0xD8,0xCA,0xED,0xE8,0x43,0x19,0x8B,0x25,0xE1,0x8B,
						0xDD,0xFD,0x3B,0x11,0xFB,0x63,0xB3,0x2B,0xD6,0x27,0x3D,0xBE,0x7B,0x5A,0xEB,0xDB,
						0x38,0x51,0x78,0xF2,0xB2,0x9E,0xD8,0xE,0xB,0xB8,0xE,0xA2,0x84,0xE8,0xE4,0x84,
						0x26,0x63,0x70,0xF2,0x43,0x8A,0x98,0xD1,0xEC,0xB0,0xF7,0x30,0xC6,0xB0,0x81,0x2B,
						0x32,0x33,0x18,0x51,0x8C,0x2E,0xD8,0xCA,0xED,0xE8,0x43,0x19,0x8B,0x25,0xE1,0x8B,
						0xDD,0xFD,0x3B,0x11,0xFB,0x63,0xB3,0x2B,0xD6,0x27,0x3D,0xBE,0x7B,0x5A,0xEB,0xDB,
						0x89,0x23,0x54,0x1C,0x64,0x82,0x71,0x21,0xCD,0xA3,0x35,0xF8,0x5E,0x3A,0x14,0xE4,
						0xEF,0x34,0xE0,0xBD,0x4C,0x1F,0x3B,0x55,0x88,0x8C,0x52,0x99,0x95,0x7C,0xBC,0x46,
						0x4E,0x88,0xFC,0x34,0xC6,0x2C,0x31,0x9D,0xD5,0x65,0x1D,0xDC,0x65,0x5F,0x92,0x95,
						0x2A,0x63,0x7B,0xFD,0x47,0x86,0x98,0xD3,0xEC,0xBB,0xF3,0x3C,0xCA,0xBD,0x8A,0x20,
						0x89,0x23,0x54,0x1C,0x64,0x82,0x71,0x21,0xCD,0xA3,0x35,0xF8,0x5E,0x3A,0x14,0xE4,
						0xEF,0x34,0xE0,0xBD,0x4C,0x1F,0x3B,0x55,0x88,0x8C,0x52,0x99,0x95,0x7C,0xBC,0x46,
						0x4E,0x88,0xFC,0x34,0xC6,0x2C,0x31,0x9D,0xD5,0x65,0x1D,0xDC,0x65,0x5F,0x92,0x95,
						0x2A,0x63,0x7B,0xFD,0x47,0x86,0x98,0xD3,0xEC,0xBB,0xF3,0x3C,0xCA,0xBD,0x8A,0x20
					};

int is_in(int *x, int n, int a)
{
	for(int i = 0; i < n; i++)
		if(x[i] == a) return 1;
	return 0;
}

double count()
{
	QuESTEnv env = createQuESTEnv();

    //Value of Pi
    const double pi = 3.14;

    const int p = 8, P = 256;
    const int n = 8, N = 256;
    Qureg qubits = createQureg((p + n), env);
    initZeroState(qubits);
    //reportQuregParams(qubits);
    //reportQuESTEnv(env);

    //Unitary matrix creation corresponding to the "oracle_e"
	ComplexMatrixN e = createComplexMatrixN(n);
    for(int i = 0; i < 256; i++)
    {
    	if(oracle_e[i] == 0xE || oracle_e[i] == 0xB) e.real[i][i] = -1;
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

//Grover's search to find "x" which satisfies "E(x) XOR E(x XOR Delta_in) = Delta_fin"
int search()
{
	QuESTEnv env = createQuESTEnv();

    Qureg qubits = createQureg(8, env);
    initZeroState(qubits);
    //reportQuregParams(qubits);
    //reportQuESTEnv(env);

    //Unitary matrix creation corresponding to the "oracle_e"
	ComplexMatrixN e = createComplexMatrixN(8);
    for(int i = 0; i < 256; i++)
    {
    	if(oracle_e[i] == 0xE || oracle_e[i] == 0xB) e.real[i][i] = -1;
    	else e.real[i][i] = 1;
    }

    //Apply Hadamard on all the n-qubits to create superposition
    //Here each element in superposition is considered to represent "x" 
    for(int i = 0; i < 8; i++)
        hadamard(qubits, i);
    
    int targs[] = { 0,1,2,3,4,5,6,7 };
    for(int gi = 0; gi < 5; gi++)
    {
    	//Marking
    	multiQubitUnitary(qubits, targs, 8, e);
    	
        //Diffusion
        for(int i = 0; i < 8; i++)
            hadamard(qubits, i);
        for(int i = 0; i < 8; i++)
            pauliX(qubits, i);
        multiControlledPhaseFlip(qubits, targs, 8);
        for(int i = 0; i < 8; i++)
            pauliX(qubits, i);
        for(int i = 0; i < 8; i++)
            hadamard(qubits, i);
    }

    qreal prob;
    int x = 0;
    for(int i = 0; i < 8; i++)
    {
    	int outcome = measureWithStats(qubits, i, &prob);
    	if(outcome) x ^= (outcome << i);
    }

    destroyQureg(qubits, env);
    destroyComplexMatrixN(e);
    destroyQuESTEnv(env);

    return x;
}

//Using Grover's search to find last round best subkey
int find_best_subkey(int *x, int x_n)
{
	//Oracle evaluation results for different subkeys
	//For each subkey (represented by index of the array), corresponding array element represents a number of "x" for which "E_R(x) XOR E_R(x XOR Delta_in) = Delta_out" holds
	//"E_R" is the encrption function only evaluated for first "R" rounds of the cipher
	//Delta_in is 0x40 and Delta_out is 0x0A
	//Each index in array represents a 4-bit "subkey"
	//The oracle elements is found by evaluating encryption function classically (can't be done as quantum simulation)
	int oracle_subkeys[16];
	for (int i = 0; i < 16; i++)
		oracle_subkeys[i] = 0;
	for (int i = 0; i < 16; i++) //Loop for each subkey
	{
		for (int j = 0; j < x_n; j += 2) //Loop for each data (ciphertext) last round backward operation
		{
			//Performing last round operation in reverse
			block c1 = oracle_enc(x[j]) ^ i, c2 = oracle_enc(x[j + 1]) ^ i; //Applying subkey "i" on ciphertext
			c1 = s_box_inv(c1);
			c2 = s_box_inv(c2);
			if ((c1 ^ c2) == 0x0A)
				oracle_subkeys[i]++;
		}
	}

	QuESTEnv env = createQuESTEnv();

    Qureg qubits = createQureg(4, env);
    initZeroState(qubits);

    //Unitary matrix creation corresponding to the "oracle_subkeys"
	ComplexMatrixN e = createComplexMatrixN(4);
    for(int i = 0; i < 16; i++)
    {
    	if(oracle_subkeys[i] == x_n / 2) e.real[i][i] = -1; //subkey that holds for all the "x"
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
int find_complete_key(block P, block C, int subkey)
{
	//Create oracle corresponding to remaining first 4 key bits
	//Each index in array represents the first 4-bit "subkey"
	//Oracle stores the entry "1" corresponding to the correct key and 0, otherwise
	int oracle_rem_bits[16];
	for (int i = 0; i < 16; i++) //Loop for each subkey
	{
			block temp_C = C;
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
    /*
    for(int i = 0; i < 4; i++)
    {
    	int outcome = measureWithStats(qubits, i, &prob);
    	if(outcome) f_subkey ^= (outcome << i);
    }
    */
    qreal max_prob = 0;
    for(int i = 0; i < 16; i++)
    {
    	prob = getProbAmp(qubits, i);
    	if(prob > max_prob)
    	{
    		max_prob = prob;
    		f_subkey = i;
    	}
    }

    destroyQureg(qubits, env);
    destroyComplexMatrixN(e);
    destroyQuESTEnv(env);

    return ((f_subkey << 4) ^ (subkey));
}

int main (int narg, char *varg[])
{
	//Count all "x" such that "E(x) XOR E(x XOR Delta_in) = Delta_fin" holds
	//This part can be SKIPPED as it is not proposed in algorithm 2
    int c = round(count());
    printf("\nThe number of x which satisfies \"E(x) XOR E(x XOR Delta_in) = { 0x0E, 0x0B }\" are %d", c);
    printf("\nFinding all such values of x!");
    
    //Find all "x" such that "E(x) XOR E(x XOR Delta_in) = Delta_fin" holds
    int *x = (int *)malloc(sizeof(int) * c);
    int x_size = 0;
    while(c)
    {
    	int a = search();
    	if(!is_in(x, x_size, a))
    	{
    		x[x_size++] = a;
    		x[x_size++] = a ^ 0x40; // a XORed with Delta_in (0x40)
    		c -= 2;
    	}
    }
    printf("\nThe possible values of x are:");
    for(int i = 0; i < x_size; i++) printf(" 0x%02X", x[i]);

    //Find the last round best subkey such that "E(x) XOR E(x XOR Delta_in) = Delta_out" holds
	int subkey = find_best_subkey(x, x_size);
	printf("\nThe best subkey for last 4-bits is 0x%02X", subkey);

	//Complete the last-round subkey using grover's search
	//To check for key correctness send plaintext and corresponding ciphertext
	int complete_key = find_complete_key(x[0], oracle_enc(x[0]), subkey);
	printf("\nThe complete key for last round is 0x%02X\n", complete_key);

    return 0;
}