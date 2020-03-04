/*
*************************************************
*   About: The code implements Algorithm 1      *
*   Usage: Run in command prompt "A1.exe"       *
*   Author: Ashwini Kumar Malviya               *
*   Email: ashwinixar@gmail.com                 *
*************************************************
*/

#include <stdio.h>
#include <stdlib.h>
#include "QuEST.h"

int main (int narg, char *varg[])
{
    //Differential used to attack is (Delta_in, Delta_out) = (0x40, 0x0A) i.e. E_R(x) XOR E_R(x XOR Delta_in) = Delta_out
    //E_R represents encryption results after R rounds of the cipher
    //For complete cipher, "Delta_out" produces "Delta_fin" possible differentials
    //Thus the oracle must search for all "Delta_fin" 
    //Oracle evaluation results of "E(x) XOR E(x XOR Delta_in) = Delta_fin"
    //Delta_in is 0x40 and Delta_fin is { 0x0E, 0x0B }
    //The differential has probability of 0.0234375. So, "6" pairs (out of 256) must hold
    //Hence, Grover's iteration must be repeated for 2^(h/2) = 2^(5.415/2) = 6.534 times or 6 times
    //Each index in array is "x" and corresponding entry is oracle evaluation result
    int oracle_e[256] = {   0x38,0x51,0x78,0xF2,0xB2,0x9E,0xD8,0xE,0xB,0xB8,0xE,0xA2,0x84,0xE8,0xE4,0x84,
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
                            0x2A,0x63,0x7B,0xFD,0x47,0x86,0x98,0xD3,0xEC,0xBB,0xF3,0x3C,0xCA,0xBD,0x8A,0x20 };

    QuESTEnv env = createQuESTEnv();

    Qureg qubits = createQureg(8, env);
    initZeroState(qubits);
    reportQuregParams(qubits);
    reportQuESTEnv(env);

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
    
    //According to Grover's search complexity
    //Ideal number of iterations for repeating grover's operator to get the correct result
    //times = (pi/4) * sqrt(2^n/m), here n = 8-qubits, m = 6 (correct results i.e. 4 entries in oracle_e has "Delta_out = 0xE")
    //So, times = (pi/4) * sqrt(2^n/m) = (3.14/4) * sqrt(256/6) = 5.127 or 5
    //One can choose either 6 iterations (according to attack) or 5 iterations according to Grover's time complexity
    int times = 5;
    printf("\nRunning Grover's operator %d times\n", times);

    int targs[] = { 0,1,2,3,4,5,6,7 };
    for(int gi = 0; gi < times; gi++)
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
    
    printf("\nCircuit output:\n");

    qreal prob;
    //Prints all the states along with its corresponding probability
    for(int i = 0; i < 256; i++)
    {
        prob = getProbAmp(qubits, i);
        printf("Probability amplitude of %02X: %f\n", i, prob);
    }

    //Measures the qubits to get the final result "x"
    int x = 0;
    for(int i = 0; i < 8; i++)
    {
        int outcome = measureWithStats(qubits, i, &prob);
        if(outcome) x = x ^ (0x1 << i);
        printf("Qubit %d collapsed to %d with probability %f\n", i, outcome, prob);
    }

    //Print the measured result
    printf("\nThe measured value of x is %02X\n", x);

    //Print the quantum differential distinguisher's result
    if(oracle_e[x] == 0xE || oracle_e[x] == 0xB) printf("Concrete\n");
    else printf("Random\n");

    destroyQureg(qubits, env);
    destroyComplexMatrixN(e);
    destroyQuESTEnv(env);

    return 0;
}