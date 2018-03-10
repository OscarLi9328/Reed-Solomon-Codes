/**
 *  COM5140 Error-Correcting Codes
 *  Project1 ImpNementation of Viterbi Decoding
 *  ID  : 101011246
 *  NAME: Ming-Ju, Li
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#define N 63
#define r 21
#define K 42 
#define ERA_NUM -1
#define BUFF_SIZE 1000 

unsigned long long SEED;
unsigned long long RANV;
int RANI = 0;
const int alpha[N] = { 1, 2, 4, 8, 16, 32, 3, 6, 12, 24, 48, 35, 5, 10, 20, 40, 19, 
						38,	15, 30, 60, 59, 53, 41, 17, 34, 7, 14, 28, 56, 51, 37, 9, 
						18, 36, 11, 22, 44, 27, 54, 47, 29, 58, 55, 45, 25, 50, 39, 13, 
						26, 52, 43, 21, 42, 23, 46, 31, 62, 63, 61, 57, 49, 33};
int xr[N] = {0};
int gen[22] = {58, 62, 59, 7, 35, 58, 63, 47, 51, 6 ,33 ,43 ,44 ,27, 7, 53, 39, 62, 52, 41, 44, 1};

int pdeg(int *a);
int multi(int a, int b);
int divd(int a, int b);
int poly_val(int x, int degree, int *polynomial);
void Euclid(int xr[], int S[], int vj[], int rj[], int u, int v);
void poly_add();
void poly_multi(int *poly0, int *poly1, int *poly2, int deg0, int deg1); 
void poly_div(int *poly0, int *poly1, int *q, int *rr);
void poly_sub(int *pola0, int *poly1, int *sub, int degree); //subtract polynomial 0 with polynomial 1
void modulo(int *poly, int xr[]); //perform modulo to x^r

int main(void){
	FILE* f; // file to be read
	char strin[BUFF_SIZE] = " "; //buffer to store the received data
	int i, j, k; 
	int u, v; //degree of vj(x) and rj(x), respectively
	double e0 = 0; // the number of the erasure errors
	int R[N], R1[N] = {0}; //received vector and the modified received vector
	int S[N + 1] = {0}; //syndrome polynomial; degree: r
	int S0[N + 1] = {0}; //modified syndrome polynomial; degree: r
	int loc[2] = {0}; //(1 - alpha^i x)
	int t2 = 0;
	int E[N] = {0}; //error pattern
	char *temp = (char*)malloc(2 * sizeof(char)); 
	int *sigma0 = (int*)calloc(N, sizeof(int)); //erasure error locator polynomial: degree: e0 + 1
	int *sigma1 = (int*)calloc(N, sizeof(int)); //error locator polynomial: degree: e1 + 1
	int *sigma = (int*)calloc(N, sizeof(int)); //error-and-erasure locator polynomial; degree: e0 * e1
	int *sigmap = (int*)calloc(N, sizeof(int)); //derivative of error-and-erasure locator polynomial.
	int temp1;
	int vj[N + 1] = {0}, rj[N + 1] = {0}; //degree of vj is less than u, rj is less than v
	int count = 0;
	int check1 = 0, check2 = 0, decode = 0;
	int cw[N + 1] = {0}; //the decoded codeword
	int degree = 1;
	
	//Read the input file
	f = fopen("in.txt","r");
    if (f == NULL) {
		fputs ("File error\n", stderr); 
		system("pause");
		return 0;
	}
	
	for(i = 0; i < BUFF_SIZE; i++){
		fscanf(f, "%c", &strin[i]);
		printf("%c", strin[i]);
	} 
	printf("\n\nR = \n");
    int eraflag = 0;
    char *pch = (char*)malloc(2 * sizeof(char));
    pch = strtok(strin, " ");
    i = 0;
    while(pch != NULL){
        if (!strcmp(pch, "*")){  // step 0: locate the erasure error, using -1 represent it.
            R[i] = ERA_NUM;
            eraflag = 1;
        }
        else
            R[i] = atoi(pch);
            
        pch = strtok(NULL, " ");
        i++;
    }
    for(i = 0; i < N; i++) printf("%d ", R[i]);
    
	fclose(f);
		
//===============================================//
//                    Decoding                   //
//===============================================//
	
	//step 1: store the locations of the erasures: "erasure locator polynomial"
	//let R1 be the modified received vector
	for(i = 0; i < N; i++){
		if(R[i] == ERA_NUM){
			R1[i] = 0;
			e0++;
		}
		else{
			R1[i] = R[i];
		}
		printf("%d ", R1[i]);
	}
	u = (r - e0) / 2;
	v = (int)ceil((r + e0) / 2) - 1;
	if(e0 / 2 >= 21){
		printf("Exceed decoding capability!\n");
		return 0;
	}
	
	//calculate the erasure locator polynomial Sigma0(x):
	sigma0[0] = 1; //otherwise, the polynomial cannot be multiplied
	for(i = 0; i < N; i++){
		if(R[i] == ERA_NUM){
			loc[0] = 1;
			loc[1] = alpha[i];
			poly_multi(sigma0, loc, sigma0, pdeg(sigma0), 1);
		}
	}
	
	//step 2: calculate the syndromes: 
	for(j = 1; j <= r; j++){
		for(i = 0; i < N; i++){
			S[j - 1] ^= multi(R1[i], alpha[(i * j) % N]);
			
		}
	}
	
	for (i = 0; i < r; i++)
       t2 += S[i];
    if (eraflag == 0 && t2 == 0){
        puts("No error.");
        printf("\ncodeword: \n");
        for(i = 0; i < N; i++){
        	cw[i] = R1[i];
        	printf("%d ", cw[i]); 
		}
        return 0;
    }   
	//step 3: compute S0(x) = sigma0(x) * S(x) mod (x^r)
	xr[r] = 1;
	poly_multi(sigma0, S, S0, e0, r);
	printf("\n\nS0 = \n");
	for(i = 0; i < e0 + r; i++) printf("%d ", S0[i]);
		
	//step 4: perform Euclid's algorithm to get sigma1(x) = vj(x) / vj(0), w(x) = rj(x) / vj(0)
	Euclid(xr, S0, vj, rj, u, v); //vj(x) is sigma1(x)
	
	//Then, sigma(x) is the product of sigma0(x) and sigma1(x), i.e., sigma(x) = sigma0(x) * sigma1(x)
	poly_multi(vj, sigma0, sigma, pdeg(vj), e0 + 1);
	
	//Time-Domain Completion
	if(sigma[0] == 0 || pdeg(rj) >= e0 + pdeg(vj)){
		printf("\nDecode fail. No. 1\n");
		system("pause");
		return 0;
	}
	else{
		count = 0;
		for(i = 0; i < N; i++){
			//sigma(alpha^-i)
			check1 = poly_val(alpha[(N - i) % N], pdeg(sigma), sigma);

			//sigma'(alpha^-i)
			for(j = 0; j <= pdeg(sigma); j++){
				if(j % 2 == 1){
					sigmap[j - 1] = sigma[j];
				}
				else sigmap[j - 1] = 0;
			}
			check2 = poly_val(alpha[(N - i) % N], pdeg(sigmap), sigmap);
			
			if(check1 == 0 && check2 != 0){
				temp1 = poly_val(alpha[(N - i) % N], pdeg(rj), rj);
				E[i] = divd(temp1, check2);
				count++;
				
			} 
			else{ 
				E[i] = 0;
			}
		}
		if(count == pdeg(sigma)){
			printf("\ncodeword is:\n");
			for(i = 0; i < N; i++){
				cw[i] = R1[i] ^ E[i];
				printf("%d ", cw[i]);
			}
		}
		else{ //DECODE FAIL
			printf("\nDecode fail. No. 2\n");
			system("pause");
			return 0;
		}
	}
}

void Euclid(int xr[], int S[], int vj[], int rj[], int u, int v){
	int *qi = (int*)calloc(N, sizeof(int));
	int *temp = (int*)calloc(N, sizeof(int));
	int *temp1 = (int*)calloc(N, sizeof(int));
	int *temp2 = (int*)calloc(N, sizeof(int));
	int *rj_0 = (int*)calloc(N, sizeof(int));
	int *rj_1 = (int*)calloc(N, sizeof(int));
	int *vj_0 = (int*)calloc(N, sizeof(int));
	int *vj_1 = (int*)calloc(N, sizeof(int));
	int i, j;
	
	vj_1[0] = 1; //v-1 = 0, v0 = 1;
	
	for(i = 0; i < N; i++){ //r-1 = x^r, r0 = S(x);
		rj_0[i] = xr[i]; 
	}
	for(i = 0; i < N; i++){ //r-1 = x^r, r0 = S(x);
		rj_1[i] = S[i];
	}
	
	if(pdeg(rj_1) > v){ //check if Euclid's is needed
		while(1){
			for(i = 0; i < N; i++){ 
				qi[i] = 0;
				temp[i] = 0;
				temp1[i] = 0;
				temp2[i] = 0;
			}
			
			poly_div(rj_0, rj_1 , qi, temp); //qi = r(i-2) / r(i-1)
			
			for(i = 0; i < N; i++){ //rj is the remainder of r(j-2) to r(j-1)
				rj_0[i] = rj_1[i];
				rj_1[i] = 0;
				rj_1[i] = temp[i];
			}	 
	
			//v(i) = v(i-2) - v(i*-1) * qi;
			poly_multi(qi, vj_1, temp1, pdeg(qi), pdeg(vj_1)); //temp1 = v(i-1) * qi;
			poly_sub(vj_0, temp1, temp2, pdeg(temp1)); //temp2 = v(i-2) - v(i-1) * qi;
			
			//shift to the next step
			for(i = 0; i < N; i++){
				vj_0[i] = vj_1[i];
				vj_1[i] = temp2[i];	
				
			}
			
			if(pdeg(vj_1) <= u && pdeg(rj_1) <= v){ //terminate when degree of vj <= u & rj <= v
				printf("\nFinish!!\n");
				break;
			}
	
		}
		
	}
	else{ //if not, ri = x^r, vi = 1;
		for (i = 0; i < N; i++)
            rj[i] = rj_1[i];
        	vj[0] = 1;
	}
	
	if(vj_1[0] != 0){
		for(i = 0; i < N; i++){
			rj[i] = divd(rj_1[i], vj_1[0]); //w(x) = rj(x) / vj(0);
			vj[i] = divd(vj_1[i], vj_1[0]); //sigma1(x) = vj(x) / vj(0);
		}
	}
		
	free(temp1);
	free(temp);
	free(rj_0);
	free(rj_1);
	free(vj_0);
	free(vj_1);
	free(qi); 
}

//poly2 = poly1 * poly0, with deg(poly0) = deg0, deg(poly1) = deg1
void poly_multi(int *poly0, int *poly1, int *poly2, int deg0, int deg1){  //work
	int i, j;
	int temp[N] = {0};
	
	for(i = 0; i <= deg0; i++){
		for(j = 0; j <= deg1; j++){
			if(i + j <= r){
				temp[i + j] ^= multi(poly0[i], poly1[j]);
			}
		}
	} 
		
	for(i = 0; i <= deg0 + deg1; i++){
		poly2[i] = temp[i];
	}
}

//poly0 = poly1 * q + rr
void poly_div(int *poly0, int *poly1, int *q, int *rr){ //work
	int deg0 = pdeg(poly0), deg1 = pdeg(poly1), deg2; 
	int *temp1, *temp;
	int q_temp[N] = {0}; //the coefficient of the highest degree of the quotient
	int i, j, k, l;
	
	if(deg1 == -1){ //meaningless to divide by a constant
		printf("\nInvalid operation!");
		return;
	}
	 
	temp1 = (int*)calloc(N, sizeof(int));
	temp = (int*)calloc(N, sizeof(int)); //literally the remainder of each division operation.
	
	for(k = 0; k < N; k++){
		temp[k] = poly0[k]; //copy the dividend to a temperary array
	}
	
	while(pdeg(temp) >= deg1){
		deg2 = pdeg(temp) - deg1; //degree of the quotient
		q_temp[deg2] = divd(temp[pdeg(temp)], poly1[deg1]); //the coefficient of the quotient
		
		poly_multi(q_temp, poly1, temp1, deg2, deg1); //multiply the coefficient to polynomial 1
		poly_sub(temp, temp1, temp, pdeg(temp1)); //subtract the polynomial to the dividend
		
		for(i = 0; i <= deg0; i++){
			temp1[i] = 0; 
			q[i] ^= q_temp[i]; //copy the temporary quotient to the output quotient
			q_temp[i] = 0;
		} 
	}
	
	for(l = 0; l < deg0; l++){
		rr[l] = temp[l]; //copy the remainder
	}
		
	free(temp);
	free(temp1);
}

//polynomial subtraction
void poly_sub(int *poly0, int *poly1, int *sub, int degree){ //work
	int i;
	int *temp0 = (int*)calloc(N,  sizeof(int)), *temp1 = (int*)calloc(N, sizeof(int)); 

	for(i = 0; i < N; i++){ //copy the two polynomials to temperary array
		temp0[i] = poly0[i];
		temp1[i] = poly1[i];
	}
	
	for(i = 0; i <= degree; i++){
		sub[i] = (temp0[i] ^ temp1[i] ); //the negative number won't be transformed into positive automatically. So I add 63 to make it positive
	}
	
	free(temp0);
}

//GF(64) multiplication 
int multi(int a, int b){ //work
	int i, j;

	if(a * b == 0) return 0;
	else{
		for(i = 0; i < N; i++)
			if(alpha[i] == a) break;
		for(j = 0; j < N; j++)
			if(alpha[j] == b) break;
	}	
	
	return alpha[(i + j) % N];
}

//GF(64) division
int divd(int a, int b){ //work
	int i, j;
	
	if(a == 0 && b != 0) return 0;
	else if(b == 0) return -999;
	
	for(i = 0; i < N; i++)
		if(alpha[i] == a) break;
	for(j = 0; j < N; j++)
		if(alpha[j] == b) break;
	
	return alpha[((i - j) + N) % N];
}

int poly_val(int x, int degree, int *polynomial){
	int i, j;
	int y = 1, z = polynomial[0];
	
	for(i = 1; i <= degree; i++){
		if(polynomial[i] != 0){
			for(j = 1; j <= i; j++){ //degree
			y = multi(y, x);	
			}
			y = multi(y, polynomial[i]);
			z ^= y;
			y = 1;
		}
		else{
			z ^= 0;
		} 
	}	
	return z;	
}

int pdeg(int *a){ //work
	int i;
	int deg;
	
	for(i = N; i >= 0; i--){
		if(a[i] != 0){
			return i;
		}
	}
	
	return -1;
}
