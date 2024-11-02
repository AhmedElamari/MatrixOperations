
// matlib24.cpp : This file contains the 'main' function and others to test matrix routines
//

// matlib24.cpp : A simple matrix library program.
// Template Written by Prof Richard Mitchell    6/12/23
// Adapted by:   << Ahmed Elamari >>

#include <string>
#include <iostream>
#include <cstdbool>
using namespace std;

const int maxM = 16;	// allows for matrices up to size 4*4

struct myMat {			// structure to store a matrix
	int numRows;		// number of rows
	int numCols;		// number of columns
	int data[maxM];		// data are stored in row order
};

myMat zeroMat(int r, int c) {
	// create a matrix with r rows and c columns, filled with zeros
	myMat m;			// define matrix
	m.numRows = r;		// set size
	m.numCols = c;
	for (int ct = 0; ct < maxM; ct++) m.data[ct] = 0;	// set elems to 0
	return m;			// return matrix
}

int getIndex(myMat m, int r, int c) {
	// returm index of element in row r, col c of matrix m
	return r * m.numCols + c;
}

myMat matError(string errstr) {
	// if this is called, an error has occured ... output errstr and return 0 size myMat
	cout << "Error : " << errstr << "\n";
	return zeroMat(0, 0);
}

int intError(string errstr) {
	// if this is called, an error has occured ... output errstr and return 0 
	cout << "Error : " << errstr << "\n";
	return 0;
}

myMat mFromStr(string s) {
	// create a matrix from string s
	// string of form "1,2;3,4;5,6"   ; separates rows and , separates columns ... No error check
	int ms;
	if (s.length() > 0) ms = 1; else ms = 0;
	myMat m = zeroMat(ms, ms);						// is s empty create 0*0 matrix, else start with 1*1 matrix
	int mndx = 0;									// used to index into array
	string sub = "";								// sub string - numbers between , or ; set empty
	for (int ct = 0; ct < s.length(); ct++) {		// each char in turn
		if ((s[ct] == ',') || (s[ct] == ';')) {	// if , or ; then sub contains number
			m.data[mndx++] = stoi(sub);				// convert number to integer, put in data array
			sub = "";								// clear sub string
			if (s[ct] == ';') m.numRows++;			// if found ; indicates an extra row
			else if (m.numRows == 1) m.numCols++;	// if , then (if in row 1) increase count of columns
		}
		else sub = sub + s[ct];						// add character to sub string
	}
	if (sub.length() > 0) m.data[mndx++] = stoi(sub);// add last sub string
	return m;
}

void printMat(const char* mess, myMat m) {
	// mess is string to be printed, followed by matrix m
	cout << mess << " = " << "\n";				// print message
	for (int r = 0; r < m.numRows; r++) {		// do each row
		for (int c = 0; c < m.numCols; c++)		// do each column
			cout << m.data[getIndex(m, r, c)] << "\t";	// outputing the element then tab
		cout << "\n";							// output new line at end of row
	}
	cout << "\n";								// and end of Matrix
}

myMat mGetRow(myMat m, int row) {
	// create a matrix from m, having one row
	myMat res = zeroMat(1, m.numCols);		// create a matrix of 1 row
	for (int col = 0; col < m.numCols; col++)		// for each element in row
		res.data[col] = m.data[getIndex(m, row, col)];		// copy col element to res
	return res;
}

myMat mGetCol(myMat m, int col) {
	// create a matrix from m, having one col
	myMat res = zeroMat(m.numRows, 1);	// adjust arguments so column vector
	for (int row = 0; row < m.numRows; row++) // for each element in the row
		res.data[row] = m.data[getIndex(m, row, col)]; // copies row element to res
	return res;
}

myMat mSetCol(myMat m, int col, myMat v) {
	// insert v into given col in m
	if (m.numRows != v.numRows) // error handling for numRows matched in vector multiplication / single row operations
		return matError("Matrix/Vector should have same number of rows");
	else {
		myMat res = m;
		for (int row = 0; row < m.numRows; row++)
			res.data[getIndex(m, row, col)] = v.data[row];
		return res;
	}
}

int dotProd(myMat v1, myMat v2) {
	if (v1.numRows != v2.numRows != v1.numCols != v2.numCols) {
		std::cout << "Error dimensions are incorrect" << std::endl; // error handling to ensure diemnsions match for dotprod
		return { 0 }; // returns 0 if dimensations are not compatible for dot product.
	}


	int res = 0;// setting counter for incrementing result of dotprod
	for (int i = 0; i < v1.numCols; i++) { //loops through vector elements
		res += v1.data[i] * v2.data[i]; //multiplies data of vector 1 and 2 and adds result to current res 
	}
	return res; // returns res of dotprod of v1 and v2
}


void testVecs(myMat A, myMat C) {
	// test vector routines ... get row from m1, col from m3, do dot product of these
	cout << "Testing Vector routines" << "\n";
	printMat("A row 0", mGetRow(A, 0));
	// display row 0 of m1
	printMat("C col 0", mGetCol(C, 0));
	printMat("C col 1", mGetCol(C, 1));
	printMat("C col 2", mGetCol(C, 2));

	// display col 1 of m2

	cout << "Dot prod of these is " << dotProd(mGetRow(A, 1), mGetCol(C, 1)) << "\n\n"; // tests dotprod
	cout << "Dot prod of A row 1 and C row 1 is " << dotProd(mGetRow(A, 0), mGetRow(C, 1)) << "\n\n"; //tests dotprod error
}


myMat mTranspose(myMat m) {
	// generate a new matrix which is the transpose of m
	myMat res = zeroMat(m.numCols, m.numRows);		// change arguments
	for (int r = 0; r < m.numRows; r++) // loops throguh row
		for (int c = 0; c < m.numCols; c++) // loops through col
			res.data[getIndex(m, c, r)] = m.data[getIndex(m, r, c)]; // transposes matrice
	return res;
}

myMat mAdd(myMat m1, myMat m2) {
	//the equiv elements in m1 and m2
	if (m1.numRows != m2.numRows || m1.numCols != m2.numCols) {	//error handling, matrix diemsnions match for addition
		std::cout << "M1 and m2 are not equivlent" << std::endl;
	}

	myMat res = zeroMat(m1.numRows, m1.numCols);		// change arguments
	for (int r = 0; r < m1.numRows; r++) // loops through row
		for (int c = 0; c < m1.numCols; c++) //loops through col
			res.data[getIndex(m1, r, c)] = m1.data[getIndex(m1, r, c)] + m2.data[getIndex(m2, r, c)]; //adds the index of current loop in m1 and m2 together and stores the current index in specific row and col in res
	return res; // returns res
}
/**
* Performs scalar multiplication or divison on a matrix
*
* @param m the input matrix to be modified via scalar
* @param s, scalar value which to multiply/divise each element of matrix
* @param isMult flag to indiciate operation type: 1 for multiplication, 0 for divison.
*/
myMat Scalar(myMat m, double s, int isMult) {
	// multiply or divide all elements in m by s
	if (s == NULL) { // check if scalar value is provided
		cout << "Input value for scalar multiplication/division" << std::endl;
		return { 0 }; // ensure there is a scalar input. // return a zero matrix as indiciation of error
	}
	// check for divison by zero to avoid undefined behaviour
	if (s == 0) {
		cout << "Cannot divide by zero" << endl;
		return zeroMat(0, 0); // return zero matrix with no dimensions
	}
	// create new matrix to store result
	myMat res = zeroMat(m.numRows, m.numCols);
	for (int i = 0; i < m.numRows; i++) { // iterate over each element in matrix
		for (int j = 0; j < m.numCols; j++) {
			if (isMult == 1) { // determines operation based on flag
				// multiply each element of matrix by scalar value
				res.data[getIndex(res, i, j)] = m.data[getIndex(m, i, j)] * s;
			}
			else {
				// divise each element of matrix by scalar value
				res.data[getIndex(res, i, j)] = m.data[getIndex(m, i, j)] / s;
			}

		}
	}
	// return result matrix after scalar operation.
	return res;
}

myMat mMult(myMat A, myMat B) {
	// generate a matrix which is the product of m1 and m2
	myMat res = zeroMat(A.numRows, B.numCols);
	if (A.numCols == B.numRows) { // if the Numcols and numrows are equivalent, or numCols of B and rows of A is equiv, than multiply each 
		for (int r = 0; r < A.numRows; r++) { // for each element in row
			for (int c = 0; c < B.numCols; c++) { // for each element in col
				res.data[getIndex(res, r, c)] = dotProd(mGetRow(A, r), mGetCol(B, c)); // multiplies the row and col and stores the result in res
			} // end of col loop
		} // end of row loop
	} 
	else {
		std::cout << "The column of matrix 1 and matrix 2 row must match." << std::endl; // error handling for matrix multiplication
		return zeroMat(0, 0); // returns zero matrix if dimensions do not match
	}
	return res; // returns res
}

void testMatOps(myMat A, myMat B, myMat C) {
	// test matrix operations m1 + m2; m1 + m3 (not poss) m1 + m3'; m1 * m3; m3 * m2
	cout << "Testing Add, Transpose, Multiply routines" << "\n";
	printMat("m1+m2", mAdd(A, B));
	printMat("m1 + m3", mAdd(A, C));
	printMat("m1 + m3'", mAdd(A, mTranspose(C)));
	printMat("m1*m2", mMult(A, B));
	printMat("m2*m1", mMult(C, A));
	printMat("m1*m3", mMult(A, B));

}


//  commented out code you will do as part of Summative Assessment later in term

myMat mSubMat(myMat m, int row, int col) {
	// return matrix m but omitting specified row and column
	// if time add code to check matrices of right size
	myMat res = zeroMat(m.numRows - 1, m.numCols - 1);		// change arguments
	// write code to do sub mat
	int Index = 0; // tracks matrix in flattened to increment through new matrix.
	for (int i = 0; i < m.numRows; i++) { // loops through row in original matrix
		if (i == row) continue; // skips the row we want to exclude
		for (int j = 0; j < m.numCols; j++) { // loops through col in original matrix
			if (j == col) continue; // skips col we want to exclude
			res.data[Index++] = m.data[getIndex(m, i, j)]; //grabs new matrix and copies the elements into res and then increments on index++
		} // end of col loop
	} // end of row loop
	return res; // returns res
}

int mDet(myMat m) {
	// compute determinant of matrix m
	if (m.numRows != m.numCols) { // error handling for determinant, matrix must be square
		cout << "Matrix must be square" << endl;
		return 0; // returns 0 if matrix is not square
	}
	if (m.numRows == 1) // if the number of rows is 1, then return the only element in the matrix
		return m.data[getIndex(m, 0, 0)]; // returns the only element in the matrix
	else if (m.numRows == 2) // However if the number of rows is 2, then return the determinant of the 2x2 matrix
		return m.data[getIndex(m, 0, 0)] * m.data[getIndex(m, 1, 1)] - m.data[getIndex(m, 0, 1)] * m.data[getIndex(m, 1, 0)];  //ad- bc logic implemented to return 2x2 determinant
	else { // recursive call to find determinant of 3x3 or 4x4 matrix
		int sum = 0; // setting sum to 0 for incrementing
		for (int col = 0; col < m.numCols; col++) { //loops through the columns
			int plusOrMinus = (col % 2 == 0) ? 1 : -1; // sets the sign in order to check whether to add or subtract by a matrix
			sum += m.data[getIndex(m, 0, col)] * plusOrMinus * mDet(mSubMat(m, 0, col)); // recursively calls mSubMat to find determinant of 3x3 or 4x4 matrix implementation theory of working out determinant
		}
		return sum;// returns the sum of the determinant
	} // end of else
}
myMat mAdj(myMat m) {
	// return adjoint of matrix m
	if (m.numRows != m.numCols) {
		return matError("Adjoint can only be found for square matrices."); // error handling for adjoint, matrix must be square
	} 
	myMat res = zeroMat(m.numRows, m.numCols); // create a matrix of same size as input matrix
	if (m.numRows == 1) { // if the number of rows is 1, then return the only element in the matrix
		res.data[0] = 1; // returns the only element in the matrix
		return res; // returns res
	}
	else if (m.numRows == 2) {
		res.data[getIndex(res, 0, 0)] = m.data[getIndex(m, 1, 1)]; // swaps the elements in the matrix to find adjoint
		res.data[getIndex(res, 1, 0)] = -m.data[getIndex(m, 1, 0)]; // swaps the elements in the matrix to find adjoint and subtracting element before swap
		res.data[getIndex(res, 0, 1)] = -m.data[getIndex(m, 0, 1)];// swaps the elements in the matrix and subtracting element before swap
		res.data[getIndex(res, 1, 1)] = m.data[getIndex(m, 0, 0)]; // swaps the elements in the matrix to find adjoint
		return res; // returns res
	}
	else { // recursive call to find adjoint of 3x3 or 4x4 matrix
		for (int row = 0; row < m.numRows; row++) {  // loops through the rows
			for (int col = 0; col < m.numCols; col ++) {  // loops through the columns
				int PlusOrMinus = ((row + col) % 2 == 0) ? 1 : -1;  // sets the sign in order to check whether to add or subtract by a matrix
				res.data[getIndex(res, col, row)] = PlusOrMinus * mDet(mSubMat(m, row, col)); // recursively calls mSubMat to find determinant of 3x3 or 4x4 matrix implementation theory of working out adjoint
		} // end of col loop
	} // end of row loop
	return res; // returns res
	} // end of else
}


void testMatEqn(myMat A, myMat b) {
	// solve Ax = b using Cramer  and using Adjoint
	if (A.numRows != A.numCols || A.numRows != b.numRows) { // error handling for matrix equation, matrix must be square and b must be a column vector
		cout << "Matrix A must be square and b must be a column vector" << endl; // error handling for matrix equation
		return; // returns if matrix is not square or b is not a column vector
	}
	
	//Adoint rule:
	myMat x = zeroMat(A.numRows, 1); // create a matrix of same size as input matrix
	x = Scalar(mMult(mAdj(A), b), mDet(A), 0); // multiply adjoint of A by b and divise by determinant of A
	if (mDet(A) == 0) { // error handling for matrix if mDet is 0 then matrix is singular
		cout << "Matrix A is singular" << endl;
		return; // returns if matrix is singular

	}
	
	   //crammers rule
	 else {
		   myMat x = zeroMat(A.numRows, 1); // create a matrix of same size as input matrix

		   if (mDet(A) == 0) {
			   cout << "Matrix A is singular" << endl;// error handling for matrix if mDet is 0 then matrix is singular
			   return;
		   }
		   for (int col = 0; col < A.numCols; col++) { // loops through the columns
			   myMat Dx = A; // create a matrix of same size as input matrix
			   Dx = mSetCol(Dx, col, b); // set the column of A to b
			   x.data[col] = mDet(Dx) / mDet(A); // divise determinant of Ax by determinant of A
		   }
		   printMat("x using Cramer's rule is = ", x); // prints the result of the matrix equation
	  }
	  
}


int main()
{
	cout << "32013680\n";	// change to your student number
	myMat A, B, C, D;						// create  matrices

	A = mFromStr("2,7,4;9,8,10;2,10,7");			// change numbers to those in A from Q1 on web page, as specified on the sheet
	B = mFromStr("99;217;144");			// ditto but B
	C = mFromStr("1,1;2,2;3,3");
	D = mFromStr("2,4;3,7");
	printMat("submatrix: ", mSubMat(A, 0, 0));
	printMat("Adj", mAdj(A)); // test adjoint
	cout << "Det is " << mDet(A) << endl; // test determinant
	testMatEqn(A, B); // test matrix equation
	

	// test the vector routines
//	testMatOps(m1, m2, m3);					// test the add, transpose and multiply
	return 0;
}
