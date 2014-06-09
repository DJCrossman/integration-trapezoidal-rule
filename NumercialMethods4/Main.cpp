/*
 *	Main.cpp
 *  
 *  A program that will numerically integrate this function from -10 to +10
 *  using multiple segment trapezoidal rule i.e evaluate the following integral.
 *  f(x) = 2 - 5x + 10x^2 + 1/2*x^3
 *
 *  Created on: April 10, 2014
 *  Author: David Crossman
 */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

using namespace std;

double method(double);
double methodIntegrated(double, double);
double calculateTrapezoidalRule(double, double, double);

int main () {
	double error;		// true percent relative error
	double n;			// number of segments
	double a = -10;		// lower bound
	double b = 10;		// upper bound

	cout << "Please enter the number of segments" << endl;
	cin >> n;

	ofstream outFile("../test/data.csv");
	
	// Trapezoid Rule
	outFile << n << " Segment Trapezoidal Rule" << endl;
	outFile << "n,|e|%" << endl;

	for (double i = 1; i <= n; i++) {
		error = calculateTrapezoidalRule(i, a, b);
		outFile << i << "," << error << endl; 
	}

	outFile.close();

	return 0;
}

double method(double x) {

	return 2 - 5*x + 10*pow(x, 2) + 0.5*pow(x, 3);
}

double methodIntegrated(double a, double b) {

	return (2*b - 5*pow(b, 2)/2 + 10*pow(b, 3)/3 + pow(b, 4)/8) - (2*a - 5*pow(a, 2)/2 + 10*pow(a, 3)/3 + pow(a, 4)/8);
}

double calculateTrapezoidalRule(double n, double a, double b) {

	double h = (b - a) / n;
	double I = method(a) + method(b);
	for (double i = 1; i <= n-1; i++) {
		I += 2*method(a + i*h);
	}
	return abs((methodIntegrated(a,b) - I*h*0.5) / (methodIntegrated(a,b))) * 100;
}