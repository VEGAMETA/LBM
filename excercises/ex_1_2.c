#define _USE_MATH_DEFINES
#include <stdio.h>
#include <math.h>

const double CELSIUS_ZERO = 273.15;
const double BOLTSMAN = 1.380649 * 10e-23;
const double AVOGADRO = 6.0221408 * 10e23;
const double HYDROGEN_MOLAR_MASS = 3.344 * 10e-27;
const double HYDROGEN_DENSITY = 0.09;

double average_molecular_speed(double temp, double mass) {
    return sqrt((8 * BOLTSMAN * temp) / (M_PI * mass));
}

float root_mean_squared_molecular_speed(double temp, double mass) {
    return sqrt((3 * BOLTSMAN * temp) / mass);
}

void exercise_1() {
    printf("Exercise 1:\nCalculate the average and rms speed of hydrogen molecule at room temperature (25C)\n\n");
    printf("Average speed: %.2f m/s\n", average_molecular_speed(CELSIUS_ZERO + 25, HYDROGEN_MOLAR_MASS));
    printf("RMS speed: %.2f m/s\n\n", root_mean_squared_molecular_speed(CELSIUS_ZERO + 25, HYDROGEN_MOLAR_MASS));
}

double maxwell_boltzmann_distribution(double temp, double mass, double vel) {
    return 4 * M_PI * pow(mass / (2 * M_PI * BOLTSMAN * temp), 3 / 2) * pow(vel, 2) * exp((-mass * vel * vel)/(2 * BOLTSMAN * temp)); 
}

double number_of_molecules(double number_density, double mass, double temp) {
    return number_density * sqrt((BOLTSMAN * temp) / (2 * M_PI * mass));
}

double pressure(double number_density, double temp) {
    return number_density * BOLTSMAN * temp;
}

double get_number_density(double density, double mass){
    return density / mass;
}

void exercise_2() {
    double number_density = get_number_density(HYDROGEN_DENSITY, HYDROGEN_MOLAR_MASS);
    printf("Exercise 2:\nFor a Maxwellian distribution function, calculate the total number of molecules striking unit area of a wall per unit time. What will be the pressure on the wall?\n\n");
    printf("Number of molecules: %e\n", number_of_molecules(number_density, HYDROGEN_MOLAR_MASS, CELSIUS_ZERO + 25));
    printf("Pressure: %.2f Pa\n\n", pressure(number_density, CELSIUS_ZERO + 25));
}
