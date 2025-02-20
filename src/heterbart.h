/*
 *  BART: Bayesian Additive Regression Trees
 *  Copyright (C) 2017 Robert McCulloch and Rodney Sparapani
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/GPL-2
 */

// Define a header guard to prevent multiple inclusions of this file
#ifndef GUARD_heterbart_h // Check if the macro GUARD_heterbart_h is not defined
#define GUARD_heterbart_h // Define GUARD_heterbart_h to prevent future inclusions

// Include necessary headers for this file
#include "bart.h"          // Include the header file "bart.h"
#include "heterbartfuns.h" // Include the header file "heterbartfuns.h"
#include "heterbd.h"       // Include the header file "heterbd.h"

// Define the class "heterbart" which inherits from "bart"
class heterbart : public bart
{
public:
    // Default constructor for heterbart, which calls the default constructor of the base class 'bart'
    heterbart() : bart() {}

    // Constructor that takes a size_t m and passes it to the 'bart' constructor
    heterbart(size_t m) : bart(m) {}

    // Member function to print or display something related to 'heterbart'
    void pr();

    // Member function to draw something, probably a graphical representation
    // It takes a pointer to a double array (sigma) and a reference to a random number generator (gen)
    void draw(double *sigma, rn &gen);
};

// End of header guard
#endif
